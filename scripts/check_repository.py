#!/usr/bin/env python3
"""Read-only checks of the shipped manuscript repository (standard library).

This checks file identity, basic citation fields and local documentation links.
It never runs a compiler, a research solver, CV, or data generation, and writes
no files. Citation checks are not YAML parsing or CFF schema validation.
"""
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlsplit
import hashlib
import json
import os
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = "SHA256SUMS.json"
SOURCE = "provenance/source-v016.json"
PDFS = ("manuscript/main.pdf", "manuscript/ESM_1.pdf", "manuscript/complete.pdf")
SHIPPED_AUX = {"manuscript/xrefs_main.aux", "manuscript/xrefs_supplement.aux"}
AUX_SUFFIXES = {".aux", ".log", ".out", ".toc", ".blg", ".fls", ".fdb_latexmk",
                ".lof", ".lot", ".nav", ".snm", ".vrb", ".idx", ".ilg", ".ind"}
SHA_RE = re.compile(r"^[0-9a-f]{64}$")


def forbidden(relative):
    parts = PurePosixPath(relative).parts
    return (any(p.lower() in {"submission", ".venv", ".library"} for p in parts)
            or any(p.lower().startswith("titlepage.") for p in parts)
            or PurePosixPath(relative).suffix.lower() == ".rds")


def ignored(relative):
    path = PurePosixPath(relative)
    if any(p in {".git", "__pycache__"} for p in path.parts):
        return True
    if path.name == ".DS_Store" or path.suffix.lower() in {".pyc", ".pyo"}:
        return True
    # Local credential/private environment files are excluded by .gitignore.
    # Never open them; a manifest entry for one is rejected separately.
    if path.name in {".env", ".Rhistory", ".RData"} or path.name.startswith(".env."):
        return True
    if path.suffix.lower() in {".key", ".pem"}:
        return True
    if relative in SHIPPED_AUX:
        return False
    if path.parts and path.parts[0] == "manuscript":
        return ("build_logs" in path.parts or path.suffix.lower() in AUX_SUFFIXES
                or path.name.endswith((".synctex.gz", ".run.xml"))
                or relative == "manuscript/supplement.pdf")
    return False


def safe_relative(value):
    if not isinstance(value, str) or not value or "\\" in value:
        return False
    path = PurePosixPath(value)
    return (not path.is_absolute() and ".." not in path.parts
            and path.as_posix() == value and value != ".")


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def scrub_markdown(text):
    """Remove code blocks/inline code before checking documentation links."""
    lines, fence = [], None
    for line in text.splitlines():
        match = re.match(r"^\s{0,3}(`{3,}|~{3,})", line)
        if match:
            marker = match.group(1)
            if fence is None:
                fence = marker
            elif marker[0] == fence[0] and len(marker) >= len(fence):
                fence = None
            lines.append("")
        elif fence is None:
            lines.append(re.sub(r"(`+).*?\1", "", line))
        else:
            lines.append("")
    return "\n".join(lines)


def markdown_targets(text):
    """Extract inline/image and reference-definition link destinations."""
    text = scrub_markdown(text)
    for match in re.finditer(r"!?\[[^\]\n]*\]\(\s*", text):
        start = match.end()
        if start < len(text) and text[start] == "<":
            end = text.find(">", start + 1)
            if end != -1:
                yield text[start + 1:end]
            continue
        index, nesting = start, 0
        while index < len(text):
            char = text[index]
            if char == "\\" and index + 1 < len(text):
                index += 2
                continue
            if char.isspace() and nesting == 0:
                break
            if char == "(":
                nesting += 1
            elif char == ")":
                if nesting == 0:
                    break
                nesting -= 1
            index += 1
        if index > start:
            yield re.sub(r"\\([() ])", r"\1", text[start:index])
    for match in re.finditer(r"^\s{0,3}\[[^\]\n]+\]:\s*(?:<([^>]+)>|(\S+))", text, re.M):
        yield match.group(1) or match.group(2)


def main():
    failures, counts = [], {}

    def fail(code, path=None):
        item = {"check": code}
        if path is not None:
            item["path"] = path
        if item not in failures:
            failures.append(item)

    def read_json(relative):
        try:
            return json.loads((ROOT / relative).read_text(encoding="utf-8"))
        except (OSError, UnicodeError, ValueError):
            fail("missing_or_invalid_json", relative)
            return None

    # Candidate inventory is the reviewed filesystem minus documented local
    # build/cache outputs. Private directories are rejected without opening them.
    candidates = set()
    for directory, directories, filenames in os.walk(ROOT, followlinks=False):
        base = Path(directory)
        keep = []
        for name in directories:
            path = base / name
            relative = path.relative_to(ROOT).as_posix()
            if ignored(relative):
                continue
            if forbidden(relative):
                fail("private_directory_in_candidate", relative)
            elif path.is_symlink():
                fail("symlink_in_candidate", relative)
            else:
                keep.append(name)
        directories[:] = keep
        for name in filenames:
            path = base / name
            relative = path.relative_to(ROOT).as_posix()
            if forbidden(relative):
                fail("private_file_in_candidate", relative)
            elif ignored(relative):
                continue
            elif path.is_symlink():
                fail("symlink_in_candidate", relative)
            elif relative != MANIFEST:
                candidates.add(relative)
    counts["candidate_files"] = len(candidates)

    manifest = read_json(MANIFEST)
    entries = manifest.get("files") if isinstance(manifest, dict) else None
    if not isinstance(entries, list) or manifest.get("algorithm") != "SHA256":
        fail("manifest_schema")
        entries = []
    expected, hash_cache = {}, {}
    for entry in entries:
        if not isinstance(entry, dict) or not safe_relative(entry.get("path")):
            fail("invalid_manifest_entry")
            continue
        relative = entry["path"]
        if relative in expected:
            fail("duplicate_manifest_path", relative)
        expected[relative] = entry
        if relative == MANIFEST or forbidden(relative) or ignored(relative):
            fail("excluded_path_in_manifest", relative)
            continue
        if not isinstance(entry.get("sha256"), str) or not SHA_RE.fullmatch(entry["sha256"]):
            fail("invalid_sha256", relative)
            continue
        size = entry.get("bytes")
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            fail("invalid_byte_count", relative)
            continue
        path = ROOT / relative
        if relative not in candidates:
            fail("manifest_file_missing_or_not_candidate", relative)
            continue
        try:
            hash_cache[relative] = digest(path)
            if hash_cache[relative] != entry["sha256"]:
                fail("manifest_sha256_mismatch", relative)
            if path.stat().st_size != size:
                fail("manifest_size_mismatch", relative)
        except OSError:
            fail("file_unreadable", relative)
    for relative in sorted(candidates - set(expected)):
        fail("candidate_missing_from_manifest", relative)
    counts["manifest_entries"] = len(entries)
    counts["manifest_files_hashed"] = len(hash_cache)

    source = read_json(SOURCE)
    copied = source.get("copied_files") if isinstance(source, dict) else None
    if not isinstance(copied, list):
        fail("source_manifest_schema")
        copied = []
    source_seen = set()
    for entry in copied:
        if not isinstance(entry, dict) or not safe_relative(entry.get("path")):
            fail("invalid_source_entry")
            continue
        relative = entry["path"]
        if relative in source_seen:
            fail("duplicate_source_path", relative)
        source_seen.add(relative)
        if forbidden(relative) or ignored(relative) or relative not in candidates:
            fail("source_file_missing_or_excluded", relative)
            continue
        expected_hash = entry.get("sha256")
        if not isinstance(expected_hash, str) or not SHA_RE.fullmatch(expected_hash):
            fail("invalid_source_sha256", relative)
            continue
        try:
            actual = hash_cache.get(relative) or digest(ROOT / relative)
            if actual != expected_hash:
                fail("copied_source_sha256_mismatch", relative)
        except OSError:
            fail("source_file_unreadable", relative)
    counts["copied_source_files"] = len(copied)
    for relative in PDFS:
        if relative not in candidates or (ROOT / relative).stat().st_size == 0:
            fail("required_pdf_missing_or_empty", relative)
    counts["required_pdfs"] = len(PDFS)

    citation_required = ("cff-version", "message", "title", "version", "authors", "preferred-citation")
    try:
        text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        fields = dict(re.findall(r"^([a-z][a-z-]*):\s*(.*?)\s*$", text, re.M))
        for field in citation_required:
            if field not in fields:
                fail("citation_basic_field_missing", field)
        for field in ("cff-version", "message", "title", "version"):
            if not fields.get(field, "").strip(" \t\"'"):
                fail("citation_basic_field_empty", field)
        if not re.search(r"^\s+-\s+family-names:\s*\S", text, re.M):
            fail("citation_author_name_missing")
        if not re.search(r"^\s+type:\s*\S", text, re.M):
            fail("citation_preferred_type_missing")
    except (OSError, UnicodeError):
        fail("citation_file_unreadable", "CITATION.cff")
    counts["citation_basic_fields"] = len(citation_required)

    local_links = markdown_files = 0
    for relative in sorted(candidates):
        if PurePosixPath(relative).suffix.lower() != ".md":
            continue
        markdown_files += 1
        try:
            text = (ROOT / relative).read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            fail("markdown_unreadable", relative)
            continue
        for destination in markdown_targets(text):
            if not destination or destination.startswith("#"):
                continue
            try:
                link = urlsplit(destination)
            except ValueError:
                fail("markdown_invalid_link", relative)
                continue
            if link.scheme or link.netloc:
                continue
            path = unquote(link.path)
            if not path:
                continue
            local_links += 1
            if path.startswith("/"):
                fail("markdown_absolute_local_link", relative)
                continue
            target = ((ROOT / relative).parent / path).resolve()
            try:
                target.relative_to(ROOT)
            except ValueError:
                fail("markdown_link_leaves_repository", relative)
                continue
            if not target.exists():
                # Report the documentation file, never a URL or embedded value.
                fail("markdown_local_target_missing", relative)
    counts["markdown_files"] = markdown_files
    counts["local_links_checked"] = local_links
    result = {"status": "FAIL" if failures else "PASS", "counts": counts,
              "citation_check": "basic fields only; no YAML/CFF schema validation",
              "files_modified": 0, "new_fit_or_cv_runs": 0, "failures": failures}
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, TypeError, KeyError):
        # Do not echo file contents or exception values into review logs.
        print(json.dumps({"status": "FAIL", "failures": [{"check": "unexpected_file_or_metadata_error"}],
                          "files_modified": 0, "new_fit_or_cv_runs": 0}))
        sys.exit(1)
