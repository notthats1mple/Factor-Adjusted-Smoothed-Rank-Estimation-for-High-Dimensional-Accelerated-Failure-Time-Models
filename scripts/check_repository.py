#!/usr/bin/env python3
"""Validate the public code tree, frozen algorithms, links, and SHA256 manifest.

Only Python's standard library and Git are required. Run --write-manifest after
all release edits; the default mode checks an existing manifest without writes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tempfile
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = "SHA256SUMS.json"
PROVENANCE = "provenance/algorithm-source.json"
FORBIDDEN_SUFFIXES = {
    ".pdf", ".tex", ".eps", ".bib", ".bbl", ".sty", ".cls", ".clo",
    ".aux", ".blg", ".rds", ".rda", ".rdata", ".zip", ".tar", ".gz",
    ".xz", ".7z", ".rar", ".dmg", ".pkg", ".pyc", ".pyo", ".pem", ".key",
}
FORBIDDEN_PARTS = {
    "manuscript", "manuscripts", "paper", "papers", "figures_source",
    "figures", "tables", "appendices", "submission", "archive", "archives",
    "reports", "data", "runtime", "migration", ".venv", ".library",
    ".cache", "__pycache__", "node_modules", ".git",
}
PRIVATE_PATHS = re.compile(
    r"/(?:Users|Volumes|home)/[^\s\"'<>]+|[A-Za-z]:[\\/](?:Users|Documents and Settings)[\\/]",
    re.IGNORECASE,
)
SECRET_PATTERNS = {
    "GitHub token": re.compile(r"\b(?:github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})\b"),
    "API secret": re.compile(r"\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{32,}\b"),
    "AWS access key": re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b"),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |DSA |OPENSSH |ENCRYPTED )?PRIVATE KEY-----"),
    "credential-bearing URL": re.compile(r"https?://[^/\s:@]+:[^/\s@]+@"),
}
HEX_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


class ValidationError(Exception):
    """A release check failed; messages never include credential values."""


def validate_relative_path(name: str) -> None:
    path = PurePosixPath(name)
    if (
        not name or path.is_absolute() or ".." in path.parts
        or "\\" in name or any(ord(char) < 32 for char in name)
        or path.as_posix() != name
    ):
        raise ValidationError("Invalid repository-relative path in release metadata")


def public_files() -> dict[str, Path]:
    command = ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, check=False)
    if result.returncode:
        raise ValidationError("Git could not enumerate the release working tree")
    try:
        names = {item.decode("utf-8") for item in result.stdout.split(b"\0") if item}
    except UnicodeDecodeError as exc:
        raise ValidationError("Release filenames must use UTF-8") from exc
    files = {}
    for name in sorted(names):
        validate_relative_path(name)
        path = ROOT / name
        # Staged deletions are absent from Git's index; unstaged deletions are
        # absent from the working tree. Neither belongs to this public snapshot.
        if not os.path.lexists(path):
            continue
        if path.is_symlink():
            raise ValidationError(f"Symlinks are outside the public release: {name}")
        if not path.is_file():
            raise ValidationError(f"Expected a regular public file: {name}")
        files[name] = path
    if not files:
        raise ValidationError("The public file set is empty")
    return files


def read_json(path: Path) -> object:
    def unique_object(pairs):
        obj = {}
        for key, value in pairs:
            if key in obj:
                raise ValidationError(f"Duplicate JSON key in {path.name}")
            obj[key] = value
        return obj

    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, OSError) as exc:
        raise ValidationError(f"Cannot read valid JSON from {path.name}") from exc


def check_content(files: dict[str, Path]) -> tuple[dict[str, str], dict[str, str]]:
    hashes = {}
    texts = {}
    for name, path in files.items():
        relative = PurePosixPath(name)
        parts = {part.lower() for part in relative.parts}
        basename = relative.name.lower()
        if (
            relative.suffix.lower() in FORBIDDEN_SUFFIXES
            or parts.intersection(FORBIDDEN_PARTS)
            or basename == ".env" or basename.startswith(".env.")
            or basename in {".rhistory", ".ds_store", ".netrc", ".git-credentials"}
        ):
            raise ValidationError(f"File is outside the code-only release: {name}")
        payload = path.read_bytes()
        if b"\0" in payload:
            raise ValidationError(f"Binary content is outside the public text release: {name}")
        try:
            text = payload.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValidationError(f"Public files must be UTF-8 text: {name}") from exc
        if PRIVATE_PATHS.search(text):
            raise ValidationError(f"Private absolute path found in {name}")
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(text):
                raise ValidationError(f"{label} pattern found in {name}; value suppressed")
        texts[name] = text
        if name != MANIFEST:
            hashes[name] = hashlib.sha256(payload).hexdigest()
    return hashes, texts


def check_frozen_sources(files: dict[str, Path], hashes: dict[str, str]) -> int:
    if PROVENANCE not in files:
        raise ValidationError(f"Required provenance file is missing: {PROVENANCE}")
    record = read_json(files[PROVENANCE])
    if not isinstance(record, dict) or not isinstance(record.get("files"), list) or not record["files"]:
        raise ValidationError("Algorithm provenance requires a nonempty files list")
    seen = set()
    for item in record["files"]:
        if not isinstance(item, dict):
            raise ValidationError("Invalid algorithm provenance entry")
        name, expected = item.get("path"), item.get("sha256")
        if not isinstance(name, str) or not isinstance(expected, str) or not HEX_SHA256.fullmatch(expected):
            raise ValidationError("Algorithm provenance requires paths and SHA256 digests")
        validate_relative_path(name)
        if name in seen:
            raise ValidationError(f"Duplicate frozen source: {name}")
        seen.add(name)
        if name not in hashes:
            raise ValidationError(f"Frozen source is missing from the public file set: {name}")
        if hashes[name] != expected:
            raise ValidationError(f"Frozen source differs from its archived SHA256: {name}")
    return len(seen)


def prose_only(markdown: str) -> str:
    lines = []
    fence_char = None
    fence_length = 0
    for line in markdown.splitlines():
        match = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if match:
            marker = match.group(1)
            if fence_char is None:
                fence_char, fence_length = marker[0], len(marker)
            elif marker[0] == fence_char and len(marker) >= fence_length:
                fence_char = None
            continue
        if fence_char is None:
            lines.append(line)
    return re.sub(r"(`+).*?\1", "", "\n".join(lines))


def markdown_targets(markdown: str):
    prose = prose_only(markdown)
    # Both ordinary/image links and reference definitions are checked. A
    # destination may use angle brackets, a title, and percent-encoded spaces.
    inline = re.compile(r"!?\[[^\]\n]*\]\(\s*(?:<([^>\n]+)>|([^\s)]+))(?:\s+[^\n]*?)?\s*\)")
    definitions = re.compile(r"^ {0,3}\[[^\]\n]+\]:\s*(?:<([^>\n]+)>|(\S+))", re.MULTILINE)
    for pattern in (inline, definitions):
        for match in pattern.finditer(prose):
            yield match.group(1) or match.group(2)


def check_markdown_links(files: dict[str, Path], texts: dict[str, str]) -> int:
    checked = 0
    for name, text in texts.items():
        if PurePosixPath(name).suffix.lower() != ".md":
            continue
        for target in markdown_targets(text):
            parsed = urlsplit(target)
            if parsed.scheme or parsed.netloc:
                if parsed.scheme.lower() in {"file", "javascript", "data"}:
                    raise ValidationError(f"Unsupported Markdown link scheme in {name}")
                continue
            decoded = unquote(parsed.path)
            if not decoded:
                continue  # A same-document anchor.
            if "\\" in decoded:
                raise ValidationError(f"Markdown local link must use forward slashes in {name}")
            base = ROOT if decoded.startswith("/") else (ROOT / name).parent
            resolved = (base / decoded.lstrip("/")).resolve()
            try:
                relative = resolved.relative_to(ROOT).as_posix()
            except ValueError as exc:
                raise ValidationError(f"Markdown link leaves the repository in {name}") from exc
            is_public_file = relative in files
            is_public_directory = resolved.is_dir() and (
                relative == "." or any(item.startswith(relative + "/") for item in files)
            )
            if not is_public_file and not is_public_directory:
                raise ValidationError(f"Markdown local target is missing from release: {name} -> {decoded}")
            checked += 1
    return checked


def check_manifest(files: dict[str, Path], hashes: dict[str, str]) -> None:
    if MANIFEST not in files:
        raise ValidationError(f"Missing {MANIFEST}; finalize edits then use --write-manifest")
    expected = read_json(files[MANIFEST])
    if not isinstance(expected, dict) or any(
        not isinstance(name, str) or not isinstance(digest, str) or not HEX_SHA256.fullmatch(digest)
        for name, digest in expected.items()
    ):
        raise ValidationError("SHA256SUMS.json must map each public filename to its SHA256")
    if set(expected) != set(hashes):
        missing = sorted(set(hashes) - set(expected))
        extra = sorted(set(expected) - set(hashes))
        detail = []
        if missing:
            detail.append("unlisted public files: " + ", ".join(missing))
        if extra:
            detail.append("manifest-only files: " + ", ".join(extra))
        raise ValidationError("Manifest file set differs: " + "; ".join(detail))
    changed = [name for name in hashes if expected[name] != hashes[name]]
    if changed:
        raise ValidationError("Manifest SHA256 differs: " + ", ".join(changed))


def write_manifest(hashes: dict[str, str]) -> None:
    payload = json.dumps(hashes, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=ROOT, prefix=".manifest-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(payload)
        os.replace(temporary, ROOT / MANIFEST)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true", help="Validate the release, then write SHA256SUMS.json")
    args = parser.parse_args()
    try:
        files = public_files()
        hashes, texts = check_content(files)
        frozen = check_frozen_sources(files, hashes)
        links = check_markdown_links(files, texts)
        if args.write_manifest:
            write_manifest(hashes)
            print(f"WROTE {MANIFEST}: {len(hashes)} public files; {frozen} frozen sources; {links} local links")
        else:
            check_manifest(files, hashes)
            print(f"PASS: {len(hashes)} public files; {frozen} frozen sources; {links} local links")
        return 0
    except (ValidationError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
