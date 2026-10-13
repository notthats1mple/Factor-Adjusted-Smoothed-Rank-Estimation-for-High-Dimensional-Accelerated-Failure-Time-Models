# Reproducibility guide

## Review the saved material

Open `manuscript/main.pdf` for the article and `manuscript/ESM_1.pdf` for proofs, implementation details, and supplementary evidence. `manuscript/complete.pdf` provides continuous reading. Table and figure values retain their v0.16 source files.

## Check this repository

From the repository root:

```sh
python3 scripts/check_repository.py
```

This check uses Python's standard library. It checks the shipped-file manifest, source identities, citation metadata, documentation links, and exclusion of private submission material. It does not require R or load fitted objects.

## Compile

Install a TeX Live or MacTeX environment providing pdfLaTeX and the packages named in `manuscript/preamble.tex`. Python 3 exports cross-references. Frozen bibliography files are included.

```sh
sh manuscript/build.sh
```

The command creates the main, supplement, and combined reading PDFs. It does not invoke a solver, cross-validation, or a simulation. Full build logs are generated locally in `manuscript/build_logs/` and are excluded from Git.

The two shipped `xrefs_*.aux` files bootstrap cross-document references; keep them when making a clean source copy. Running the build regenerates them from the newly compiled documents. Rebuilding can change PDF metadata and file hashes. Use a separate checkout to preserve the shipped PDF snapshot for checksum comparison.

## Optional formatting tools

`python3 scripts/format_references.py` refreshes the frozen bibliography from the supplied reference source. `scripts/reformat_existing_figures.py` is a separate figure-formatting utility requiring PyMuPDF, Ghostscript, and Poppler; the repository already includes the formatted figures. A manuscript build does not need that utility.

## Research execution

The current repository candidate contains the manuscript and its publication-build material. Complete fitting/CV source, saved training inputs, and large trajectories remain in the local study archive. Including a research-code release requires a separate content inventory and environment instructions; the current build does not claim to reproduce the six-cohort computation.

## Local preparation verification

A separate clean source copy compiled all three reading documents on 2026-10-14: 46 pages for the main text, 31 pages for Online Resource 1, and 76 pages for the continuous version. The continuous version has its own bibliography and pagination. All page texts match the shipped PDFs; the complete comparison is recorded in `provenance/clean-build.json`. The repository retains the original PDFs, including their original hashes.
