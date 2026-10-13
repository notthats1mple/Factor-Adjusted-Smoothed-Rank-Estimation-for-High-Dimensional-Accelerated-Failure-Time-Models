# Factor-Adjusted Smoothed Rank Estimation for High-Dimensional Accelerated Failure Time Models

This repository provides the v0.16 manuscript, supplementary proofs, figures, and LaTeX build materials for peer review.

The paper develops sparse coefficient estimation for right-censored accelerated failure time models with factor-structured covariates. It connects an estimated-factor representation to the original sparse slope, establishes finite-sample bounds under stated model and design conditions, and links numerical accuracy to the coefficient analysis. Computational diagnostics and training-only cross-validation use six fixed synthetic cohorts. The research algorithm source, training inputs, and large fitted objects remain in the study archive; this repository contains the paper and its compilation materials.

**Authors:** Weineng ZHU and Yang JU. **Status:** unpublished author review draft, version 0.16. The manuscript uses the *Lifetime Data Analysis* submission format. Author declarations and the corresponding-author designation await confirmation; the supplied PDFs retain that draft status.

## Read the paper

| Document | Purpose |
|---|---|
| [Main manuscript](manuscript/main.pdf) | Model, estimator, theory, computation, and numerical studies |
| [Online Resource 1](manuscript/ESM_1.pdf) | Proofs, numerical evidence notes, and implementation details |
| [Complete reading version](manuscript/complete.pdf) | Main text and appendices in one document |

These PDFs are copied unchanged from the v0.16 author package. Their source identity and the repository adaptations are documented in [Provenance](docs/PROVENANCE.md).

The [review guide](docs/REVIEW_GUIDE.md) maps the paper's arguments to the relevant sections and suggests how to report comments. [Reproducibility](docs/REPRODUCIBILITY.md) describes compilation, reference checks, and the available verification records.

## Build

From the repository root:

```sh
sh manuscript/build.sh
```

The build uses pdfLaTeX from TeX Live or MacTeX and Python 3. Python exports cross-references using its standard library. The included bibliography files support the normal build. Consult [Provenance](docs/PROVENANCE.md) for the Springer class files and the separate local typography override.

The command compiles `main.tex`, `supplement.tex`, and `complete.tex`, then copies `supplement.pdf` to `ESM_1.pdf`. It produces document files only. Additional tools used for PDF checks or figure formatting are listed in [Reproducibility](docs/REPRODUCIBILITY.md).

## Cite and discuss

Use the unpublished-manuscript entry in [CITATION.cff](CITATION.cff) when citing this version. For a review comment, identify the document, section or equation label, and the proposed correction or question. Numerical comments should identify the scenario, method, and tuning value so they can be matched to the reported evidence.


## Reuse terms

The authors have not yet assigned a reuse license to their manuscript and project scripts. The third-party Springer template retains its upstream conditions, recorded in [Third-party material](docs/THIRD_PARTY.md).
