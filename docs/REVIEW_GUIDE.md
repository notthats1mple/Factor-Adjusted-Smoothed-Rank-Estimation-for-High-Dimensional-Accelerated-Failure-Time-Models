# Review guide

The central question is how estimated factor coordinates support sparse coefficient estimation in a high-dimensional, right-censored AFT model. The manuscript develops the target representation, statistical bounds, and numerical accuracy interface in that order.

## Reading map

| Question | Main manuscript | Supporting material |
|---|---|---|
| What coefficient is being estimated, and what does factor adjustment change? | Sections 2–3: ordinary and additional-factor models, partial penalty, reference representation, and frozen prediction map | Appendices A–B |
| How do censoring, pairwise fluctuations, and estimated design enter the coefficient bounds? | Section 4: conditional analysis and primitive factor-design conditions | Appendices C–D |
| How is the objective solved and numerical accuracy interpreted? | Section 5: reference proximal gradient, complete KKT residual, acceleration, and qualified active-face refinement | Appendices E and G |
| What do the saved computations and training selections show? | Section 6: active curvature, refinement, validation paths, and selected model summaries | Appendix F and the supplementary trajectory figure |
| What follows from these results? | Section 7: interpretation and the next empirical questions | Reproducibility notes and source documentation |

Read [the main paper](../manuscript/main.pdf) with [Online Resource 1](../manuscript/ESM_1.pdf), or use [the complete version](../manuscript/complete.pdf) for continuous reading.

## Points for technical review

1. **Target and representation.** Check how the ordinary observed-covariate slope is represented in estimated-factor coordinates, and how the additional-factor model introduces its separate approximation term. Distinguish the sparse residual coefficient from the deployed prediction slope.
2. **Conditions and proof dependencies.** Follow the common-error and conditional-censoring assumptions into score centering and local curvature, then the factor-design conditions into the coefficient bounds. The fixed-tuning conditions are stated with the theoretical results.
3. **Numerical interface.** Check the full-coordinate KKT criterion and its statistical tolerance implication. The reference proximal-gradient theorem, executed accelerated recursion, and fixed-face Newton step have separately specified roles.
4. **Evidence and interpretation.** Relate spectral diagnostics to residual directions and work at a common accuracy. In the training-selection results, retain the failed candidate and boundary selections when interpreting the chosen models.

## Reporting a comment

A useful comment identifies:

- the document and stable section, theorem, equation, figure, or table label;
- the statement under review;
- the mathematical issue, evidence discrepancy, or wording improvement;
- the proposed correction, with any supporting derivation or source.

For numerical comments, include the scenario, method, and penalty value. For a proof issue, identify the first step whose justification needs attention and the assumptions it uses. This keeps scientific questions distinct from presentation changes.

## Repository and document status

This is the unpublished v0.16 author review draft. The supplied PDFs preserve the author package, including declarations that await confirmation. Compilation instructions, reference checks, and available verification records are in [Reproducibility](REPRODUCIBILITY.md); the origin of the manuscript, template files, and repository adaptations is recorded in [Provenance](PROVENANCE.md).

