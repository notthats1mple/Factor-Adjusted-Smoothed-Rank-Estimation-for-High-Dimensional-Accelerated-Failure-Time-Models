# Contributing

Use an issue to describe a reproducible defect, an interface improvement, or a documentation correction. Include the commit, environment output, command, a small input, expected behavior, and observed behavior. Fixed arrays are especially useful for numerical findings.

For a code change:

1. Explain the computational contract or user behavior being changed.
2. Add or update a focused test that distinguishes the corrected behavior.
3. Run `python3 scripts/check_repository.py` and `Rscript --vanilla scripts/run_checks.R`.
4. Update interface documentation and provenance when a preserved implementation changes.

Numerical changes should specify the objective, coefficient coordinates, penalty partition, tolerance, and stopping rule they affect. Preserve explicit failure outcomes and return fields needed to audit a computation. Changes to a numerical contract belong in a clearly identified implementation/version with corresponding tests and configuration.

The repository keeps archived algorithm identities in `provenance/algorithm-source.json`. A change to a tracked original must be described as a modified implementation, with its source and change recorded accurately.

Share inputs that you are authorized to publish. Keep manuscript drafts, private datasets, fitted research objects, credentials, and machine-specific environments out of pull requests.

The authors currently reserve all rights to the code. Proposing or accepting a contribution does not by itself establish a public reuse license; discuss licensing with the authors before distributing derivative software.
