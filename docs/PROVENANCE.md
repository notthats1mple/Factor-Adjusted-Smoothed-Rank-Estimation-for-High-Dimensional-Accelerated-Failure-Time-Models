# Source provenance

The estimation modules and `renv.lock` are copied byte for byte from the archived research project. Their source paths and SHA256 values are recorded in [`provenance/algorithm-source.json`](../provenance/algorithm-source.json). This file also identifies the complete local archive by name and SHA256.

The public release adds a loader, fixed-array examples, regression checks, dependency restoration, repository checks, and documentation. Its solver configuration is a subset of the archived v0.12 configuration, retaining the original numerical controls. The extraction record is stored alongside the public configuration.

`SHA256SUMS.json` covers the published files. Run `python3 scripts/check_repository.py` to check those hashes, frozen algorithm identity, local documentation links, and publication contents. Maintainers can regenerate the release manifest after reviewing an intentional change:

```sh
python3 scripts/check_repository.py --write-manifest
python3 scripts/check_repository.py
```

The complete study inputs, fitted objects, trajectories, failure records, manuscript, and historical reports remain in the authors' local archive. Public tests use explicit small arrays and numerical reference fixtures. Their results validate the distributed code and its interfaces.
