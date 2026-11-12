# Reproduction

From any working directory, run Python3 on `reproduce.py` using its actual path. For example from the repository root:

```sh
python3 05-knowledge/results/session287_padic_cores_20261112/reproduce.py
```

The standard-library-only harness runs four engines normally and with `-O`, parses their JSON and compares both to the stored receipt. All test predicates use explicit errors rather than removable Python assertions. Each lane's proof describes its finite universe, inherited inputs and all-order argument.

The reset lane includes reviewed Session286 `core.py` and `transfer.py` copies for self-contained execution, not as new algorithms. The independent audit engine uses separate literal/skeleton implementations rather than importing lane results.

`REPLAY.json` records the final full harness replay. `audit/HASHES.json` binds reviewed sources. `RELEASE.md` records publication and aggregate documentation-check scope.
