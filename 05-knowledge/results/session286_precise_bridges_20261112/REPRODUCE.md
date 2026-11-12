# Reproduction and scope

From any working directory:

```
python3 /path/to/session286_precise_bridges_20261112/reproduce.py
```

The root runner executes eleven standalone standard-library verifiers both normally and with Python optimization enabled, compares parsed output with each stored receipt, and prints a final receipt. It uses explicit checks, not removable assertions. `REPLAY.json` records the final run. Allow several minutes and a few hundred MB for the mixed-alphabet finite base. No network inputs or external packages are required.

The automata verifier deterministically regenerates its attachment table and baseline graph table. All other proof-control receipts are read-only inputs to the runner. The additive padded-bank control directly checks all six advertised witnesses. Wildcard search scripts are exploratory extras, not part of the proof replay. The wildcard search scripts use the inherited Session285 producer; their bounded search records are not all-order claims.

`audit/REVIEW.md` separates mathematical review from finite controls and inherited inputs. Independent verifiers rebuild literal packet data and selected path counts without importing the producer implementation. Exact finite matrix identities support the written induction arguments; replay alone is not a proof of the global Chan–Pak conjecture.
