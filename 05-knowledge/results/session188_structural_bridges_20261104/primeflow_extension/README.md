# Session 188 prime-flow extension

* `PROOF.md`: finite-lattice interval universality, explicit nonmodular
  pentagon, sound terminating multi-source propagation, a coupled arithmetic
  contradiction, and a parity obstruction beyond interval propagation.
* `verify.py`: independent standard-library implementation and checks.
* `normal.json`, `optimized.json`: identical PASS receipts.

Run `python verify.py > normal.json` and
`python -O verify.py > optimized.json`, then compare the files.

The checks cover 51 naturally ordered lattice relations through six elements,
1,625 order/meet/join comparisons, 1,383 binary interval patterns, a fully
explicit prime-arithmetic N5 witness, and 400 allocation-propagation controls
including 176 feasible cases and 206 preserved witnesses. Both contradiction
examples have exact finite row domains recorded in the receipts.

The universality theorem is an existence theorem using Bang–Zsigmondy. Only the
explicit pentagon's prime realization is computed; the 51 finite lattices are
symbolic encoding checks. All OPN implications are conditional partial-branch
obstructions. No support-wide or unrestricted nonexistence claim is made.
