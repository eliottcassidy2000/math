# Quadratic cusp phase research packet

2026-10-19. Frozen initial checkpoint; later refinements belong in a separate extension.

## Files

- `parabolic_escape_precision_20261019.md`: standalone proofs, literature scope, base/digit limitations, and one explicitly open next-order hypothesis.
- `verify_parabolic_escape_precision_20261019.py`: directed-interval phase certificate plus explicitly numerical finite controls.
- `parabolic_escape_precision_20261019.out`: captured successful verifier output.
- `MANIFEST.sha256`: hashes of the three initial research files and this README.

## Reproduce

Run `python verify_parabolic_escape_precision_20261019.py` from this directory. Requires Python 3 and mpmath; tested with Python 3.12.14 and mpmath 1.3.0. Runtime is a few seconds here. No external network access, repository writes, or generated data inputs are needed. The script prints to standard output and modifies no files.

The exact theorem proofs do not depend on finite experiments. The interval computation uses `mpmath.iv` at 50 digits and an explicit analytic tail inequality. Other printed data use 70-digit arithmetic and are intentionally not represented as interval-certified escape integers or as proofs of asymptotic rates.

## Scope and review

The leading pi-law and additive O(1) refinement are established prior literature. The exact-defect proof, endpoint-sum phase, adjacent spacing, staircase interval, and Cauchy occupation law are independently derived new-to-this-session results; no worldwide novelty or priority claim is made. The next-order pi/3 coefficient remains OPEN in this frozen note.

The proof and strict-threshold conventions have been mathematically reviewed; the companion verifier passed. No separate independent audit is claimed in this initial checkpoint.

These results do not imply digit normality, primality transfer, a global analytic embedding, any LRC statement, or a solution to a Collatz problem.
