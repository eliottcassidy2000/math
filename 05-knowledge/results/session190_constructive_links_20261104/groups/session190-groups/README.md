# Symmetry, Maschke, and metagraph edges

**Status:** elementary proofs with independently reproducible exact controls; peer review pending. No claim of literature novelty.

The most actionable result is that every labelled path-complement metagraph edge has trivial ordinary stabilizer. Under relabelling plus global converse, its stabilizer is1 orC2, with a one-reflection test. Exact all-order edge-orbit formulas are proved in `PROOF.md`.

Other results distinguish planar O(2) from the five finite SO(3) families, give a literal stabilizer diagram for1,C2,7,21 inside AGL(1,7), and explain why Paley7's group needs six real dimensions but only three binary dimensions. Maschke's characteristic condition is retained throughout. Existing Fermat-prime rigidity and the raw7,21 spectrum exclusions are credited, not reclaimed.

## Reproduce

    python3 verify.py --check
    python3 -O verify.py --check

Python standard library only. About two seconds per run on the session machine. `results.json` is deterministic. Every mathematical check uses an explicit exception rather than an assertion.

Universes:

- All labelled tournaments and actual metagraph edges at n3,4,5, with Hamilton paths checked independently by permutations and subset dynamic programming.
- All normalized off-path edge cubes at n3 through7. The n8 row is formula-only.
- All5040 Paley7 vertex permutations; all189 Hamilton paths and neighbours.
- All128 binary seven-coordinate vectors and all168 invertible binary3x3 matrices.
- All circulant choices at n3,5,7.

Files: `PROOF.md`, `SOURCES.md`, `verify.py`, `results.json`, normal and optimized check receipts, and `SHA256SUMS`.
