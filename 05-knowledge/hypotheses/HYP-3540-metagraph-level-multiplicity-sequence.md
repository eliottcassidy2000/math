---
id: HYP-3540
title: The iso-class arc-flip metagraph eigenvalue d-2k has multiplicity = dim of S_n-invariants at level k of Q_{C(n,2)}; this level-multiplicity profile is a Burnside/orbit-counting sequence with a closed form
status: OPEN (klein-2026-06-29-S1) — multiplicities computed exactly n=3..6; closed form / OEIS match not yet identified
source: klein-2026-06-29-S1
depends_on:
  - THM-584
related:
  - HYP-3538
  - HYP-2080   # self-complementary / Burnside
results:
  - 04-computation/r-block-spectra-antipodal.py
  - 05-knowledge/results/r-block-spectra-antipodal.out
---

# HYP-3540 — The metagraph level-multiplicity sequence

## Claim

By THM-584 the iso-class arc-flip metagraph adjacency `A` (the S_n-quotient of the hypercube `Q_d`,
`d=C(n,2)`) has eigenvalues `d - 2k` with multiplicity

```
mult(d-2k) = dim of S_n-invariant functions at level k of Q_d
           = number of S_n-orbits of k-subsets of the d arcs (pairs of vertices)
           = number of iso classes of "k-edge digraph patterns on the arc set", a Burnside count.
```

Conjecture: this level-multiplicity profile has a clean closed form (a polynomial-in-`n` /
Burnside-over-S_n expression per `k`), and the **R-even / R-odd** block dimensions
`(A000568±SC)/2` are its even-`k` / odd-`k` partial sums.

## Evidence (exact, n=3..6)

Block spectra (`r-block-spectra-antipodal.py`), with multiplicities:

- **n=4** (`d=6`): R-even `{-2, 2, 6}` (mult 1,1,1); R-odd `{0}` (mult 1). Levels k=0..6 contribute
  invariants 1,?,?,… summing to 4 classes.
- **n=5** (`d=10`): R-even `6,2(×4),-2(×3),-6,10` → eigenvalue/mult: 10:1, 6:1, 2:4, −2:3, −6:1
  (sum 10 = dim V_+); R-odd `0,4` (mult 1,1; sum 2 = dim V_-).
- **n=6** (`d=15`): R-even `15:1, 11:1, 7:5, 3:10, -1:13, -5:4, -9:1` (sum 34); R-odd
  `9:1, 5:5, 1:8, -3:6, -7:2` (sum 22).

Note `mult(d) = 1` always (Perron, level 0), and `mult(d-2)` = number of S_n-orbits of single arcs =
1 (all arcs equivalent), explaining the single second-from-top eigenvalue.

## Why it matters

The multiplicities ARE the per-level orbit counts — the same Burnside machinery that gives A000568
itself (`Fix(σ)` over S_n), now refined by edge-count `k`. A closed form would give:
- the full metagraph spectrum at every `n` without enumeration (the eigenvalues `d-2k` are known; only
  the multiplicities are missing);
- a level-graded refinement of A000568 (`= Σ_k mult(d-2k)`, its total) and of SC
  (`= Σ_k (-1)^k mult(d-2k) = dim V_+ − dim V_-`, the signed/antipodal Euler characteristic — this one
  is automatic from the block dimensions, verified 2,2,8,12 for n=3..6, not open).

## Next steps

1. Compute `mult(d-2k)` = #S_n-orbits of k-subsets of `C(n,2)` arcs via Burnside, n up to ~10, and
   search OEIS (this is "number of graphs/digraph-edge-patterns on n labeled-then-quotiented nodes
   with k arcs" — likely a known 2-variable table).
2. Verify the signed identity `SC(n) = Σ_k (-1)^k mult(d-2k)` (antipodal Lefschetz/Euler relation).
3. Connect to THM-584's mod-4 law and to the merged-metagraph invariants reflection.
