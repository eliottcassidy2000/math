# AK certificate rung ladder — frozen witnesses (klein-S691, 2026-07-28)

All FINITE-EXACT (exact rational forcing closure, re-verified fresh at
record time; engine `04-computation/ak_forcing_engine.py`, superset game
`ak_mode3_v2.py`). Game hierarchy and caveats in
`06-writeups/AK-FORCING-WORKBENCH-klein-S691.md` §4b. NO AK(α) claim is
made anywhere — scores are certificate arithmetic; the certificate ⟹ AK
theorem is the external benchmark's.

## Ladder state (2026-07-28 evening)

| game | record | witness size | status |
|---|---|---|---|
| loose (literal spec text) | → 1 (unsound) | 14/9 at n=9 | spec bug, report upstream |
| strict (= verifiable, fixed) | **13/7 ≈ 1.857** | n=8, two witnesses | robust attractor |
| per-suffix labels, merge-free | **7/4 = 1.750 (= KT99 exponent)** | n=8 | this file |
| + legal-slot merges (mode3-v2) | **12/7 ≈ 1.714** | n=7 eff. | annealing → 5/3 |

## The 7/4 witness (per-suffix game, [2,2,2], m=10, r=4, n=8, t=0)

Coordinates (a,b,c) ∈ {1,2}³; slot ((prefix),(suffix)) → label; a-axis =
stage 1 (per-suffix labels — intuitive-style step gluing), then b, c.

- a-edges (1,b,c)−(2,b,c): (b,c)=(1,1)→(1,2); (1,2)→(2,1); (2,2)→(1,1);
  (2,1)→(1,2).
- b-edges (a,1,c)−(a,2,c): only at c=2: a=1→(2,1); a=2→(0,1).
- c-edges (a,b,1)−(a,b,2): (2,2)→(2,1); (2,1)→(1,2); (1,1)→(1,1);
  (1,2)→(0,1).
- Seeds: (1,0)@(1,1,2); (0,1)@(1,1,1); (1,1)@(1,2,1); any-of-4@(2,1,1)
  (e.g. (1,0)) — the last seed is the "finisher": without it exactly
  (2,1,1) is stuck at 7/8.
- X = {(0,0),(1,0),(0,1),(1,1),(1,2),(2,1)}; all nonzero a+b ≠ 0 ✓.

Derivation path (recorded for mechanism value): the annealer's 12/7
mode3-v2 witness had ONE merge (1,1,1)~(2,1,1); replacing that merge by
label (1,2) leaves exactly one vertex unfired (7/8), and any single seed
finishes ⟹ 7/4. So on this witness the identification is worth exactly
1/4 of score: merge(12/7) → edge+seed(7/4). The merge = free hub; the
compilation tax here is one seed.

## Verifiable-format expressibility: OBSTRUCTED for this witness

Checked all 6 axis orders: the verifiable format demands non-final axes
have suffix-UNIFORM labels and all-suffix edge presence (f_i domain =
prefixes; cost d_{i+1}···d_k each). The witness needs suffix-dependent
labels on ≥2 axes ⟹ not expressible. The strict/verifiable record hence
remains 13/7. Open compilation lemmas (= the benchmark's implicit
equivalence claims):
1. per-suffix step labels → prefix-varying deeper products (plausible via
   variation-by-depth);
2. identifications → verifiable structure (unknown; worth trying to
   REFUTE score-preservation instead: if merges are strictly stronger,
   the two published formats are inequivalent).

## Next rungs

5/3 ≈ 1.6667 < 1.675 (the benchmark target, modulo soundness+format);
integer-feasible shapes at small n: (m+r, n−t) = (10,6), (15,9), (20,12).
12/7 → 5/3 needs one more unit of savings at n−t = 21 scale or a
6-effective-vertex design at m+r = 10.
