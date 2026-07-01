# The minimal-flip iso-covering subcube: a clean gauge-fixing rule exact for n≤5 that breaks at n=6

*kind-pasteur-2026-07-01-S10. The owner's observation — at n=4 all 4 iso classes are reachable by flipping **2** arcs (not the naive 3), provided the 4 fixed arcs are in "a certain configuration" — turns out to be the first case of a sharp gauge-fixing problem. The clean rule is exact for n=3,4,5 and **breaks at n=6** (proven), a coverage obstruction that arrives before the arithmetic one at n=8.*

## The problem (HYP-3803)

Tournaments on `n` vertices = the hypercube `Q_m`, `m = C(n,2)` (one bit per arc-pair = its orientation). `S_n` permutes the coordinates; **iso classes = `S_n`-orbits**, counted by `A000568(n) = 1,1,2,4,12,56,456,6880,…`. An **iso-covering subcube of dimension `k`** = a choice of `k` FREE arcs `D` together with a fixed orientation of the other `m−k` FIXED arcs, whose `2^k` tournaments hit **every** iso class. Let `k_min(n)` = the least such `k`.

- **Counting bound:** `k_min(n) ≥ ⌈log₂ A000568(n)⌉`.
- **The naive Hamiltonian-path tiling** fixes only the `n−1` base-path arcs and frees the `m−(n−1)` tiles — a *weak* gauge (for n=4 it frees 3 arcs where 2 suffice).

## The result

Computed exactly (exhaustive S_n canonical forms n≤5; union-find over all `2^{15}` tournaments + exhaustive subcube scan at n=6):

| `n` | `m` | `#classes` | counting bound | **`k_min`** | tight? | optimal free-arc shape |
|---|---|---|---|---|---|---|
| 3 | 3 | 2 | 1 | **1** | ✓ | 1 arc (partition `(2,1)`) |
| 4 | 6 | 4 | 2 | **2** | ✓ | **perfect matching** (partition `(2,2)`) |
| 5 | 10 | 12 | 4 | **4** | ✓ | **triangle + edge** (partition `(3,2)`) |
| 6 | 15 | 56 | 6 | **7** | ✗ **+1** | irregular (block-within *fails*) |

For `n ≤ 5` the minimizers are **exactly** the within-block arc sets of the *balanced vertex partition*, and nothing else — the exhaustive n=5 search returned only the 10 `(3,2)`-block sets. So the owner's "certain configuration rule" is:

> **Partition the vertices into two balanced blocks. FREE the within-block arcs; FIX the between-block (complete bipartite) arcs, oriented to order the blocks. `k_min =` the within-block arc count.**

At n=4 the two blocks are `{a,b},{c,d}`: the 2 within-arcs `(a,b),(c,d)` are a **perfect matching** (free), and the 4 fixed arcs are the complementary `K_{2,2}` 4-cycle. That is precisely the owner's n=4 picture.

## The gauge reframe

`k_min` is the dimension of the **smallest axis-aligned subcube meeting every `S_n`-orbit**. Dually, the `f = m − k` fixed arcs **gauge-fix the relabeling symmetry**: fixing the between-block bipartite orientation makes the block structure rigid, so the free within-block flips land in distinct classes. Minimizing free arcs = maximizing the fixed between-block **cut** = the **balanced 2-partition** (Turán's max-cut `⌊n²/4⌋`). The natural ceiling is the **symmetry entropy**: `f ≤ log₂|S_n| = log₂(n!)`, because there are only `≈ 2^m/n!` orbits to separate. Fixing arcs spends the `log₂(n!)` bits of `S_n`; the naive Ham-path spends only `n−1` of them.

## Two obstructions, at two different `n`

There are two independent reasons the clean rule can fail, and they arrive at **different** thresholds:

1. **Coverage obstruction — arrives at n=6 (the actual break).** At n=6 the balanced `(3,3)` block-within subcube has `2^6 = 64` tournaments and the counting bound (6) *permits* covering — yet it reaches only **42/56** classes, and no k=6 subcube whatsoever exceeds **47/56** (exhaustive over all `5005×512`). So `k_min(6) = 7` (proven; a k=7 covering exists), **one over** the counting bound, and the optimal free-set is no longer block-within. The rigid bipartite gauge *over-identifies* classes: distinct within-block fillings collapse under the residual block-swap and internal symmetry.

2. **Counting obstruction — arrives at n=8.** The balanced within-count is `C(n,2) − ⌊n²/4⌋`; it equals the counting bound for `n ≤ 7` but at **n=8** drops to `12 < ⌈log₂ 6880⌉ = 13`: the max-cut `⌊n²/4⌋` (quadratic) finally **overtakes the symmetry entropy `log₂(n!)`** (`~ n log n`), so the balanced partition removes *more* arcs than the symmetry can afford and cannot even count-cover.

The crossover `⌊n²/4⌋ ⋛ log₂(n!)` sits exactly between n=7 (`12 < 12.30`) and n=8 (`16 > 15.30`) — the project's canonical critical window. But the *constructive* break (obstruction 1) already kills tightness at **n=6**: geometry defeats counting before counting defeats itself.

## Creative directions (be free)

- **`k_min(7)` and the gap.** Does `k_min(n) − ⌈log₂ A000568⌉` grow? n=6 gives `+1`. If it grows like the "collapse defect" of the rigid gauge, that defect is a new tournament invariant. (n=7 needs canonical forms over `2^{21}`; nauty-style refinement or orbit sampling.)
- **The break's structure.** The n=6 covering k=7 free-sets are irregular — characterize them. Are they "block-within + one repair arc"? The repair arc = the extra bit paying for the collapse.
- **Even-graph dual.** Repeat on the even-graph metagraph `E_n` (`V = 2,3,7,16,54`). Its denser, chordal-then-not structure should give a *different* `k_min` break point — a dual diagnostic.
- **Gauge capacity = `m − k_min`** as a sequence: `2,4,6,8,…?` (n=3..6). Is it `2(n−2)`? (n=3→2, n=4→4, n=5→6, n=6→8 ✓ so far — a clean conjecture, `f_max = 2(n−2)`, i.e. `k_min = C(n,2) − 2(n−2)`.) Check n=7: predicts `k_min = 21 − 10 = 11` (vs counting bound 9, gap `+2`).
- **Gray code on the metagraph.** Can the covering subcube be ordered so consecutive tournaments are adjacent in `G_n`? A Hamiltonian path through the orbit transversal.
- **Covering-code view.** `k_min` is a covering radius of the orbit space under the subcube family — connect to the domination number of `G_n`.

## The `f_max = 2(n−2)` conjecture

The fixed-arc capacities `f_max = m − k_min` are `2, 4, 6, 8` for `n = 3,4,5,6` — exactly `2(n−2)`. If it holds, `k_min(n) = C(n,2) − 2(n−2) = (n²−5n+8)/2`, giving `k_min(7)=11, k_min(8)=16`, with the gap over the counting bound widening (`0,0,0,1,2,3,…`). This says the gauge can fix exactly `2(n−2)` arcs — twice the base-path length `n−1`, minus 2. A target for n=7.

## Connection to the LRC finish (why this sits next to the census)

The same **entropy-vs-geometry threshold** governs the LRC(14) census: the "clean structured minimizer" (dilated-AP / polygon `(Z/N)*`) is sharp at small size and the pattern strains exactly where a rigid symmetry stops separating cases. In both, a `Z₂` symmetry is the gauge — here tournament **complementation** `T ↦ Tᵒᵖ` (= the dihedral inversion `ι` on the vertex-heptagon), there the LRC scale/antipode. This turn's sharpest-lever check confirms the census side: over 16376 eleven-cores the **pentagon `(Z/10)*` binds** (`313/9702 = 0.032261`); the **sporadic two-clash `(Z/19)` minimum** is `389/12012 = 0.032384`, just **0.38% above** the pentagon — and **both clear `1/36`** (`1.166×`). The binding minimizer is the full-orbit polygon, not the sporadic clash; the tight locus is the balanced polygon, exactly as the iso-covering optimum is the balanced partition.

## Honest status

- **Verified:** `k_min(n)` and the block-within characterization for n=3,4,5 (exhaustive); `k_min(6)=7` (union-find + exhaustive k=6 scan, global best 47/56); the crossover arithmetic; the sharpest-lever census min.
- **Conjectural:** `f_max = 2(n−2)` (fits n≤6, untested n≥7); `k_min(7)`; the gap growth; the even-graph dual.
- **Reframe, not theorem:** the gauge/entropy picture is a lens; the coverage obstruction at n=6 is the one hard, proven fact that a pure counting argument misses.

— Related: `07-reflections/the-isomorphism-class-graph.md`, `merged-metagraph-invariants.md`, `even-graphs-as-first-class.md` (the metagraph objects), the LRC census reflections `the-finish-is-a-recursive-tight-singular-series-…`, `moment-relaxation-reduces-multifar-…`, HYP-3793 (the atom/units census), Burnside/A000568 (`everything-is-the-triangle.md`), [[triangle_foundation]]. Scripts: `04-computation/tournament_iso_covering_{subcube,crossover,n6,n6_kmin,n6_exhaustive_k6}_kps.py`, `lrc14_sporadic_twoclash_min_kps.py`. HYP-3803.
