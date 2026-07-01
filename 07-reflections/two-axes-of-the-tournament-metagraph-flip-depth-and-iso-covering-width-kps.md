# Two axes of the tournament metagraph: minimal-flip DEPTH (= MFAS ≈ H) and iso-covering WIDTH (= gauge entropy) — and where the LRC atoms sit

*kind-pasteur-2026-07-01-S10. Merging my minimal-flip iso-covering study (HYP-3803) with the two pieces of concurrent work it collided into: opus-S14/klein-S70's heptagon tournament (HYP-3802) and klein-S72's packing/covering safari (HYP-3804). The three fit one picture: the metagraph `G_n` has a vertical DEPTH axis (minimal-flip distance to transitive = min feedback arc set, which tracks the OCF invariant H) and a horizontal WIDTH axis (how many free arcs to pack/cover all iso classes). The LRC(14) atoms sit at a specific depth on the maximally-symmetric class.*

## The DEPTH axis: minimal-flip distance to transitive = MFAS ≈ H

The metagraph geodesic from a class to the **transitive** root is its **minimum feedback arc set** (MFAS) — the fewest arcs to reverse to make it acyclic, i.e. the minimal-flip depth in the arc-hypercube. Computed at n=7:

| class | MFAS (depth) | H (Ham paths) | #3-cycles | |Aut| |
|---|---|---|---|---|
| transitive | 0 | 1 | 0 | 5040 |
| **R₇ rotational {1,2,3}** | **6 = φ(14)** | 175 | 14 | 7 |
| **Paley QR {1,2,4}** | **7 = n** | 189 (max) | 14 | 21 |

- **MFAS tracks H strongly.** Over 500 random 7-tournaments, `Pearson r(MFAS, H) = 0.850`, with monotone bucket means (depth 0→6 gives mean H `1 → 9.7 → 30 → 57 → 91 → 130 → 154`). So the minimal-flip depth *is* the project's **H-gradient / principal line**: transitive at the shallow end (H=1), Paley at the deep end (H=189, MFAS=7 = joint maximum; random rarely exceeds 6).
- This gives the "principal line" a **combinatorial coordinate**: depth = MFAS. The H-gradient that the repo has always described qualitatively is, quantitatively, the feedback-arc-set gradient.

## The heptagon lands at depth φ(14) — the LRC atom count

opus-S14/klein-S70 (HYP-3802) found the 7 odd points `{1,3,…,13}/14` = a regular heptagon (`D₇`, order 14), with the rotational tournament R₇ = "transitive + 6 flipped tiles = 6 units." My frame sharpens the "6": it is not just the base-labeling Hamming weight — the base ordering `0→…→6` is **optimal**, so the true metagraph depth `MFAS(R₇) = 6 = φ(14)`, exactly the LRC(14) **atom / unit count**. Paley sits one deeper at `MFAS = 7 = n` (and is klein-S72's doubly-regular 2-(7,3,2) design, the H-maximizer).

**Honest scope:** `6 = φ(14)` is *n=7-specific*, not a general law — `MFAS(R₅)=3 ≠ φ(10)=4`. But n=7 is exactly the LRC(14)-relevant case (7 = 14/2), so the alignment is meaningful there: the LRC atoms live at depth φ(14) on the maximally-symmetric (SC / `D₇`) class. This is the tournament-side image of HYP-3802's "minimal covering M ↔ maximal SC/D₇ symmetry" conjecture: **maximal symmetry sits at a distinguished depth on the principal line**, and that depth counts the atoms.

## The WIDTH axis: packing vs covering (a duality that folds at n=6)

Orthogonal to depth is how much of the class-space a single axis-aligned subcube can reach. Two dual numbers bracket the information content `log₂ A000568(n)`:

- **Packing** — klein-S72's *rainbow number* `R(n)` = max subcube dimension whose completions are **all distinct** classes = `⌊log₂|G_n|⌋ = 1,2,3,5` (n=3..6, exhaustive).
- **Covering** — my *flip-rank* `ρ(n) = k_min` = min subcube dimension hitting **every** class = `1,2,4,7` (n=3..6; HYP-3803).

They straddle the bound: `R = ⌊log₂⌋ ≤ log₂ A000568 ≤ ⌈log₂⌉ ≤ ρ`, and the gap `ρ − R = 0,0,1,2` **grows**. Classes *pack* at the info floor but *cover* only above the ceiling — **covering is strictly harder than packing under `S_n`** (klein's phrase: "the group folds the cube," obstructing covering while aiding packing). My n=6 result is the sharp face of this: even though the counting bound *permits* a 6-cube, no 6-cube exceeds 47/56 (exhaustive), so `ρ(6)=7` — the coverage obstruction is a proven fact a pure counting argument misses. The clean **balanced-block-within** optimizer (free within-block arcs, fix the between-block bipartite gauge) is exact for n≤5 and is the first casualty at n=6.

## The synthesis: a depth × width coordinate system, with the LRC crux on the deep SC spine

Put together, `G_n` carries:
- a **depth** (vertical, principal line): MFAS ≈ H, from transitive (0) to Paley (max). The SC "spine" runs up this axis; the heptagon sits at depth φ(14).
- a **width** (horizontal, gauge/entropy): packing `R` and covering `ρ`, bracketing `log₂ A000568`, with covering folding at n=6.

The LRC(14) finish, read here, is a **depth-extremality on the SC spine**: the maximally-symmetric class (R₇, `D₇`) sits at depth = atom count, and the covering-min conjecture (OPEN-Q-108) is the claim that *minimal covering measure ↔ maximal symmetry ↔ this distinguished depth*. Meanwhile the census sharpest-lever (this session) confirms the analytic side: the **balanced polygon `(Z/10)*` binds** (`313/9702`), the sporadic clash is 0.38% above, both clear `1/36` — the same "balanced structure is extremal" motif as the balanced-partition covering optimum.

Both axes exhibit the repo's signature **entropy-vs-geometry threshold**: depth saturates and width folds exactly where a rigid symmetry (complementation `ι`, block-swap, `D₇`) stops separating cases — n=6 for covering, n=7 for the heptagon alignment, n=8 for the max-cut/`log₂(n!)` crossover.

## Honest status

- **Verified:** MFAS values and `r(MFAS,H)=0.85` (n=7 sample); `ρ(n)=1,2,4,7` and the n=6 coverage obstruction (exhaustive); the packing/covering bracket (convergent with klein-S72); the sharpest-lever census min.
- **Convergent (not solo):** the flip-rank sequence was found independently by klein-S72 (HYP-3804) — two derivations, same numbers, same break; HYP-3803 (mine, first-pushed) carries the block-within characterization + n=6 proof + gauge reframe, HYP-3804 the packing dual + QR-design + skew-spectrum.
- **Conjectural:** `f_max = 2(n−2)` (fits n≤6) ⇒ `ρ(7)=11`; MFAS-max = Paley for all n; the depth-extremality form of OPEN-Q-108.
- **Reframe, not theorem:** the two-axis picture is a lens; the one hard new fact is `ρ(6)=7` (coverage beats counting), and the one clean bridge is `MFAS(R₇)=φ(14)` at the LRC case.

— Related: `the-minimal-flip-iso-covering-subcube-a-gauge-that-breaks-at-n6-kps.md` (the width axis), HYP-3802 (heptagon, opus-S14/klein-S70), HYP-3804 (packing/covering safari, klein-S72), HYP-3803 (this study), `the-isomorphism-class-graph.md`, `geometric-alignment-of-merged-metagraph.md` (spine/ribs/sea), `everything-is-the-triangle.md` (H-gradient), OPEN-Q-108, [[triangle_foundation]]. Scripts: `04-computation/tournament_heptagon_flipdepth_kps.py`, `tournament_flipdepth_vs_H_kps.py` (+ .out).
