# THM-1033 — The alignment dichotomy: aligned cores die by stability, non-aligned by misalignment (death-star-2026-07-18-S56)

**Status:** a clean **decomposition** of the inverse-theorem crux into two mechanisms, with one of them
**rigorous** on compact cores. Verified: 0 counterexamples across 2866 valid non-AP cores. **Does not
close LRC(14)** — the non-aligned mechanism and the large-`max` aligned case remain — but it is the
sharpest structural handle on "only the AP is 182-aligned," and it yields a genuinely new unconditional
bound. Source HYP-7305/7362. Scripts: `04-computation/lrc_align13_deathstar_S56.py`,
`lrc_align_split_deathstar_S56.py`.

Setting: `V = W ∪ {v_max}` primitive covering, `M(V) < 1/13`, `W` = valid non-AP core (covers 2..12,
misses 13,14). `D_W` = a denominator at which `M(W)` is attained; by the wall lattice `D_W ≤ 2·max(W)`.
The far element must satisfy `182 = lcm(13,14) ∣ v_max` (boxeph THM-1017) and the stability window
`v_max ≤ max(W)/(13δ)`, `δ = M(W)−1/13` (THM-1028).

---

## The dichotomy (verified, 0 counterexamples / 2866 cores)

Call `W` **aligned** if `13 ∣ D_W` for some optimal denominator, else **non-aligned**. Then a valid
non-AP core is eliminated in one of two disjoint ways:

| type | count | mechanism |
|---|---|---|
| **aligned** (`13∣D_W`) | 96 | **empty candidate window** — stability alone kills it |
| **non-aligned** (`13∤D_W`) | 16 (nonempty window) | **`M(V) = M(W)`** — the far element cannot cover `W`'s good set |

Not one core had a candidate that lowered `M` below `M(W)`. The AP `{1..12}` alone escapes both — it is
aligned *and* tight (`M = 1/13` at `D_W = 13`), so its (measure-zero) good set is covered by `182 = 14·13`
(the deep well).

## Theorem (aligned cores, RIGOROUS for `max(W) ≤ 34`)

**Every aligned non-AP core with `max(W) ≤ 34` has an empty candidate window** (hence `M(V) = M(W) ≥
1/13`).

*Proof.* Aligned: `13 ∣ D_W`, write `D_W = 13m`. `M(W) = p/D_W` (reduced, so `13 ∤ p`); non-AP ⟹
`M(W) > 1/13 = m/(13m)` ⟹ `p > m` ⟹ `p ≥ m+1`, so `δ = (p−m)/(13m) ≥ 1/(13m)`. The window
`v_max ∈ [182, max(W)/(13δ)]` is nonempty only if `182 ≤ max(W)/(13δ) = max(W)·m/(p−m) ≤ max(W)·m`, i.e.
`m ≥ 182/max(W)`. But `D_W = 13m ≤ 2·max(W)` gives `m ≤ 2·max(W)/13`. Composing,
`182/max(W) ≤ 2·max(W)/13`, i.e. `max(W)² ≥ 182·13/2 = 1183`, i.e. `max(W) ≥ 35`. So `max(W) ≤ 34`
forces the window empty. ∎

This is unconditional — it uses only the wall-lattice bound `D_W ≤ 2·max(W)`, stability, and covering.
It **closes the aligned case for every compact near-AP core** (the deep-well cores have `max = 12`).

## The non-aligned mechanism (verified, not yet rigorous)

For `13 ∤ D_W`, the far element's arcs (denominator `182k`, at multiples of `1/(182k)`) do **not** hit
`W`'s good-set components (at denominators tied to `D_W`), so — although the arcs are wide enough
inside the stability window — they are **mispositioned** and cannot cover, giving `M(V) = M(W)`. Made
rigorous this needs the equidistribution/density step (the covering-width lens gives width, not
position). It is the "position half" of the alignment; the aligned theorem above is the "denominator
half," now proved.

## What this buys, and the residual

- **A new unconditional bound:** aligned non-AP cores with `max ≤ 34` are eliminated with no computation
  (Theorem above) — the first rigorous elimination of a positive-measure family of non-AP cores.
- **The crux is now two clean sub-problems:** (i) aligned large-`max` cores (`max ≥ 35`: finite `m`
  window per `max`, or recurse — these have a far element within `W`); (ii) the non-aligned
  misalignment lemma (`13∤D_W ⟹ M(V)=M(W)`), an equidistribution statement.
- **Only the AP is 182-aligned-and-tight** — every other core fails the denominator test (Theorem) or
  the position test (misalignment). The residual is to make the position test rigorous and to handle
  `max ≥ 35` (renormalization, HYP-3901), which is where the last of the wall lives.

→ THM-1028, THM-1029, THM-1017, THM-1002, THM-724/726, the alignment/OCF/Freiman reflections, HYP-3901.
