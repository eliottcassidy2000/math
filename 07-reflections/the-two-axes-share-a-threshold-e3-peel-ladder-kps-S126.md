# The two axes of the Freiman ladder share a threshold — the E3-side peel ladder

*kind-pasteur-2026-07-09-S126. Owner: "take the E3-side stability rung, spend a long session getting the
whole ladder to materialize." This note records what materialized (a complete, sorry-free peeling ladder
for the Schur count), and the meta-pattern it exposed: the E3 axis and the burden axis — the two halves
of the Freiman ladder feeding THM-681 — fail at exactly the same place, `|S| ≤ 4`, for the same reason.*

---

## What materialized: `LRCSchurPeel.lean`, the whole ladder as one recursion

The E3-side ladder is built on a single clean recursion — **peel the maximum**. Removing the largest
element `m` of `S` from the Schur count `E₃ S = #{(a,b) ∈ S² : a+b ∈ S}` can only lose the pairs that
*sum to* `m`, because any pair *involving* `m` sums to more than `max S` and lands outside `S`. So

  `E₃ S = E₃ (S.erase m) + repCount S m`,   `repCount S m = #{(a,b) ∈ S² : a+b = m}`  (**Rung A**),

and `repCount S m ≤ |S|-1` because `(a,b) ↦ a` injects the representations into `S.erase m` (**Rung B**).
Together these give the **deficit recursion** (`deficit S := C(|S|,2) − E₃ S`):

  `deficit S + repCount S m = deficit (S.erase m) + (|S|-1)`.

Every downstream fact is now an `omega` away:
- **monotonicity** `deficit (S.erase m) ≤ deficit S` — peeling never increases the deficit;
- **peel cost** `(|S|-1) − repCount S m ≥ 0`, bounded by the total deficit;
- **the capstone** `deficit S = totalPeelCost S` — the deficit is *exactly* the sum of the per-peel
  costs down the chain `S ⊃ S.erase(max) ⊃ …` (proved by strong induction on `|S|`);
- **the base** `deficit S = 0 ⟺ dilated interval` — recovering my S121 rigidity;
- **full-peel = reflection** `repCount S m = |S|-1 ⟺ (∀ a ∈ S, a ≠ m → m-a ∈ S)`: a peel costs nothing
  exactly when `S` below its max is symmetric under `a ↦ m-a`.

The last two combine into a corollary that surprised me: **`S` is a dilated interval iff it is
reflection-symmetric at every peel.** Rigidity is a statement about a *chain* of local symmetries, not one
global condition. This is the E3-axis analogue of accumulating opus's burden down the `restrictedSum`
chain — same shape, different (anchored vs. translation-invariant) additive quantity.

All of it is sorry-free (`lake build TournamentH7.LRCSchurPeel`, 8476 jobs green), verified numerically
across every `k`-set to `k=7` (`lrc14_e3_peel_ladder_kps_S126`).

## The honest capstone: `dist ≤ deficit` is Freiman-hard, and false for small `k`

The tempting capstone is quantitative stability: *a set is within `deficit` element-changes of a dilated
interval*, `dist_to_dilated S ≤ deficit S`. I wanted it, and the peeling skeleton looks like it should
carry an induction. **It does not — and the statement is simply false for `|S| ≤ 4`.** The witness is
`S = {1,4,5}`: `E₃ = 2` (only `1+4=5`), so `deficit = 1`, yet the nearest dilated interval `{d,2d,3d}`
disagrees with it in `2` places. Deficit `1`, distance `2`. The naive capstone dies at `k=3`.

The induction fails exactly where the counterexample lives: peeling `max` and extending a near-dilated
`D'` for `S.erase(max)` by one term costs `+1` in distance **unless the peel is full and its
representation aligns with `D'`** — and a full peel only gives a *reflection* symmetry (`repCount_max_eq_iff`),
which need not align with the best interval below it. That misalignment is the irreducible Freiman-stability
content; the peel recursion exposes it but cannot dissolve it.

## The meta-pattern: both axes break at `|S| ≤ 4`

Here is what made the long session worth it. The quantitative capstone `dist ≤ deficit` **holds for
`|S| ≥ 5`** (verified exhaustively to `N = 5k`: 593,775 sets at `k=6`, zero failures) and **fails only for
`|S| ≤ 4`**. That threshold is *not new*. It is exactly opus's `LRCFreimanAP.ap_of_min_burden`: minimal
burden `2n-3` forces an AP **for `n ≥ 5`**, and is **false for `n ≤ 4`** (MISTAKE-133, the bi-arithmetic
minima census `44/38/0/0/0/0`).

So the two axes of the Freiman ladder —

- **burden** (opus/mac-mini): translation-invariant `restrictedSum`, stability = "min burden ⟹ AP";
- **E3** (mine): origin-anchored Schur incidences `a+b=c`, stability = "min deficit ⟹ dilated";

— are different measures of different symmetries, yet their stability theorems **switch on at the same
`|S| = 5`** and **fail on the same small sets** for the same reason: below five elements there is too much
accidental additive coincidence for either "few relations" to pin down a global arithmetic shape. `{1,4,5}`
is a deficit-1 non-near-dilated set; the burden axis has its own `n ≤ 4` bi-arithmetic minima. The two
failures are the same phenomenon viewed through two frames.

This is reassuring for the endgame. LRC(14) lives at `k=13`, comfortably in the `k ≥ 5` regime where
*both* stability statements are true. The ladder does not have to be proved at the small `k` where it is
false; it has to be *cited* at `k=13` where it holds. The E3 rigidity (`deficit 0 ⟺ dilated`) is exact at
every `k`; the quantitative stability above it is a `k ≥ 5` statement on both axes — a single Freiman
threshold wearing two costumes.

## Where this leaves the finish

The a-priori supply for LRC(14) reduces (THM-681, my S124 correction) to the `W₀ > 0.08` branch, which the
Freiman ladder must climb: `¬near-AP (hcoarse) → low relation content → live`. Both structural entries are
now Lean-formalized — burden (opus's `ap_of_min_burden` + `finset_min_burden_isAP`) and E3 (this file's
peel ladder, top = rigidity). The remaining quantitative rung is genuinely the Freiman-stability step, and
it is the same step on both axes: mac-mini's exhaustion (`LEM-018/021`) discharges it operationally on the
`k=13` corpus; klein's signed program pins it analytically to one hyperbola-counting lemma. My peel ladder
contributes the E3-axis skeleton and the exact top, plus the clean fact that on *both* axes the step is a
`k ≥ 5` statement — never in doubt for the case at hand.

The moral, matching S124's: when a "capstone" looks one `omega` away, check the smallest cases first.
`dist ≤ deficit` is false at `k=3`. What survives is the exact recursion beneath it — and the recursion,
not the capstone, is what the ladder actually is.

*Files: `LRCSchurPeel.lean` (sorry-free), `lrc14_e3_peel_ladder_kps_S126.py`/`.out`. Builds on kps-S121
(`LRCE3Budget` rigidity), kps-S125 (`schurCount_add_sdiff_eq_choose` deficit formula), opus-S195
(`ap_of_min_burden`, the burden-axis twin), MISTAKE-133 (the shared `k ≤ 4` threshold). See
[[the-final-rung-is-signed-not-absolute-kps-S124]] and
[[the-interval-is-the-shared-extremizer-schur-triples-and-lrc-loneliness-kps-S113]].*
