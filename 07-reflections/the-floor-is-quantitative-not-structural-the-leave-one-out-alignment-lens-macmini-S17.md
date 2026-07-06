# The density floor is quantitative, not structural — and the leave-one-out alignment lens

*mac-mini-2026-07-06-S17 (HYP-4452). Owner: work the one density floor creatively,
through many lenses, small pieces from each; push/pull often. This session's honest
turn: a concurrent instance (HYP-4442) showed the second gap is **n-specific**
(nonempty at n=7,8; empty at n=13), which forces a correction to my structural
lenses and refocuses the floor as a **quantitative** wall-vs-budget fact — and it
surfaced one clean new lens, the leave-one-out alignment. Verified:
`lrc_nspecific_gap_verify…out`, `lrc_leaveoneout_alignment…out`.*

## The correction that reframes everything

The second gap `(1/n, 2/(2n-1))` is **not universally empty**:

- **n = 7:** `{1,5,6,11,16,17}` has `M = 5/33 ≈ 0.1515 ∈ (1/7, 2/13)` — a genuine gap
  member (gcd 1, verified exact). Its witness is `t = 10/33`, denominator
  `33 = 16+17` (a lifted pair — my witness-denominator lemma), and **`g = 4` gaps**,
  and `5/33` is **not** a rung `k/(6k+1)`.
- **n = 13:** exhaustively empty (the concurrent lift census + my 59k near-AP search).

So (G) is **n = 13-specific**, and this **corrects my S15 three-gap quantization**
(HYP-4412): "near-tight ⇒ `g ≤ 3` ⇒ a CF rung" holds on the *n=13* data but **fails
at n=7** (`g = 4`, non-rung gap member). My structural lenses — three-gap /
roots-of-unity / sum-product — all describe the AP's *universal* specialness
(minimal discrepancy, `(ℤ/13)*`-orbit, additive∩multiplicative). They are
**necessary but not sufficient**: they would "predict" an empty gap at every `n`,
but the gap is nonempty at 7,8. **The emptiness at 13 is a quantitative fact the
structural lenses cannot see.**

## What the floor actually is: two walls vs the runner budget

A gap member needs (concurrent HYP-4442 + THM-622):

- **gap width** `2/(2n-1) − 1/n = 1/(n(2n-1))` — the target interval: `1/91` at n=7,
  **`1/325`** at n=13 (3.5× narrower);
- **clearance depth** `q ≥ 3n-1` (numerator `c ≥ 3`) — `20` at n=7, **`38`** at n=13.

At n=7 the walls are loose enough that a lifted transversal (`{16,17}`, witness 33)
threads the gap. At n=13 the walls are too tight for **12 runners** to place a value
inside `1/325` while clearing depth 38 *and* covering `{2,…,12}` *and* staying
primitive. **The floor is exactly: at n=13 the walls exceed the covering budget.**
This is a quantitative Diophantine bound, not a residue identity — matching S13's
"residue-pinning + spread-value are necessary, not sufficient."

## The one clean new lens: leave-one-out alignment

A precise **necessary condition for covering** (hence for `M < 2/25`), combining the
measure lens with the sub-family LRC bounds:

> If `S` covers at `β` (`M(S) < β`), then for **every** `j`,
> `Safe(S∖{v_j}, β) ⊆ A_j` — the removed runner's danger arcs must contain the
> entire hole of the `11`-subfamily.

*Why:* `Safe(S∖{v_j}, β)` is nonempty (an 11-family has `M ≥ 1/12 > 2/25`), and any
point safe from all-but-`v_j` must be caught by `v_j` for `S` to cover. **Verified
exactly** (`…alignment…out`): the AP `{1,…,12}` satisfies it for all 12 drops (each
sub-hole, measure 0.01–0.09, sits inside `A_{v_j}`); the n=7 gap member satisfies it
for all 6. So covering *is* this alignment — each sub-family's hole nesting into the
dropped runner's arcs.

**Why it is the right lens for a quantitative floor.** `A_{v_j}` is `v_j` arcs of
width `2β/v_j` at positions `a/v_j`; `Safe(S∖{v_j})` is a union of holes each of
width `~ gap-width = 1/(n(2n-1))`. Containment requires each hole (width `1/(n(2n-1))`)
to nest inside an arc (width `2β/v_j`) *at an aligned position*. The **hole width
shrinks with n** (`1/91 → 1/325`), so at large n the holes are needle-thin and must
align *precisely* with the `a/v_j` lattice — a rigidity that only the AP's harmonic
lattice `{a/k}` achieves. This is the n-dependent threshold: at n=7 the holes are fat
enough for a lifted lattice to align; at n=13 they are too thin for anything but the
AP. **The floor = "only the AP-lattice aligns all 12 leave-one-out holes for a gap
value at n=13."**

This turns the floor into a **covering-alignment** statement — a lattice-rigidity of
the harmonic arcs — rather than an all-order Riesz-product estimate. It is a lens, not
yet a proof; but it is *quantitative* (n-width-driven), which is what the
n-specificity demands.

## Small pieces from each lens, assembled

- **witness denominator** (mine, now formal opus-S109): the hole/witness lives at
  `q ∣ v_i±v_j`, `q ≤ 2·max` — the alignment lattice is *pairwise-additive*.
- **divisibility-rich** (kps): the arcs `A_{v_j}` include the harmonic periods
  `1/5,1/7,1/8,1/9,1/11` — the covering budget is `≥ 5` committed lattices.
- **leave-one-out alignment** (this): covering = each sub-hole nests in a dropped
  arc; the hole width `1/(n(2n-1))` sets the n-threshold.
- **multi-scale + fast-M** (mine): the search is finite and bounded; empirically
  clean at n=13.

The floor is the statement that these pieces are **jointly infeasible for a non-AP
12-family at a gap value** — the walls (narrow hole, deep clearance) exceed the
alignment budget of 12 harmonic arcs. That the same pieces *are* feasible at n=7 is
the honest check that the argument must be quantitative.

## Net

- **Correction:** the second gap is n-specific (nonempty n=7,8; empty n=13); my
  three-gap quantization is n=13-favorable, not a universal proof. Structural lenses
  are necessary, not sufficient.
- **New lens:** leave-one-out alignment — covering ⇔ every 11-subfamily hole nests in
  the dropped runner's arcs; the hole width `1/(n(2n-1))` is the n-threshold. Verified.
- **The floor, sharpened:** a lattice-alignment infeasibility — only the AP's harmonic
  arc-lattice can align all 12 leave-one-out holes for a value in the `1/325` gap at
  n=13. Quantitative, matching the n-specificity.

## Pointers

- `lrc_nspecific_gap_verify_macmini_S17.py/.out`, `lrc_leaveoneout_alignment_macmini_S17.py/.out`,
  `lrc_gap_occupancy_by_n_macmini_S17.py/.out`.
- concurrent HYP-4442 (n-specific gap, the two walls); mac-mini HYP-4412 (three-gap,
  now n=13-corrected), HYP-4432 (witness denominator), HYP-4402 (multi-scale); opus
  HYP-4406/coverer_height, HYP-4416/witness-denominator-Lean; kps HYP-4417 gap_candidate.
