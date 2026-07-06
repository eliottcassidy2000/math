# The Cluster-GCD Ladder: an absolute-height mechanism for common-factor clusters

**kind-pasteur-2026-07-05-S17 (HYP-4217, part 2).  Status: PROVED (elementary; full proof below), verified 24/24 on planted adversarial families.  Lean: not yet (this draft is the handoff).**

## Statement

**THEOREM (cluster-gcd ladder).** Let `W` be a 12-set of nonzero integer speeds
with no `2/25`-margin point (in particular any spectral-gap violator, since the
citation gives `M ≥ 1/13`).  Let `S ⊂ W` with `1 ≤ |S| ≤ 6`, let `T := W ∖ S`,
and let `d := gcd(T)` (the gcd of the absolute values).  Then

    (25 − 4·|S|) · d  ≤  50 · Σ_{w ∈ S} |w|.

Rung table (`d ≤ 50·ΣS/(25−4|S|)`):

| |S| | 1 | 2 | 3 | 4 | 5 | 6 |
|-----|---|---|---|---|---|---|
| d ≤ | 50·ΣS/21 ≈ 2.4·ΣS | 50ΣS/17 ≈ 2.9ΣS | 50ΣS/13 ≈ 3.8ΣS | 50ΣS/9 ≈ 5.6ΣS | 10·ΣS | 50·ΣS |

The pole at `|S| = 25/4` is the same six-element wall as the fee ladder
(HYP-4115/4116) — this is its dual, seen from the periodicity side.

## Why this matters

mac-mini-S55 proved the profile filters are CRT-ray-periodic: the floating
cluster `{1..5} ∪ {20,21,24,25,45,46,66}·S` passes every residue/ratio filter at
every frozen scale, so **no filter bounds absolute height**.  The ladder is the
mechanism the filters cannot supply: for that shape (`T` = the 7-cluster,
`S = {1..5}`, `ΣS = 15`, `|S| = 5`) it gives `gcd ≤ 10·15 = 150` — a gap
violator cannot carry the cluster to unbounded scale.  Every common-factor
escape route (the S11 counterexample family `{c,…,11c,12c+1}` included:
`S = {12c+1}`… there `T = c·{1..11}`, `|S| = 1`: `d = c ≤ 50(12c+1)/21` —
consistent, since that family is loose; the point is the ladder BINDS whenever
the complement is genuinely small) is closed at explicit height.  What the
ladder does NOT close: clusters with small gcd (`d = 1` tops) — the confinement
regime, which is the descent/census lane.

## Proof

Write `C := S`, `c := |C| ≤ 6`, `T` the other `12 − c ∈ [6, 11]` runners, all
divisible by `d`.

**Step 1 (citation).**  `|T| ≤ 11`, so LRC(≤13) gives `t₀` with
`‖w·t₀‖ ≥ 1/(|T|+1) ≥ 1/12 > 2/25` for every `w ∈ T`.

**Step 2 (periodicity).**  Every `w ∈ T` is `d·w′`, so `w·(t₀ + j/d) =
w·t₀ + w′·j` differs from `w·t₀` by an integer: the `T`-margin is invariant
under `t ↦ t + j/d`.  Hence ALL `d` points `t_j := t₀ + j/d`, `j = 0..d−1`,
have `T`-margin `> 2/25`.

**Step 3 (the cover obligation).**  Since `W` has no `2/25`-point, at each `t_j`
some runner is within `2/25` — by Step 2 it must be a `C`-runner: each `t_j`
lies in a tooth (radius-`2/25` arc) of some `w ∈ C`.

**Step 4 (the sharp count).**  Fix `w ∈ C`; count `J_w := #{j : ‖w t_j‖ < 2/25}`.
Let `g := gcd(w, d)`, `w = g w′`, `d = g d′` with `gcd(w′, d′) = 1`.  Then
`w·t_j mod 1 = w t₀ + (w′ j mod d′)/d′ mod 1`, and as `j` runs over `0..d−1`,
`w′ j mod d′` takes each value of `ℤ/d′` exactly `g` times (unit multiplication
is a bijection on `ℤ/d′`).  So

    J_w = g · #{k ∈ [0, d′) : ‖w t₀ + k/d′‖ < 2/25}.

The set `{x ∈ [0,1) : ‖w t₀ + x‖ < 2/25}` is one circular arc of length `4/25`
(at most two intervals after unwrapping), so the `1/d′`-grid meets it in at most
`(4/25)d′ + 2` points.  Hence

    J_w ≤ g·((4/25)d′ + 2) = (4/25)·d + 2g ≤ (4/25)·d + 2|w|.

**Step 5 (pigeonhole).**  Every `j` is counted:

    d ≤ Σ_{w∈C} J_w ≤ (4c/25)·d + 2·Σ_{w∈C} |w|,

and `c ≤ 6` gives `25 − 4c ≥ 1 > 0`:

    (25 − 4c)·d ≤ 50·Σ_{w∈C}|w|.  ∎

## Verification

`04-computation/lrc_cluster_gcd_kps_S17.py` (this session): 24/24 planted
families with `d` exceeding the rung bound are loose, with the `2/25`-witness
found at one of the periodic copies `t₀ + j/d` — the mechanism, not just the
inequality.  (Consistency: the dilated AP `c·{1..12}` with `S` = any single
element `cj`: `d = c ≤ 50cj/21` ✓ — the tight locus survives the ladder, as it
must.)

## Lean handoff (any instance; all tools exist)

1. Periodicity: elementary (integer-shift, cf. `margin_add_int` in
   `LRCGridAttainment`).
2. Step 4 in Lean: (a) `ZMod d′` unit-multiplication bijection for the
   `g`-to-1 count (`ZMod.unitOfCoprime`, `Function.Bijective.sum_comp`-style);
   (b) the arc/grid count: at most two integer intervals, each `≤ len + 1`
   (cf. the interval-counting in `LRCGapLadder`).
3. Step 5: `Finset.card_eq_sum_card_fiberwise` over the assignment
   `j ↦ (its C-runner)`.
4. Consumable shape: mirror `gap_ladder_rung` (integer-form conclusion
   `(25 − 4·S.card)·d ≤ 50·Σ`).

Suggested name: `LRCClusterGcd.lean`, `gap_gcd_rung`.
