# The Q50 conjecture and the periodicity reduction: what the absolute-height mechanism has to be

**mac-mini-2026-07-05-S55 (HYP-4119).** The owner's mod-575 CRT directive, executed;
data in `05-knowledge/results/lrc_gap_spectroscopy48_macmini_S55.out` and
`lrc_pole_necessity_macmini_S55.out`.

## Half 1: no residue filter can anchor height (the pole-necessity corollary)

THE FLOATING 7-CLUSTER: bottom {1,2,3,4,5}, cluster {20,21,24,25,45,46,66}·S,
with S ≡ 1 (mod L), L = lcm(2..25) = 26,771,144,400.  On this arithmetic
progression of scales, EVERY residue mod every q ≤ 25 is frozen, so the full
formal gap-violator profile (covering, spread, 24-compression, pair-38,
pinning at every q ≤ 25) passes at EVERY frozen scale — verified exactly to
S ≈ 2.7·10^14 — while w_(7)/w_(8) = 4S → ∞.

Consequences:
  (a) **No ladder constant C_7 exists** derivable from the profile: the
      l ≤ 6 pole of the C_l/T_l compression ladder is NECESSARY.  This is
      the ladder-side twin of opus-S81's fee-mean ceiling (2ρl ≥ 1 kills
      uniform fees at l ≥ 7): the ceiling says the ACCOUNTING dies; the
      cluster says the STATEMENT dies.
  (b) **The profile is periodic on CRT-frozen rays.**  Any profile instance
      generates infinitely many at all frozen scales.  Therefore the
      absolute-height mechanism CANNOT come from residue filters, ratios,
      or any combination — it must come from the WITNESS side: showing the
      2/25-margin point itself exists template-uniformly.

## Half 2: the witnesses ARE template-uniform (the Q50 conjecture)

WITNESS SPECTROSCOPY on the full [25,48] census (511,947 full-profile
survivors — every 12-subset of [1,48] passing all filters):
  - **100.000% have a margin-2/25 witness at a modulus q ≤ 45** (first-q
    histogram supported on [26,45]; mass concentrated on 26..33).
  - **100.000% have a witness at a COMPOSITE-CRT modulus** — q ≥ 26 dividing
    p·p' for pinned p,p' ∈ [13,25] (first composite-q supported on [26,50]).
  - Witnesses NEVER live at q ≤ 25: the pinning conditions are precisely
    "margin at every q ≤ 25 is ≤ 1/q" — the profile bans small witnesses,
    and the first free moduli (26..50) always deliver.

**THE Q50 CONJECTURE.**  Every primitive 12-family passing the formal
gap-violator profile has a rational-point witness of margin ≥ 2/25 with
modulus q ≤ 50.

If true, the gap closes ABSOLUTELY: the witness condition at q depends only
on residues mod q; the profile depends only on residues mod q' ≤ 25 (plus
scale-free ratio facts); everything in sight is determined by the residue
data mod L' = lcm{q ≤ 50 with all prime-power factors ≤ 50} — a FIXED
modulus.  The gap statement becomes a FINITE template verification: no
scale, no census-per-height, no MISTAKE-102 tail risk.

## Why q ∈ [26,50] should suffice (the mechanism sketch)

A witness at q needs the 12 residues mod q to avoid a band of ~4q/25 ≈ 4–8
values around 0 under some unit dilation.  For q ≤ 25 the pinning FORCES
band-hitting (that is what pinning is).  For q ∈ [26,50], the constraints on
W mod q come only through the pinned divisors of q (e.g. q = 46: through the
mod-23 near-pair-system and mod-2), and the remaining freedom — the lift of
the mod-23 shape to mod-46 — always contains a clearing dilation in the data.
The proof shape: per q, a finite check over "q-templates" = (pinned-divisor
shape) × (lift pattern); the shapes are few because the pinning at 23/25 is
near-rigid (11 of 12 / 10 of 12 distinct pairs).  The check couples across
q only through the joint template — the engineering is the enumeration, and
the enumeration is bounded.

## The two-piece mechanism (how this composes with the fleet)

  GAP EMPTINESS  =  (Q50 template verification)      [single-scale, finite]
                  + (cluster descent between scales)  [opus-S81's gap descent:
                     spread tops dodged at any count; multi-scale families
                     peel to single-scale clusters]

The bottom-6 alignment data (same run): w_min ≤ 4 always (spread-forced),
bottom-6 ⊄ [1,12] in 99.45% of survivors, bottom-6-max concentrated in
[19,25] — the profile survivors are "two-band" families (a low band 1..4 and
a mid band ~19..25 + top ~25..48), i.e. already visibly two-cluster shaped:
the descent's input structure is what the census actually produces.

## Status

Q50: CONJECTURE (511,947/511,947 evidence; heights ≤ 48; B=52 census also
HARD=0 at 1.97M survivors).  Next: (i) the per-q template check at q = 26..33
(the mass carriers) — enumerate mod-q templates consistent with pinned
divisors, verify clearing-dilation existence; (ii) state the composition
theorem (Q50 + descent ⟹ dichotomy loose branch) in klein's assembly
language; (iii) push the census to B=64 for more evidence cheaply.
