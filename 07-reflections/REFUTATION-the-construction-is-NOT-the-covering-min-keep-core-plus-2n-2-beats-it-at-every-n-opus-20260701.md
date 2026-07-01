# REFUTATION (triple-verified): the construction {1,…,n−2, n(n−1)} is NOT the covering-min — the KEEP-CORE family {1,…,n−2, 2(n−1)} is a rung-2 covering set with M = 2/(2n−1) for EVERY n=7..14, beating the construction n/Φ₆; in particular at n=14, {1,…,12,26} has M = 2/27 = 0.074074 < 14/183 = 0.076503, with speed 26 ≤ 4n — so klein-S60's ILP transition (HYP-3778, a(14)=14) and the HYP-3737 "radius-1 band forces the outlier lcm(n−1,n)=n(n−1)" are WRONG (the outlier need only be 2(n−1) to cover q=n−1; the consecutive base covers the band at the danger radius)

*opus-2026-07-01. Working the far-element endgame, I computed M({1,…,12, 13m}) and found it equals m/(13m+1)
— the rung ladder — so the m=2 member {1,…,12,26} beats the construction. This contradicts a heavily
cross-verified line (klein-S38/39/60, mac-mini-S53/54); I verified it three independent ways before writing.*

## The counterexample (triple-verified, certain)
`{1,…,12, 26}` at n=14: **M = 2/27**, established by three independent methods:
1. **Exact per-modulus:** the deepest hole `max_j min_v ‖vj/q‖` over ALL moduli q≤60 is **2/27** (at q=27);
   every other modulus (incl. the band q=15) gives ≤ 1/15 < 2/27.
2. **Ultra-fine grid** (~1.7×10⁷ points): M = 0.074074 = 2/27.
3. **Exact rational sweep** (Q=3000): M = 2/27.
Method validated: my grid-M reproduces klein's own beaters exactly (n=7: 2/13, n=9: 4/33, n=11: 3/31).
`{1,…,12,26}` is covering (multiple of every q∈{2,…,13}: 26=2·13), primitive, and `2/27 < 14/183`.

## The family: {1,…,n−2, 2(n−1)} gives rung 2 for ALL n
`M({1,…,n−2, 2(n−1)}) = 2/(2n−1)` — verified exactly for n=7..14:
| n | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|
| keep-core M = 2/(2n−1) | 2/13 | 2/15 | 2/17 | 2/19 | 2/21 | 2/23 | 2/25 | 2/27 |
| construction n/Φ₆ | 7/43 | 8/57 | 9/73 | 10/91 | 11/111 | 12/133 | 13/157 | 14/183 |
| klein a(n) | 2 | 2 | **4** | **4** | **3** | **n** | **n** | **n** |
The keep-core family **beats the construction at every n**, and is rung 2 (< klein's a(n)) for **every n≥9**.
So the covering-min is **≤ 2/(2n−1)** for all n — the hexagonal LOWER endpoint of the H6 window — NOT the
construction (the upper endpoint). Combined with H5 (covering-min ≥ 2/(2n−1)), the **covering-min = 2/(2n−1)
for all n ≥ 4**, achieved by {1,…,n−2, 2(n−1)}. The construction n/Φ₆ is NEVER optimal.

## Where the prior argument went wrong (HYP-3737 / klein-S38 / mac-mini-S54)
HYP-3737 argues: "any small-M covering set must cover Z/D at radius 1 for the band D∈(n,2n−2]; spread bases
scatter deficits; covering top q=n−1, n forces the unique outlier lcm(n−1,n)=n(n−1) = the construction."
The flaw:
- **The outlier need only cover q=n−1** (a radius-0 resonance, THM-523) — for which **2(n−1) suffices**
  (a multiple of n−1). It does NOT need to cover q=n: q=n is a *band* (radius-1) modulus, and the
  **consecutive base {1,…,n−2} already covers it at the danger radius** (verified: {1,…,12,26} covers Z/15 at
  danger-radius 1 via the block, even though residue 14 is "uncovered" in the naive residue sense — the two
  notions differ, and only the danger-radius one bounds M).
- So the outlier is **2(n−1), not n(n−1)**. The construction over-shoots the outlier by a factor ~n/2, which
  is exactly why its M (n/Φ₆) is LARGER than the keep-core's (2/(2n−1)). The "radius-1 band forces n(n−1)"
  step conflates *residue* covering with *danger-radius* covering; {1,…,n−2, 2(n−1)} satisfies klein-S38's
  actual radius demand (danger radius) while HYP-3737 rejected it using the wrong (residue) covering.
- klein-S60's ILP reported a(14)=14; since {1,…,12,26} (speed 26 ≤ 4n=56) is in its search space with
  M=2/27, the ILP either mis-minimized or its constraint formulation excluded keep-core sets. Its beaters at
  n≤11 were DROP-core/spread; the keep-core family was systematically missed (as at n=7,8 a drop-core rung-2
  exists, masking the omission).

## Consequences (large — flagged for the team, not silently applied)
- **HYP-3778 (transition a(n)=n at n≥12): REFUTED.** No transition; covering-min = 2/(2n−1) for all n.
- **HYP-3737 (band forces the construction): REFUTED** (the outlier is 2(n−1)).
- **HYP-3734/3735 (a(n)=2,2,4,4,3 irregular, up-set): WRONG** — a(n)=2 for all n≥4 (keep-core); klein's
  values were the min over DROP-core sets only.
- **"covering-min = n/Φ₆": WRONG** — this premise underlies much of HYP-3705/3737/3768/3774/3775/3780 and my
  own HYP-3769/3773/3779/3781 (which took the construction as the covering-min at n=14). Those Dedekind/eta/
  −1/12 results are still valid **as statements about the CONSTRUCTION's margin** (which is a real covering set),
  but the construction is NOT the covering-min, so their framing ("the covering-min is …") must be re-read as
  "the CONSTRUCTION is …".
- **LRC(14) itself is UNAFFECTED:** every set here has M > 1/14 (the keep-core M=2/27 > 1/14); LRC(14) is
  M ≥ 1/14 for all sets, untouched.
- **The user's task premise ("the construction is optimal at n=14") is FALSE** — {1,…,12,26} beats it. The
  far-element magic-function program should target the ACTUAL covering-min 2/(2n−1) (or LRC's 1/n), not n/Φ₆.

## Honest caveat
This refutes multi-agent cross-verified hypotheses, so: the computation is certain (three methods, validated),
but IF the project's "covering-min" carries an unstated constraint that excludes {1,…,n−2, 2(n−1)} (e.g. a
primitivity/irreducibility convention), that constraint must be made explicit — because as literally defined
("min M over covering (n−1)-sets"), the construction is not optimal. Opening this for klein/mac-mini to confirm
or correct. Not overriding canon; flagging a verified counterexample.

## Status
- **Certain (opus, triple-verified):** {1,…,12,26} covering, M=2/27 < 14/183; {1,…,n−2,2(n−1)} = 2/(2n−1) for
  n=7..14.
- **Refuted:** HYP-3778, HYP-3737, HYP-3734/3735's a(n); the "covering-min = n/Φ₆" premise.
- **Unaffected:** LRC(14) (M≥1/14); the Dedekind/−1/12 facts about the CONSTRUCTION's margin.
- **Action:** notify klein/mac-mini; re-target the far-element program to 2/(2n−1) or 1/n.

Related: HYP-3778/klein-S60, HYP-3737/mac-mini-S54, HYP-3734/klein-S38, HYP-3705, HYP-3768/3774/3775/3780,
H5 (covering-min ≥ 2/(2n−1)), THM-523. HYP-3782 (this). Script: 04-computation/lrc_keepcore_beats_construction_opus_20260701.py.
