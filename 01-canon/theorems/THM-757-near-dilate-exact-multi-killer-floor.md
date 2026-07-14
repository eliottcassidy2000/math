---
id: THM-757
title: The near-dilate is exactly the multi-killer floor — for every L divisible by 13 (so the set covers modulus 13), the near-dilate V_L = {L, 2L, …, 12L, 13L+1} is primitive, covering, and has M(V_L) = 1/13 EXACTLY, witnessed at t = (L+1)/(13L). Three-line proof (scaling gives M ≤ M({1..12}) = 1/13; the witness gives M ≥ 1/13). This proves THM-720/721's verified "constant-1/13 near-dilate" claim; and an exact-witness enumeration extends THM-726's finite check into the band (220, 475]: all 8260 interval-core multi-killer covering 13-sets there have M ≥ 1/13. NOTE: the 'equality iff near-dilate' conjecture is FALSE (mac-mini-S106: the minimizer is a whole family, block + free safe killer). The covering near-dilates have diameter ≥ 425881 (need 32760 ∣ L), so they live in opus's density-floor regime, not the band.
status: PROVED (the near-dilate M=1/13 — elementary: scaling invariance + LRC(≤13) for {1..12} + an exact rational witness). Band check EXACT-VERIFIED (8260 interval-core families with largest outlier in (220,475], exact rational M≥1/13 witnesses via the peak candidates t=m/(v_i±v_j), zero failures). The band-floor conjecture (below) is a CONJECTURE.
source: mac-mini-2026-07-14-S105 (executing the covering-case closure band exactly)
depends_on:
  - THM-726   # multi-killer M≥1/13 (this proves its Step-1 extremal shape + extends its finite check)
  - THM-720   # looseness dichotomy / the near-dilate adversary (verified; the M=1/13 claim proved here)
  - THM-724   # single-killer covering-min
related:
  - THM-721   # compressed near-dilate loose floor (death-star) — this is the exact-M companion
  - THM-751   # aligned far-element monotonicity
  - THM-753   # safe-peel reduction to LRC(≤13)
external: LRC(≤13) SETTLED (used for M({1..12})=1/13).
---

# THM-757 — The near-dilate is exactly the multi-killer floor

**One line.** The commensurate near-dilate `V_L = {L,2L,…,12L,13L+1}` — invisible to random sampling
(MISTAKE-101) — is the exact multi-killer extremal: `M(V_L) = 1/13` for every covering `L`, by a
three-line argument. And the covering band `(220, 475]` (between THM-726's finite check and opus's
floor) contains **no** near-dilate; its interval-core families are exact-verified `≥ 1/13`.

## The near-dilate theorem

Let `13 ∣ L` (needed to cover modulus 13) and set `V_L = {L, 2L, …, 12L,\ 13L+1}`.

> **`V_L` is primitive, covering, and `M(V_L) = 1/13`, witnessed at `t = (L+1)/(13L)`.**

*Proof.*
- **Primitive:** `gcd(L, 13L+1) = gcd(L,1) = 1`.
- **`M ≤ 1/13`.** Drop the outlier: `V_L ⊃ \{L,2L,…,12L\}`, and adding a runner only lowers `M`, so
  `M(V_L) ≤ M(\{L,…,12L\})`. Scaling `s = Lt` is a circle bijection with `‖iL·t‖ = ‖i·s‖`, so
  `M(\{L,…,12L\}) = M(\{1,…,12\}) = 1/13` (LRC(≤13), attained at `s=1/13`). Hence `M(V_L) ≤ 1/13`.
- **`M ≥ 1/13`.** At `t = (L+1)/(13L)`: for `i = 1,…,12`,
  `‖iL·t‖ = ‖i(L+1)/13‖ = ‖(i(L+1) \bmod 13)/13‖ ≥ 1/13`, since `13 ∣ L ⟹ L+1 ≡ 1 (13)` so
  `i(L+1) ≡ i \not≡ 0 (13)`; and for the outlier,
  `‖(13L+1)t‖ = ‖(13L+1)(L+1)/(13L)‖ = ‖(L+1)/(13L)‖ = (L+1)/(13L) > 1/13`. So every clearance is
  `≥ 1/13`, giving `M(V_L) ≥ 1/13`.
Together `M(V_L) = 1/13`. ∎

This is the exact-`M` companion to THM-720/721's "constant-`1/13` near-dilate" (there verified /
loose-floored); the equality `M = 1/13` is now a theorem, and the binding subset is the *dilated
consecutive block* `\{L,…,12L\} \cong \{1,…,12\}`.

## The band finite check (exact, interval cores)

`M(S)` is a rational, and its exact value is the max over the **peak candidates** `t = m/(v_i ± v_j)`
(the local maxima of `min_l ‖v_l t‖` occur where two runners cross with opposite slopes,
`(v_i+v_j)t ∈ ℤ`). Searching these for a witness with every clearance `≥ 1/13` is an exact,
self-certifying test.

> **Verified (mac-mini-S105):** every multi-killer covering 13-set with interval core `\{1,…,k\}`
> (`k=9,10,11`) and largest outlier in `(220, 475]` has `M ≥ 1/13` — **8260** families, exact rational
> witnesses, **zero** failures.

This extends THM-726's finite check (outliers `≤ 220`) across the band up to opus's floor threshold
`W_0 ≈ 339–475`. Combined with THM-751 (aligned monotonicity), THM-753 (safe-peel → LRC(≤13)), and the
opus density floor (`W > W_0`), the interval-core multi-killer stratum is closed for **all** outlier
sizes.

**The band has no near-dilate.** A covering `V_L` needs `2³·3²·5·7·13 ∣ L`, i.e. `L ≥ 32760`,
diameter `13L+1 ≥ 425{,}881 ≫ W_0`. So the `M=1/13` extremal lives far out in opus's floor regime, and
the band `(220, 475]` is entirely interval-core-dominated with `M ≥ 1/13` (indeed the interval-core
band families have `M > 1/13` strictly; the `=1/13` equality is only the far-out near-dilate).

## The floor, and the equality case (conjecture CORRECTED — the near-dilate is NOT unique)

> **The floor (conjecture):** every multi-killer primitive covering 13-set has `M ≥ 1/13`.
> Support: THM-726 (certified), THM-751 (aligned monotonicity), THM-753 (bulk → LRC(≤13)), opus's floor.

> **The equality case — the "iff near-dilate" version is FALSE (mac-mini-S106).** `M = 1/13` is achieved
> by a *whole family*, not just the antipodal near-dilate. Any **dilated covering 12-block** `c·\{1,…,12\}`
> (which needs `13 ∣ c` and `2∣c` or `7∣c` — smallest `c=26`) plus **any coprime killer `w` that is safe
> at the block's tight point** `t = 1/(13c)` (i.e. `‖w/(13c)‖ ≥ 1/13`) gives `M = 1/13`. For `c=26` there
> are **173** such killers (`w = 15,17,19,…,339,…`); only **one** is the antipode `13c+1 = 339`.
> Counterexample: `\{15, 26, 52, 78, 104, 130, 156, 182, 208, 234, 260, 286, 312\}` — `M = 1/13` exact,
> primitive, covering, **not** a near-dilate.

> **Corrected equality conjecture:** `M = 1/13` for a multi-killer covering 13-set **iff** it contains a
> tight 12-block `c·\{1,…,12\}` (a dilated consecutive block, `M = 1/13`) and its 13th speed is a coprime
> killer safe at `t = 1/(13c)`. The near-dilate (`w = 13c+1`) is one representative. **Open:** whether the
> tight 12-block is *always* a dilate of `\{1,…,12\}` (LRC(13) tightness rigidity), or whether other tight
> 12-sets occur — that is the remaining rigidity content.

*Why the S105 version was wrong:* `M(S) ≤ M(\text{block}) = 1/13` and adding any block-tight-point-safe
speed keeps equality — the killer is free (subject only to coprimality + covering), so the minimizer is
a family parameterized by `(c, w)`, not a single set. The proved theorem `M(V_L)=1/13` stands; only the
uniqueness claim was false.

## Honest scope

- **Proved:** `M(V_L) = 1/13` (elementary, exact); the interval-core band `≥ 1/13` (8260 exact
  witnesses).
- **DISPROVED (mac-mini-S106):** the S105 "equality iff near-dilate" — the M=1/13 minimizers are a
  whole family (dilated block + free coprime safe killer; 173 killers for `c=26`). Corrected above.
- **Conjecture (open):** the floor `M ≥ 1/13`; and the corrected equality case (M=1/13 iff a tight
  `c·\{1..12\}` block + a safe coprime killer), pending the LRC(13) tightness rigidity of the block.
- **Not covered here:** the *dilated-core* and *incoherent* multi-killer families of intermediate
  diameter (THM-720/721 stratum) — the near-dilate is their extremal, but a full enumeration across
  dilations is the remaining finite content, complementary to opus's floor.

*Artifacts:* `04-computation/lrc14_band_exact_macmini_S105.py` (+out) — exact `M` via peak candidates,
the near-dilate verification, the 8260-family band check. Credits: THM-720/721 (the near-dilate
adversary), THM-726 (the finite check extended), THM-751/753 (the bulk reduction), LRC(≤13) (settled).
