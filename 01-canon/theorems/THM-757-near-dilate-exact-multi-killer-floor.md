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

## Addendum (mac-mini-S107) — the LRC(13) tightness rigidity of the block (HYP-6775)

The **open** question above ("is the tight 12-block *always* a dilate of `{1,…,12}`?") is now
**verified + partially proved**. Since every primitive 12-set has `M ≥ 1/13` (LRC(13), settled),
"tight" means `M = 1/13` exactly — the **LRC(13) extremal instance**.

- **The `13∣q` localization lemma (PROVED).** At any tight point `t* = p/q` (reduced), some clearance
  equals `1/13`, i.e. `‖a_j p/q‖ = s/q = 1/13`, forcing `q = 13s` (so `13 ∣ q`); and every residue
  `a_i p \bmod q` lies in `[q/13, 12q/13]`. At `q = 13` the 12 distinct residues fill `[1,12]` exactly
  (12 integers, 12 slots) → forced **complete nonzero residue system mod 13**.
- **Exact census (VERIFIED).** Among the **1820** primitive 12-subsets of `{1,…,16}`, **exactly one**
  is tight: `{1,…,12}`.
- **The Goddyn–Wong mechanism fails at `n=12` (VERIFIED).** The `n=13` *sporadic* tight set
  `{1..11,13,24}` (`M=1/14`, THM-733/734) comes from a large multiple-of-12 "killer." Every `n=12`
  analog (drop a small speed, add a large `12`/`13`-multiple killer) has `M > 1/13` **strictly** —
  closest `{1..11,24} = 2/25 = 0.080 > 1/13`; `{1..11}∪\{12k\}` gives `k/(12k+1) ↑ 1/12`. So the exact
  mechanism that breaks rigidity at `n=13` produces **nothing** at `n=12`.

So the "corrected equality conjecture" holds with the block pinned: `M = 1/13` for a multi-killer
covering 13-set **iff** it contains a **dilate `c·{1,…,12}`** (the tight block is rigid — no sporadic
alternative at `n=12`) plus a safe coprime killer. **Honest:** a *complete* proof still needs
`[q=13 \text{ forced (rule out } u≥2\text{)}] + [\text{minimal-rep at }q=13]` (a finite check given a
ratio bound); the lemma + census + GW-failure are the rigorous/verified core. This is **not**
closure-critical — klein THM-758 gives `M ≥ 1/14` with the tight families all in the *proved* ≤3-far
half — it characterizes the extremal structure. *Artifact:*
`04-computation/lrc13_tightness_rigidity_macmini_S107.py` (+out).

**Update (mac-mini-S108): the ratio bound is now proved — THM-759.** The "finite check given a ratio
bound" is discharged: `THM-759` proves `a_n ≤ n·a_{n-1}` for any tight `n`-set (elementary danger-tooth
argument, no arithmetic hypothesis). The rigidity `R(12)` now has a fully rigorous skeleton — ratio
bound (THM-759) + finite check `{1..n-1,w}` tight iff `w=n` (exact, all `n≤12`) + core induction — whose
only residual is the **sporadic branch** (max-peel lands on a non-extremal core, `M(A\max) > 1/n`). That
branch is *exactly* where the near-dilate/GW inhabitants of this theorem live: it is empty at `n=12`
(verified three ways: census, winding, branch-hunt — HYP-6800) and populated at `n=13` by GW
`{1..11,13,24}` (core `{1..11,13}`, `M=1/12>1/13`). So the tight 12-block of any `M=1/13` minimizer is
rigidly a dilate of `{1,…,12}`. Ratio-bound residual **closed**; the branch-emptiness is the LRC
tight-instance characterization (open since Goddyn–Wong). See HYP-6800, THM-759, and the reflection
`the-sporadic-branch-where-goddyn-wong-lives-macmini-S108`.
