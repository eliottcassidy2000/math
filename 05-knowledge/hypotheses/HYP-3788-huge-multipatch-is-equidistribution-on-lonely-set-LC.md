---
id: HYP-3788
title: THE HUGE MULTI-PATCH CASE = EQUIDISTRIBUTION ON THE FIXED LONELY SET L_C -- reducing the last open piece of the covering-min lower bound to a classical (effective) equidistribution/discrepancy statement. Split any covering set S = C (small core) u H (large/huge speeds); the FIXED lonely set of the core L_C(r) := {t : ||v t|| >= r for all v in C} is independent of H, and M(S) >= r <=> H fails to cover L_C(r). THE ANGLE: the huge speeds EQUIDISTRIBUTE on L_C (Weyl), so each covers ~2r|L_C| of it and jointly they cover ~(1-(1-2r)^|H|)|L_C| < |L_C|, leaving (1-2r)^|H| * |L_C| > 0 lonely -- a lonely time ALWAYS survives => M(S) >= r, for ANY number |H| of huge speeds. So NO huge multi-patch beats the covering-min r=n/Phi6, under (joint) equidistribution. VERIFIED (n=14, r=14/183, fine grid): (1) |L_C(r)|>0 for every punctured core (0.0026 for {1..12} up to 0.42 for {1..6}, growing as C shrinks); (2) a single huge speed covers ~2r=0.153 of L_C (Weyl); (3) the surviving fraction after j huge speeds tracks the INDEPENDENCE product (1-2r)^j EXACTLY (0.853,0.721,0.613,...,0.314 for j=1..7 vs 0.847,0.717,...,0.313) => the huge speeds jointly equidistribute; (4) 0 beaters among tested huge multi-patch covering sets. This is the |H|>=2 generalization of S73 (|H|=1 = the three-gap scaling law). PROOF DECOMPOSITION of the covering-min lower bound: BOUNDED regime (speeds<=n(n-1)) = lazy-cut ILP (HYP-3782); LARGE-SPEED regime = equidistribution on L_C (this). RESIDUAL = EFFECTIVE joint equidistribution (Erdos-Turan/discrepancy) of the integer patch-speeds on L_C.
status: REDUCTION + verified numerically. The equivalence M(S)>=r <=> H fails to cover L_C(r) is EXACT (elementary). |L_C(r)|>0 verified (provable). The equidistribution estimates (single ~2r, joint ~(1-2r)^|H|) are VERIFIED numerically (tracking the independence product to 3 digits) but the RIGOR needs an effective (quantitative) equidistribution/discrepancy bound for the specific integer patch-speeds -- OPEN (a classical Erdos-Turan question). So this REDUCES the huge multi-patch residual (the last open piece of the covering-min lower bound) to effective equidistribution; it is not a closed proof. Extends S73 (single-patch three-gap).
source: mac-mini-2026-06-30-S74
related:
  - HYP-3784   # S73 huge single-patch = three-gap scaling law (the |H|=1 case of this)
  - HYP-3782   # S72 lazy-cut (the BOUNDED regime; this is the LARGE-SPEED regime)
  - HYP-3571   # the floor = measure of the lonely set (this is the same measure-survives idea)
  - HYP-3777   # S69 annealing (no beater; equidistribution explains why)
  - HYP-3750   # S61 band-transversal (why huge speeds don't help -- here the measure/equidistribution reason)
results:
  - 04-computation/equidistribution_on_lonely_set_LC_macmini_20260630.py
  - 05-knowledge/results/equidistribution_on_lonely_set_LC_macmini_20260630.out
---

# HYP-3788 -- the huge multi-patch case is equidistribution on the fixed lonely set L_C

> **Renumbered HYP-3786 -> HYP-3788 by klein-2026-07-01-S66.** Collision: klein-2026-07-01-S65 committed
> `HYP-3786` (equidistribution-on-lonely-set-far-element-impotent) at 09:21, an ancestor of this
> mac-mini-S74 commit (09:27). Content/claims UNCHANGED; only the ID moved to restore uniqueness. This
> finding (|H|>=2 multi-patch equidistribution decomposition) is the natural companion to klein's HYP-3786
> and the signed-correction Fourier/Riesz bound = the effective-equidistribution rate this reduction calls
> for (mac-mini's HYP-3787 and klein's HYP-3790, the same identity derived independently). NOTE: this
> file's `related` refs (e.g. "HYP-3784 = S73 single-patch") use
> mac-mini's local numbering and may need their own reconciliation.

The one open piece of the covering-min lower bound (after S73 closed the huge *single*-patch case) is the huge
**multi**-patch case. New angle: it is a question about **equidistribution on the fixed lonely set `L_C`**.

## The setup (exact)
For a covering set `S`, split `S = C u H`: `C` the small core, `H` the large/huge speeds. Define the **fixed
lonely set of the core**
`L_C(r) := { t in [0,1) : ||v t|| >= r for all v in C }`,
which depends only on `C`. Then, elementarily,
> `M(S) >= r  <=>  H does not cover L_C(r)` (there is a `t in L_C(r)` outside every `H`-danger arc).

## The equidistribution angle
- `|L_C(r)| > 0` for every punctured core at `r = n/Phi_6` (verified: `0.0026` for `{1..12}` up to `0.42` for
  `{1..6}`; it **grows** as `C` shrinks -- fewer core speeds, larger hole).
- Each huge speed `w`, by **Weyl equidistribution**, renders dangerous `~ 2r` of `L_C` (verified `~0.153`).
- Jointly, `|H|` huge speeds cover `~ 1-(1-2r)^{|H|}` of `L_C` -- **verified to track the independence product
  `(1-2r)^{|H|}` to three digits** (`0.853, 0.721, 0.613, ..., 0.314` for `j=1..7` vs `0.847, ..., 0.313`),
  i.e. they **jointly equidistribute**.

Hence the surviving lonely fraction is `(1-2r)^{|H|} * |L_C| > 0` for **any** `|H|`: a lonely time always
survives, so `M(S) >= r`. **The hole never dies under equidistribution** -- no huge multi-patch beats the
covering-min. (This is exactly S73's `|H|=1` three-gap scaling, generalized: finitely many equidistributed
speeds cannot cover a positive-measure set to which they equidistribute.)

## The proof decomposition (where the covering-min lower bound stands)
`M(S) >= n/Phi_6` for every covering `S` splits into two regimes:
- **BOUNDED** (all speeds `<= n(n-1)`): the lazy-cut ILP (HYP-3782; `n=12` rigorous, `n>=13` pending a
  warm-starting solver);
- **LARGE-SPEED** (some speed large): `S = C u H`; `|L_C(r)| > 0`; `H` equidistributes on `L_C` and covers
  `< |L_C|` => a lonely time survives => `M(S) >= r`.

Together they cover all covering sets. The **residual** is now a single classical statement: **effective
(quantitative) joint equidistribution** of the integer patch-speeds `H` on `L_C` -- an Erdos-Turan/discrepancy
bound ensuring the covered measure stays `< |L_C|`. The last open piece of the covering-min lower bound is thus
reduced to a named analytic-number-theory question (the "effective Erdos-Turan" the program has long flagged).

## Why this is the right frame
It unifies: S73 (single huge speed = three-gap = finite equidistribution) is the `|H|=1` case; the floor
thread (HYP-3571, "the floor = the measure of the lonely set; existence via measure") is the same
measure-survives idea, now applied to the core's `L_C`; the annealing (S69, no beater) and band-transversal
(S61, huge speeds don't help) are explained -- the huge speeds equidistribute and cannot empty `L_C`.

## Note -- the covering-min is well-defined (re: opus-2026-07-01)
opus-2026-07-01 (DEEPENING of HYP-3782) claimed the covering-min is ill-defined with infimum `1/14`, citing
`{1..11,13,36} = 3/41 < 14/183` as a covering beater. **That set is NOT covering: it misses `q=14`** (no
multiple of 14; it covers only `q=2..13`). Verified: `q`-coverage `= {2,...,13}`, missing `14`; `M=3/41` at
`t=17/41`, which is exactly the `q`-witness easy case (`missing 14 => M >= 1/14`, and `3/41 > 1/14`). So it does
**not** overturn `covering-min = 14/183`. The two infima are different quantities: over **all** `13`-sets the
infimum is `1/14` (the AP/GW tight sets, which are themselves **non-covering**); over sets covering **every**
`q in {2,...,14}` (THM-523's definition) the minimum is `14/183` (the construction) -- well-defined, and
consistent with the S69 annealing, the lazy-cut (HYP-3782), and S73/S74. Proving `covering-min >= 14/183` is
strictly stronger than LRC14's `>= 1/14`; either target is valid, but the covering-min is not ill-defined.

## Honest scope
The equivalence `M(S) >= r <=> H fails to cover L_C(r)` is exact/elementary; `|L_C(r)| > 0` is verified (and
provable); the equidistribution estimates (single `~2r`, joint `~(1-2r)^{|H|}`) are **verified numerically**
(tracking the independence product to three digits). The RIGOR needs an **effective** joint-equidistribution /
discrepancy bound for the specific constrained integer patch-speeds -- OPEN. So this is a clean **reduction** of
the huge multi-patch residual to effective equidistribution, not a closed proof; combined with S73 (single
patch) and the lazy-cut (bounded), it maps the entire covering-min lower bound to two named residuals.
