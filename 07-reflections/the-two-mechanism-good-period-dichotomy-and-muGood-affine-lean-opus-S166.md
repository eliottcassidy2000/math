---
source: opus-2026-07-09-S166
status: (math) the good-period capstone splits into TWO OPPOSITE-mechanism branches on the
  longest-AP axis -- near-AP (cluster-the-AP, ELEMENTARY, klein LEM-012, DONE) vs dissociated
  (phases-already-spread, needs the r_N/W-hat a-priori, kps-S91) -- explaining why they tile all L
  with no gap and why the dissociated last-mile is analytic not elementary. (lean) muGood_affine
  (dilation ∘ translation invariance of the covering good-set) PROVED, kernel-pure -- the
  WLOG-normalize primitive, completing the S155-S164 dilation-invariance formalization.
tags:
  - lrc14
  - good-period
  - capstone
  - longest-AP
  - lean
  - dilation-invariance
---

# The two-mechanism good-period dichotomy, and muGood_affine in Lean

**opus-2026-07-09-S166.** Owner: finish the math, then formalize. The fleet closed the good-period
capstone to one branch overnight (klein-S196 LEM-012 proved near-AP; mac-mini formalized the S165
spine). This note (a) explains the branch structure — why near-AP is elementary and dissociated is
analytic — and (b) formalizes the affine-invariance primitive the whole reduction rests on.

## (math) The capstone's two branches use OPPOSITE mechanisms

The good period `j*` (THM-527-A) is governed by the longest AP `L` (opus-S164, dilation-invariant).
The two branches that tile all `L` use **opposite** gap-creation mechanisms:

- **Near-AP (`L >= k-5`): CLUSTER the AP.** klein-S196 LEM-012 (elementary): Dirichlet-cluster the
  `L`-term sub-AP into a tiny arc (`‖jd/V‖` small, `j <= ceil(7(L-1)/(L-k+6))`); the complement is a
  single gap `> (m+1)V/7` (`m = k-L <= 5`), which the `m` strays cut into `<= m+1` pieces, the largest
  `> V/7`.  The gap is MADE by concentrating structure.
- **Dissociated (`L <= k-6`): the phases are ALREADY spread.** Here `m = k-L >= 6`, so
  `S = 1-(m+1)/7 <= 0`: no sub-AP can be clustered small enough, and `m >= 6` strays CAN fill any gap
  (`m+1 >= 7` pieces, each possibly `<= V/7`).  The elementary gap-split provably fails.  Instead the
  good period comes from EQUIDISTRIBUTION: a dissociated set has few small additive relations, so the
  dilated phases `{e_i j/V}` decorrelate fast in `j`, and `W(j/V) ~ E[W] > 0` at small `j` -- the
  `r_N < 1` route (opus-S165), whose a-priori target is the near-resonance `W-hat` sum (mac-mini-S61),
  bounded by LEM-011.  The gap APPEARS from spreading, not concentrating.

So the dichotomy is intrinsic: **concentrate (near-AP) vs spread (dissociated)**, the two ways to make
a `>1/7` gap, split exactly at the pigeonhole threshold `m = 5` (`L = k-5`).  This is why LEM-012 is
elementary (it only needs the concentrate side) and the dissociated last mile is genuinely analytic
(it is the equidistribution/`W-hat` side) -- there is no elementary escape for `m >= 6` (the `>= 6`
strays defeat any single-cluster pigeonhole).  The remaining work is exactly kps-S91's `r_N` a-priori
bound for the dissociated branch; opus-S165 reduced it to `|Corr_N| < N(6/7)^k` with a 16% margin.

## (lean) muGood_affine — the WLOG-normalize primitive

`TournamentH7/LRCGoodDilation.lean` (builds, kernel-pure `[propext, Classical.choice, Quot.sound]`):

> **`muGood_affine`**: `muGood θ (E.image (fun e => c * e + m)) = muGood θ E`  for `0 < c`.

So `muGood` (equivalently `D3`) depends ONLY on the affine-dilation class of the speed set -- every
family reduces to its PRIMITIVE representative `(E − min E)/gcd`.  It composes this session's
`muGood_dilate` (dilation) with `LRCTailDiameter.muGood_translate` (translation).  This is the
reduction the whole k=11 covering tail (opus-S155–S164) and the dilation-invariant good period `j*`
rest on -- now a formal, kernel-pure primitive.  Full file: `emptyArc_dilate`, `good_dilate`,
`good_add_one`, `good_add_natCast`, `muGood_fold`, `muGood_dilate`, `muGood_affine`.

(The `#print axioms` on `muGood_dilate` and `muGood_affine` both report `[propext, Classical.choice,
Quot.sound]` — no `sorryAx` — retroactively confirming the S163 `muGood_dilate` purity too.)

## Ledger / next

- MATH: the good-period capstone is two opposite-mechanism branches (concentrate/near-AP DONE
  elementarily; spread/dissociated = the analytic `r_N` last mile, kps-S91) split at `m = 5`; no
  elementary escape for the dissociated side.
- LEAN: `muGood_affine` done (WLOG-normalize, kernel-pure), completing the covering good-set
  dilation-invariance formalization.
- NEXT: kps-S91's `r_N < 1` a-priori for dissociated (`|Corr_N| < N(6/7)^k`, opus-S165 target); on the
  Lean side, the good-period assembly (LEM-010 branches → good period) building on mac-mini's nodes.
- Files: `LRCGoodDilation.lean`; builds on klein-S196 LEM-012, mac-mini-S61, opus-S164/S165.
