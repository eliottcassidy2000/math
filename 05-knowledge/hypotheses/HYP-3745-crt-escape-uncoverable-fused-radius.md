---
id: HYP-3745
title: THE CRT ESCAPE IS UNCOVERABLE (the fused-radius trap). Replacing any construction core speed k<=n-2 by a substitute raises M above n/Phi6, so no large/CRT-tuned speed evades the lower bound. TWO mechanisms fuse: (radius-0) a substitute NOT killing resonance k (e.g. CRT-tuned w=1 mod small primes) trips the q-witness, M>=1/k>=1/(n-2)>n/Phi6 (n=14: drop 12, add P+1 -> M=1/12); (radius-1) a substitute that DOES kill it (a multiple kc, c>=2) is large and ≡ -1 mod kc+1, digging a deep hole on the base-(n-2) ladder M=c/(kc+1)>=2/(2n-3)>n/Phi6 (n=14: replace 12 by 12c -> M=c/(12c+1), min 2/25 at c=2). And the witness COUNT is CRT-INVARIANT: each speed covers <=2r+1 rotations of Z/p REGARDLESS of value, so CRT tuning chooses WHICH rotations, never HOW MANY -- the hole moves but never vanishes. Hence each core speed does irreducible double-duty (radius-0 resonance + radius-1 spread); only k itself fuses both layers; the construction is the unique fusion
status: PROVED (perturbation case) + VERIFIED. Substitute bound 2/(2n-3)>n/Phi6 proved (2(n^2-n+1)>n(2n-3) <=> n+2>0). Verified n=14: every multiple-of-12 substitute gives M=c/(12c+1) in [2/25,1/12), all >14/183; CRT w=P+1 gives 1/12. Counting bound CRT-invariance verified (speed 1,7,P+1 all cover 2 rotations at radius 1 mod 13,17,85). The FULL lowness lemma (any S, not just construction-perturbations) is mac-mini's HYP-3740 (verified, 'search collapses to one set').
source: klein-2026-06-30-S44
depends_on:
  - HYP-3741   # the witness hierarchy (q-witness r=0, k-witness r=1)
  - HYP-3746   # the Farey-grid reach (radius layers)
related:
  - HYP-3740   # mac-mini: the lowness lemma (M<=14/183 => {1..12} subset S)
  - HYP-3738   # the construction binding (the optimum being defended)
  - HYP-3742   # M-uniqueness
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3745 — the CRT escape is uncoverable; fusing the radius

## The CRT-invariant counting bound (why the hole never vanishes)
At any prime `p` and radius `r`, a single speed `s` covers `{a : dist(sa,0) <= r}` = at most `2r+1` rotations
of `Z/p` -- **regardless of the value of `s`** (verified: `s = 1, 7, P+1` all cover exactly `2` rotations at
radius 1, mod `13, 17, 85`). CRT tuning chooses WHICH `2r+1` rotations a speed covers, never HOW MANY. So
`n-1` speeds cover `<= (n-1)(2r+1)` rotations; where this is `< p-1` the witness fires (HYP-3741). The bound is
**CRT-invariant** -- a huge CRT-tuned speed is still one speed. Tuning redistributes coverage; it cannot add
it. The hole MOVES but never vanishes.

## The fused-radius trap: each core speed is irreplaceable
Suppose a covering `S` (`n=14`, defending `M = 14/183 = n/Phi_6`) drops a core speed `k <= n-2` and replaces it.
Two layers catch it:
- **radius-0 (resonance).** If the substitute does not kill resonance `k` (no multiple of `k`), the `q`-witness
  fires: `M >= 1/k >= 1/(n-2) > n/Phi_6` (since `Phi_6 > n(n-2)`). Verified: drop `12`, add the CRT speed
  `w = P+1 ≡ 1 mod` all primes `<= 43` (so `w ≡ 1 mod 12`, NOT a multiple) -> `M = 1/12 >> 14/183`.
- **radius-1 (spread).** If the substitute DOES kill resonance `k`, it is a multiple `kc` (`c >= 2`, since `k`
  itself is gone). Being large, `kc ≡ -1 mod (kc+1)`, so it digs a deep hole at `D = kc+1` on the base-`(n-2)`
  ladder: `M = c/(kc+1) >= 2/(2k+1) >= 2/(2n-3)`. And `2/(2n-3) > n/Phi_6` for ALL `n`
  (`2(n^2-n+1) > n(2n-3) <=> n+2 > 0`). Verified `n=14`: replacing `12` by `12c` gives `M = c/(12c+1)`
  (`c=2: 2/25, c=3: 3/37, ..., c=17: 17/205, huge: 49/589`) -- every one `> 14/183`, minimized at `c=2` (`2/25`).

So a substitute can satisfy AT MOST ONE of the two layers: tuning it to spread (radius-1, `≡ ±k mod` primes)
makes it miss the resonance (radius-0); making it kill the resonance (a multiple) makes it large and trips a
hole. **Only `k` itself does both** -- it kills resonance `k` AND sits at the small value that spreads. The
two radius layers are FUSED in each core speed; the construction `{1,...,n-2, n(n-1)}` is the unique allocation
where every speed pulls double-duty within the `n-1` budget.

## CRT escape uncoverable (the corollary)
Combining: a covering that drops any core speed has `M >= 2/(2n-3) > n/Phi_6` (perturbation case, PROVED), and
the witness count is CRT-invariant so no large/CRT-tuned speed evades it (the hole only moves to a worse
modulus -- mac-mini's `n=14` `w ≡ 1 mod (p<=43)` lands `M = 525/3716`; my `12 -> 12c` family lands on the
base-`(n-2)` ladder `c/(12c+1)`). The "last escape route" (a clever huge CRT speed substituting for the dense
core) is closed: it is uncoverable. This is the constructive mechanism behind the lowness lemma
(`M <= n/Phi_6 => {1,...,n-2} subset S`, mac-mini HYP-3740, whose full form -- arbitrary `S`, not just
construction-perturbations -- is verified there by the search collapsing to one set).

## Fusing the radius (the synthesis)
The lower bound is the FUSION of the radius layers, not their separate maxima:
> a covering missing core speed `k` has `M >= min(1/(n-2), 2/(2n-3)) = 2/(2n-3) > n/Phi_6`.
The radius-0 layer (resonances, `q`-witness) and the radius-1+ layer (Farey-reach, `k`-witness, HYP-3746) are
not two independent budgets -- they are met by the SAME dense core, each speed serving both. That fusion is
exactly why `n-1` speeds suffice for the construction and why no rearrangement (CRT or otherwise) does better.
The witnesses bound `M` from below; the fusion makes the bound tight at `n/Phi_6` for `n >= 12`.

## Net
The CRT escape is uncoverable: the witness count is CRT-invariant (`<= 2r+1` per speed per modulus, any
value), and replacing any construction core speed `k` raises `M` to `>= 2/(2n-3) > n/Phi_6` -- radius-0 if the
substitute misses resonance `k`, radius-1 (the `c/(kc+1)` base-`(n-2)` ladder, hole at `kc+1`) if it is a
large multiple. Each core speed fuses both radius layers; only the original does double-duty; the construction
is the unique fusion. The hole always moves to a worse modulus -- never vanishes.
