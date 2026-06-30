---
id: HYP-3701
title: ATTACKING the owner's reframe (LRC covering-min = n/Phi_6(n) = 14/183) -- the inequality n/Phi_6(n) >= 1/n is TRIVIAL (n^2 >= n^2-n+1), so the content is the EXTREMALITY, and the honest answer is that n/Phi_6(n) is NOT the covering-min in general: REFUTED for n=4,5,6, where a 'drop-2 + tuned large speed' covering set beats the construction and reaches 2/(2n-1) < n/Phi_6(n) (n=4 {1,3,4}=2/7; n=5 {1,3,4,5}=2/9; n=6 {1,3,4,5,18}=2/11; 2n-1 = the signed-LRC modulus C, =27=3^3 at n=14). BUT these beaters DO NOT SCALE: at n=14, {1,3,..,13,14m} bottoms at 9/83~0.108 >> 14/183, so the construction (which KEEPS 2 and uses n(n-1)=182 to cover q=13,14 at once) STANDS as the tightest known. So the covering-min EXTREMAL FAMILY is n-DEPENDENT (drop-2-split for small n; the n(n-1)-construction at n=14) -- there is NO single uniform optimal construction; that instability is the exact nature/pattern of the razor-thin M-edge. Both candidate floors are prime-3 flavored (n/Phi_6 = Eisenstein/hexagonal Q(sqrt-3); 2/(2n-1)=2/C signed-LRC). The margin (M - 1/n) ~ 1/n^2 -> 0 (the M-edge thins with n)
status: COMPUTED + verified exact (construction M=n/Phi_6(n) n=4..14; covering-min 2/(2n-1) for n=4,5,6 exhaustive-search-confirmed; n=14 drop-2+14m best 9/83). Honest correction of the owner's premise: n/Phi_6(n) is NOT the covering-min (refuted n<=6), and the extremal family transitions; n=14 construction unrefuted but the uniform-extremality picture is false.
source: mac-mini-2026-06-30-S42
related:
  - HYP-3706  # klein-S24: Phi_6 = Eisenstein/hexagonal (p6m), Kershner thinnest-covering analogy; LRC covering-min = hexagonal covering (open)
  - HYP-3700  # the disproof edge is isolated (the gap-edge); this is the M-VALUE edge (covering-min vs 1/n)
  - HYP-3548  # the two razor-thin lines (this refines Line 1, the gap-M line: covering-min n-dependent, margin ~1/n^2)
  - THM-523   # the covering reduction (LRC <=> M(S)>=1/n for covering sets)
results:
  - 04-computation/covering_min_extremal_family_macmini_20260630.py
  - 05-knowledge/results/covering_min_extremal_family_macmini_20260630.out
---

# HYP-3701 -- the covering-min extremal family transitions with n

The owner reframed the open problem: the off-cusp covering-min is `n/Phi_6(n) = 14/183` (positive measure, M
tightest), the inequality `n/Phi_6(n) >= 1/n` is trivial, so the real content is whether `n/Phi_6(n)` is
genuinely the covering-min (no covering set beats the construction `{1,..,n-2, n(n-1)}`). I attacked it.

## The trivial inequality and the construction
`M({1,..,n-2, n(n-1)}) = n/Phi_6(n) = n/(n^2-n+1)` (verified exact, n=4..14). `Phi_6(n) = N(n-zeta_6)` is the
Eisenstein/hexagonal norm (`Q(sqrt-3)`, prime-3; klein-S24 HYP-3706). And `n/Phi_6(n) >= 1/n <=> n^2 >=
n^2-n+1 <=> n >= 1` -- trivially true. So everything is in the EXTREMALITY.

## The construction is NOT the covering-min (refuted, n=4,5,6)
Exhaustive search (covering, primitive, bounded) finds the covering-min BELOW `n/Phi_6(n)`:
| n | covering-min | minimizer | `= 2/(2n-1)`? | `n/Phi_6(n)` |
|---|---|---|---|---|
| 4 | `2/7` | `{1,3,4}` | yes | `4/13` |
| 5 | `2/9` | `{1,3,4,5}` | yes | `5/21` |
| 6 | `2/11` | `{1,3,4,5,18}` | yes | `6/31` |
The covering-min is `2/(2n-1)`, achieved by DROPPING the speed `2` (using `4` to cover `q=2`) and adding a
TUNED LARGE speed (e.g. `18=3*6` for n=6). It is strictly tighter than the construction. Note `2n-1` is the
SIGNED-LRC modulus `C` (THM-407/413; `C=27=3^3` at n=14) -- so the small-n covering-min is `2/C`.

## ...but the beaters do NOT scale to n=14 (the family transitions)
The small-n winner (drop-2 + a moderate large speed) does NOT generalize. At n=14, the family
`{1,3,..,13, 14m}` bottoms at `M = 9/83 ~ 0.108`, FAR above the construction `14/183 ~ 0.0765`. The
construction's trick -- KEEP `2` (use the consecutive `{1,..,12}`) and let the single large speed
`n(n-1)=182=13*14` cover BOTH top `q`'s -- is tighter at n=14. So:
> **The covering-min EXTREMAL FAMILY is n-DEPENDENT: a drop-2 split wins for small n (`<=6`), the
> `n(n-1)`-construction wins at n=14. There is no single uniform optimal construction.**
This instability IS the exact nature of the razor-thin `M`-edge the owner asked about: the extremal config
is not stable across n, so "is the construction the covering-min" has no clean uniform answer -- it is
false for small n and (so-far) true at n=14, with the optimal family switching in between.

## The patterns of the M-edge
1. **Two candidate floors, both prime-3**: `n/Phi_6(n)` (Eisenstein/hexagonal, `Q(sqrt-3)`) and `2/(2n-1) =
   2/C` (signed-LRC modulus). The small-n covering-min is `2/C`; at n=14 the tightest known is `n/Phi_6(n)`.
2. **The extremal family transitions** (drop-2 -> keep-2-construction around n=6-7) -- no uniform extremal.
3. **The margin thins**: `n/Phi_6(n)-1/n = (n-1)/(n Phi_6) ~ 1/n^2`; `2/(2n-1)-1/n = 1/(n(2n-1)) ~ 1/n^2`.
   Both `-> 0`: the covering-min approaches `1/n` as `n` grows (relative margin `~1/n`), so the `M`-edge
   genuinely thins -- and the open question is whether the covering-min STAYS `>= 1/n` (it does for every
   construction checked, trivially, but the EXTREMAL config keeps changing).
4. **This is the M-VALUE edge** (covering-min vs `1/n`), distinct from the GAP edge (HYP-3700, isolated by
   the doublet `0.198`). The M-edge is razor-thin (margin `~1/n^2`) AND has an unstable extremal family;
   the gap-edge is isolated. Two different "razor-thin" phenomena (refines HYP-3548's two lines).

## Honest status (correcting the premise)
`n/Phi_6(n)` is NOT the covering-min: REFUTED for n=4,5,6 (covering-min `= 2/(2n-1)`). The refutation does
NOT extend to n=14 (the small-n beaters give `9/83 >> 14/183`), so at n=14 the construction stands as the
tightest known -- but the clean claim "the construction is the covering-min" is false as a uniform
statement, because the extremal family transitions with n. The single sharpest TRUE statement is weaker than
the owner's: *for each n the covering-min is `>= 1/n` (trivially, for the known extremal family), but the
extremal family is n-dependent and the margin is `~1/n^2`* -- so a proof must handle the transition, not a
single construction. The Kershner/hexagonal optimality (klein-S24) is the right asymptotic frame for the
`n(n-1)`-construction regime, but it is not the small-n covering-min.
