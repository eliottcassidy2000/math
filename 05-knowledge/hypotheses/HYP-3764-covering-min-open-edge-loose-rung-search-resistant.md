---
id: HYP-3764
title: THE COVERING-MIN OPEN EDGE -- the construction n/Phi6 is the LOOSE rung (rung n, near 1/(n-1)), the true covering-min is a SMALL-rung SPREAD set beating it by ~1/n^2, and the claimed "n>=12 transition" is UNVERIFIED (every bounded search misses known beaters). REFRAME: M of a covering set = r/(r(n-1)+1) at Farey rung r (binding pair), INCREASING in r from 1/n (r=1, the floor, UNREACHABLE by coverings since M>1/n strictly, THM-523) to 1/(n-1) (r->inf). The construction {1..n-2,n(n-1)} = rung n = the LOOSE end near 1/(n-1); the covering-min = the SMALLEST realizable rung a(n)>=2. VERIFIED beaters (construction NOT the min): n=7 rung2 (2/13), n=8 rung2 (2/15), n=9 rung4 (4/33), n=10 rung4 (4/37), n=11 rung3 (3/31) -- all < construction. METHOD-FAILURE (the honest edge): random sampling misses ALL structured beaters; hillclimb finds n=7,8,10 but MISSES n=9,11; targeted drop<=2 finds n=9,10 but MISSES n=11 -- so NO search I have certifies "no beater at n=12,13,14", and the SAME searches that "support" the transition demonstrably fail at n=9/11. LRC-ORTHOGONAL: covering-min in [2/(2n-1), n/Phi6], both ~(1+O(1/n))/n, so LRC (rung>=2) holds with margin >=1/(n(2n-1)) REGARDLESS of the exact rung. LEADING HYPOTHESIS (bold): the transition is a search MIRAGE; low-rung beaters persist for all n and the construction is never exactly the covering-min. This challenges the evidential basis of HYP-3737/3747 (n>=12 construction) without a counterexample -- flagged, not refuted
status: OPEN / hypotheses. VERIFIED: rung monotonicity (exact); beaters n=7..11 (construction not the covering-min there, exact M); the search-failure pattern (random/hillclimb/drop<=2 each miss >=1 known beater, exact). CONJECTURAL: whether beaters persist n>=12 (H3 leading, bold) vs the real transition (H4, repo's HYP-3737). LRC-orthogonality is PROVED (rung>=2 => covering-min>=2/(2n-1)>1/n). No beater FOUND at n=12,13,14 -- but the searches are demonstrably unreliable, so this is NOT evidence for the transition.
source: klein-2026-06-30-S53
depends_on:
  - THM-523    # covering-set reduction + M>1/n strictly (=> rung>=2)
  - HYP-3732   # the Farey rung: covering-min = Farey neighbor of 1/(n-1) at rung a(n)
related:
  - HYP-3737   # "radius-1 band over-constraint forces construction n>=12" -- this challenges its basis
  - HYP-3747   # the lowness lemma (n>=12) -- rests on the transition
  - HYP-3763   # my large-multiples-forced (built in the construction regime)
  - HYP-3701   # covering-min extremal family transitions with n (MISTAKE-087)
  - HYP-3726   # margin 1/(n(2n-1)) = 1/hexagonal
  - HYP-3551   # "densest core + minimal killer" -- SUPERSEDED (that is the LOOSE rung, not the min)
results:
  - 04-computation/covering_min_open_edge_klein.py
  - 04-computation/covering_min_probe_klein.py
  - 04-computation/covering_min_targeted_klein.py
  - 05-knowledge/results/covering_min_probe_validation_klein.out
---

# HYP-3764 — the covering-min open edge

## The reframe: M is a Farey rung, and the construction is the LOOSE one
For a covering set the gap `M(S) = j/D` at a binding pair (`D` = pair sum/difference). The covering-min
sits on the Farey ladder above `1/n` (HYP-3732): `M = [0; n-1, r] = r/(r(n-1)+1)` at **rung `r`**. This
is **increasing in `r`**:
```
 rung r:   1     2      3      4     ...    n (=construction)   ...  ->inf
 M     :  1/n  2/(2n-1) 3/(3n-2) ...        n/Phi6(n)           ...  1/(n-1)
 n=14  : .0714  .0741   .0750  .0755 ...    .0765 (14/183)      ...  .0769
```
- **rung 1 = `1/n` = the LRC floor**, UNREACHABLE by covering sets (`M>1/n` strictly, THM-523).
- **rung `n` = the construction** `{1..n-2, n(n-1)}` = `n/Phi6`, sitting near the LOOSE end `1/(n-1)`.
- **the covering-min = the SMALLEST realizable rung `a(n) >= 2`** -- a *spread* set, tighter than the
  construction by `~1/n^2`.

So the long-standing picture "construction = densest core + minimal killer = the tightest covering"
(HYP-3551) is BACKWARDS: the dense construction is the **loosest** rung (`M -> 1/(n-1)`); tightness comes
from SPREAD sets at small rungs. (Superseded by MISTAKE-087; sharpened here into the rung monotonicity.)

## Verified: the construction is NOT the covering-min (n=7..11)
Exact covering-set beaters, all `< n/Phi6`:
```
 n=7 : {1,2,5,6,7,8}          M=2/13 (rung 2)  < 7/43
 n=8 : {1,4,5,6,7,11,16}      M=2/15 (rung 2)  < 8/57
 n=9 : {1,3,4,5,7,11,18,32}   M=4/33 (rung 4)  < 9/73
 n=10: {1,2,3,5,6,7,8,9,30}   M=4/37 (rung 4)  < 10/91
 n=11: (repo) 3/31 (rung 3)                    < 11/111
```
`a(n) = 2,2,4,4,3` for `n=7..11` -- a small, irregular rung, far below `n`.

## The honest open edge: every bounded search MISSES known beaters
Testing whether beaters persist to `n=12,13,14` is defeated by search unreliability (this is the crux):
```
 method               finds beater at         MISSES known beater at
 random covering       (none structured)       all of n=7..11
 hillclimb (capped-D)  n=7,8,10                n=9, n=11
 targeted drop<=2      n=9,10                  n=11  (true 3/31 not found)
```
Every method fails at some `n<=11` where a beater is KNOWN to exist. Therefore "0 beaters found at
`n=12,13,14`" (which all three methods report) is **not evidence** that the construction is the
covering-min there -- the same methods report "0 beaters" at `n=11`, wrongly. The `n=9` beater is a
**drop-3 spread** set; the `n=11` beater eludes drop-`<=2`; the deep search space is where the
covering-min lives.

## The hypotheses (free, ranked by current evidence)
- **H3 (LEADING, bold): the "n>=12 transition" is a search mirage.** Low-rung spread beaters exist for
  ALL `n`; the construction is never exactly the covering-min; `a(n)` stays a small irregular rung. There
  is no structural reason for the `2,2,4,4,3` ladder to jump to `a(n)=n` at `n=12`; the apparent jump is
  the search horizon. This challenges the *evidential basis* of HYP-3737 ("band over-constraint forces
  construction `n>=12`") and hence HYP-3747 (lowness lemma) -- **without a counterexample**, so it is a
  flagged doubt, not a refutation. Resolution: a targeted drop-`>=3` / large-spread search at `n=12,13,14`,
  or a proof that the radius-1 band genuinely forbids rung `< n` for large `n`.
- **H4 (ALTERNATIVE, repo's): the transition is real.** HYP-3737's radius-1 band `[n,2n-2]`
  over-constraint genuinely forbids low-rung spread coverings for `n >= n0 (~12)`, so `a(n)=n` there. If
  true, it needs a PROOF (the current support is search, which is unreliable).
- **H5 (PROVED, LRC-orthogonal): the exact rung does not gate LRC.** For any covering set `M >= 2/(2n-1)`
  (rung `>= 2`, since `M>1/n` strictly forces rung `>= 2`), so `covering-min in [2/(2n-1), n/Phi6]` and
  `LRC(n)` holds with margin `>= 1/(n(2n-1))` (HYP-3726) whatever `a(n)` is. The covering-min open edge is
  a *refinement* question, not an LRC gate.
- **H6 (thin-edge universal form): every margin here is `1/(quadratic)` pinned between cyclotomic /
  hexagonal reciprocals.** floor->mediant `= 1/(n(2n-1))` (hexagonal, HYP-3726); construction->floor
  `= (n-1)/(n*Phi6)` (`Phi6` = 6th cyclotomic); construction->covering-min `~ c/n^2`. The edge is
  `Theta(1/n^2)`: it never closes (rung `>=2`, LRC) and never opens past `1/(n-1)-1/n = 1/(n(n-1))`.
- **H7 (creative): `a(n)` is a covering-system / Jacobsthal realizability function.** Rung `r` is
  realizable iff the resonances `{2..n}` admit a covering set whose binding pair sums to `r(n-1)+1`; this
  is a CRT/covering-system condition on the factorization of `r(n-1)+1`, kin to the Jacobsthal function
  (HYP-2893) and the irregular Farey rung (HYP-3732). Prediction: `a(n)` has no closed form -- another
  genuine number-theoretic sequence, like `width(G_n)` and `A000568`.

## Net
The covering-min's open edge is a `~1/n^2` razor between the construction (loose rung `n`) and a
search-resistant small-rung spread set. The construction is NOT the covering-min for `n=7..11` (verified),
and the belief that it becomes the covering-min at `n>=12` rests on searches that demonstrably miss real
beaters at `n=9,11` -- so it is UNVERIFIED. My leading (bold) hypothesis is that the transition is a
mirage and beaters persist for all `n`; the safe fact is that none of this gates LRC (rung `>= 2` always).
