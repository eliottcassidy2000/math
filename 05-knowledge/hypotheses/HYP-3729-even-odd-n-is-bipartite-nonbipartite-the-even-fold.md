---
id: HYP-3729
title: The even-n / odd-n LRC split IS the BIPARTITE / NON-BIPARTITE (even/odd cycle) split -- EVEN n: the equally-spaced {1..2p-1} descends (2-adic, Mode B) to the equally-spaced {1..p-1} (n=14 -> {1..6} = 6 runners = the PROVEN Barajas-Serra range), the even cycle C_n is BIPARTITE so the apex gap lambda_min(Q(C_n))=0 (the equally-spaced = the degenerate measure-0 cusp, M=1/n), and the even-fold reduces even-speed coverings to the smaller case; ODD n: C_n is NON-BIPARTITE so the apex gap >0, and the covering-min is the irregular SPREAD frontier (n=7 {1,2,5,6,7,8} M=2/13, n=9 {1,3,4,5,7,11,18,32} M=4/33 -- speeds*a mod D land irregularly, gaps {1,1,2,2,3,4}, NOT a clean structure). This is the same even/odd = bipartite/non-bipartite = sigma-even/sigma-odd = orderable/Condorcet split running through the whole repo
status: VERIFIED (the even-fold descent {1..2p-1}->{1..p-1} for n=6,8,10,14; the beaters M=2/13, 4/33 + their irregular images; lambda_min(Q(C_n))=0 for even n, >0 for odd). The even-fold gives the TIGHTNESS (the equally-spaced achieves 1/n) and reduces even-speed coverings; the full even-n LOWER bound (no covering beats 1/n) is the hard conjecture direction (the descent's gap relation M(S)<=M(E/2)/2 is the upper-bound direction). The odd-n covering-min frontier (the irregular spread) is open (mac-mini owns it).
source: klein-2026-06-29-S35
depends_on:
  - HYP-3604   # apex gap = lambda_min(Q(C_p)) >0 iff p odd (non-bipartite)
  - THM-580    # the 2-adic descent (Mode B, the even-fold)
related:
  - HYP-3728   # Ihara/Ramanujan: bipartite = the K_3,3 case; the even/odd cycle
  - HYP-3602   # the odd 3-cycle = Condorcet = the tournament intransitivity atom
  - HYP-3615   # the small-measure regime / the equally-spaced extremal = the phi(n) units
results:
  - 04-computation/even_fold_odd_search_klein.py
  - 05-knowledge/results/even_fold_odd_search_klein.out
---

# HYP-3729 — even-n/odd-n is bipartite/non-bipartite; the even-fold

## The even-fold: the equally-spaced descends to the equally-spaced
The 2-adic descent (THM-580, Mode B) of the equally-spaced worst-case `{1..n-1}` for EVEN `n = 2p`:
> odd part `O = {1,3,...,2p-1}`, even part `E = {2,4,...,2p-2}`, `E/2 = {1,...,p-1}` -- the equally-spaced
> for `p`.
So `{1..2p-1}` folds to `{1..p-1}` (verified `n=6->{1,2}`, `8->{1,2,3}`, `10->{1,2,3,4}`, `14->{1..6}`). The
worst-case STRUCTURE is preserved by the descent, and `M` halves (`M(2p)=1/(2p)`, `M(p)=1/p`). For
**`n=14` this reaches `{1..6}` = 6 runners = exactly the PROVEN Barajas-Serra range** (LRC proven for `<=6`
moving runners). So the even-fold connects the `n=14` equally-spaced worst-case to a proven case, and reduces
all EVEN-speed coverings (`S = 2 S'`, `M(S) = M(S')/2`) to the smaller LRC. HONEST: the fold gives the
TIGHTNESS (the equally-spaced achieves `1/n`) and the even-speed reduction; the full LOWER bound (no covering
beats `1/n`) is the hard conjecture direction -- the descent's `M(S) <= M(E/2)/2` is the upper-bound side.

## The even/odd split IS bipartite/non-bipartite (the apex gap)
The cycle `C_n` controls it (HYP-3604): `lambda_min(Q(C_n)) = lambda_min(2I+A(C_n)) = 2-2cos(pi/n)` if `n`
odd, `0` if `n` even (`C_n` bipartite, eigenvalue `-2` at `k=n/2`).
- **EVEN `n`:** `C_n` BIPARTITE, **apex gap `= 0`** -- the equally-spaced is the DEGENERATE measure-0 cusp
  (`M=1/n`, the `phi(n)` units, HYP-3615); the problem folds (Mode B) to the proven `p`-case.
- **ODD `n`:** `C_n` NON-BIPARTITE, **apex gap `> 0`** -- the covering-min is the irregular SPREAD frontier.
This is the SAME split as `sigma`-even/`sigma`-odd (the lonely measure vs the witness), the
orderable/Condorcet (the even cycle is bipartite = 2-colorable = orderable; the odd cycle is the
intransitivity / Condorcet 3-cycle atom, HYP-3602), and the Ihara/Ramanujan bipartite vs Ramanujan
(HYP-3728, `K_3,3` bipartite). Even/odd `n` -> bipartite/non-bipartite `C_n` -> apex gap `0`/`>0`.

## The odd-n beaters are irregular spreads (no clean structure)
The covering-min beaters (mac-mini, MISTAKE-087):
```
n=7: {1,2,5,6,7,8}            M=2/13  (binding D=13, rot 3; speeds*3 mod 13 = {2,3,5,6,8,11}, gaps {1,1,2,2,3,4})
n=9: {1,3,4,5,7,11,18,32}     M=4/33  (binding D=33, rot 4; speeds*4 mod 33, gaps {1,1,2,4,4,5,8,8})
```
At the binding rotation the speeds' images avoid the `M`-neighborhood of `0` (covering radius `= M.numerator`)
but land IRREGULARLY (the gaps are not equal, not a difference set, not an interval) -- which is why naive
perturbation / low-speed search misses them (the large "spreader" speeds, `32 ~ 3.5n` at `n=9`). The smarter
search target is the binding modulus `D` (`13, 33` -- not a clean family) with a min-covering-radius spread;
this remains the open odd-`n` realizability frontier (mac-mini's domain).

## Convergence with mac-mini-S49 (primitivity + the Paley/Ramanujan frame) -- a three-way meet
mac-mini-S49 (HYP-3727) independently resolved the SAME back-and-forth via PRIMITIVITY, and it dovetails:
- **Parity chooses the scale `g`.** The full covering-min `= 1/n` for all `n`, achieved by the NON-PRIMITIVE
  `g.{1..n-1}` with `g = smallest prime factor of n` (the `q=n` witness, EASY): even `n -> g=2`, odd prime
  `-> g=n`, odd composite `-> g=p_min`. So my "even-fold (`g=2`)" is the even case of mac-mini's
  "parity chooses `g`"; the bipartite/non-bipartite apex-gap (`0`/`>0`) is the geometric face of the same
  split. THM-523 (canon) reduces LRC to PRIMITIVE coverings, where `M > 1/n` (the hard margin) -- the
  irregular spread frontier (my odd-`n` beaters).
- **The Paley/Ramanujan frame meets my metazeta (HYP-3728).** The primitive covering-min lives on a circulant
  mod `2n-1`, and `2n-1` is a PALEY vertex count: `n=7 -> 13` (`1 mod 4`, Paley GRAPH, Ramanujan); `n=14 ->
  27 = GF(3^3)` (`3 mod 4`, Paley TOURNAMENT). The "Ramanujan iff Ihara-RH" criterion = the Weil
  `sqrt`-bound on the speed character sums -- exactly the Ihara/Ramanujan machinery of HYP-3728, now on the
  Paley tournament on `2n-1` (the TOURNAMENT side of the metazeta). So my Bass=cut(+)cycle / Ramanujan
  (HYP-3728), my even/odd = bipartite/non-bipartite (this), and mac-mini's primitivity + Paley (HYP-3727)
  are three views of one object: the `2n-1` Paley graph/tournament and its zeta.

## Net
The even-`n`/odd-`n` LRC dichotomy is the bipartite/non-bipartite (`C_n`) dichotomy, the project's master
even/odd split. EVEN `n`: the equally-spaced worst-case folds (2-adic) to the equally-spaced one level down
(`n=14 -> {1..6}`, the proven range), the bipartite cycle has a degenerate (zero) apex gap, and even-speed
coverings reduce -- the structural even-fold is clean (the lower bound is still the hard conjecture
direction). ODD `n`: the non-bipartite cycle has a positive apex gap and an irregular spread covering-min --
the open realizability frontier. The odd cycle is the tournament's Condorcet atom; the even cycle is the
orderable/bipartite one -- the same dichotomy, from tournaments to Ihara zetas to the lonely runner.
