---
id: HYP-3546
title: EVEN PERFECT NUMBERS are exactly the arc-hypercube dimensions d=C(2^p,2)=2^{p-1}(2^p-1) at the FORCADE orders n=2^p with Mersenne-prime apex (Euclid-Euler restated through tournaments); arc-flip = Q_d-EDGE invariance realizes Forcade's per-type parity as the DISCRETE Ky Fan count (grounds OPEN-Q-059, builds on klein THM-584/587); the Mersenne-DRT doubling B_0(T_7)=T_3 is the descending symmetry gauge and the doubly-regular tournament is the Ham-Sandwich score-bisection fixed locus; LRC(14) apex 7 -> Forcade order 8 -> perfect 28 = 2*14
status: (A) THEOREM (elementary, Euclid-Euler); (B) VERIFIED arc-flip parity-invariance n=5 exhaustive (10240 pairs, 0 violations) -- Forcade realized as Q_d-edge-invariance; (C) VERIFIED B_0(T_7)=T_3 via QRs, DRT regular. LRC resonance = structural observation, not a proof.
source: mac-mini-2026-06-29-S13
related:
  - HYP-3545   # mac-mini: oriented-path parity = graded Ky Fan (the object that goes all-odd at n=2^p)
  - THM-587    # klein: the signed cycle index P_n(x) (the GF); closed HYP-3540
  - THM-584    # klein: complement = antipodal map of the arc-hypercube Q_d (d=C(n,2))
  - THM-448    # Mersenne doubling tower B_0(T_{2m-1})=T_{m-1}; skew-Hadamard DRTs
  - THM-585    # mac-mini: vertex-transitive => n|H; A038375 maximizer symmetry
  - THM-586    # mac-mini: Paley H-arithmetic (T_7 = the apex DRT)
  - HYP-3544   # klein: equivariant homology / Betti numbers (the topological refinement)
  - OPEN-Q-059 # tournament Ky Fan
external: Euclid-Euler perfect-number theorem; Forcade 1973; arXiv:2512.09332; Havet-Thomasse N=8
results:
  - 04-computation/perfect_number_forcade_gauge_macmini_20260629.py
  - 05-knowledge/results/perfect_number_forcade_gauge_macmini_20260629.out
---

# HYP-3546 -- perfect numbers are Forcade arc-hypercube dimensions

## (A) THEOREM: even perfect numbers = arc-counts of the Forcade/Mersenne tournaments
The arc-hypercube `Q_d` of order-`n` tournaments has dimension `d = C(n,2)` = the number of arcs.
At the **Forcade orders** `n = 2^p` (where Forcade 1973 / HYP-3545 give: EVERY oriented type has ODD
count -- the maximally degenerate Ky Fan count),
```
d = C(2^p, 2) = 2^p(2^p-1)/2 = 2^{p-1}(2^p-1).
```
By Euclid-Euler, an even number is **perfect** iff it equals `2^{q-1}(2^q-1)` with `2^q-1` prime. Hence
> `d = C(2^p,2)` is a perfect number  <=>  `2^p-1` is a Mersenne prime.
VERIFIED: n=4->d=6 (perfect, 2^2-1=3 prime), n=8->28 (perfect, 7 prime), n=16->120 (not, 15=3*5),
n=32->496 (perfect, 31), n=64->2016 (not, 63), n=128->8128 (perfect, 127). So the **even perfect
numbers are exactly the arc-hypercube dimensions of the Forcade orders whose apex `2^p-1` is a Mersenne
prime** -- Euclid-Euler read through tournaments: a perfect number is the arc-set of a tournament on
`2^p` vertices precisely when the all-odd Ky-Fan degeneracy (n=2^p) coincides with a prime apex.

## (B) Forcade = Q_d-EDGE invariance = the DISCRETE Ky Fan (grounds OPEN-Q-059)
The open piece of OPEN-Q-059 is to realize the per-type parity as a Ky-Fan/Z_2 invariant for ARBITRARY
`T`. The structural realization: **flipping one arc (= one edge of klein's Q_d) preserves `N_tau(T) mod 2`
for every type** -- VERIFIED exhaustively at n=5 (all 10240 (tournament, arc-flip) pairs, ZERO
parity-changing flips). Since `Q_d` is connected by arc-flips, `N_tau mod 2` is CONSTANT on `Q_d`, equal
to its value at the transitive vertex = `#{perms with descent set S}` = Fan's alternating count in the
transitive (`=` Fan's linear-order) gauge. So **Forcade's theorem IS the Q_d-edge-invariance of the
graded Ky Fan count**, and the complement-antipodal `Z_2` of THM-584 is its equivariant structure; the
signed cycle index `P_n` (THM-587) is its generating function. Remaining for a fully elementary proof:
the explicit single-arc-flip involution (the bistellar/path-following move of the Prescott-Su Ky-Fan
proof); the invariance is verified but its local involution is not yet exhibited in closed form.

## (C) The Mersenne-DRT doubling is the descending symmetry GAUGE
THM-448's Mersenne doubling `B_0(T_{2m-1}) = T_{m-1}` (out-neighborhood of vertex 0 = previous level)
is the gauge descent. VERIFIED at the apex prime: Paley `T_7` (the DRT on `7 = 2^3-1`) has
out-neighborhood of vertex 0 = `{1,2,4}` (the quadratic residues mod 7), and the induced sub-tournament
is exactly the 3-cycle `T_3 = B_0(T_7)`. The **doubly-regular tournament is the Ham-Sandwich /
score-bisection fixed locus**: every score `= (n-1)/2` (T_7: all scores 3), the antipodally-balanced
point of the score measure -- the canonical gauge in which the cut-space is bisected (klein's Ham
Sandwich reading of THM-587). At Mersenne PRIME orders the DRT is Paley (vertex-transitive, `n|H`,
THM-585/586); the doubling tower `T_7 -> T_3` shadows the perfect-number descent `28 -> 6` (consecutive
even perfect numbers, both `= C(2^p,2)`).

## (D) The LRC(14) resonance (structural observation)
LRC(14): `14 = 2*7`, apex prime `7 = 2^3-1` (Mersenne). Its Forcade order is `8`, with perfect
arc-count `d = 28 = 2*14 = C(8,2)`. Three coincidences worth cataloging:
- `n = 8` is simultaneously the Havet-Thomasse/Rosenfeld threshold (every oriented Ham path appears,
  N=8) AND the arXiv:2512.09332 arc-deletion threshold (`n>=8`) -- the threshold IS the perfect-number /
  Mersenne-apex / Forcade order.
- The divisors of the perfect number `28` are `{1,2,4} U {7,14,28} = {2^j} U {2^j * 7}` -- exactly the
  tower of orders the LRC apex generates (the apex `7`, the runner count `14`, the powers of 2). The
  perfect-number identity `sum of proper divisors = 28` ties this tower together.
- The Mersenne-DRT doubling `7 -> 3` is the LRC-apex descent `LRC(14) -> LRC(6)` (apex `3 = 2^2-1`,
  perfect `6 = C(4,2) = 2*3`). So the two smallest even perfect numbers `6, 28` are the arc-counts of the
  two smallest LRC apex-DRT orders.

## Leverage of the unused tools (next)
- **Ham Sandwich**: the DRT/regular tournament is the score-bisection fixed locus (verified). For the LRC,
  bisect the cap obstruction `M_odd` (HYP-3538) by the spectral eigenvector -- the antipodally-symmetric
  signed measure cannot be moved to zero by a symmetric perturbation (Borsuk-Ulam/Ham-Sandwich), an
  existence-without-construction certificate complementary to the descent route.
- **Tucker**: a `Z_2`-labeling of `Q_d` under the antipodal complement with no complementary edge is
  impossible (discrete Borsuk-Ulam); a Tucker labeling encoding "which runner binds" would be the finite
  checkable form of the saddle-index obstruction.
