---
id: THM-4737
title: "The Kepler-polygon reformulation of the lonely runner conjecture. Put the n runners (one stationary observer at a vertex) on the regular n-gon of perimeter 1, so that one side is the loneliness threshold 1/n. (K1) Kepler's second law holds exactly for constant-speed motion on the polygon: about the centre every runner sweeps area at the constant rate a v/2 (a = apothem), and about the observer's vertex the rate on side s is d_s v/2 with d_s = R(cos(pi/n) - cos((2s+1)pi/n)), which vanishes exactly on the two sides incident to the observer; per lap the area swept about either point is the polygon area P. (K2) Kepler criterion: the observer is lonely at time t iff every runner has nonzero angular momentum about the observer at time t, i.e. iff every observer-centred area clock Phi_i(t) is strictly increasing; each clock is flat for exactly the fraction 2/n of its period (equal sides, equal times) and the lonely times are the common strict-increase set, of measure m(D) with m(D) > 0 iff kappa(D) > 1/n. (K3) Polygon vertex times: at t = j/N the runners sit on vertices of the regular N-gon and the observer is lonely iff every d j mod N lies in [N/n, N - N/n]; for N = n this says n divides no speed, so the n-gon alone proves the conjecture for every speed set with no multiple of n (at t = 1/n); the polygonal order N_min(D) (least N with a lonely vertex time) is a new invariant of the speed set, equal to n on the tight families {1..k} and {1,3,4,7}, {1,3,4,5,9}, and unbounded in general. (K4) For real speeds the Kepler flow t -> tD on the k-torus is Haar-equidistributed on its closure, the annihilator torus of the relation lattice L(D); the conjecture for D is equivalent to the zonotope M [1/n, 1-1/n]^k (M a basis matrix of L(D)) containing an integer point, which holds whenever rank L(D) <= 1 (the linear form of a relation has range of length >= 3(1 - 2/n) >= 1 on the box), while integer speeds are exactly the maximal-rank case in which the closure is the Kepler curve itself. Elementary throughout; the content is the dictionary and the invariant, not a proof of the conjecture"
status: >
  PROVED (pre-audit; elementary): (K1) by the angular-momentum formula (1/2)|r x v| on each side (distance from the
  point to the line of the side), (K2) from (K1) and the side geometry, (K3) by reading t = j/N on the N-gon, (K4) by
  Kronecker-Weyl on the annihilator torus and the range of a linear form on a box (relations among distinct positive
  reals have coefficient sum >= 3).
  FINITE-EXACT (kepler_polygon_lrc.out, polygon_orders.out): (K1) rates and per-lap areas numerically exact at n = 3..14;
  (K2) 20,000 random (D,t) pairs with no mismatch between ||v t|| >= 1/n and 'no runner on an incident side'; (K3) the
  n-gon law 'time 1/n is lonely iff n divides no speed' exhaustive for all primitive k-sets of speeds <= M (k = 2..6);
  exact kappa(D) over the same ranges (no violation of kappa >= 1/(k+1); tight sets exactly {1,2}, {1,2,3}, {1,2,3,4},
  {1,3,4,7}, {1,2,3,4,5}, {1,3,4,5,9} in range, each with lonely measure 0 and N_min = n); N_min histograms and
  maxima (15, 18, 17 for k = 2, 3, 4 in range; records 3, 5, 8, 11, 15, 16, 17, 19 on the pairs (1, m), m <= 30000).
  CITED: THM-639 (the runner/tournament dictionary and the Haar law on the annihilator loop), THM-1820 (relation
  lattice pairing), THM-527 (lonely-density reformulation), the oracle-S517 / codex-S529 / codex-S431 polygon and
  tournament lenses (the tight witness {1..k} at t = 1/n is the regular n-gon minus the observer, the half-turn chord
  tournament), Kronecker-Weyl, Cusick's view-obstruction formulation, Barajas-Serra (seven runners),
  Sungkawichai-Trakulthongchai (CORE-PAPERS: LRC for k <= 12 nonzero speeds, preprint), Kravitz / Giri-Kravitz (spectra).
  NOT PROVED: the lonely runner conjecture for any new k; any bound on N_min(D); the fractional-chromatic bridge is
  recorded in HYP-9350 as CITED-from-memory pending verification.
source: mac-mini-2026-10-20-kepler (owner's prompt: lonely runner in the general sense, each runner count on a convex regular polygon loop, Kepler's second law used creatively, tournaments and Hadwiger-Nelson in view)
depends_on:
  - 01-canon/theorems/THM-639-hamiltonian-path-frame-for-runner-families.md
  - 01-canon/theorems/THM-1820-lrc-is-a-moment-nullcone-problem-relation-lattice-pairing.md
related:
  - 01-canon/theorems/THM-527-lrc-lonely-density-reformulation-and-the-bounded-spread-compact-reduction.md
  - 05-knowledge/results/lrc_polygon_tournament_wheel_s517.out
  - 05-knowledge/results/lrc_polygon_outside_inside_s529.out
  - 05-knowledge/results/lrc_neighbor_tournament_lift_s431.out
  - 05-knowledge/hypotheses/HYP-9350-polygonal-order-kepler-flow-and-the-hadwiger-nelson-bridge-of-the-lonely-runner-conjecture.md
scripts:
  - 04-computation/experiments/kepler_lrc_20261020/kepler_polygon_lrc.py
  - 04-computation/experiments/kepler_lrc_20261020/polygon_orders.py
  - 04-computation/experiments/kepler_lrc_20261020/polygonal_order_bound.py
  - 04-computation/experiments/kepler_lrc_20261020/pair_one_m_growth.py
outputs:
  - 04-computation/experiments/kepler_lrc_20261020/polygonal_order_bound.out
  - 04-computation/experiments/kepler_lrc_20261020/pair_one_m_growth.out
  - 04-computation/experiments/kepler_lrc_20261020/kepler_polygon_lrc.out
  - 04-computation/experiments/kepler_lrc_20261020/polygon_orders.out
---

# THM-4737 -- lonely means every runner is sweeping area about you

**PROVED (pre-audit, elementary) + FINITE-EXACT** as itemized in the status.

## 0. The lonely runner conjecture and the polygon

For `D = {v_1 < ... < v_k}` positive integers (or reals), `kappa(D) = sup_t min_i ||v_i t||`
(`||.||` = distance to the nearest integer). The conjecture (Wills 1967, Cusick 1973) is
`kappa(D) >= 1/(k+1)`: with `n = k+1` runners including the stationary observer, at some time
every runner is at circular distance at least `1/n` from the observer. Proved for `k <= 6`
(Barajas-Serra) and, by a computer-assisted sieve, for `k <= 12` (CORE-PAPERS, preprint); the
repo's LRC(14) programme is the case `k = 13`.

Take the loop to be the **regular n-gon of perimeter 1** (side `1/n`, circumradius
`R = 1/(2n sin(pi/n))`, apothem `a = R cos(pi/n)`, area `P = a/2`), with the observer at
vertex `0` and runner `i` at arc-length position `v_i t mod 1`, counter-clockwise, at constant
speed `v_i` (laps per unit time). Side `s` runs from vertex `s` to vertex `s+1`; a runner is
on side `s` when `s/n <= v_i t mod 1 < (s+1)/n`.

## 1. (K1) Kepler's second law on the polygon

For a point moving at constant speed along a straight line, the rate at which it sweeps area
about a fixed point `Q` is `(1/2)|r x v| = (1/2) dist(Q, line) * speed` -- constant on the
line (angular momentum about `Q` is conserved in free motion). On side `s` of the regular
`n`-gon:

- about the **centre**, `dist = a` for every side, so the rate is `a v_i / 2`, **constant in
  time**: the regular polygon is the loop on which constant arc speed IS Kepler motion about
  the centre (equal areas in equal times), and the n equal triangles over the sides are the
  equal-time cells `[s/n, (s+1)/n)`;
- about the **observer's vertex** `O = (R, 0)`, the distance to the line of side `s` is

```text
d_s = a - R cos((2s+1) pi / n) = R [cos(pi/n) - cos((2s+1) pi/n)],      d_0 = d_(n-1) = 0,
```

  so the rate is `d_s v_i / 2`, **zero exactly on the two sides incident to the observer**
  (the runner moves along a line through `O`) and positive on the other `n - 2` sides.

Per lap the swept area about either point is the polygon area `P` (the triangle fan from any
interior or boundary point of a convex polygon covers it once): `sum_s d_s / (2n) = P`.
Checked numerically to machine precision for `n = 3..14` (kepler_polygon_lrc.out (A)). `#`

## 2. (K2) The Kepler criterion

Runner `i` is within circular distance `< 1/n` of the observer iff it is on an open incident
side (positions in `(-1/n, 1/n)` mod 1 are exactly the two incident sides minus their far
vertices). Hence, with `L_i(t) = d_{s_i(t)} v_i / 2` the angular momentum of runner `i`
about the observer:

```text
the observer is lonely at t   <=>   L_i(t) > 0 for every i
                             <=>   every observer area clock Phi_i(t) = int_0^t L_i is strictly increasing at t,
```

up to the measure-zero boundary (a runner exactly at the far vertex of an incident side is
at distance exactly `1/n`, lonely-admissible, with `L_i = 0` at that instant). Each `Phi_i`
is continuous, piecewise linear, non-decreasing, with slope `0` for exactly the fraction `2/n`
of each lap (two of the `n` equal-time sides) and mean slope `v_i P`. The lonely set
`Lambda(D) = {t : all L_i(t) > 0}` has measure `m(D) = |{t in [0,1) : all ||v_i t|| >= 1/n}|`
(for integer speeds the picture is 1-periodic), and `m(D) > 0 <=> kappa(D) > 1/n`, while
`kappa(D) = 1/n` exactly (the **tight** case) means the lonely set is finite: the runners
touch the far vertices simultaneously. 20,000 random `(D, t)` tests show no mismatch between
the arithmetic and the geometric criterion (kepler_polygon_lrc.out (B)). `#`

So the conjecture reads: *for every speed set, at some time every runner is sweeping area
about the observer.* The observer is "lonely" in Kepler's sense -- it is a focus about which
nobody is momentarily radial.

## 3. (K3) Polygon vertex times and the polygonal order

At `t = j/N` every runner sits on a vertex of the regular `N`-gon inscribed in the same loop
(position `v_i j / N mod 1`). The observer is lonely at `t = j/N` iff

```text
v_i j mod N  in  [N/n, N - N/n]   for every i.                                     (V)
```

For `N = n` the window `[1, n-1]` is "any vertex but the observer's", so **time `1/n` is lonely
iff `n` divides no speed** (exhaustively confirmed for all primitive sets of speeds `<= M`,
`k = 2..6`). Thus the conjecture is immediate on the `n`-gon for every speed set without a
multiple of `n`; the whole difficulty is the speeds divisible by `n`, which sit on the observer
at every `n`-gon vertex time.

Define the **polygonal order** `N_min(D)` = the least `N >= 2` for which some vertex time
`j/N` is lonely. It is a well-defined invariant of the speed set (every integer speed set has
a lonely rational time, because `kappa` is attained at a rational time `m/(2v_i)` or
`m/(v_i +- v_j)`), it equals `n` on the tight families `{1, ..., k}` (time `1/n`, the regular
`n`-gon minus the observer of oracle-S517) and on `{1,3,4,7}`, `{1,3,4,5,9}`, it is `2` when all
speeds are odd (everyone antipodal at `t = 1/2`), and it grows without bound, but slowly:
on the pairs `(1, m)` the records are `N_min = 3, 5, 8, 11, 15, 16, 17, 19` at
`m = 2, 6, 12, 18, 72, 3960, 8640, 16560` (pair_one_m_growth.py), the record `m` being
highly composite; the maxima over all primitive sets of speeds `<= M` are `15` (`k = 2`,
`M = 200`), `18` (`k = 3`, `M = 45`), `17` (`k = 4`, `M = 24`) (polygonal_order_bound.out).
Histograms over all primitive speed sets are in the outputs; the growth law is HYP-9350 (a).

## 4. (K4) Real speeds: the Kepler flow, its closure, and the zonotope form

For real `D` the trajectory `t -> t D mod 1` in the torus `T^k` is Kepler motion in the sense
of Liouville: it preserves Haar measure on its closure `K_D`, the annihilator torus
`{x : m . x in Z for all m in L(D)}` of the relation lattice `L(D) = {m in Z^k : m . D = 0}`
(Kronecker-Weyl; this is THM-639 (C) and THM-1820's dictionary). Let `M` be an `r x k` basis
matrix of `L(D)` (`r = rank L(D)`, saturated, so `K_D` is connected of dimension `k - r`).
Then

```text
kappa(D) >= 1/n   <=>   K_D meets the closed box B = [1/n, 1-1/n]^k   <=>   M B contains a point of Z^r,
```

the second equivalence because `x in B` with `M x in Z^r` is exactly `x in B cap K_D`, and
`M B` is the zonotope `M c + sum_i (1 - 2/n) M e_i [-1/2, 1/2]` (`c` the centre of the box).
Consequences:

- `r = 0` (Q-independent speeds): `K_D = T^k`, lonely times are free, `m(D) = (1 - 2/n)^k`
  (the classical trivial case, also THM-1820).
- `r = 1`, `L(D) = Z m`: the linear form `m . x` ranges over an interval of length
  `sum_i |m_i| (1 - 2/n)` on `B`; a relation among distinct positive reals has
  `sum_i |m_i| >= 3` (two terms force `a v_1 = b v_2` with `a + b >= 3`), so the length is
  `>= 3(1 - 2/n) >= 1` for `n >= 3` and the interval contains an integer: **the conjecture holds
  for every real speed set with at most one independent rational relation.**
- Integer (or commensurable) speeds are the opposite extreme, `r = k - 1`: `K_D` is the Kepler
  curve itself, and the zonotope statement is the classical view-obstruction form in the
  quotient. The Kepler/Haar reduction is exactly the standard reduction of the real conjecture
  to the integer one; the integer case carries all the difficulty.

## 5. Tournaments at the vertices

At a vertex time `j/N` the occupied vertices carry the half-turn chord orientation
(`i -> i'` when `i'` is clockwise within half a turn of `i`), the vertex-time tournament of
codex-S431; for the tight family at `t = 1/n` it is the rotational tournament on the regular
`n`-gon minus the observer (oracle-S517), and THM-639's dictionary (steps = Hamiltonian path,
balanced lattice = tiles, step reversal = the grid mirror) is the runner-world reading of the
tiling model. The Kepler clock adds one thing to that dictionary: the observer's loneliness is
the condition that its two neighbours in the circular order both lie at least one side away,
i.e. the observer is a vertex of the vertex-time tournament with no neighbour on an incident
side; nothing beyond this identification is claimed here.

## 6. Files

`04-computation/experiments/kepler_lrc_20261020/kepler_polygon_lrc.py` (`.out`: (A) rates,
(B) criterion, (C) exact `kappa`, lonely measure, `N_min` over all primitive `k`-sets of speeds
`<= M`, (D) tight family) and `polygon_orders.py` (`.out`: the `n`-gon law, tight sets with
their `N_min`, extremal `N_min` sets, the fraction of sets settled by polygons of size `<= n`
and by sizes `n, 2n, 3n, 4n`, growth of `max N_min` with `M`).
