---
id: THM-4737
title: "Kepler polygon clocks with closed one-sided boundary guards, conditional finite polygonal order, and the real-orbit closure formulation"
status: >
  PROVED elementary, with K2/K3/K4 scope corrected 2026-10-20: (K1) by the swept-area formula (1/2)|r x v| on each side;
  (K2) by the per-runner maximum of the two one-sided area rates, including tight vertices; (K3) the vertex criterion,
  finiteness of N_min iff the loneliness threshold is attained, and unbounded finite examples for every fixed k>=2;
  (K4) the closed-orbit-closure/supremum equivalence, with actual-time attainment supplied separately for rank<=1.
  VERIFIED numerical (kepler_polygon_lrc.out): K1 rates at n=3..14 and 20,000 random interior K2 probes; these did not
  verify the missing boundary. FINITE-EXACT (kepler_polygon_lrc.out, polygon_orders.out): (K3) the
  n-gon law 'time 1/n is lonely iff n divides no speed' exhaustive for all primitive k-sets of speeds <= M (k = 2..6);
  exact kappa(D) over the same ranges (no violation of kappa >= 1/(k+1); tight sets exactly {1,2}, {1,2,3}, {1,2,3,4},
  {1,3,4,7}, {1,2,3,4,5}, {1,3,4,5,9} in range, each with lonely measure 0 and N_min = n); N_min histograms and
  maxima (15, 18, 17 for k = 2, 3, 4 in range; records 3, 5, 8, 11, 15, 16, 17, 19 on the pairs (1, m), m <= 30000).
  CITED: THM-639 (the runner/tournament dictionary and the Haar law on the annihilator loop), THM-1820 (relation
  lattice pairing), THM-527 (lonely-density reformulation), the oracle-S517 / codex-S529 / codex-S431 polygon and
  tournament lenses (the tight witness {1..k} at t = 1/n is the regular n-gon minus the observer, the half-turn chord
  tournament), Kronecker-Weyl, Cusick's view-obstruction formulation, Barajas-Serra (seven runners),
  Sungkawichai-Trakulthongchai (CORE-PAPERS: LRC for k <= 12 nonzero speeds, preprint), Kravitz / Giri-Kravitz (spectra).
  NOT PROVED: the lonely runner conjecture for any new k or an all-source finite upper bound on N_min(D).
  REFUTED: HYP-9350's recalled fractional-chromatic equality; D={1,3,4,7} has chi_f=4 and 1/kappa=5.
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

# THM-4737 -- closed loneliness through one-sided area rates

**Correction, 2026-10-20.** Strict area-flow positivity loses boundary-only
witnesses; rational maximization does not prove that the threshold is reached;
and closure intersection alone is a supremum statement. The repaired forms
below retain those coordinates. [Proof, explicit hostiles and exact controls](../../05-knowledge/results/runner_area_clock_20261020j/boundary_and_colouring.md).
The fractional-colouring equality recalled in HYP-9350 is also refuted there.

## 0. The lonely runner conjecture and the polygon

For `D = {v_1 < ... < v_k}` positive integers (or reals), `kappa(D) = sup_t min_i ||v_i t||`
(`||.||` = distance to the nearest integer). The supremum formulation is
`kappa(D) >= 1/(k+1)`. The actual-time formulation asks, with `n = k+1`
runners including the stationary observer, for a time when every runner
is at circular distance at least `1/n`. For integer speeds attainment
is automatic; the real-orbit qualification appears in K4 below.
Literature status and runner-count conventions are maintained in
[CORE-PAPERS](../../05-knowledge/reference/CORE-PAPERS.md); the repository's
LRC(14) programme has `k=13` and is not itself closed by an external preprint.

For n>=3 take the loop to be the **regular n-gon of perimeter 1** (side `1/n`, circumradius
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

## 2. (K2) The closed Kepler criterion

Let `L_i^-(t),L_i^+(t)` be the two one-sided derivatives of the observer's
swept-area clock. On a side interior both equal `d_s v_i/2`; this is the
area rate, half the unsigned angular momentum. Then exactly, including
all polygon vertices,

```text
the observer is lonely at t  <=>  for every i, max(L_i^-(t),L_i^+(t))>0.
```

At a far endpoint of an incident side one rate is zero and the other is
positive; at the observer both are zero. This proves the criterion. The
maximum is taken separately for each runner at the same time. For D={1,2}
at t=1/3, the two runners need opposite sides; no single-sided positivity
test retains this tight witness. Discarding measure-zero endpoints would
discard every lonely time of this example.

Each area clock is continuous, piecewise linear, non-decreasing, with
zero rate on fraction2/n of a lap and mean rate v_i P. For integer speeds
the closed lonely set has the same measure m(D) as its strict interior,
and m(D)>0 iff kappa(D)>1/n. Equality kappa(D)=1/n gives a nonempty finite
lonely set modulo1; kappa(D)<1/n would give the empty set. Zero measure
alone does not distinguish these. The original20,000 random interior
probes did not test the closed endpoints; the linked correction packet
does so exactly. `#`

## 3. (K3) Polygon vertex times and the polygonal order

For integer speeds, at `t = j/N` every runner has phase at a vertex of a
regular `N`-gon representing the phase circle (position `v_i j / N mod 1`).
The observer is lonely at `t = j/N` iff

```text
v_i j mod N  in  [N/n, N - N/n]   for every i.                                     (V)
```

For `N = n` the window `[1, n-1]` is "any vertex but the observer's", so **time `1/n` is lonely
iff `n` divides no speed** (exhaustively confirmed for all primitive sets of speeds `<= M`,
`k = 2..6`). Thus the conjecture is immediate on the `n`-gon for every speed set without a
multiple of `n`; the whole difficulty is the speeds divisible by `n`, which sit on the observer
at every `n`-gon vertex time.

Define the **polygonal order** `N_min(D)` = the least `N >= 2` for which some vertex time
`j/N` is lonely, or infinity if none exists. The finite piecewise-linear maximum
`kappa` is attained at a rational candidate `m/(2v_i)` or `m/(v_i +- v_j)`.
Consequently N_min is finite **iff kappa reaches1/n**; rational attainment
alone does not prove this threshold. It equals `n` on the tight families `{1, ..., k}` (time `1/n`, the regular
`n`-gon minus the observer of oracle-S517) and on `{1,3,4,7}`, `{1,3,4,5,9}`, it is `2` when all
speeds are odd (everyone antipodal at `t = 1/2`). Unboundedness is proved
below; the observed growth is slow in the stated finite ranges:
on the pairs `(1, m)` the records are `N_min = 3, 5, 8, 11, 15, 16, 17, 19` at
`m = 2, 6, 12, 18, 72, 3960, 8640, 16560` (pair_one_m_growth.py), the record `m` being
highly composite; the maxima over all primitive sets of speeds `<= M` are `15` (`k = 2`,
`M = 200`), `18` (`k = 3`, `M = 45`), `17` (`k = 4`, `M = 24`) (polygonal_order_bound.out).
Histograms over all primitive speed sets are in the outputs; the growth law is HYP-9350 (a).

Unbounded finite values do have an elementary unconditional proof. Fix k>=2
and L>=3, put m=lcm(2,...,L), and take D={1,m,2m,...,(k-1)m}. Every
denominator N<=L fails because the speed m is at the observer at every
j/N. But t=1/2+1/(km) is lonely: m is even, the multiple-speed phases
are j/k for1<=j<k, and the speed1 distance is at least5/12>1/(k+1).
Thus L<N_min(D)<=km for this explicit family. This proves the
unboundedness part, not the proposed upper growth bound.

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
For a nonclosed real orbit this equivalence concerns the supremum kappa.
Actual hitting follows from an interior safe point by density, or from
closed-orbit attainment. It does not follow from density and a boundary
intersection alone.
Consequences:

- `r = 0` (Q-independent speeds): `K_D = T^k`, lonely times are free, `m(D) = (1 - 2/n)^k`
  (the classical trivial case, also THM-1820).
- `r = 1`, `L(D) = Z m`: the linear form `m . x` ranges over an interval of length
  `sum_i |m_i| (1 - 2/n)` on `B`; a relation among distinct positive reals has
  `sum_i |m_i| >= 3` (two terms force `a v_1 = b v_2` with `a + b >= 3`), so the length is
  `>= 3(1 - 2/n) >= 1` for `n >= 3`. For n>=4 the length is strictly
  greater than1, so an integer lies in the interval's interior, giving
  an interior safe point and hence an actual time by density. At n=3,
  rank one means the two speeds are commensurable and the orbit is closed.
  Therefore **the conjecture holds for every real speed set with at most
  one independent rational relation.**
- Integer (or commensurable) speeds are the opposite extreme, `r = k - 1`: `K_D` is the Kepler
  curve itself, and the zonotope statement is the classical view-obstruction form in the
  quotient. The Kepler/Haar reduction is exactly the standard reduction of the real conjecture
  to the integer one; the integer case carries all the difficulty.

## 5. Tournaments at the vertices

At a vertex time `j/N`, orient a pair from i to i' when its clockwise
phase difference lies strictly in(0,1/2). Coincident and antipodal pairs
are retained as ties; this is a tournament only when neither kind of tie
occurs. Even a tight regular configuration can have antipodal pairs.
The orientation alone does not determine loneliness: retain the phase
gaps, observer, runner labels and common time. With those data, the
observer is lonely iff both adjacent circular gaps are at least1/n.
Far endpoints of the incident polygon sides are allowed. The exact
same-orientation/different-loneliness hostile is recorded in the
[runner packet](../../05-knowledge/results/runner_area_clock_20261020j/README.md).

## 6. Files

`04-computation/experiments/kepler_lrc_20261020/kepler_polygon_lrc.py` (`.out`: (A) rates,
(B) criterion, (C) exact `kappa`, lonely measure, `N_min` over all primitive `k`-sets of speeds
`<= M`, (D) tight family) and `polygon_orders.py` (`.out`: the `n`-gon law, tight sets with
their `N_min`, extremal `N_min` sets, the fraction of sets settled by polygons of size `<= n`
and by sizes `n, 2n, 3n, 4n`, growth of `max N_min` with `M`).
