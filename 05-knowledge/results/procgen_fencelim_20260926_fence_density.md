# Friedman's fences: the density constant lambda, the Corner Lemma audit, angle-potential certificates, and constructions

Session collatz-procgen-20260922, lane "fencelim", 2026-09-26.  Builds on
`procgen_smallgraph_20260926_small_graphs_arithmetic.md` section 2 (Propositions F1-F5, bound `U(n)`).

**Status header.** All finite facts are re-checked by
`04-computation/experiments/procgen_fencelim_20260926_run.py`, whose output
`05-knowledge/results/procgen_fencelim_20260926.out` has 19 `[OK]` lines and 0 failures.

| item | status |
|---|---|
| (1) junction inequality `e_j + t_j + r_j >= 3` | PROVED. Tight types: T, **X**, Y, L (the orchestrator's list missed X) |
| (2) face angle identity; global form `2C = sum_j (k_j - alpha_j) + 2H + 2c_o` | PROVED |
| (3) Corner Lemma `sum_fields (kappa_f - 3) = sum_j (3 alpha_j - k_j)/2 - kappa_o - 3H - 3c_o <= n - 3` | PROVED; equality for the unit square and the n = 5, 6 records |
| (4)+(5) `A + 2 mu sqrt(pi A) <= 2(mu+rho) n - 6 rho`, so `lambda <= (6-P_5)/(8-P_5) = 0.5224525` | PROVED |
| checks: 26 exact configurations + the reproduced n = 6 record; all 48 records satisfy the bound | FINITE-EXACT / VERIFIED |
| B1: every angle-potential certificate has `2(mu+rho) >= (4 - 12^(1/4))/(6 - 12^(1/4)) = 0.5167670 > 1/2` | PROVED. The family cannot prove `lambda = 1/2` |
| B2: the torus/convex potential LP attains 0.5167670; extremal primal = regular hexagons "with roundabouts" + unit squares | VERIFIED (LP) |
| B3: if every field has a convex outer boundary, then `A + 2 mu sqrt(pi A) <= 0.5168084 n` | PROVED (computer-assisted, float64 margin >= 2e-9) |
| B4/B6: non-convex fields bounded through their hull give back the kappa-LP value (crude, water-filling, and exact hull-angle bounds); the obstruction is the notched regular pentagon | EMPIRICAL; a pocket-geometry bound is OPEN |
| C1: if every field has at most 4 convex corners then `A < n/2`; if every field has at least 5 then `A <= (n-3)/2` | PROVED |
| C2-C5: Cairo geometry, unit-side bounds, mean field 0.5105, Cairo-balanced typed pentagons never beat the square | EMPIRICAL / VERIFIED |
| a periodic pattern with density `> 1/2` | **none found**; `lambda = 1/2` is OPEN |

**Answer in one line.**
- `1/2 <= lambda <= 0.5224525` (PROVED).
- For configurations whose fields are convex, `lambda <= 0.5168084` (PROVED).
- The angle-potential method provably cannot reach 1/2. Its extremal relaxation violates fence lengths.
- No construction beats 1/2.

## A. Audit of the orchestrator's bound (1)-(5)

### A.0 Conventions

A configuration is a finite set of `n` closed unit segments (fences) such that two fences meet in at
most one point, which is an end of at least one of them, and every end lies on another fence.  A
*junction* is a point where some fence ends.  At a junction `j` there are `e_j >= 1` ends and
`t_j in {0,1}` fences passing through `j` in their interior (two such fences would cross).  Every end
lies on another fence, so `t_j = 0` forces `e_j >= 2`.  The rays from `j` along the fences are
`d_j = e_j + 2 t_j >= 2` pairwise distinct directions (two equal directions would be an overlap).
They cut the plane around `j` into `d_j` sectors of angles `theta > 0` summing to `2 pi`.
`r_j = #{theta >= pi}`, `k_j = #{theta < pi}`, `alpha_j = 2 - r_j`.  Since the angles sum to `2pi`,
`r_j <= 2`, with `r_j = 2` only for a straight join (`d_j = 2`, two sectors `pi`).
Note `pi alpha_j = sum_{theta<pi} theta + sum_{theta>=pi} (theta - pi)`.

The plane graph has the junctions as vertices and the fence pieces between consecutive junctions
as edges.  Every vertex has degree `d_j >= 2`, so every component contains a cycle.  Faces are the
components of the complement.  Fields are the bounded faces; `C` is their number.  A *corner* of a
face is a sector (at some junction) lying in that face; a face may own several corners at one
junction (pinch).  `kappa_f` = number of corners of `f` with angle `< pi` (convex corners), `h_f` = number
of boundary walks of `f` minus 1 (holes), `H = sum h_f`, `c_o` = number of boundary walks of the
unbounded face (= components not enclosed by a field), `kappa_o` = its number of convex corners.

### A.1 Junction inequality (PROVED)

`3 alpha_j - k_j = 6 - 3 r_j - (d_j - r_j) = 6 - 2 r_j - e_j - 2 t_j`, so
`3 alpha_j - k_j <= e_j  <=>  e_j + t_j + r_j >= 3`.

*Proof.* If `t_j = 1`, the through fence splits the plane into two half-planes of angle `pi` each.
If one side carries no end, that side is a single sector of angle exactly `pi`, so `r_j >= 1`; since
`e_j >= 1` we get `e_j + t_j + r_j >= 3`.  Otherwise both sides carry ends and `e_j >= 2`.
If `t_j = 0`, then `e_j >= 2`; for `e_j = 2` the two sectors sum to `2pi`, so one is `>= pi`.  ∎

**Equality cases** (`e + t + r = 3`): `(e,t,r) = (1,1,1)` **T** (one stem on a through fence);
`(2,1,0)` **X** (two stems on opposite sides of a through fence, meeting at the same point: allowed,
the two stems meet in a common end); `(2,0,1)` **L** (two ends, not collinear); `(3,0,0)` **Y** (three
ends, all sectors `< pi`).  The orchestrator's statement "tight exactly for T, Y, L" omits X.  It is
harmless for the bound: X is equivalent in all counts to two T's.
Checked on every junction of the configurations of A.6 (junction types are classified exactly;
X occurs in the `1 x 2` rectangle example).

### A.2 Face angle identity (PROVED)

For a bounded face `f`, push each boundary walk slightly into `f`.  This gives a smooth domain with
`1 + h_f` boundary curves, which is homotopy equivalent to `f`.  By Gauss–Bonnet (Hopf's
Umlaufsatz), the total turning of these curves is `2 pi chi(f) = 2 pi (1 - h_f)`.  A corner of angle
`theta` contributes turning `pi - theta`; this holds also at pinches and along bridges.  Every vertex
has degree `>= 2`, so no corner has angle `2 pi`.  Hence

`sum_{corners of f} theta = pi (m_f - 2 + 2 h_f)`   (`m_f` = number of corners),

which is `sum_{theta < pi} theta + sum_{theta >= pi} (theta - pi) = pi (kappa_f - 2 + 2 h_f)`.

For the unbounded face, cut it off by a large circle: the Euler characteristic is `1 - c_o` and the
circle turns by `2pi`.  So `sum (pi - theta) = -2 pi c_o`, i.e.
`sum_{theta<pi} theta + sum_{theta >= pi}(theta - pi) = pi (kappa_o + 2 c_o)`.
Every sector is a corner of exactly one face.  Summing over all faces:

`pi sum_j alpha_j = pi (K - 2C + 2H + 2 c_o)`, with `K = sum_j k_j`, i.e. `2C = sum_j (k_j - alpha_j) + 2H + 2c_o`.

### A.3 Corner Lemma (PROVED)

`sum_f kappa_f = K - kappa_o`.  Substituting `C` from A.2:

`sum_fields (kappa_f - 3) = sum_j (3 alpha_j - k_j)/2 - kappa_o - 3H - 3 c_o <= (1/2) sum_j e_j - 3 = n - 3`,

using A.1, `sum_j e_j = 2n` (every end is at a junction), `kappa_o, H >= 0`, and `c_o >= 1`.
Equality needs every junction to be T, X, Y or L, no convex corner on the outer face, no hole, and
one component.  The unit square, the `n = 5` record and the reproduced `n = 6` record all attain it
(A.6).

### A.4 Per-face inequality (PROVED)

Let `P_k = 2 sqrt(k tan(pi/k))` be the perimeter of the regular `k`-gon of unit area
(`P_3 = 4.559014`, `P_4 = 4`, `P_5 = 3.8119353`, `P_6 = 3.7224194`, decreasing to `2 sqrt(pi)`).

* *Hull vertices are convex corners.* Let `v` be an extreme point of `conv(cl f)`.  It is a polygon
  vertex.  Every sector of `f` at `v` lies in a supporting half-plane.  A sector of angle `>= pi`
  would be exactly that half-plane, and then `v` would lie inside a segment of `cl f`, so it would not
  be extreme.  Hence every sector of `f` at `v` is a convex corner, and
  `#vertices(conv f) <= kappa_f`.  Also `kappa_f >= 3`.
* *Perimeter.* Let `P_f` be the boundary length of `f` counted with multiplicity: both sides of a
  bridge inside `f`, and the hole boundaries.  Then `P_f >= perimeter(conv f)`.  By the discrete
  isoperimetric inequality for convex polygons with at most `kappa_f` vertices (Zenodorus/L'Huilier;
  CITED, classical) and monotonicity of `P_k`:
  `P_f >= P_{kappa_f} sqrt(area(conv f)) >= P_{kappa_f} sqrt(a_f)`.
* *The scalar inequality.* With `mu = 1/(8 - P_5) = 0.2387738` and
  `rho = (4 - P_5)/(2(8 - P_5)) = 0.0224525`, `s = sqrt(a) in [0,1]`:
  `s^2 - mu P_k s` is convex in `s`, so its maximum on `[0,1]` is `max(0, 1 - mu P_k)`.  One needs
  `1 - mu P_k <= 2 rho (k - 3)`:
  * `k = 3`: `mu P_3 = 1.08857 > 1` (slack `0.0886`);
  * `k = 4`: `1 - 4 mu = (4 - P_5)/(8 - P_5) = 2 rho`, **equality**;
  * `k = 5`: `1 - mu P_5 = (8 - 2P_5)/(8 - P_5) = 4 rho`, **equality**;
  * `k = 6`: `0.11119 <= 0.13471` (slack `0.0235`);
  * `k >= 7`: `1 - mu P_k < 1 - 2 sqrt(pi) mu = 0.15357 < 8 rho = 0.17962`.

  (A4 also scans `k <= 400`.)  Hence `a_f <= mu P_f + 2 rho (kappa_f - 3)` for every field,
  with equality only for the unit square and the unit-area regular pentagon.

### A.5 The finite bound (PROVED)

Summing A.4 over the fields and using A.3 (note `rho > 0` and `kappa_f - 3 >= 0`):
`A <= mu sum_f P_f + 2 rho (n - 3)`.  Every fence has two sides, each in one face, so
`sum_f P_f + P_o = 2n`, where `P_o` is the boundary length of the unbounded face (checked exactly in
A.6).  Let `U` be the union of the closed fields.  A non-junction point of `dU` lies on a fence that
has a field on at most one side, so `H^1(dU) <= P_o`.  The isoperimetric inequality gives
`H^1(dU) >= 2 sqrt(pi A)`.  Therefore

**`A + 2 mu sqrt(pi) sqrt(A) <= 2(mu + rho) n - 6 rho`,  and  `lambda <= 2(mu + rho) = (6 - P_5)/(8 - P_5) = 0.52245246`.**

The equality analysis of A.3–A.4 makes this the optimum of the "kappa-LP" relaxation.  In that
relaxation a field is described only by its number of convex corners, subject to
`sum (kappa_f - 3) <= n` and `sum P_f <= 2n`.  Per fence, the optimum uses `x_5 = 1/(8-P_5) = 0.4775`
regular pentagons of area 1 and `x_4 = 2 rho = 0.0449` unit squares.  The dual certificate is the
Corner Lemma with prices `(mu, rho)`.

### A.6 Exact checks (FINITE-EXACT; runner lines A1-A9)

Engine `procgen_fencelim_20260926_geom.py`: exact arithmetic over `Q` or `Q(sqrt 3)` (exact signs),
or 50-digit arithmetic for the numerically optimised `n = 6` record.  For each configuration it
checks the rules (no crossing, no overlap, no loose end, unit lengths) and builds the plane graph by
half-edge traversal.  It assigns holes by component, classifies every sector exactly, and verifies:
- A.1 at every junction;
- the per-face angle identity (numerically, to `10^-30`);
- the global identity of A.2 and the Corner Lemma identity of A.3, both exactly;
- `sum (kappa_f - 3) <= n - 3`, `sum P_f + P_o = 2n`, A.4 on every field, and A.5.

| configuration | n | C | H | c_o | kappa_o | sum(kappa_f-3) | A | junction types |
|---|---|---|---|---|---|---|---|---|
| grids 1x1 ... 4x4 (7 grids) | 4..40 | 1..16 | 0 | 1 | 0 | = C | = C | L, e3t0r1, e4t0r0 |
| quasi-square k-ominoes k=3,5,7,10,13 | 10..34 | k | 0 | 1 | 1 | k | k | L, e3t0r1, e4t0r0 |
| L-tromino, S-tetromino, 8-cell ring (enclosed centre) | 10,13,24 | 3,4,9 | 0 | 1 | 1,2,0 | 3,4,9 | 3,4,9 | |
| n=5 record: unit square + chord (0,3/5)-(4/5,0) | 5 | 2 | 0 | 1 | 0 | **2 = n-3** | 1 | L4 T2 |
| square + two chords / square + mid chord | 6 / 5 | 3 / 2 | 0 | 1 | 0 | **n-3** | 1 | L, T |
| 1x2 rectangle with an X junction | 9 | 4 | 0 | 1 | 0 | 4 | 2 | L, e3t0r1, T, **X** |
| brick wall (running bond, half bricks) | 25 | 10 | 0 | 1 | 0 | 10 | 9 | L, e3t0r1, T |
| regular hexagon + 3 spokes (Y), Q(sqrt3) | 9 | 3 | 0 | 1 | 0 | 3 | 2.598 | Y, L, e3t0r1 |
| unit triangle nested in unit square, Q(sqrt3) | 7 | 2 | **1** | 1 | 0 | 1 | 1 | L |
| two squares + bridge fence | 9 | 2 | 0 | 1 | 4 | 2 | 2 | L, T |
| two squares sharing a corner (pinched outer face) | 8 | 2 | 0 | 1 | 2 | 2 | 2 | L, e4t0r0 |
| figure-8 component in a 4x4 frame (pinched hole; fields > 1) | 24 | 3 | 1 | 1 | 0 | 5 | 16 | L, straight joins, e4t0r0 |
| two disjoint squares | 8 | 2 | 0 | **2** | 0 | 2 | 2 | L |
| n=6 record, reproduced (unit-sided pentagon + unit chord) | 6 | 2 | 0 | 1 | 0 | **3 = n-3** | 1.4758530 | L5 T2 |

The `n = 6` record was re-optimised (SLSQP, then Newton at 50 digits).  Its fields are a pentagon of
area exactly 1 and a quadrilateral of area 0.4758530, with total `1.47585296`; this matches the page's
`1.47585+`.

**Pythagorean pinwheel tilings cannot be made finite by restriction (PROVED).**  The Pythagorean
tiling by squares of sides `a` and `1 - a` is a fence tiling.  Its maximal segments have length 1,
and every junction is a T, so its density is `(a^2 + (1-a)^2)/2 < 1/2`.  Take any nonempty finite
subset in which every end is supported.  An end of a fence in the subset is supported only by the
fence it met in the tiling, so every junction of the subset is again a T.  But a finite
configuration has at least 3 end-only junctions (the extreme points of the hull, Proposition F1).
Hence the support-closure of any finite patch is empty.  Pruning disc patches with `a = 3/5, 4/5`
indeed ends empty.  A finite pinwheel therefore needs a frame with L/Y corners.  The torus version
is checked in part C.

**Records (A8, A9).**  All 48 best-known values `3 <= n <= 50` (Friedman's page, fetched
2026-09-26, sha256 `78504fff...fa27`, same page as the smallgraph lane) satisfy
`A(n) <= B(n)`, where `B(n)` is the largest root of A.5.  The smallest gap `B(n) - record` is
`0.0768`, at `n = 4`.  `B(n)` improves the previous lane's `U(n)` by 5–10% (`B/U` from 0.897 to
0.946).

| n | record | B(n) | U(n) | record/B |
|---|---|---|---|---|
| 3 | 0.43301 | 0.71628 | 0.79881 | 0.605 |
| 4 | 1 | 1.07677 | 1.17348 | 0.929 |
| 5 | 1 | 1.45615 | 1.56854 | 0.687 |
| 6 | 1.47585 | 1.84903 | 1.97853 | 0.798 |
| 7 | 2 | 2.25219 | 2.40010 | 0.888 |
| 8 | 2.10306 | 2.66351 | 2.83097 | 0.790 |
| 10 | 3.04687 | 3.50512 | 3.71457 | 0.869 |
| 12 | 4 | 4.36608 | 4.62070 | 0.916 |
| 13 | 4.16199 | 4.80229 | 5.08047 | 0.867 |
| 17 | 6.01086 | 6.57636 | 6.95415 | 0.914 |
| 24 | 9.02394 | 9.75983 | 10.32699 | 0.925 |
| 33 | 13.01887 | 13.94535 | 14.77450 | 0.934 |
| 42 | 17.02275 | 18.19754 | 19.30250 | 0.935 |
| 50 | 20.55829 | 22.01632 | 23.37474 | 0.934 |

(The full table for `n = 3..50` is printed by the runner.)

## B. Angle-potential certificates

### B.0 The family

A certificate is a triple `(mu, rho, Phi)`:
- `mu >= 0` is the price of perimeter and `rho >= 0` the price of a fence end;
- `Phi: (0, 2 pi) -> R` is a corner potential with `Phi(pi) = 0`.  This is no loss: a field may carry
  any number of straight corners, so `Phi(pi) >= 0` is forced, and every face constraint below is
  implied by the same constraint with the straight corners removed.

The requirements:

- **(J)** At every junction, `sum over its sectors of Phi <= rho e_j`.  It suffices to check four
  families:
  - *fan(m)*: `m >= 1` stems on one side of a through fence, i.e. `m+1` convex sectors summing to
    `pi`, with bound `rho m`.  T is `m = 1`.  X and every junction with stems on both sides are sums
    of two fans, using `Phi(pi) = 0`.
  - *conv(e)*: `e` ends, all sectors convex, summing to `2 pi`.  Y is `e = 3`.
  - *refl(e)*: `e` ends with one reflex sector.  L is `e = 2`.
  - Ends-only junctions with a straight sector are implied by `fan(e-2)`.
- **(F)** Every field satisfies `a_f <= mu P_f + sum over its corners of Phi`.
  - For a *convex* field with corner angles `theta_i`, L'Huilier's inequality (CITED; sanity check
    B5) gives `P^2 >= 4 a g`, where `g = sum cot(theta_i/2)`.  Since `a - 2 mu sqrt(a g)` is convex
    in `sqrt a`, (F) becomes `sum Phi(theta_i) >= max(0, 1 - 2 mu sqrt(g))`.
  - Non-convex fields go through the convex hull (B.4).
- **(O)** `Phi(theta) >= a + c theta` for convex `theta`, and `Phi(psi) >= c(psi - pi)` for reflex
  `psi`, with `c >= 0` and `a + c pi >= 0`.  By the angle identity A.2, every outer-face walk and every
  hole walk then has `sum Phi >= (a + c pi) kappa + 2 pi c >= 0`.

Summing (F) over the fields and using (J) and (O):
`A <= mu (2n - P_o) + 2 rho n`, hence `A + 2 mu sqrt(pi A) <= 2(mu + rho) n` and `lambda <= 2(mu+rho)`.

The kappa-certificate of part A is the member `Phi_0(theta) = rho (3 theta/pi - 1)` for convex
`theta` and `Phi_0(psi) = 3 rho (psi - pi)/pi` for reflex `psi`.  On any junction its sum is
`rho(3 alpha_j - k_j)`, and on any field it is `2 rho (kappa_f - 3)`.

### B.1 Theorem B1: the family stalls above 1/2 (PROVED)

**Theorem B1.** Every certificate in the family satisfies

`2(mu + rho) >= lambda_Phi := (8 - P_6)/(12 - P_6) = (4 - 12^(1/4))/(6 - 12^(1/4)) = 0.51676701...`

In particular, no angle-potential certificate proves `lambda = 1/2`.

*Proof.* Five of the requirements already force it:
- `T(60,120)`: `Phi(60) + Phi(120) <= rho`;
- `T(90,90)`: `2 Phi(90) <= rho`;
- the unit square: `4 Phi(90) >= 1 - 4 mu`;
- the regular hexagon of area 1 (`g = 2 sqrt 3`, `2 sqrt g = P_6 = 2 * 12^(1/4)`):
  `6 Phi(120) >= 1 - mu P_6`;
- the tiny equilateral triangle (`a -> 0`): `3 Phi(60) >= 0`.

Adding six times the first to the hexagon and triangle constraints gives `6 rho >= 1 - mu P_6`.  The
square and the second give `2 rho >= 1 - 4 mu`.  Minimising `mu + rho` over these two half-planes
puts the optimum at their intersection, `mu = 2/(12 - P_6)`, `rho = (1 - 4mu)/2`, with value
`2(mu+rho) = (8 - P_6)/(12 - P_6)`.  ∎  (runner B1)

**The obstruction is a fractional tiling.** It uses:
- regular hexagons of area 1, each of whose `120°` corners sits at a T junction;
- the `60°` partner corners, absorbed by tiny triangles, so that three hexagons meet around a small
  pinwheel "roundabout" instead of at a Y;
- unit squares at `T(90,90)` junctions.

It satisfies every counting and angle constraint the family can see: the Corner Lemma, end counts,
angle matching at junctions, and L'Huilier.  It is **not** realisable.  Look at a hexagon corner at a
roundabout T.  Along the stem, the side starts at a fence end.  Along the through fence, the fence
runs on to the next roundabout, where it ends.  So every hexagon side is a whole fence minus
roundabout-sized bits, of length about 1, and the hexagons would have area about 2.6.  The
potentials never see fence lengths.

### B.2 The LP attains the stall value (VERIFIED)

`procgen_fencelim_20260926_potlp.SampledLP` restricts to angles on a grid (exact angle sums), convex
fields, `fan/conv` junctions and the torus.  It gives
- `0.51676716` at 6° (runner B2);
- `0.516767` at 3° and 2° (exploration runs).

This equals `lambda_Phi`, and its dual solution is exactly the fractional tiling above: weights
1.4497 `T(60,120)` + 0.5503 `T(90,90)` per fence; faces 0.2416 hexagons, 0.2752 squares, 0.4832 tiny
triangles.  So `lambda_Phi` is not only a lower bound on the family.  On the convex/torus relaxation
it is the family's optimum (up to discretisation; see B3 for the rigorous version).

### B.3 Theorem B3: a rigorous bound for convex fields (PROVED, computer-assisted)

**Theorem B3.** Suppose every field of a configuration has a convex outer boundary (holes and the
unbounded face are unrestricted).  Then `A + 2 mu sqrt(pi A) <= 0.5168084 n` with `mu = 0.2416`.
In particular, the density of such configurations is at most `0.5168084`.

*Certificate.*
- `Phi` is piecewise linear with nodes at multiples of 2° on `[0, pi)`.  There is a separate node
  value at `pi-`, and `Phi(pi) = 0`.  The reflex potential is piecewise linear on `(pi, 2pi)`.
- It is found by a cutting-plane LP (`PLLP`, `nmode='none'`, 1612 rows) and re-verified after a
  safety bump (`certify_convex`).  Every family has slack `>= 2e-9`, far above the float64 error
  (`< 1e-13`) of the sums involved.

*Why finitely many checks suffice.*
- **Junctions.** Fix the cells `[k_i d, (k_i+1) d]` of the sectors of a junction family.  The set
  `{theta in cells : sum theta = S}`, with `S` a multiple of `d`, is a polytope whose vertices are
  grid points: `k - 1` coordinates sit at cell ends, and the last is then an integer multiple of `d`.
  `sum Phi` is linear on the cell, so its maximum is at a vertex.  The junction constraints for
  **all** real angles therefore reduce exactly to node tuples; dynamic programming checks every one.
- **Fields.** On a cell, `g = sum cot(theta_i/2)` is convex, so it is at least its tangent plane `l`
  at the cell centre.  `sqrt g >= Lt(g)`, where `Lt` is the piecewise-linear interpolant of `sqrt` on
  `[1, 5]`, capped at `sqrt 5`.  `Lt` is concave and nondecreasing, so
  `sum Phi + 2 mu Lt(l(theta))` is concave on the cell and attains its minimum at a vertex.  At a
  vertex, `l >= g - E` with the Taylor term `E = (1/2)(d/2)^2 max (cot(t/2))''`.
- **Small angles.** A field with an angle below 22° has `g > cot(11°) >= 5`, so only `sum Phi >= 0`
  is needed there (and `2 mu sqrt 5 >= 1`).
- **Tails.** `kappa > 10`, `m > 10`, `e > 9` are closed by the linear minorant (O) and a linear
  majorant `Phi <= A + B theta`.

(runner B3; slacks printed there.)

### B.4 Non-convex fields: where the method stops (EMPIRICAL / OPEN)

For a non-convex field only the convex hull is available.  The polygon itself obeys no
L'Huilier-type inequality: a narrow notch adds convex corners of tiny angle but almost no perimeter.

- **Crude hull bound (`P >= P_kappa sqrt a`).** With this bound for non-convex fields, the LP returns
  to the kappa-LP value: 0.5225297 at 4° (runner B4).  The binding primal consists of "star" fields
  (five 72° tips and reflex dents), whose hull the crude bound treats as a regular pentagon.
- **Water-filling hull bound** (`g_hull >= min{sum tan(x_v/2) : sum x_v = 2 pi, 0 <= x_v <= pi - theta_v}`):
  still too weak.
- **Exact hull-angle formulation (runner B6).** This is the strongest hull-based treatment tried.
  Let `phi_v` be the hull angle at a convex corner `v`.  Then:
  - `phi_v = theta_v` at hull vertices that are not pocket ends;
  - `phi_v in [theta_v, pi)` at pocket ends, and there are at most `2r` of them, since every pocket
    between consecutive hull vertices contains a reflex corner;
  - `phi_v = pi` for corners inside pockets.

  Moreover `sum phi_v = (kappa-2) pi`, `g_hull = sum cot(phi_v/2)`, and the reflex excess is
  `E = sum (phi_v - theta_v)`.  The reflex potentials are bounded below by linear minorants
  `gamma_r + beta_r E` (separately for `r = 1, 2, >= 3`).  The constraint matrix (differences
  `phi - theta`, one sum row, boxes) is totally unimodular, so node tuples suffice.

  Result: **0.5227722 at 6°, 0.5225764 at 4° — no gain.**
  - An earlier run gave 0.5213.  It was an artefact: a cut-key collision stopped the cutting planes
    early.  Fixed; the final re-separation now reports worst violation `< 1e-16`.
  - The binding primal is regular pentagons plus **notched regular pentagons**.  Such a field has a
    regular-pentagon hull, four 72° tips as pocket ends, one exact 108° corner, and two 252° dents
    supplied by L junctions whose convex side is 108°.  The junctions are `T(72,108)`, and there are
    unit squares at `T(90,90)`.
  - The hull relaxation scores such a field as a regular pentagon (area 1, `P = P_5`).  But each
    pocket is a 36-36-108 triangle standing on a full hull edge.  At unit hull area each costs area
    0.105 and adds perimeter 0.179, so the real field has `a/P` about 0.19.

**So for general fields the method is blocked by the hull relaxation itself, not by its constants.**
Notches decouple a field's corner angles from its roundness, at no cost in the Corner Lemma, because
L is tight.  **OPEN:** a pocket-geometry inequality.  A pocket spanning a hull edge of length `b` with
end angles `d_v, d_w` loses area of order `b^2 sin d_v sin d_w / sin(d_v + d_w)` and gains perimeter
of order `b`.  With such an inequality B3 should extend to all configurations.

**Verdict on B.**
- The angle-potential family **stalls at exactly `lambda_Phi = 0.5167670 > 1/2`** (PROVED lower bound,
  LP-attained).  So it can never prove `lambda = 1/2`.
- It does give a proved improvement, 0.52245 → 0.51681, for fields with convex outer boundaries.
- For general fields it is currently stuck at the kappa-LP value.
- The information it lacks is fence length: every fence is a unit segment, and its pieces on each side
  sum to exactly 1.

## C. Constructions and obstructions

### C.1 Proved obstructions

**Proposition C1 (PROVED).**
- (a) A field with at most 4 convex corners has `a <= P/4`, because `P >= P_kappa sqrt a >= 4 sqrt a`
  and `a <= 1`.  Hence if every field has at most 4 convex corners, `A <= (2n - P_o)/4 < n/2`.
- (b) If every field has at least 5 convex corners, the Corner Lemma gives `2C <= n - 3`, so
  `A <= C <= (n-3)/2`.

So a configuration with `A > n/2` must mix round fields (at least 5 corners, `P < 4` at area near 1)
with fields of at most 4 corners that relieve the corner budget.  This is exactly the pentagon/square
mix of the kappa-LP optimum.  (runner C1)

**Proposition C6 (PROVED).**  Take the Pythagorean fence tilings (squares of sides `a` and `1-a`,
maximal segments of length 1, all junctions T).  Per period they have 2 fences, 2 fields and
4 T junctions, so `sum (kappa_f - 3) = 2 = n` (Corner Lemma equality on the torus).  Their density is
`(a^2 + (1-a)^2)/2 < 1/2`.  No finite patch is support-closed (A.6).

### C.2 The ladder of relaxations (each extremal primal fails the next constraint)

| relaxation | value | extremal primal | killed by |
|---|---|---|---|
| kappa-LP (Corner Lemma + isoperimetry) | 0.5224525 | 0.4775 regular pentagons (area 1) + 0.0449 unit squares per fence | angles: regular-pentagon corners cannot pair at T/Y junctions |
| angle potentials (B1, B2) | 0.5167670 | regular hexagons with roundabout corners + unit squares + tiny triangles | fence lengths: hexagon sides are whole fences |
| unit-side mean field (C4, angles ignored) | 0.5104818 | pentagons with 3 Y + 2 T corners (3 whole-fence sides, `P >= 3.92288`) + P(2Y,3T) + 2.1% unit squares | angles (C5) |

- *Unit sides (C3; cyclic-polygon theorem, CITED).*  At a Y (or any ends-only) corner both sides
  are fences ending there, so a face side between two such corners is a whole fence: length exactly 1
  in tight configurations, at least 1 in general.  A convex face with `y` such corners has at least
  `y` unit sides.  The least perimeter of an area-1 pentagon with `u` unit sides is
  `3.92288, 3.87306, 3.83819, 3.81194` for `u = 3, 2, 1, 0`.
- *Mean field.*  In a tight torus configuration, `n = sum_f (y_f/2 + t_f/4)` (`y_f` end-end corners,
  `t_f` through-end corners), and the angle balance `sum_f (2 y_f + 3 t_f - 12) = 0` is the Corner
  Lemma.  With the unit-side perimeters, the resulting LP gives 0.51048 (runner C4).
- *Angles kill it in the Cairo case (C5, EMPIRICAL).*  Take a pentagon with 3 Y + 2 T corners whose
  own corners close up its junctions: Y angles summing to 360° and T angles to 180°, as in the Cairo
  tiling.  With its 3 forced whole-fence sides and area at most 1, the maximum of `a/P` is exactly
  `1/4`, over all four stem patterns and 100 starts each.  It is attained by a unit square with a
  clipped corner; it never exceeds the square.
- *The Cairo tiling itself (C2).*  4-valent vertices at `(0,0), (1,1)` (period 2), 3-valent vertices
  at `(a, a-1)` with `a = (3 + sqrt 3)/6`.  Its pentagons have area exactly 1, angles 90°/120° and
  perimeter 3.8637.  But a straight line Y-X-Y through a 4-valent vertex has length 1.633 > 1, so no
  fence decomposition exists.  Any Cairo-type pattern with fences must shrink the through lines to
  length 1, and then the all-pentagon count (C1(b)) caps it at 1/2.

### C.3 What a construction would need (OPEN)

- By C1 and the Corner Lemma, a pattern with density `> 1/2` must, per fence:
  - be almost tight: `sum (kappa_f - 3) = n` minus the slack;
  - be pentagon-heavy;
  - carry a few percent of 4-corner fields;
  - have most pentagons of area almost exactly 1 with `P < 4`.
- Pentagons with 3 end-end corners have 3 whole-fence sides, so `P >= 3.9229` and
  `a - P/4 <= 0.019`.  Rounder pentagons need more through-end corners.
- The remaining constraint that none of the relaxations sees is **fence bin-packing**.  Each of the
  `2n` fence sides is a bin of length exactly 1, partitioned by the stems on that side into face-sides:
  - at most one face-side longer than 1/2 fits in a bin;
  - whole-fence sides fill a bin alone.
- Honeycomb, Cairo and regular-pentagon relaxations all violate this bin structure.
- **Result: no periodic pattern with density `> 1/2` was found, and none of the relaxations yields
  one.**
- The margins are small: 0.52245, 0.51677, 0.5105, and 1/4 exactly in the Cairo-balanced class.
  Together this is **evidence (EMPIRICAL, not proof) for `lambda = 1/2`**.
- A proof would need fence-level (bin) accounting.
- A disproof would need a large-period pentagon/square pattern that satisfies the bin structure with
  pentagon area within about 1% of 1.

## D. Typing to the repository

| relation | type | comment |
|---|---|---|
| Corner Lemma ↔ Euler/discharging counting (Kuratowski lane: `E <= 3V - 6`; smallgraph F1: `#fields = n + c - V0`) | **REAL** (same mechanism) | Euler/Gauss–Bonnet plus a local inequality at vertices (`e + t + r >= 3`, cf. face degree `>= 3`) gives a global linear bound.  The Corner Lemma is a discharging identity with charges `3 alpha_j - k_j` at junctions and `kappa_f - 3` on faces |
| `(mu, rho, Phi)` ↔ potential certificates of THM-4486 (mean-payoff game values certified by potentials) | **ANALOGY** (shared method: LP duality with potentials) | both are dual certificates checked exactly on finite families, and both have a primal (a strategy there, a fractional tiling here) that shows where the method stalls.  No mathematical object is shared |
| tight junctions {T, Y, L} ↔ the owner's triples {Petersen, K_{3,3}, K_5} | **ANALOGY at most, and the list is a quadruple** | the equality cases of `e + t + r = 3` are **T, X, Y, L**.  They are *optimal* local types, not obstructions, and X is two T's merged.  Any "triple" reading is NUMEROLOGY |
| discrete half (junction/end counts, Corner Lemma) + continuous half (Zenodorus/L'Huilier, planar isoperimetry) | **REAL** | the bound is literally an LP combining the two.  Each refinement (angles, unit sides) adds a discrete-geometric coupling between them |
| the stall value contains the honeycomb constant `12^(1/4)` (`P_6 = 2 * 12^(1/4)`) ↔ Hales's honeycomb theorem | **REAL identity / ANALOGY in mechanism** | the regular hexagon is the roundest polygon whose corners can all be 120° T-corners, with the 60° partners in tiny roundabouts.  Hales's theorem (arXiv math/9906042, CITED as context, not used) is about partitions into equal-area cells.  It would give only the weaker `12^(-1/4) = 0.537` |

## E. What is OPEN

1. `lambda = 1/2`?  (Friedman's question.)  Known: `1/2 <= lambda <= 0.5224525`.  For configurations
   whose fields are convex, `lambda <= 0.5168084`.
2. Extend B3 to non-convex fields.  Every hull-based bound fails on notched regular pentagons, so
   this needs a pocket-geometry inequality (area and perimeter lost to a notch on a hull edge).  If it
   is found, the expected result is `lambda <= 0.51681` in general.
3. A fence-level (bin-packing) certificate.  This is the only ingredient the relaxation ladder of C.2
   has not used.  A certificate that reaches 1/2 must use it.
4. A construction beating 1/2 must be a large-period, almost tight pentagon/square pattern (C.3).
   A natural test class is the all-T sub-family (every end a stem, `F = n` on the torus, Y junctions
   replaced by pinwheel roundabouts).  It admits pinwheel pentagons with no whole-fence side.

## F. Reproduction

`python3 -u 04-computation/experiments/procgen_fencelim_20260926_run.py` takes about 20 s, with peak
RSS 169 MB.  Output: `05-knowledge/results/procgen_fencelim_20260926.out` (19 `[OK]`, 0 failures).

Helpers, all `04-computation/experiments/procgen_fencelim_20260926_*`:
- `geom.py`: exact planar engine over `Q`, `Q(sqrt 3)` or 50-digit arithmetic;
- `audit.py`: the configurations, including the `n = 6` record re-optimisation;
- `records.py`: Friedman's table;
- `potlp.py`: the rigorous `PLLP` + `certify_convex`, the exploration `SampledLP`, and the non-convex
  (crude, water-filling, pocket) variants;
- `typed.py`: typed pentagons, cyclic unit-side bounds, Cairo data, mean field.

Scratch notes: `scratch/procgen_fencelim/NOTES.md`.
