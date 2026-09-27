# Friedman's fences, per-fence accounting: the pinwheel lemma, a typed piece-potential LP, and pinwheel constructions

Session collatz-procgen-20260922, lane "fence2", 2026-09-26.  Builds on
`01-canon/theorems/THM-4509-fence-corner-lemma-and-density-bound.md` and its note
`05-knowledge/results/procgen_fencelim_20260926_fence_density.md` (B.1: stall value 0.5167670 and its
hexagon-roundabout obstruction; B.4: the notched pentagon; C.1-C.5).  The setting is the same: `n` unit
fences, no crossings, every end on another fence, fields of area `<= 1`, `lambda = lim A(n)/n`.  Known:
`1/2 <= lambda <= 0.5224525` (PROVED), and `<= 0.5168084` for convex fields (computer-assisted).

**Status header.**  Every finite claim is re-checked by
`04-computation/experiments/procgen_fence2_20260926_run.py`.  Its output
`05-knowledge/results/procgen_fence2_20260926.out` has 15 `[OK]` lines and ends ALL CHECKS PASSED.

| item | status |
|---|---|
| A.1 sector lemma = orchestrator's (P1), extended: every non-straight corner has an end ray; reflex corners have two and occur only at `t = 0` junctions | PROVED |
| A.2 corner types O / I / E / R; landing balance `#O = #I = Lambda` | PROVED |
| A.3-A.4 side lemma (whole side = exactly `j + 1`) and counting identity `#whole - #e0 = #E + #R` per walk (= (P2), generalized to all walks) | PROVED |
| A.5 pinwheel lemma: no whole side iff pure O or pure I; every field with all sides `< 1` is a convex pinwheel; overhang `1 - s`; partner angle `pi - theta` at T/X, less at fans | PROVED ((P2) correct up to two amendments) |
| A.6 fence-side pieces `[1a][0]...[0][1b]` or `[2]`; pinwheel pairing `s + s' = 1`; chain version | PROVED |
| A.7 THM-4509's two extremal objects (hexagon roundabouts, notched regular pentagons) are impossible | PROVED |
| A1-A9 checks on 26 plane configurations, the `n = 6` record and 5 tori | FINITE-EXACT |
| B.1 typed piece-potential certificates: weak duality `density <= 2(beta + rho)` | PROVED (model: torus, no straight join, fields without holes; convex, or with reflex corners and their junction constraints) |
| B.3 LP values (coarse grid): separable model 0.509626; non-separable model **0.501160** (repaired bound 0.5025); with one reflex corner 0.501736 (repaired 0.5046); pure-pinwheel model 0.500072 | EMPIRICAL (grid potentials, heuristic field pricing; torus, no straight joins) |
| C.1 torus identities: density `= det/n`, `C = n - T0`; all-T/X gives density = mean field area | PROVED |
| C.2 Pythagorean tilings `(a^2+b^2)/2`; flipped squares (area 1, perimeter 4) | PROVED / FINITE-EXACT |
| C.2 pinwheel pentagons with roundabouts: `a_P + a_T <= 1 - c/2` (density `< 1/2`) | EMPIRICAL |
| a periodic pattern with density `> 1/2` | **none found** |
| D Conjecture D (a typed certificate with value exactly 1/2) | OPEN; fails at the coarse grid (0.50116); the finer-grid run is unconverged |
| `lambda = 1/2` | OPEN |

**Answer.**
- **(P1)** is correct, and extends to reflex corners (both fences end there).  **(P2)** is correct after
  two amendments (A.5), and holds in a stronger form for every boundary walk of every field:
  `#whole sides - #(e = 0) sides = #E + #R` (A.4).  With the landing balance `#O = #I`, the fence-side
  structure `[1a][0]...[0][1b]` / `[2]` and the pinwheel pairing `s + s' = 1` (A.6), these facts make
  both extremal objects of THM-4509 impossible (A.7).  All PROVED and checked exactly.
- Put into an LP (typed piece potentials on fence pieces; weak duality PROVED), they lower the
  relaxation from 0.5168 to **0.50116** for periodic configurations with convex fields (0.50174 with one
  reflex corner allowed), on a coarse 20-degree grid.  Status EMPIRICAL (heuristic field pricing; no
  straight joins; torus); repaired bounds 0.5025 / 0.5046.
- What remains is a Cairo-like fractional mixture: area-1 pentagons with two whole sides and two
  Y-type corners (perimeter 3.9335), unit squares, and hexagons with nearly flat corners (B.4).
- No periodic pattern beating 1/2 was found.  All pinwheel families stay below it (C.2), and the
  pure-pinwheel model's LP sits at 0.500072, on a degenerate limit of the Pythagorean tilings.
- `lambda = 1/2` stays OPEN.  The smallest statement that would settle the convex periodic case is
  Conjecture D (a typed certificate of value exactly 1/2).

## A. Per-fence structure: the orchestrator's (P1), (P2) and their extensions (PROVED)

### A.0 Conventions

As in THM-4509 / fencelim A.0: fences are closed unit segments; two fences meet in at most one point,
which is an end of at least one of them; every end lies on another fence.  At a junction `j` there
are `e_j` ends and `t_j in {0,1}` through-fences; the rays cut the plane into sectors.  A face is
traversed along its boundary walks with the face on the left.  At a corner (a sector of the face at a
junction) the walk arrives along the *incoming* ray and leaves along the *outgoing* ray.  Put

- `iota(c) = 1` if the fence of the incoming ray ends at the junction, else 0;
- `o(c) = 1` if the fence of the outgoing ray ends at the junction, else 0.

A *side* of a walk is the maximal straight run between two consecutive non-straight corners
`c, c'`; its *end count* is `e(side) = o(c) + iota(c')`.  A *landing* is a junction side of a
through fence (one of the two open half-planes) that carries at least one stem; `Lambda` is the
number of landings.

### A.1 Sector lemma (this is (P1), sharpened)

Let a sector at a junction be bounded by the consecutive rays `r1, r2`.

1. If its angle is not `pi`, at least one of `r1, r2` is an end ray.
2. If it is reflex (`> pi`), then `t_j = 0` and both rays are end rays.
3. If it is straight (`= pi`), either both rays belong to the through fence (type **S0**: `o = iota = 0`)
   or both are end rays of two collinear fences (a *straight join*, type **S2**: `o = iota = 1`).

*Proof.*  Two through rays belong to the same fence (two fences through one point would cross), so
they are opposite and the sector between them, containing no other ray, has angle exactly `pi`.  This
gives 1.  If `t_j = 1`, the two through rays split the plane into two closed half-planes and every
sector lies in one of them, so no sector is reflex; this gives 2.  An end ray is never collinear with
the through fence (that would be an overlap), so a sector bounded by an end ray and a through ray has
angle in `(0, pi)`; this gives 3.  ∎

So **(P1) is correct**, for reflex corners as well (both fences end there).  It fails only for S0
corners (the flat side of a T), which are not corners of the polygon.

### A.2 Corner types and the landing balance

By A.1 every non-straight corner has one of the types

| type | `(o, iota)` | where |
|---|---|---|
| **O** | (1, 0) | extreme sector next to a through ray on a stem side; the walk arrives along the through fence and leaves along the stem |
| **I** | (0, 1) | the other extreme sector of a stem side; arrives along the stem, leaves along the through fence |
| **E** | (1, 1), convex | between two stems of a fan, or at a junction with `t = 0` (Y, conv(e)) |
| **R** | (1, 1), reflex | only at `t = 0` junctions (L, refl(e)) |

**Landing balance (PROVED).** Every landing contributes exactly one O corner and one I corner, and
every O or I corner arises this way.  Hence, over all faces (the outer face included),
`#O = #I = Lambda`.  In a tight configuration (all junctions T, X, Y, L) the O and I corners of a
landing have angles `theta` and `pi - theta`; at a fan (`m >= 2` stems on one side) they sum to less
than `pi`, the rest being E corners.

*Proof.* With the through fence along the x-axis and `m >= 1` stems at angles
`0 < phi_1 < ... < phi_m < pi` on the upper side, the sector `(0, phi_1)` is entered along the first
stem and left along the positive x-ray (type I), the sector `(phi_m, pi)` is entered along the
negative x-ray and left along the last stem (type O), and the middle sectors are bounded by two stems
(type E).  An O or I corner has a through ray, hence lies at a `t = 1` junction next to the through
fence.  ∎

### A.3 Side lemma (the length part of (P2))

Let a side run from `c` to `c'` and contain `j` straight joins (S2 corners) in its interior (S0
corners do not interrupt it).  It is covered by `j + 1` collinear fences `F_0, ..., F_j`, consecutive
ones meeting at the joins; `F_1, ..., F_{j-1}` end at both of their joins, so each is a whole fence of
length 1 on the side.  `F_0` ends at `c` iff `o(c) = 1`, in which case it contributes length exactly
1; otherwise it passes through `c` and contributes a length in `(0, 1)`.  The same holds for `F_j`
and `iota(c')`.  Hence

| `e = o(c) + iota(c')` | length of the side |
|---|---|
| 2 (**whole side**) | exactly `j + 1` |
| 1 | in `(j, j + 1)` |
| 0 | in `(max(0, j-1), j + 1)` |

In particular a side with `e = 2` and no straight join is a whole fence, and every side shorter than 1
has `j = 0` and `e <= 1`.  A whole side has no landing on its own side of the fence (a landing would
be an O or I corner inside it).

### A.4 The counting identity (this is (P2), generalized)

For every boundary walk with non-straight corners `c_1, ..., c_m` and sides `s_i = [c_i, c_{i+1}]`:

`sum_i e(s_i) = sum_i (o(c_i) + iota(c_i)) = #O + #I + 2 #E + 2 #R`,  so
**`sum_i (e(s_i) - 1) = #E + #R`  and  `#{whole sides} - #{e = 0 sides} = #E + #R >= 0`.**

For a convex field this is the orchestrator's `sum e_i >= kappa` (equality iff there is no E corner).
A side is whole iff it goes from a corner of type O/E/R to a corner of type I/E/R.  Everything here is
per boundary walk and per corner, so pinches (a walk visiting a junction twice, with one corner each
time), bridges (a fence with the same face on both sides: the walk runs along both sides) and holes
(extra walks) need no separate treatment; runner A1 checks all of them.

### A.5 Pinwheel lemma (the structural part of (P2))

**(a)** A walk has no whole side iff all its non-straight corners are O, or all are I.  *Proof.*  If
some corner has `o = 1` (types O, E, R), the next non-straight corner must have `iota = 0`, i.e. be O;
inductively every corner is O, and then no E/R corner (which has `iota = 1`) can occur.  Otherwise no
corner has `o = 1`, so all are I.  ∎

**(b)** Such a walk has no reflex corner, so all its turnings are `>= 0`; its total turning is `+2pi`,
so it is the outer walk of a bounded face and a convex polygon traversed once.  A face all of whose
walks have no whole side therefore has a single walk: it is a **convex pinwheel without holes**.

**(c)** Every field all of whose sides are shorter than 1 is a convex pinwheel (whole sides have length
`>= 1`).  The hypothesis "convex" in (P2) is therefore not needed: it is a consequence.

**(d) Structure.**  Let `f` be a pure-O pinwheel with corners `c_1, ..., c_kappa` (counterclockwise),
angles `theta_i` and sides `s_i = [c_i, c_{i+1}]` with no straight join (automatic if `|s_i| < 1`).
Then `s_i` is carried by one fence `F_i` that ends at `c_i` and passes through `c_{i+1}`; it overhangs
`c_{i+1}` by exactly `1 - |s_i|`.  At `c_{i+1}` the stem `F_{i+1}` lands on `F_i`, on the side of
`f`.  The face containing the I corner of that landing (the *partner* `W_{i+1}`):
- has the opposite orientation there (an I corner);
- has angle `pi - theta_{i+1}` if that side of `F_i` carries only this stem at `c_{i+1}` (T or X
  junction), and less otherwise (fan);
- has, along the overhang, a piece of length `<= 1 - |s_i|` (the part of `F_i` from `c_{i+1}` to the
  next landing on that side, or to the end of `F_i`), with equality iff `F_i` carries no further
  landing on that side.

So (P2) is correct with two amendments: "every corner is a T or X junction" should read "every corner
is at a junction with a through fence, as the extreme sector beside it" (T or X in tight
configurations, a fan otherwise), and "the partner corner has angle `pi - theta`" holds exactly when
the landing has one stem, and as `<=` in general.  Pure-I pinwheels are the mirror images.

### A.6 Fence sides and the pinwheel pairing

Orient a fence side so that its faces are on the left.  If it carries `k` landings, its pieces are, in
order, **1a** (from the fence end to an O corner), `k - 1` pieces **0** (from an I corner to an O
corner) and **1b** (from an I corner to the fence end); if `k = 0` it is one whole piece **2**.  The
piece lengths add up to exactly 1.  Without straight joins every face side is exactly one piece.

**Pinwheel pairing (PROVED).**  If at a landing the O corner belongs to a pure-O pinwheel `f` and the I
corner to a pure-I pinwheel `W` (their sides at the landing without straight joins), then that fence
side carries exactly one landing and `|s_f| + |s_W| = 1`.  *Proof.*  `f`'s side ending at the landing
starts at a corner with `o = 1`, i.e. at the fence's end; `W`'s side starting there ends at a corner
with `iota = 1`, the other end of the fence.  ∎

**Chain version (PROVED).**  On every fence side the first piece is a side of the O-field of the
first landing and the last piece a side of the I-field of the last landing.  So if a field has a
0-side `sigma` whose start landing has a pure-O pinwheel as O-field and whose end landing has a pure-I
pinwheel as I-field, then that fence side is exactly `[tau_1][sigma][tau_2]` (two landings) and
`|tau_1| + |sigma| + |tau_2| = 1`.  In particular a 0-side between two *roundabout* corners (corners
whose partners are tiny pinwheels) has length `1 - O(size of the roundabouts)`.

### A.7 Non-convex fields, and what the lemmas kill

For a non-convex field A.4 gives, walk by walk, `#whole sides >= #E + #R >= #R`: a walk with `r`
reflex corners has at least `r` whole sides (each of length `>= 1`), and a hole walk (total turning
`-2pi`) has at least three reflex corners, hence at least three whole sides.

Consequences for the two extremal objects of THM-4509:
- **B1's hexagon roundabouts are impossible.**  A regular hexagon of area 1 has sides `0.6204 < 1`, so it
  is a pinwheel (A.5c), say pure O.  At each of its corners the partner (the I corner of the landing) is
  a corner of a tiny roundabout triangle, whose sides are `< 1`, so it is a pure-I pinwheel.  The
  pinwheel pairing (A.6) gives hexagon side + triangle side `= 1`, so the triangle side is `0.3796`:
  not tiny.
- **B.4's notched regular pentagon is impossible.**  Its two 252-degree dents are R corners, so it has
  at least two whole sides of length `>= 1`; but all its sides are hull edges (`0.7624`) or pocket sides
  (`0.4712`).

## B. The typed piece-potential LP (a relaxation with fence-level constraints)

### B.1 The certificate family (weak duality PROVED; model assumptions stated)

Model: a periodic configuration (torus) in which every field is convex and no face side contains a
straight join.  Then every face side is exactly one fence-side piece (A.6), every corner is O, I or E
(A.1-A.2, no R because every face is a convex field), and every junction is either a landing
junction (`t = 1`) or a `t = 0` junction with only convex sectors (conv(e), `e >= 3`).

A certificate consists of

| symbol | meaning |
|---|---|
| `rho` | price of a fence end |
| `beta` | price of a fence side |
| `Psi1(l, th)` | potential of a **1a** piece of length `l` whose end corner is an O corner of angle `th` (by mirror symmetry the same function for a **1b** piece of length `l` whose start corner is an I corner of angle `th`) |
| `Psi0(l, th_I, th_O)` | potential of a **0** piece of length `l` from an I corner of angle `th_I` to an O corner of angle `th_O` (TypedLP2; the first model used the separable form `P0(l) + T0(th_I) + T0(th_O)`) |
| `Psi2` | potential of a whole piece (length exactly 1) |
| `PhiE(th)` | potential of an E corner |

subject to

- **(F)** every field (every convex polygon with typed corners whose whole sides have length 1 and
  other sides length `< 1`, area `<= 1`): `a <= sum over its sides of the piece potentials + sum over
  its E corners of PhiE`;
- **(S)** every fence side, with any number `k >= 0` of landings at any positions, any fans
  (`m >= 1` stems per landing, the middle sectors being E corners) and any angles:
  `sum of the piece potentials + sum of PhiE over the fan sectors <= beta + rho * (number of stems)`;
- **(J)** every `t = 0` junction: `sum of PhiE over its sectors <= rho e`.

**Weak duality (PROVED).**  Every O corner is the end corner of exactly one piece (its incoming side,
on the through fence), every I corner is the start corner of exactly one piece (its outgoing side),
every E corner is a fan sector or a sector of a `t = 0` junction, and every fence end is a stem of a
landing or an end at a `t = 0` junction.  Summing (F) over the fields and regrouping by fence sides and
`t = 0` junctions gives `A <= 2n beta + 2n rho`, so the density is at most `2(beta + rho)`.  The same
argument works for fields with reflex corners (type R, potential `PhiR(psi)`, charged in (F) and in (J)
at L and refl(e) junctions); fields with holes stay excluded.

The kappa-certificate and every convex angle-potential certificate of THM-4509 embed into the family
(`Psi1(l, th) = mu l + Phi(th)`, `Psi0 = mu l + Phi(th_I) + Phi(th_O)`, `Psi2 = beta = mu`,
`PhiE = Phi`: (F) becomes THM-4509's face constraint, (S) its fan constraints), so with fine enough
grids the family can only improve on 0.5168084 (Theorem B3).  The new ingredients are exactly the
per-fence facts of Part A: whole sides have length 1, a fence side's pieces add up to 1, the O corner
at a landing belongs to the piece before it and the I corner to the piece after it, and their angles add
up to `pi` minus the fan sectors.

### B.2 Computation

Potentials are grid functions (length step `1/NL`, angle step `pi/NA`) with bilinear/trilinear
interpolation.
- (S) and (J) are imposed **exactly** for all real data: on each grid cell the constraints are
  multilinear, so their maxima over the polytopes `{sum l = 1}`, `{sum th = pi}` are attained at grid
  points (vertex argument, as in THM-4509 B3), and all partitions with any number of landings are
  encoded by an exact dynamic program with auxiliary LP variables (`Z[s][jO]`, `Y[s][jI]`, fan and
  junction tables).  An independent DP (`grid_check`) re-verifies them.
- (F) is imposed by cutting planes: for each cyclic type word over {O, I, E} (up to rotation and
  mirror; `kappa <= kmax`), a multistart SLSQP search over angles and side lengths (analytic
  gradients; closure, convexity, `a <= 1`, whole sides `= 1`) looks for the most violated field; plus a
  library of random typed polygons; plus an L1 centering step that stabilises the cutting planes.
- Mirror symmetry: reflecting a configuration swaps O and I and reverses every walk, and the mirror
  image of a certificate is again a certificate.  Averaging, one may take `Psi1` common to 1a and 1b
  pieces and `Psi0(l, a, b) = Psi0(l, b, a)`; then a type word and its mirror image give the same field
  constraint, and words are enumerated up to rotation and mirror.  (Without the symmetry of `Psi0`
  the mirror words would have to be priced separately: an early run violated this and produced an
  invalid certificate; imposing the symmetry left the LP value unchanged.)
- Status: **EMPIRICAL**.  The field pricing is heuristic, and tails are not covered
  (`kappa > kmax`, fans with more than 4 stems on a side, conv(e) with `e > 6`; the runner re-checks
  conv(e) and fans up to 12 exactly and fields up to `kappa = 8` by pricing).  A measured violation
  `delta` is repaired by adding `delta/3` to every piece potential, to `rho` and to `beta` (every
  field has at least 3 sides; a fence side with `k` landings has `k + 1` pieces and at least `k`
  stems), at cost `4 delta/3` in the bound.

### B.3 Results (EMPIRICAL)

All runs: torus; no straight joins; fans with at most 4 stems on a side and conv(e) with `e <= 6`
inside the LP (the runner re-checks larger fans and conv(e) up to 12 exactly, and fields up to
`kappa = 8` by pricing).  "value" is `2(beta + rho)` of the final certificate; "delta" is the largest
violation found by a fresh verification (multistart pricing of every type word, tails included);
"bound" is the repaired value `value + 4 delta/3`.

| model | grid (length, angle) | fields | value | delta | bound |
|---|---|---|---|---|---|
| TypedLP (separable `Psi0`), convex | 0.1, 20deg | `kappa <= 6` | 0.509626 | 2.9e-03 | 0.513541 |
| same, O/I corners only (all junctions T/X) | 0.1, 20deg | `kappa <= 6` | 0.509636 | - | - |
| **TypedLP2 (non-separable `Psi0`), convex** | 0.1, 20deg | `kappa <= 7` | **0.501160** | 9.8e-04 | **0.502468** |
| TypedLP2 + one reflex corner | 0.1, 20deg | `kappa <= 6` (`<= 5` with R) | 0.501736 | 2.2e-03 | 0.504621 |
| TypedLP2, finer angles (unconverged) | 0.1, 10deg | `kappa <= 6` | 0.5007 after two rounds (violations 0.04) | - | - |
| pure-pinwheel model (D.2) | 0.1, 20deg | `kappa <= 7` | 0.500072 | 7.3e-04 | 0.501042 |

A fresh verification with other random starts keeps finding violations of order `1e-3` that the
LP's own pricing missed; the repaired bounds absorb what is found, but a more thorough search could
find larger ones.  This is why every number in this table is EMPIRICAL.

Reference points: kappa-LP 0.5224525; angle potentials 0.5167670 (stall), 0.5168084 (THM-4509 B3,
convex fields, rigorous).  The all-T/X (O/I only) separable model gives 0.509636, the same as with
E corners (0.509626): in the separable model the excess did not need Y junctions at all.

### B.4 The extremal fractional solutions (what the relaxations still allow)

**Separable model (first LP).**  Its optimum (coarse grid, 0.509626; all-T/X model 0.509636, same
structure, `beta = 1/4`, `rho = 0.0048`) is a fractional mixture, per fence:
- 0.47-0.49 tiny pinwheel triangles ("roundabouts", angles such as (60,60,60), (80,50,50), (73,80,27));
- 0.23 near-unit squares of type IOIO (two whole sides, two 0-sides of length `1 - O(eps)`: rungs of
  length 1 between two rails that stop an `eps` beyond the square);
- 0.18-0.21 hexagons of area 1 and perimeter 3.855 (types IIOIIO and IIOOIO, two whole sides, e.g.
  angles (98, 144.4, 117.6) twice and sides (0.568, 0.36, 1) twice), with `x = a - P/4 = +0.036`;
- a few pentagons (EOIOI, IIIOO: area 1, perimeter 3.88);
- fence sides: about half whole, half chains `[1a][0][1b]` with two landings.

**The exploit, and why it is not geometry.**  In a real configuration the hexagon's 0-side of length
0.36 runs between two corners whose partners are roundabouts; by the chain version of the pinwheel
pairing (A.6) the fence carrying it would have length `0.36 + O(eps)`, not 1.  The separable potential
`P0(l) + T0(th_I) + T0(th_O)` cannot see this: it lets the LP attach the hexagon's corner angles to one
fence side and its 0-side length to another.  Making the 0-piece potential a genuine function
`Psi0(l, th_I, th_O)` (TypedLP2) removes the exploit.

**Non-separable model (TypedLP2).**  Linking each 0-piece's length to its two end angles removes the
roundabout-hexagon exploit: the value drops from 0.5096 to 0.50116 (coarse grid; `kappa <= 7`), and
the extremal mixture changes completely.  Per fence, at the final certificate:
- 0.13 unit squares with four E corners (the grid) and 0.03 near-unit squares of type EEIO;
- 0.17 pentagons of type EIEOI (and its variant EOIOI) with area 1 and perimeter 3.9335
  (`x = a - P/4 = +0.0166`): angles about (120, 86, 145, 94, 95), sides (1, 0.55, 0.74, 1, 0.64):
  two whole sides, two E corners (at `t = 0` junctions), one O and two I corners (near-right-angle
  T junctions);
- 0.16 hexagons of types EIEIOI / EIEIEI with area 1 and perimeter 4.035-4.05 (`x` about `-0.01`):
  three whole sides and two nearly flat corners (140-167 deg);
- 0.06 tiny roundabout triangles; 0.44 conv(e) junctions and 1.46 whole fence sides per fence.
The whole excess `0.00116 = sum x_f` comes from the pentagons (`+0.0029`) minus the hexagons
(`-0.0017`).  So what the typed family still allows is a *Cairo-like* mixture: pentagons with two
whole sides and two Y-type corners, fed by unit squares and by hexagons whose extra corners are almost
flat.  What it does not see is global: the angles around each `t = 0` junction must be supplied by
actual neighbouring fields in a consistent cyclic order, and the two sides of a fence are coupled at
its ends (a stem end carries the O corner on one side and the I corner on the other).

**Non-convex extension (one reflex corner, TypedLP2).**  Allowing fields with one reflex corner (at
L or refl(e) junctions; word length up to 6 convex corners, 5 when a reflex corner is present) raises
the value to 0.50174.  The good pentagons EIEOI (perimeter 3.9328) stay; the near-flat hexagons are
replaced by fields of type EIROOI / EIEIRO (area 1, perimeter 4.035-4.05) with a reflex corner of
about 211-217 deg at an L junction, and by limits of EIEOI with a zero-length side at a reflex corner;
conv(e) junctions drop to 0.12 per fence and reflex junctions appear (0.38 per fence).  So the reflex
corners are used only as cheap partners, never to make a field "rounder" than its hull, which is what
the notched pentagon of THM-4509 B.4 did (A.7 forbids it).


## C. Constructions

### C.1 What a periodic pattern beating 1/2 must look like (PROVED)

On the torus every face is a field, so

- **density = det(lattice)/n** (fields tile the torus; the only constraint is `a_f <= 1` for all fields);
- **Euler:** `C = n - T0`, where `T0` is the number of junctions with `t = 0` (a `t = 1` junction adds an
  edge to its through fence, a `t = 0` junction does not).  So a configuration all of whose junctions
  are T or X has `C = n` and **density = mean field area**; each Y, L or grid-type junction removes a
  field.
- `A - n/2 = sum_f x_f` with `x_f = a_f - P_f/4`.  A field with `x_f > 0` has at least 5 convex
  corners, `P_f < 4`, `a_f > pi/4`, at most 3 whole sides (each has length `>= 1`), and at most
  3 E/R corners (A.4).  If it has no whole side it is a convex pinwheel (A.5).

### C.2 Pinwheel families (PROVED where stated)

- **Pythagorean tilings** (verified exactly on tori, runner A3): one pure-O and one pure-I square,
  every landing a pinwheel pair, density `(a^2 + b^2)/2 < 1/2`.
- **Pure pinwheel tilings** (every field a pinwheel; equivalently every junction is T or X and every
  fence side carries exactly one stem, with no straight join): `C = n`, every fence side is a pinwheel
  pair `(s, 1 - s)` (A.6), and the landing balance gives `1/mean(kappa_O) + 1/mean(kappa_I) = 1/2`
  (PROVED).  So if the O-fields are pentagons, the I-fields have mean `kappa = 10/3`, and their sides
  are the complements `1 - s_i` of the pentagons' sides.  Example: regular pentagons of area 1 as
  O-fields force every I-side to be `1 - 0.7624 = 0.2376`, and then (0.4 pentagons and 0.6 small
  I-fields per fence) the density is below 0.43.  The pure-pinwheel LP is in B.3.
- **Pinwheel pentagons with roundabouts** (runner C2, EMPIRICAL).  An O-pinwheel pentagon needs cheap
  partners; the cheapest are tiny I-triangles ("roundabouts").  The pinwheel pairing then forces the
  pentagon's sides before those corners to have length `1 - t_k`.  With three roundabout corners and
  two corners paired with mirror pentagons (sides `s`, `1 - s`), the only arrangement that closed in
  the search is `L s L (1-s) L`, and `max (a_P + a_T) = 0.8065, 0.9014, 0.9503` at triangle diameter
  `c = 0.2, 0.1, 0.05` (runner; an exploratory run gave 0.9801 at `c = 0.02`), i.e. about `1 - c`:
  the gain `a_T ~ c^2` never pays for the linear loss, and the density `(a_P + a_T)/2 < 1/2`.
- **The limit object** (PROVED, exact, runner C1): at `c -> 0` the maximiser is the *flipped square*:
  the unit square with its corner triangle `(1,0),(1,1),(1-s,1)` replaced by the congruent triangle
  with exchanged legs on the other side of the hypotenuse.  It is a convex pentagon with sides
  `(1, s, 1, 1-s, 1)`, area exactly 1, perimeter exactly 4 for every `s in (0,1)`: density-neutral
  (`x = 0`), with three whole sides.  This is fencelim C5's "unit square with a clipped corner".

## D. What would prove lambda = 1/2

### D.1 The smallest statement

**Conjecture D (typed certificate at 1/2).**  There are `rho >= 0` and potentials
`Psi1(l, th)`, `Psi0(l, th_I, th_O)`, `PhiE(th)`, with `Psi2 = beta = 1/4 - rho`, such that (F), (S),
(J) of B.1 hold.

- By B.1 (weak duality, PROVED), Conjecture D implies `density <= 1/2` for every periodic
  configuration with convex fields and no straight join.  Together with the non-convex, hole,
  straight-join and outer-face extensions listed in E, it would give `lambda = 1/2`.
- Conjecture D is equivalent, by LP duality, to "the typed relaxation has value exactly 1/2".
- Values forced by the tight configurations (PROVED): the grid forces `Psi2 = beta` and
  `PhiE(90deg) = rho` (unit square: `4 Psi2 + 4 PhiE(90) >= 1`, `Psi2 <= beta`, conv(4):
  `4 PhiE(90) <= 4 rho`); the Pythagorean tilings with `a -> 1` force `Psi1(1, 90deg) = 1/4` and
  `Psi1(0, 90deg) = 0`.
- So the whole difficulty is concentrated in how `Psi1` and `Psi0` bend away from the perimeter
  price `l/4` near the corner angles of round fields: every field with `x = a - P/4 > 0` must be paid
  for by a negative deviation at its partners, and the pinwheel pairing / chain constraints (S) say
  exactly which partners exist.

### D.2 The orchestrator's pinwheel-pentagon form

For pure pinwheel configurations (every field a pinwheel; then every fence side has exactly one
landing and all junctions are T or X), (S) reads `Psi1(l, th) + Psi1(1 - l, pi - th) <= 1/4`.
Writing `Psi1 = l/4 + g(l, th)`, the statement "every pinwheel field pays its excess to its
partners" becomes

  **(P)**  `g(l, th) + g(1 - l, pi - th) <= 0` and, for every convex polygon with sides `s_i < 1`,
  angles `th_i` and area `<= 1`:  `a - P/4 <= sum_i g(s_i, th_{i+1})`.

(P) implies density `<= 1/2` for pure pinwheel configurations (PROVED: sum over fields, regroup by
fence sides).
- The one-parameter choice `g = -k cos(th) l (1 - l)` makes the pairing an identity, but it fails:
  the best `k` (about 0.25) leaves a pentagon with `a - P/4 - sum g = +0.0028` (EMPIRICAL).
- The LP over grid functions `g` (pure-pinwheel model: every field a pinwheel, every fence side one
  pinwheel pair; coarse grid, `kappa <= 7`) converges to **0.500072** (EMPIRICAL; the LP's own pricing
  ends with violations below `1e-8`, a fresh verification with other random starts finds `7e-4`, so
  the repaired value is 0.50104).  Its extremal mixture is degenerate: 0.5 tiny roundabout triangles per
  fence and 0.49 pinwheel "pentagons" of area 1 that are unit squares with one corner cut by a side of
  length 0.1 (angles about 172, 90, 92, 97, 89; perimeter 3.998, `x = +0.0005`).  So for pure
  pinwheel configurations the per-fence accounting stays within `7e-5` of 1/2, and the residue sits on
  a degenerate limit of the Pythagorean tilings; (P) with value exactly 1/2 is plausible but OPEN.

### D.3 What the optimal certificates look like (EMPIRICAL)

The final certificate (TypedLP2, coarse grid) is the perimeter price with a signed angle correction:
`Psi1(l, th) - l/4` is about `-0.02 .. -0.013` for acute corners (`th <= 60deg`) and
`+0.01 .. +0.02` for obtuse ones (`th >= 100deg`) at mid lengths, and within `0.003` of 0 at `l = 0`
and `l = 1`, so the pinwheel pairing `Psi1(l, th) + Psi1(1 - l, pi - th)` never exceeds
`beta + rho = 0.25058`; `Psi2 = beta = 0.23859`, `rho = 0.01199`, and `PhiE` is about `rho` for
`th >= 80deg` and 0 below.  The smooth one-parameter version `g = -k cos(th) l (1 - l)` (D.2) is the
simplest function of this shape; it fails by 0.0028.

## E. What is OPEN

1. `lambda = 1/2`?  Known: `1/2 <= lambda <= 0.5224525` (PROVED); `<= 0.5168084` for convex fields
   (THM-4509 B3, computer-assisted).  This lane: the per-fence relaxations of B.3 (EMPIRICAL).
2. Make B rigorous: an exact (interval / cell) verification of the field constraints (F) for typed
   polygons, as THM-4509 B3 did with L'Huilier for untyped convex fields, plus the tails
   (`kappa > kmax`; fans and conv(e) beyond the checked sizes).
3. Extend the model: fields with holes or several reflex corners, straight joins (a side made of
   several pieces), and finite configurations (the outer face), as THM-4509 did with its (O)
   constraints.  Until then B.3 bounds periodic densities only.
4. Conjecture D (a typed certificate of value exactly 1/2).  At the coarse grid it fails by 0.00116.
   A run with 10-degree angles, started from the coarse cuts, stood at 0.5007 after two rounds with
   violations still 0.04 and was stopped (8 minutes per round), so the effect of refinement is
   unsettled; it may hold in the limit.  A stronger family would add the linkage the present one
   ignores: the two sides of a fence at each of its ends (a stem end carries the O corner on one side
   and the I corner on the other), i.e. a fence-level LP.
5. A construction beating 1/2 would have to beat every relaxation above; none of the pinwheel families
   does (C.2), and the remaining fractional optimum (B.4) is a Cairo-like mixture whose realisability
   is open.

## F. Reproduction

`python3 -u 04-computation/experiments/procgen_fence2_20260926_run.py > 05-knowledge/results/procgen_fence2_20260926.out`
(about 297 s; peak RSS 66 MB; a second run gave identical output up to the timing line).

Files, all `04-computation/experiments/procgen_fence2_20260926_*`:
- `torus.py`: exact engine for periodic (torus) and plane configurations (Fractions, `Q(sqrt 3)` or
  50-digit numbers); corner typing O / I / E / R / S0 / S2; sides with end counts;
- `struct.py`: the checks A1-A9 of Part A;
- `constr.py`: Pythagorean tori, flipped squares, the pentagon-roundabout family;
- `lp.py`: the typed piece-potential LPs (`TypedLP`: separable 0 pieces; `TypedLP2`: non-separable),
  exact grid encodings, field pricing, certificate I/O and verification;
- `lpdriver.py`: the cutting-plane driver that produced the certificates (run from a scratch
  directory; its docstring lists the run sequence);
- `cert_sep.json`, `cert_main.json`, `cert_nonconvex.json`, `pincert.json`: the certificates (grid
  potentials and the extremal polygons with their LP weights).
The runner also uses the fencelim lane's `procgen_fencelim_20260926_geom.py` and `..._audit.py`
(read only).  Logs, pickles and notes of the LP runs are in `scratch/procgen_fence2/`.
