# Session 19: exact polygon/areal-clock transport for Lonely Runner

Status: PROVED elementary transport statements and exact hostile controls;
no general Lonely Runner proof. Staged outside the repository, no Git writes.
Source check: 2026-10-21. Repo baseline inspected: `de0463bc8ddf`.

## Immediate result and literature correction

The polygon/Kepler idea has a precise positive core. A regular polygon with a
number of sides divisible by the runner count transports the **sharp loneliness
threshold exactly** from perimeter/area phase to polar angle. Taking an even
multiple also preserves the half-turn tournament, including its ties. The
transport changes the dynamics but does not prove that the target is reached.

Do not retain the repository's dated external headline that fourteen total
runners is the first unknown case without qualification. Jaan Allikvere's
[arXiv:2609.02604v2](https://arxiv.org/html/2609.02604v2), dated 2026-09-24,
states Theorem 1.1 for fourteen and fifteen total runners. This is a
computer-assisted preprint, not independently audited here. Its section 6
explicitly says the full computation has not been independently reimplemented
or formally verified. It links code and certificates at
[fourteen-runner archive](https://doi.org/10.5281/zenodo.22066772) and
[fifteen-runner archive](https://doi.org/10.5281/zenodo.22667683).
General LRC is not proved by that paper.

The inherited baseline is Sungkawichai–Trakulthongchai,
[arXiv:2604.23906v2, Theorem 1.3](https://arxiv.org/html/2604.23906v2),
which proves the cases through twelve nonzero speeds, or thirteen total
runners. Their notation `LRC(k)` counts nonzero speeds; the repository's
`LRC(14)` counts total runners. Allikvere's `LRC(14)` therefore denotes the
fifteen-total-runner case, not the repository's fourteen-total-runner case.

## 1. Inheritance and live board

Closest inherited mechanism: `HYP-2009-lrc-polygon-outside-inside-arcs.md`
retains the observer clasp, outside cyclic gaps, half-turn chord orientations,
and endpoint/diameter ties. It is OPEN, not a proof. Its predecessor
`HYP-1995-lrc-twin-roots-of-unity-bridge.md` distinguishes geometric circulants
from Paley tournaments and records even regular polygons as tie walls.

Canonical hostile: the tight arithmetic progression has boundary-only lonely
times; positive-measure or strict-inequality substitutes lose these. The Oct 5
correction in `01-canon/MISTAKES.md` further warns that for speeds
`1,...,N-1`, the lonely times mod one are the primitive fractions `a/N`, not
all polygon vertices. The exact N=6 control below returns only 1/6 and 5/6.

Corrected near miss: regular-score or Hamiltonian-path statistics do not
determine the geometric target. `THM-1820-intransitivity-is-schur-convexity-not-covariant-vanishing.md`
itself contains a correction banner retracting its claim that the Hamiltonian
path count is a function of the score sequence. No statement below depends on
its uncorrected prose.

Least-used relevant sidecar: the distinction between the fixed-center lattice
hit and all translated lattice cosets in the zonotope quotient. It is exactly
the common-origin versus arbitrary-start distinction.

Board:

1. Anchor: N-runner orbit, closed threshold 1/N, marked observer and common time.
2. Geometric lens: symmetric polygon, normalized swept area, polar angle.
3. Tournament: half-turn orientation, collision/antipode ties, metric slack.
4. Niche: fixed-coset zonotope hit versus all-coset covering.
5. Wildcard: simultaneous preservation of rational thresholds by LCM symmetry.

## 2. Definitions

Write T=R/Z and d(x,y)=min_{k in Z}|x-y-k|. There are N >= 2 total runners, with distinct original velocities,
with positions x_i(t)=v_i t mod 1, all starting together. For a specified
observer i, loneliness means d(x_i(t),x_j(t)) >= 1/N for every j != i.

A half-turn relation has i -> j when 0 < (x_j-x_i) mod 1 < 1/2.
It is a tournament only when no pair collides or is antipodal. At the walls
retain the tie type; do not manufacture an orientation and call it equivalent.

Let H:R->R be a strictly increasing continuous lift with H(x+1)=H(x)+1.
Its induced circle map is an orientation-preserving homeomorphism.

## 3. Proposition A: exact threshold-preserving homeomorphisms

For 0 < delta <= 1/2, the following are equivalent:

(a) For every x,y, d(x,y) >= delta iff d(H(x),H(y)) >= delta.

(b) For every x, H(x+delta)=H(x)+delta.

Proof for delta<1/2. Fix x. The forbidden open arc is (x-delta,x+delta)
mod one. An increasing homeomorphism preserving its complement must send
its two boundary points to H(x)-delta and H(x)+delta, in their cyclic order.
The lift interval x < x+delta < x+1 forces the latter equality in R rather
than merely modulo one. Conversely, (b), its backward-shift counterpart,
and monotonicity carry the forbidden arc exactly to the forbidden arc.
For delta=1/2 the allowed set is the singleton antipode; its preservation is
exactly (b).

Consequences:

- If delta=p/q in lowest terms, (b) is equivalent to commutation with 1/q,
  by Bezout and commutation with integer translations.
- Preserving thresholds 1/N and 1/M is equivalent to commuting with
  1/lcm(N,M). This concerns all pairs, not one observed configuration.
- Preserving thresholds with unbounded denominators forces H(x)=x+c:
  H(x)-x is continuous and has positive periods tending to zero.
- The half-turn relation and both its wall types are preserved iff H
  commutes with 1/2. Necessity already follows from antipodal ties.

This is a sharp criterion, not just a sufficient symmetry construction.

## 4. Proposition B: swept-area coordinates and rotational symmetry

Let C be a continuous positive radial graph about O,
q(theta)=r(theta)(cos(theta),sin(theta)), of total enclosed area
A=(1/2) integral_0^{2pi} r(theta)^2 dtheta. Convex curves containing O
strictly in their interior are examples; convexity is not necessary.

The normalized swept-area coordinate is

F(y) = (1/(2A)) integral_0^{2pi y} r(theta)^2 dtheta,

extended by F(y+1)=F(y)+1. It is a strictly increasing homeomorphism.
Let H=F^{-1}. Equal areal sweep x(t)=v t means normalized polar angle
y(t)=H(vt), with the appropriate lift.

For an integer L >= 2, if C has L-fold rotational symmetry, then r(theta+2pi/L)=r(theta), so
F(y+1/L)=F(y)+1/L and H(x+1/L)=H(x)+1/L. Thus every threshold whose
denominator divides L is transported exactly by Proposition A.

Conversely, if the area-to-angle map preserves threshold 1/L for every
pair, Proposition A gives the shift identity for F. Differentiating that
identity (r is continuous) yields r(theta+2pi/L)^2=r(theta)^2, hence the
L-fold symmetry. So for this class the symmetry requirement is also necessary.

This shows precisely what Kepler's equal-area law contributes: a uniform
phase coordinate. It does not on its own force simultaneous separation.

## 5. Proposition C: regular polygons, perimeter speed, and tournament lift

Let P be a convex polygon containing O in its interior, parametrized
counterclockwise by arclength s. On an edge at distance a_e from O,

dA/ds = det(q(s),q'(s))/2 = a_e/2.

Therefore normalized swept-area phase equals normalized perimeter arclength
for every boundary point iff all edge support distances a_e are equal.
Equivalently, P is tangential about a circle centered at O. In particular,
this holds for every regular m-gon with integer m >= 3 about its center. The formula holds on
each open edge and the integrated identity passes through vertices.

For a regular m-gon, let P_m(x) be its unit-perimeter phase parametrization.
Rotational symmetry gives

P_m(x+1/m) = R_(2pi/m) P_m(x).

Its lifted polar angle H_m(x) hence obeys H_m(x+1/m)=H_m(x)+1/m.
For N dividing m, putting p_i(t)=P_m(v_i t) gives the exact equivalence

d(v_i t,v_j t) >= 1/N
iff d(H_m(v_i t),H_m(v_j t)) >= 1/N.

If m is even, det(p_i,p_j) has the same sign as
sin(2pi(v_j-v_i)t), and is zero for exactly the same collision/antipodal
pairs. The radius factors are positive and H_m preserves half-turns.
Thus m=2N works for all N. The mapped runners move at constant perimeter
speed and constant signed areal speed, but not constant polar-angular speed.

The general LRC for N is therefore equivalent to loneliness for this
specific equal-area flow on a regular 2N-gon when distance is measured by
polar angle (or perimeter phase). This is an exact reformulation only.

Within a face whose outward-normal direction is theta_j, put beta=pi/m
and psi=theta-theta_j in [-beta,beta]. For inradius a,
r=a sec(psi), and the normalized area fraction within that face is

u = (tan(psi)+tan(beta))/(2m tan(beta)),  0 <= u <= 1/m.

This gives an explicit inverse chart and identifies the nonlinearity:
psi=arctan((2m u-1)tan(beta)). No extra dynamical theorem has been inserted.

### Honest mechanical realization: a circular billiard

The regular polygon is a periodic orbit of a particle undergoing specular
reflection at the boundary of its circumdisk, not a smooth inverse-square
Kepler orbit. This realization makes equal-area conservation exact even at
the corners. For unit circumradius, let w_j be consecutive polygon vertices,
phi=2pi/m, ell=2sin(phi/2), and

u^-=(w_j-w_(j-1))/ell, u^+=(w_(j+1)-w_j)/ell.

Then u^+ - u^- = -2sin(phi/2) w_j and u^- dot w_j=sin(phi/2).
Thus u^+=u^- -2(u^- dot w_j)w_j, the specular-reflection law. Its impulse
is radial, so det(w_j,u^+)=det(w_j,u^-)=cos(pi/m), the apothem.
Along each straight flight det(q,q') is constant as well. At speed c,
signed areal speed is c cos(pi/m)/2 on every segment and unchanged at every
bounce. Scaling the speed changes energy and angular momentum but not the
same polygonal orbit. This supplies an exact conservation-law analogy,
without claiming Newtonian gravity or a new loneliness existence mechanism.

## 6. Three-predicate LCM corollary (root's refinement)

Choose L=lcm(N,6), use the regular L-gon's H_L, then radially project to
the unit circle z_i=exp(2pi i H_L(x_i)). This preserves simultaneously:

1. the LRC threshold 1/N;
2. the half-turn orientation and collision/antipode ties;
3. Euclidean unit-distance edges **on that unit circle**, since
   |z_i-z_j|=1 iff the phase difference is +/-1/6.

The third assertion concerns the radially projected circle, not chords of
the physical polygon. The unit-distance graph of the whole unit circle is a
disjoint union of six-cycles, hence bipartite. This construction supplies no
high-chromatic plane unit-distance graph and no six-color forcing theorem.

## 7. Sharp hostiles and lost data

### 7.1 An odd polygon can reverse the tournament

For an equilateral triangle with vertices (1,0), (-1/2,sqrt(3)/2),
(-1/2,-sqrt(3)/2), use perimeter phases x=1/12 and y=9/16.
Then y-x=23/48<1/2, so the phase edge points x->y. But

P_3(x)=(5/8,sqrt(3)/8), P_3(y)=(-1/2,-3sqrt(3)/16),
det(P_3(x),P_3(y))=-7sqrt(3)/128<0.

The polar-angle edge is reversed. The 1/3 loneliness threshold still
survives: threshold preservation and half-turn preservation are different
requirements. Replacing the triangle by the regular hexagon repairs both.

### 7.2 Polygon chords are not the phase metric

On a unit-side square, a vertex pair one side apart and the pair of midpoints
on two adjacent sides both have normalized perimeter separation 1/4.
Their squared Euclidean chord lengths are respectively 1 and 1/2.
Angle or perimeter predicates cannot silently be replaced by a fixed
Euclidean chord-length predicate. Retain radii/edge location if using chords.

### 7.3 Kepler ellipses do not provide arbitrary threshold transport

An eccentric ellipse about a focus has radial function
r(theta)=p/(1+e cos(theta)), 0<e<1. It has no nontrivial rotational symmetry
about that focus. By Proposition B its equal-area-to-angle map preserves no
threshold 1/N universally for N>=2. The area metric remains faithful, but
true-anomaly separation is not the original LRC separation. A centered
ellipse has twofold symmetry, not the focused Kepler law; it preserves the
half-turn relation but not 1/N for N>=3. The circle is the common exception.

These statements invoke only the equal-area parametrization, not an assertion
that polygonal paths solve a smooth inverse-square-force equation.

An explicit focus-ellipse orientation reversal works for every 0<e<1.
At eccentric anomalies E_1=pi/2 and E_2=3pi/2, mean anomalies are
M_1=pi/2-e and M_2=3pi/2+e. Their forward difference is pi+2e>pi.
The focus-centered position vectors are (-ae,b) and (-ae,-b), where
b=a sqrt(1-e^2). Their true anomalies are nu_1=acos(-e) and
nu_2=2pi-acos(-e), with forward difference 2acos(e)<pi. The mean-anomaly
half-turn edge is reversed by the true-anomaly picture. This is exact and
does not require a numerical eccentricity or an approximate Kepler solver.

### 7.4 Actual regular runner configurations are too strong

For N>=3 take velocities {0,1,...,N-2,N}. If their N positions ever formed
a regular N-gon, the stationary position 0 would be one vertex and speed1
would satisfy t=j/N. Then speed N also sits at 0, a collision. Thus no
regular N-gon configuration occurs. For N=3, speeds (0,1,3) make the
stationary runner lonely at t=1/2, with separation1/2, despite never forming
a regular triangle. A regular reference **track** is legitimate; requiring
the runners themselves to form a regular polygon is not equivalent to LRC.

### 7.5 Scalar areal conservation loses synchronization

The primitive three-speed tuples (1,2,6) and (1,3,5) have identical total
speed 9 and hence identical aggregate signed areal sweep on any fixed
area-parametrized curve. Their maximum stationary-observer loneliness is
respectively 2/7 and 1/2 (exactly checked below). A scalar sum of area rates
therefore does not recover the target or its slack. Preserve individual
phases and the common time, not just conserved total area.

## 8. Common-origin guardrail, with a current primary hostile

Blanco–Criado–Santos,
[arXiv:2603.24784v3, Proposition 5.11](https://arxiv.org/html/2603.24784v3),
gives a shifted counterexample with speeds (1,2,3,4,5) and starts
(0,46,38,47,72)/94. Its maximum over time is 15/94<1/6.
Our independent rational calculation verifies this **fixed starting vector**
by both the piecewise-linear lower envelope and a closed bad-interval cover.
It does not independently prove the paper's optimization over all starts.
The four maximizing times are 63/188, 39/94, 45/94, and 55/94.
With the same speeds and all starts zero, the maximum is1/6.

Therefore independent phase shifts are not a harmless repair to a geometric
map. Common initial phase is a mathematical hypothesis, not presentation data.

## 9. Exact torus-to-zonotope map, not a planar-polygon shortcut

Take k=N-1 positive integer relative speeds v and delta=1/N. In T^k,
LRC asks whether the orbit t v meets [delta,1-delta]^k modulo Z^k.
Let pi:R^k->R^k/Rv and Lambda=pi(Z^k). For rational/integer v this
is a lattice in the (k-1)-dimensional quotient. Put Z=pi([0,1]^k) and
c=pi((1/2,...,1/2)). Directly from ker pi=Rv,

exists t with t v mod Z^k in [delta,1-delta]^k
iff (c+(1-2delta)(Z-c)) intersect Lambda is nonempty.

Proof: a cube point q has pi(q)=pi(z) for some z in Z^k iff q-z=t v.
Projection of the smaller cube is precisely the displayed central homothet.

This is a **fixed-center lattice hit**. Requiring all starting shifts s
replaces Lambda by every coset pi(s)+Lambda, a covering problem; the stronger
statement is false by the hostile in section8. The quotient keeps the
arithmetic lattice and dimension k-1=N-2. For N=14 that dimension is12,
not2. A picture of a regular polygon is not automatically this zonotope.

Primary background: Henze–Malikiosis,
[arXiv:1609.01939](https://arxiv.org/abs/1609.01939), and Beck–Hoşten–Schymura,
[Lonely Runner Polyhedra](https://arxiv.org/abs/1606.01783).
The explicit cube quotient above is proved here, so no covering-radius
strengthening is being imported tacitly.

## 10. Bridge ledger and exact remaining obligation

| Source -> target | Map | Preserved predicate | Lost information / sidecar | Decisive test |
|---|---|---|---|---|
| Circle phase -> symmetric boundary | x -> P_L(x) | area phase, cyclic order | Euclidean metric; retain phase/radius | square chord hostile |
| Area phase -> angle | H=F^-1 | threshold1/N iff N-fold symmetry | generic gaps/slack; retain H^-1 | Proposition A/B |
| Half-turn tournament -> polygon determinants | i->j by det sign | edges and ties if L even | threshold slack still absent | odd-triangle reversal |
| Polygon angle -> unit circle | radial projection | threshold/tournament; unit edges if6 divides L | physical polygon chords | C6 unit-distance components |
| Runner orbit -> regular point configuration | no general map | only a sufficient special witness | actual reachability destroyed | {0,1,...,N-2,N} |
| Orbit -> lattice zonotope | quotient by Rv | fixed-center lattice hit iff LRC | time lift and lattice needed | shifted counterexample |
| Individual area clocks -> sum | summation | total areal rate | synchronization and metric slack | (1,2,6) vs (1,3,5) |

The next valid proof obligation is still: for each primitive relative-speed
vector, show that the common equal-area trajectory reaches the intersection
of the closed angle-safe sets around the marked observer. The geometric
transport neither proves this intersection is hit nor supplies a monotone
quantity forcing it. A useful future tournament certificate must retain
the observer, phase gaps or equivalent threshold slacks, common time,
endpoint ownership, and collision/antipode ties. An unweighted orientation
alone is insufficient.

The last insufficiency already has a three-vertex exact witness: phases
(0,2/5,3/5) and (0,1/5,3/5) have the same labelled directed three-cycle,
without ties. The marked observer0 is lonely at threshold1/3 in the first
configuration and not in the second. Thus even the labelled geometric
tournament, not just its isomorphism class, loses the needed metric predicate.

The positive result is a faithful change of representation with a sharp
symmetry criterion. The stopped shortcuts are actual regular orbit forcing,
arbitrary focused Kepler reparametrization, Euclidean polygon-chord
substitution, and arbitrary-phase covering. None is needed for the valid
polygon/area-clock equivalence.

## 11. Reproduction and scope

Run `python3 verify_geometry.py > verification.json` and
`python3 -O verify_geometry.py > verification.optimized.json` in this directory.
Both modes retain explicit checks and must give identical receipts. Rational
starting phases are normalized modulo one, including negative or large inputs.
The script uses only Python's standard library. Rational lower-envelope
calculations and hostile determinants are exact. The general propositions
are proved above, not inferred from tests. Additional floating checks for
N=3,...,20 and m=N,2N test 3,636 threshold boundaries, 18,180 off-boundary
threshold pairs, and 13,635 even-polygon orientations. They are labeled
diagnostics, not finite-exact proofs or a runner census.

Repository startup read was bounded. `agents/processor.py --status --peek
--limit 8` was unavailable because this checkout has no `.machine-id`;
no identity was invented. No repository file was changed by this worker.
