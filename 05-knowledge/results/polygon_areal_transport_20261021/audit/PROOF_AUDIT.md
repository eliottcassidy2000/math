# Session 19 geometry: independent audit

Date: 2026-10-21 UTC. Verdict: **PASS, within the scope below.**

The frozen producer report and checker are identified in `producer_freeze.json`.
No producer or repository file was edited by this auditor. The root made the
explicit hypothesis and checker-hardening repairs listed below before freezing.

## 1. Permitted conclusion and excluded conclusions

- **PROOF-REVIEWED:** universal threshold transport for orientation-preserving
  circle homeomorphisms; the positive-continuous-radial-graph symmetry criterion;
  normalized area/perimeter equivalence for centered tangential convex polygons;
  regular-polygon transport, the radial LCM bridge, and circular-billiard realization.
- **FINITE-EXACT, independently replayed:** the specified shifted start vector
  has maximum loneliness 15/94, attained at precisely four stated times; the
  unshifted AP and other literal hostile controls agree with the report.
- **SOURCE-VERIFIED, NOT MATHEMATICALLY AUDITED:** the pinned 2026 preprints make
  the cited claims. Their large computations and full proofs were not audited.
- **NOT ESTABLISHED:** a general LRC proof, an orbit-hitting mechanism, independent
  certification of the fourteen/fifteen-runner preprint, optimization over all
  shifted starting positions, a high-chromatic plane unit-distance construction,
  or a theorem converting polygon chords into the phase metric.

The positive result is an exact representation change, not a new existence
proof. In particular, retaining a dated repository headline as the current
external literature status would be misleading, but the external preprint is
not thereby promoted to independently audited repository canon.

## 2. Proposition A: homeomorphism threshold criterion

For 0 < delta < 1/2, fixing x makes the allowed set the closed oriented arc
[x+delta,x+1-delta]. Its image under an increasing degree-one lift is the
corresponding arc based at H(x). Equality of these sets determines both ordered
endpoints, so H(x+delta)=H(x)+delta. Conversely, the shift identity and its
backward version carry the open forbidden arc exactly onto the open forbidden
arc. Taking complements retains all equality endpoints. This proves both
implications, not merely sufficiency.

For delta=1/2 the allowed set is the single antipode. Preservation gives exactly
H(x+1/2)=H(x)+1/2. The report explicitly handles this case; the ordinary two-ended
arc argument is not used on a collapsed allowed arc. Delta=0 or delta>1/2 would
make the predicate trivial and are excluded.

For reduced delta=p/q, Bezout combines shifts p/q and 1 into shift 1/q; negative
integer combinations follow by reversing the shift identity. Conversely, q
successive small shifts imply all p/q shifts. Thus simultaneous 1/N and 1/M
preservation is equivalent to 1/lcm(N,M) commutation. Unbounded denominators
give periods tending to zero for the continuous function H(x)-x; approximation
of any real displacement by multiples of these periods makes this function
constant. No differentiability is needed here.

The half-turn relation is preserved because each open forward semicircle maps
onto the corresponding one. Homeomorphism injectivity preserves collisions, and
the shift identity preserves antipodes. Conversely, preserving antipodal ties
already forces that identity. Discarding ties would change the claimed object.

## 3. Proposition B: areal clock and exact shape hypotheses

A positive continuous radial graph gives a simple star-shaped boundary about O,
a finite positive area A, and a C1 strictly increasing swept-area lift F with
F'(y)=pi*r(2*pi*y)^2/A. Convexity is sufficient but not necessary. Positivity and
a complete once-around radial graph are important stated assumptions; arbitrary
curves, self-intersections, or a center outside the enclosed region are not
silently included.

L-fold symmetry gives a periodic derivative and equal area in each 1/L sector,
so F and its inverse H commute with 1/L. Conversely, Proposition A gives this
shift for H, hence for F. Differentiation gives the same shift for r squared;
positivity then gives the shift for r itself. Therefore the shape criterion is
an iff for integer L >= 2. It is not a claim that area conservation guarantees
a simultaneous lonely time.

## 4. Proposition C: perimeter, determinants, and mechanics

For a counterclockwise convex polygon with O strictly inside, the determinant
of the boundary point with unit edge tangent is the positive edge support
distance a_e. Thus dA/ds=a_e/2 on each open edge. Equality of the normalized area
and arclength coordinates at every point forces this derivative to be constant
on every edge, and conversely equal a_e gives the integrated equality. Equal
support distances say precisely that every edge is tangent to the same circle
centered at O. Continuity passes the integrated identity through vertices;
no pointwise tangent is falsely asserted at a corner.

A regular nondegenerate m-gon (integer m >= 3) has the required support-distance
and rotational identities. If N divides m, Proposition A gives exact closed
threshold transport. If m is even, the additional half-turn identity preserves
determinant signs and both zero types. For N >= 2, m=2N is nondegenerate and
handles the delta=1/2 case. Original LRC velocities are now explicitly distinct.

The face formula integrates a squared secant and has the correct factor:
u=(tan(psi)+tan(beta))/(2*m*tan(beta)). It takes values 0 and 1/m at the two
vertices, and inversion gives the stated arctangent expression.

The mechanical realization is a billiard in the circumdisk, not a gravitational
polygon orbit. The incoming velocity's radial component is sin(pi/m), so its
reflection changes velocity by -2*sin(pi/m)*w_j, exactly the difference of the
adjacent unit chord velocities. The impulse is radial and hence has zero
angular torque about O. The angular momentum determinant equals cos(pi/m)
on every unit-speed segment and across every bounce. At speed c the area rate
is c*cos(pi/m)/2. Changing c changes the energy and angular momentum, not the
polygonal geometric trajectory. The proof is valid at all m >= 3.

## 5. LCM three-predicate bridge

For L=lcm(N,6), N, 2, and 6 all divide L. The resulting angle map preserves the
LRC threshold and the complete half-turn relation with ties. Radial projection
onto the **unit circle** is essential for the third predicate: a unit chord
there is equivalent to circular phase difference +/-1/6. This predicate need
not describe a chord of the original physical polygon.

The full unit-circle unit-distance graph has edges given by addition of +/-1/6.
Every component is exactly a six-cycle. Hence the construction is bipartite and
does not force six colors or settle the plane chromatic-number problem.

## 6. Literal hostiles and endpoint audit

The independent checker does not import or call the producer's mathematics.
For each claimed maximum R it intersects the complete closed individual safe
interval unions at level R, retaining singleton intersections. Separately, it
produces a closed bad-interval cover at radius R. The latter proves maximum <= R;
the safe intersection gives precisely all equality times. This is a different
verification path from enumerating lower-envelope crossings.

- Shifted v=(1,2,3,4,5), s=(0,46,38,47,72)/94: maximum 15/94 < 1/6,
  equality times exactly 63/188, 39/94, 45/94, 55/94. Safe set at 1/6 is empty.
- All-zero starts with the same speeds: maximum 1/6, equality only at 1/6 and
  5/6. The AP boundary set is also exactly checked for N=2,...,30 against the
  primitive fractions a/N. These are finite controls, not their general proof.
- (1,2,6) and (1,3,5) have equal speed sum 9 but maxima 2/7 and 1/2. Thus the
  common scalar areal rate loses individual synchronization and target slack.
- The odd-triangle witness is recomputed by exact perimeter interpolation in
  Q(sqrt(3)): phase gap 23/48, determinant -7*sqrt(3)/128.
- Square chords with perimeter separation 1/4 have squared lengths 1 and 1/2.
- Labelled three-cycle word (+,-,+) occurs in both configurations; observer
  slack is 1/15 for (0,2/5,3/5) and -2/15 for (0,1/5,3/5).
- 2,304 ordered rational phase pairs on a regular hexagon are checked with exact
  Q(sqrt(3)) coordinates for determinant orientation, radial unit-distance,
  the 1/6 and 1/3 closed thresholds, and the singleton antipodal 1/2 threshold. These
  are exact finite controls, not substitutes for the general proof.

For the focused ellipse hostile, 0<e<1 makes the mean-anomaly gap pi+2e lie
strictly between pi and 2*pi, while the true-anomaly gap 2*acos(e) lies strictly
between 0 and pi. The focus vectors have positive determinant 2*a*e*b. The
literal choice a=5, e=3/5, b=4 gives determinant 24. A focused eccentric ellipse
has a unique radial minimum direction, excluding any nontrivial rotational
symmetry. A centered eccentric ellipse has exactly twofold rotational symmetry.
The circle is expressly excluded from these failures.

The non-regular-configuration witness uses N>=3, so speed 1 is present. If the
N positions were regular and included the stationary vertex 0, speed 1 would
force t=j/N; speed N would then collide with 0. The N=3 tuple (0,1,3) nevertheless
has an observer-lonely time t=1/2. A regular reference track and an attained
regular configuration are therefore not equivalent requirements.

## 7. Quotient and common-origin guardrail

For integer nonzero relative vector v and k=N-1, Rv is rational and pi(Z^k)
is a rank k-1 lattice in the quotient. Projection of the closed smaller cube
is exactly c+(1-2*delta)(Z-c). A projected cube point is a lattice point exactly
when a lift differs from an integer vector by t*v. This proves the displayed
fixed-center equivalence, including its endpoints. At N=2 the quotient has
dimension zero and the formula remains valid.

With initial vector s, the lattice condition uses the coset pi(s)+Lambda.
Requiring all starts is therefore all-coset coverage, which the explicit
shifted hostile refutes. Irrational directions would need extra care about
whether the projected group is a lattice; the quotient section is explicitly
restricted to integer speeds. The quotient dimension for N=14 is 12, not 2.

## 8. Freeze, repairs, and reproduction

Before this audit's final freeze, the root explicitly added N>=2, distinct
original velocities, integer L>=2, and m>=3; normalized arbitrary rational
starts modulo one; added input guards; changed assertions into explicit
require calls; and explicitly tested all four shifted maximizing times,
integer-shift invariance, and the two scalar-area maxima. The report's stated
conclusions did not change.

The final producer checker has zero Assert AST nodes. Its normal and -O
outputs are byte-identical to both saved producer receipts. The independent
checker also uses explicit failures, tests its failure path, and produces
byte-identical normal and -O receipts. The producer's floating diagnostics
remain labelled as diagnostics; this audit does not promote them to proof.

Reproduce from the audit directory:

    python3 verify_geometry_independent.py > independent-normal.json
    python3 -O verify_geometry_independent.py > independent-optimized.json
    cmp independent-normal.json independent-optimized.json

Reproduce the producer with its `verify_geometry.py`, with and without `-O`,
and compare against its saved `verification.json`. Exact frozen hashes are
listed in `producer_freeze.json`; publication files are enumerated by the
manifest. Exclude interpreter caches and any temporary files.

## 9. Literature and inherited-source check

`SOURCE_PINS.json` records the versioned primary URLs and repository paths.
Allikvere's arXiv:2609.02604v2 (2026-09-24), Theorem 1.1, states the fourteen-
and fifteen-total-runner result; section 6 reports the absence of independent
reimplementation/formal verification. Its convention counts nonzero speeds,
so LRC(14) there means fifteen total runners. The two archive links match the
paper. These are checked source statements, not an audit of its certificates.

Sungkawichai–Trakulthongchai arXiv:2604.23906v2, Theorem 1.3, gives k<=12 under
the same nonzero-speed convention. Blanco–Criado–Santos arXiv:2603.24784v3,
Proposition 5.11, gives the cited shifted vector and value. This audit verifies
the fixed vector, not the optimization over all starting vectors in that paper.

At repository baseline de0463bc8ddfe46adc61120fa98b9747055fe34f, the exact cited
HYP-2009 and HYP-1995 files are OPEN. The particular THM-1820 file has the
correction banner described by the producer. The October 5 AP correction
retains primitive fractions, not all polygon vertices. Those inheritance
statements were checked directly without adopting uncorrected historical prose.
