# The rhombus carrier, a recovered five-colour ceiling, and genuine field escape

2026-10-11. **PROVED, scoped and self-audited:** coordinate-ring colourings,
the field reduction, and the geometric/arithmetic implications below.
**FINITE-EXACT:** explicit finite colour certificates and assemblies.
**VERIFIED by unconditional PARI descent:** the stated elliptic ranks and
torsion groups. **OPEN:** a forcing six-colour context in Q(zeta55).
No independent peer audit or new plane-wide chromatic bound is claimed.

## Inheritance and an important stronger result recovered during the session

The [unit-component packet](cyclotomic_unit_component_20261011.md) supplied
the rhombus, the marked Cartesian no-gain theorem, and a square-class test
for a deliberate cross-edge. The [boundary packet](boundary_palette_compression_20261011.md)
requires a full relation and its lift. The [pentagon context theorem](pentagon_context_number_fields_20261011.md)
already proves that a five-colourable external context cannot force the
required bad lists under its stated private-port interface.

The hostile is an authentic unit graph with no palette gain; the corrected
near miss is equating an additional unit edge with a stronger boundary
relation. The sidecar is the **actual coordinate ring**, including permitted
denominators, rather than the name of a large ambient field.

During the inheritance search, [THM-4558, residue colourings](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md)
pointed to an existing five-colouring of the q=19 residue graph. This is
stronger than the initial ring-specific five-colour screen here. We recovered
and directly verified a small colour table, then applied it to the entire
degree-eight carrier. The finite-field idea is credited to
[decalion89, local colourings, section 13](https://github.com/decalion89/chromatic-number-of-the-plane/blob/cb660f322dce47a773f347d39818ffc46e6ca72c/notes/local_colourings.md).
That repository was read at the pinned commit; its unrelated claims are not
imported. Our short initial SAT searches on 361 vertices returned UNKNOWN;
the explicit table below resolves the needed upper bound independently.

## G1. All dyadic assemblies of the new rhombus and pentagon are three-colourable

Set

    alpha=zeta5, beta=sqrt(-11), w=sqrt(-55)=sqrt(5)*beta,
    eta=(-3+w)/8,
    R=(0,-eta,3/4,-bar(eta)),
    K=Q(alpha,beta) < L=Q(zeta55).

The rhombus R is the inherited induced C4. Its diagonals have lengths
3/4 and sqrt(55)/4. Both it and the unit pentagon lie in the ring

    A0=Z[alpha,w,1/2].

The map

    r0:A0 -> F5,       alpha -> 1, w -> 0, 1/2 -> 3      (1)

is a ring homomorphism and commutes with conjugation, which acts trivially
on its target. One can check this directly in the free basis
alpha^i w^j, 0<=i<4, 0<=j<2: the defining relations Phi5(alpha)=0 and
w^2=-55 both vanish after substitution. The quadratic extension has degree
two because Q(alpha)'s only quadratic subfield is real Q(sqrt5).

For any unit difference u in A0,

    r0(u)^2=r0(u*bar(u))=1,

so the unit graph maps to the cycle on F5 with steps +/-1. The residue
colour table

    residue: 0 1 2 3 4
    colour:  0 1 0 1 2                                   (2)

is proper. Since A0 contains a unit C5, its complete unit graph has
chromatic number **exactly 3**.

This covers arbitrary translations, sums, products, conjugations and dyadic
divisions that stay in A0, including deliberate point identifications and
every additional unit edge in the completed geometry. It is stronger than
checking a fixed number of rhombus layers. Inverting 5, or taking roots not
in this ring, is not among these closure operations.

## G2. A useful half-angle, then the stronger full-field ceiling

The unit

    rho=(sqrt(5)+beta)/4,       rho^2=eta                 (3)

lies in A1=Z[alpha,beta,1/2], but does not lie in A0. Otherwise (1) would
give r0(rho)^2=1 from its norm and r0(rho)^2=r0(eta)=-1 simultaneously.
This is a real change
of representation: beta=w/sqrt5 introduces a denominator 5 in the old basis.
On A1 the two maps r_+,r_- to F5 send alpha to 1 and beta to +/-2.
Conjugation exchanges them. A unit u satisfies

    r_+(u)r_-(u)=1,

so r_+ itself gives a five-colouring. This also holds after inverting any
rational integers coprime to 5. It is a sufficient ring bound, not an exact
chromatic count.

A natural first attempted escape was

    v=(2+beta)/(2-beta)=(-7+4beta)/15.

It has norm 1 and a pole at one prime over 5. Its induced unit path is

    0, 1, 1+v, 1/15.                                   (4)

The nonedge squared distances are 16/15, 196/225 and 1/225, all different
from 1. Thus it genuinely escapes the preceding integral-at-5 test.

**But it does not escape the stronger result:**

    3 <= chi(G(K)) <= 5,       K=Q(zeta5,sqrt(-11)).      (5)

At 19, Frobenius sends alpha to alpha^-1 and beta to -beta, because
19=-1 mod5 and -11 is a nonsquare modulo19. It is precisely complex
conjugation, of order two. Thus each prime above19 has residue field F_(19^2)
and is conjugation-fixed. Every unit direction is a local unit there;
the inherited residue-colouring lemma applies to all components of G(K),
without a restriction on the coordinates' denominators.

Write that residue field as F19[i], i^2=-1. For x^2+y^2=1, the possible
values of x+y are exactly

    D={1,2,4,6,9,10,13,15,17,18}.

The following table c on 0,...,18 colours every pair with difference in D:

    (4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0).           (6)

Hence (x,y)->c(x+y mod19) is a five-colouring of the entire residue graph.
All twenty circle directions and all nineteen translates are checked
exactly. The lower bound in (5) is the inherited pentagon. We do not claim
K actually needs five colours.

Therefore **no change of generators, roots or denominators staying in K
can produce the required six-chromatic context**. The pole unit (4) remains
a useful control, but ceases to be a viable six-colour route once (5) is
recovered. This is why prior results must be checked before elaborating a
weaker local escape.

There is a further **odd-degree construction barrier**. The extension L/K
has degree5. Suppose two distinct centres a,b in K have a common unit-circle
point z in L. Set t=(z-a)/(b-a). The two norm equations give

    t+bar(t)=1,      t*bar(t)=1/N(b-a),
    t^2-t+1/N(b-a)=0.

Its degree over K is at most2 and divides5, so t, and hence z, is already
in K. Thus arbitrarily many rational operations, conjugations and valid
two-circle intersections **inside L** cannot escape K from these seeds.
This applies the elementary degree law, not a new constructibility theorem.

An explicit replacement generator is gamma=zeta11. With beta defined by
the quadratic Gauss sum, it satisfies the irreducible quintic over K

    X^5+(1-beta)X^4/2-X^3+X^2-(1+beta)X/2-1=0.          (7a)

Its five roots are gamma^r for r=1,3,4,5,9. The program verifies the
identity in the degree40 power basis. This generator supplies the missing
degree-five extension; calling a quadratic construction in a larger ambient
field would not. It also leaves the old q=19 residue carrier: a primitive
eleventh root cannot lie in F_(19^2), whose multiplicative group has order360.

## G3. Overlaps were tested, including a genuinely larger carrier

First keep the original pentagon P fixed and successively take Minkowski
sums with eta R, R, eta^2 R, identifying equal coordinates. Every pair is
tested for unit distance; no abstract edge list is substituted for geometry.

| Assembly stage | Vertices | Unit edges | Chromatic number | Interior vertices left by relative degree-5 peeling |
|---|---:|---:|---:|---:|
| P+eta R | 18 | 35 | 3 | 0 |
| previous + R | 72 | 212 | 3 | 0 |
| previous + eta^2 R | 216 | 780 | 3 | 211 |

The first product has actual collisions and is not the old injective
Cartesian product. Nevertheless, (2) colours the entire family. The last
row is a useful hostile: a nonempty relative degree-six core does not imply
six-colour forcing. Degree peeling is only a sufficient extension test.

To leave K, adjoin a genuine eleventh root gamma=zeta11. Let P11 be the
eleven partial sums of 1,gamma,...,gamma^10. P+P11 already generates L as
a field: its difference set contains both alpha and gamma. Put B=P+P11.
The exact controls then give:

| Full-field assembly | Vertices | All unit edges | Exact chi |
|---|---:|---:|---:|
| B | 54 | 108 | 3 |
| B+R | 216 | 648 | 3 |
| B+eta R | 188 | 550 | 3 |
| B+eta gamma R | 206 | 613 | 3 |
| B+R+eta gamma R | 824 | 3276 | 3 |
| Five K-Galois conjugates of B+eta gamma R, united | 994 | 3063 | 3 |
| Five K-Galois conjugates of B+R+eta gamma R, united | 3976 | 16228 | 3 |

Every row retains the pentagon's original coordinates. Rows containing +R
also retain all four original rhombus ports, without rescaling them. The
last graphs lie beyond the field ceiling (5), have genuine overlaps, and
still have verified three-colourings. The saved output contains
every colouring and a digest of the deterministic coordinates and completed
edge set. An exact evaluation modulo331 rejects impossible unit pairs;
survivors receive full power-basis norm tests. The small control also tests
all rejected pairs directly. There is no floating-distance tolerance.

The Galois operation is especially relevant to fixed context. Every
automorphism over K fixes the original pentagon and rhombus pointwise and
commutes with conjugation, so preserves all unit edges and nonedges. Each
copy has the same labelled boundary relation. Gluing only at those ports,
with independent interiors and no new cross-edges, would intersect identical
relations and give no gain. The actual unions test the stronger operation:
coordinate coincidences and **all cross-copy unit edges** are retained. They
remain three-coloured in the declared cases. This does not assert that every
Galois union is harmless.

These are finite failed forcing candidates, not evidence that L has
chromatic number three. They locate the next obligation: a full-field
overlap must actually delete allowed boundary colour partitions. Merely
using more of the field, adding a small-scale connector, or creating a
dense core is insufficient.

## G4. The golden cross-edge gate has an elliptic quotient

For a golden chord of length phi=(1+sqrt5)/2 and a rational spacing d>0,
the question |phi+d u|=1 with u*bar(u)=1 in L requires

    Delta=4 phi^2 d^2-(1-phi^2-d^2)^2
          =A+B phi,
    D=d^2, A=4D-1-D^2, B=2D-1,
    Delta/55 is a square in F=L intersect R.             (7)

The inherited square-class criterion is sufficient as well, before any
vertex-collision or full-incidence requirements are imposed. Taking the
norm to Q gives the useful exact factorization

    N(Delta)=(D^2-7D+1)(D^2-3D+1).                     (8)

Thus the golden trace recursion produces both 3=phi^2+phi^-2 and
7=phi^4+phi^-4 here. This is a mechanism, not just a repeated numeral.

Gal(F/Q) is cyclic of order20. Its quadratic subfield is Q(sqrt5), and
its unique degree-four subfield is obtained by adjoining
sqrt(11(5+sqrt5)/2). Consequently a Q(sqrt5) element square in F is either
a square already or 11(5+sqrt5)/2 times one. Its rational norm must be a
square or 5 times a square. Applied to (7), this is a necessary condition
on (8).

Set t=d-1/d. After dividing (8) by d^4, that necessary condition becomes

    (t^2-1)(t^2-5)=epsilon v^2,       epsilon in {1,5}.   (9)

For t!=0 the map X=epsilon t^2, Y=epsilon^2 t v gives

    E_epsilon: Y^2=X(X-epsilon)(X-5epsilon).             (10)

**VERIFIED global arithmetic input:** PARI/GP's unconditional 2-descent
returns rank bounds [0,0] for E_1 and E_5, and in each case torsion
Z/2 x Z/2, with exactly the three displayed roots as affine torsion points.
The [official algorithm contract](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Elliptic_curves.html#ellrank)
distinguishes these bounds from an analytic-rank/BSD inference; no analytic
rank is used. The reproduction retains the complete outputs.

It follows that t!=0 would have t^2=1 or5. A rational t cannot have square5;
t=+/-1 gives no rational solution of d-1/d=t. The t=0 case gives d=1.
Thus, **using that verified descent input**, the only positive rational
spacing admitting (7) is d=1. It is realizable by fifth-root directions.
The earlier denominator<=120 scan found only d=1, but does not supply this
all-range conclusion.

This gate concerns a rational spacing against a golden chord. It does not
exclude arbitrary algebraic spacings or all rhombus arrangements. For a
spacing sqrt55*h, h rational, the corresponding twists include

    E_11: Y^2=X(X-11)(X-55),

which has verified rank one and point (180,1950). This is a positive new
arithmetic carrier. It is not yet a rhombus: that point gives r=6/11 and
the remaining lift needs r^2+4/55=224/605 to be a rational square, which it
is not. A productive helper is to enforce the extra lift and field-square
conditions on this rank-one family, then authenticate collisions and the
whole colour relation. A rational point on a quotient curve alone loses
exactly those conditions.

## G5. What is retained from the pasted tiling idea

The useful proposal is a three-label system with a binary side bit and
two geometric scales. The repo already has the exact version in
[the pair-bit and seam packet](plane_cycles_pants_and_six_seven_20261011b.md):
same-label edges must form bipartite graphs, with relative bit parity kept
when fragments are glued. A rhombus alone supplies no odd-cycle obstruction.

Incoming commit **9d30772597** arrived during this session and was read before
filing. Its [square-cell conflict theorem](square_cell_conflicts_20261011g.md)
is stronger than the preliminary rectangle checker considered here: exact
half-open ownership, and an optimum of **eight** for the specified forty-cell
motif with colour periods8a in both directions, over all scales. A bounded
diagonal seam cannot repair that periodic bulk. Larger colour periods and
position-dependent refinements are outside its scope.

Its [diagonal theorem](diagonal_palette_holonomy_20261011g.md) separately
proves a sharp nine-colour requirement under fixed commuting reflection and
diagonal-translation palette actions with every phase orbit of length3.
Changing the seam rule or action leaves those hypotheses. Its signed
component equations retain all crossing edges, including parallel constraints.
These results replace a blanket objection to periodicity.

The concurrent [multiport hostile](spindle_tournament_interfaces_20261011g.md)
has a three-tip relation excluding all-equal colours although every pair
projection allows both equal and unequal colours. Its triangle geometry is
not imported into our triangle-free field L. The transferable lesson is that
Galois gluing must retain the **joint** port relation; pairwise summaries
cannot certify forcing.

Nor does the square motif's quotient eight-clique give an eight-vertex unit
graph: different pairs can use different translated cells and different
interior points. Their witnesses do not share one simultaneous point choice.
This source distinction prevents an apparent contradiction with the exact
three-coloured point assemblies above.

The pasted text is inspiration, not a proved obstruction: periodicity alone
does not force a unit conflict, nonmeasurability supplies no automatic
protection, and reflection across a line moves a point by twice its distance
to the line. A point at distance1/2, not1/sqrt2, has a unit-separated mirror.
Swapping two distinct colours within a pair is not monochromaticity.

## Reproduction and next target

[Program](../../04-computation/experiments/changing_controllers_rhombus_20261011.py),
[output](changing_controllers_rhombus_20261011.json). Dependencies are Python,
SymPy, PySAT and PARI/GP. Ordinary and optimized runs are compared. The proofs
are self-audited; the elliptic rank outputs are independently executable
algorithmic inputs, not a claimed handwritten descent.

The remaining target is an authenticated full-field context with six-colour
boundary support strong enough for the marked pentagon. The search must
leave K, retain the original ports, and measure the full extension relation.
The current 3976-point example satisfies the field and geometry requirements
but fails the palette requirement. Its explicit colouring is a stopping
reason for that candidate, not for the full-field route.
