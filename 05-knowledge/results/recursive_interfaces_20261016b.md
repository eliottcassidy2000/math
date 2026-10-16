# Marked recursive interfaces: golden hinges, Sturm certificates and Lawvere

2026-10-16, clock_duality_oct15. **PROVED, scoped elementary deductions:**
the rational-radius classification in Q(zeta55), infinitely many mixed-family
contacts (using classical elliptic-curve facts), the finite rotation-event
reduction, the marked boundary composition, and the exact-clock obstruction.
**FINITE-EXACT:** the accompanying rational and cyclotomic controls.
**CITED:** the attached papers and classical Lawvere/Hermite mechanisms;
the papers' complete geometric and formal proofs are not audited here.
**OPEN:** six-colour forcing, uniform paired-reset cancellation, universal
Collatz entry and planar Jacobian. No new theorem ID or priority claim.

[Program](../../04-computation/experiments/recursive_interfaces_20261016b.py),
[exact receipt](recursive_interfaces_20261016b.json).

## 1. Inheritance, source intake and live board

The anchor is a useful exact bridge between recursive geometry and arithmetic
interfaces. The niche is real-root certification; the wildcard is whether
the golden and Paley quadratic fields classify usable rotation events.

The closest proved mechanisms are:

- [THM-4551, Fibonacci monotile boundary words are leaper flights](../../01-canon/theorems/THM-4551-fibonacci-monodromy-monotile-boundary-words-are-leaper-flights.md).
  It already identifies the **ordered** boundary words with centred lattice
  cutting sequences. MISTAKE-564 repairs the earlier intercept argument.
- [THM-4592, eleven-square torus and period-five Collatz cycle set](../../01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md).
  Its period-five labels and equivariance are proved; its general golden
  reader bijectivity is labelled FINITE-EXACT, not an all-period theorem.
- [THM-1960, regular-seed spectral substitution](../../01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md)
  and [THM-1975, the path-cover compositional invariant](../../01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md): cyclic substitution needs a path-cover
  profile, not just the scalar Hamiltonian-path count.
- [Finite algebraic witnesses](finite_algebraic_witness_20261007.md) already
  implements selected-root Sturm checks; [carry-corrected Hermite forms](trinomial_trace_sign_empty_core_certificate_sep06.md)
  already owns the trace-form operation. Neither is introduced as new here.
- [Hinge rotations and forced radii](hn_type_m_hinge_rotations_and_forced_radii_20261014.md)
  already gives the equal-radius rotation-field law. Section 2 specializes
  it completely to Q(zeta55), then derives a mixed-family genus-one condition.

The canonical hostile is a quotient that retains all counts but loses the
selected boundary assignment. The corrected near miss is MISTAKE-529:
neither the numeral 2 in an integer tournament identity nor a self-similar
diagram supplies Lawvere's weakly point-surjective evaluation map. The
least-used coordinate here is the selected real embedding at an isolated
geometric contact event.

The six-concept board is **marked geometry / ordered recursion / real
embedding / field of rotations / boundary relation / integer clock**.
Every connection below states what survives and what is lost.

### Attached sources

1. PINGYOU LTD, *The Fibonacci Geometry of Eleven Squares*, local
   `/Users/e/Downloads/eleven-squares.pdf`: sections 1--2, 6, 8 and Appendix C.
   Its feasibility/minimum statement concerns its specified two-parameter
   contact family. Its torus needs a mark to recover the packing; the
   unmarked torus is insufficient. This distinction was already recovered
   in [the 2026-10-07 routing note](eleven_squares_routing_bridges_20261007.md).
   Figure 1 was inspected; the octic and its Sturm table are checked below.
2. PINGYOU LTD, *A Chiral Aperiodic Polygon with Fibonacci Monodromy*, local
   `/Users/e/Downloads/fibonacci-monotile.pdf`: Theorems B--D and the signed
   boundary substitution on page 2, inspected visually because overbars
   disappear in text extraction. The ordered free-group map and its
   abelianization -Q^3 are different levels of information. We inherit
   THM-4551 rather than re-proving a weaker word-count statement.
3. Ioannis Tsiokos, *A Strongly Aperiodic Monotile in Three Dimensions*,
   local `/Users/e/Downloads/2609.19214v1.pdf`, sections 1.2--1.3 and the
   displayed geometric encoding. The author explicitly calls this a proof
   submission. Its useful mechanism is **every admitted tiling decodes,
   and decoded parents obey the same rule**. This session has not built its
   Lean repository or audited the named compiler hooks.
   [Goodman-Strauss's later account](https://arxiv.org/html/2609.24779v1),
   sections 1--2, independently presents the forced-parent mechanism and
   distinguishes it from a rep-tile that merely permits hierarchical tilings.

Incoming commit `dc233fcf45` preregisters the Collatz dense census from
2^40 to 2^42. Its Mellin fits and aliasing observations are finite predictive
models, not a uniform energy result. In particular its D(s) is an exponential
polynomial, so the rational-polynomial Sturm checker below does not directly
certify its complex zeros. No pending preregistered test is treated as passed.

## 2. Two complete rational-radius families in Q(zeta55)

Put K=Q(zeta55), with its actual complex embedding and ordinary conjugation.
A **unit chord at squared radius s** consists of a,b in K satisfying

    a bar(a) = b bar(b) = s,       (a-b) bar(a-b) = 1.

The centre is 0; a different common centre can be translated to 0, provided
the centre is retained as part of the marked data.

**R1 (PROVED).** For rational s>1/4, such a chord exists if and only if

    s=(1+11q^2)/4  or  s=(1+55q^2)/4,    q in Q, q != 0.       (1)

These two nondegenerate families are disjoint. At s=1/4 the antipodal chord
is possible; below 1/4 it is impossible even over C.

**Proof.** Divide the chord by the unit u=a-b. Then alpha=a/u has norm s
and alpha-1=b/u has the same norm. Subtracting the two norm equations gives
alpha+bar(alpha)=1. Consequently

    alpha=(1+sqrt(1-4s))/2.                                  (2)

Conversely (2) and alpha-1 give a chord whenever that square root lies in K.
The quadratic subfields of K are exactly Q(sqrt5), Q(sqrt(-11)) and
Q(sqrt(-55)): its Galois group is (Z/5)^* x (Z/11)^*=C4 x C10,
with exactly three nontrivial quadratic characters. The three displayed
subfields account for them (the quadratic Gauss sums supply the generators).
Since 1-4s<0 is rational, its square root lies in K exactly in the square
classes -11 and -55. This proves (1). If both forms agreed, nonzero rational
q,r would satisfy q^2=5r^2, impossible. The diameter gives the lower bound.

This is a classification of **rational equal-radius chords**, not of all
rotations or all finite unit graphs in K. Golden irrational radii and
unequal-radius contacts are still available.

For d=11 or 55, an explicit family is

    alpha_q=(1+q sqrt(-d))/2,
    s_q=(1+d q^2)/4,
    gamma_q=(d q^2-1+2q sqrt(-d))/(1+d q^2).

It obeys

    |gamma_q|=1,       gamma_q alpha_q=alpha_q-1.             (3)

Thus the user's difference-family idea has a literal geometric realization:
the common difference is exactly 1; norm and the imaginary quadratic field
record the circle and its rotation. For odd integral q, alpha_q is an
algebraic integer, with polynomial X^2-X+s_q. Rational q is allowed for
the geometric classification without claiming integrality.

At q=1 the squared radii are 3 and 14. The first rotation is
(5+sqrt(-11))/6, the familiar Moser hinge angle. It lies in K although a
Moser spindle cannot: K has no sqrt(-3), hence no equilateral unit triangle.
The availability of an angle does not supply the needed seed graph.

**R2 (PROVED stopping boundary).** For either family and q nonzero rational,
the orbit {gamma_q^j alpha_q : j in Z}, with its induced unit-distance graph,
is a doubly infinite path. In particular it is two-colourable.

Indeed gamma_q lies in Q(sqrt(-d)), whose only roots of unity are +/-1;
the formula excludes both. All orbit points are distinct. A fixed circle
has precisely the two points at unit distance from one point, obtained by
gamma_q and its inverse. Thus two orbit points are adjacent exactly when
their exponents differ by 1. The 12-point finite controls agree in Q[z]/Phi55.
An indefinitely iterable legal generator therefore need not produce a
forcing context. A second circle or centre must contribute cross-edges.

### A mixed-family cross-edge and a genus-one helper problem

The q=1 circles of squared radii 3 and 14 cannot have a common-centre unit
cross-edge: sqrt14-sqrt3>1. Changing the second rational parameter repairs
this obstruction. Let j=sqrt(-11), k=sqrt(-55), and put

    a=(1+j)/2,       b=(3+j)/5,
    gamma=(5+j)/6,   delta=(3+k)/8.

Then, exactly,

    |a|^2=3,  |b|^2=4/5,  |a-b|=1,
    |a-gamma*a|=1,   |b-delta*b|=1.                          (4)

The radius 4/5 is the d=55 family at q=1/5. Thus the two rotation fields
really meet along a unit cross-edge. Equation (4) is a small usable gluing
gadget, not a four- or six-chromatic construction.

Fix the first squared radius 3 and let the second be s=(1+55r^2)/4.
The quadratic contact discriminant is Delta=s^2-8s+4. The branch
Delta=-11 v^2 becomes

    y^2=-275r^4+150r^2-3,       y=4v.                       (5)

The quartic is squarefree. Its smooth projective double cover has genus one,
and (r,y)=(1/5,8/5) is a rational point. This is a concrete route from
unit-distance gluing to an elliptic-curve problem, rather than an analogy
between unrelated genus labels. The other branch is the -55 twist and is
not discarded. The exact census over reduced 0<r<1 with denominator at
most 100 finds only r=1/5 across these two square-class branches. It is not
a classification of the rational points.

**R2b (PROVED): infinitely many mixed-family contacts.** The smooth
projective curve C of (5) has infinitely many rational points, giving
infinitely many distinct second radii with a unit contact to radius squared
3 in K. This does not assert their exact Mordell--Weil rank or chromatic force.

Here is a short independently checkable arithmetic proof. The map

    f:C -> E, (r,y) -> (X,Y)=(-275r^2,-275ry),
    E: Y^2=X^3+150X^2+825X                              (5a)

extends from the affine equation to the smooth projective curves. It is
nonconstant, and sends Q0=(1/5,8/5) to P=(-11,-88), and Q1=(1/5,-8/5)
to -P. Good reductions at 7 and 17 have exactly 6 and 20 points,
respectively. Prime-to-p torsion injects under good reduction: the 7-part
is excluded by the 17-count and the 17-part by the 7-count; all remaining
torsion has order dividing gcd(6,20)=2. Since P has nonzero Y, it is not
2-torsion and has infinite order. The counts and nonsingularity are exact
finite checks, not numerical estimates of rank.

Give C the elliptic group law with origin Q0. The nonconstant morphism
f-P sends origin to origin and hence is a group homomorphism. It sends Q1
to -2P, so Q1 has infinite order on C. Therefore C(Q) is infinite, and the
finite-to-one r-coordinate takes infinitely many values. This uses only
the standard [good-reduction torsion injection](https://public.websites.umich.edu/~asnowden/teaching/2013/679/L08.html#reduction-of-torsion-points)
and [origin-preserving morphism theorem](https://math.mit.edu/classes/18.783/2025/LectureNotes4.pdf).
The quartic's leading coefficient is negative and its value at r=0 is -3,
so its rational points are affine with r nonzero; there are no rational
zeros of the quartic. Thus y is nonzero and Delta=-11(y/4)^2<0, exactly
the physical cross-contact condition. The corresponding second radius
always belongs to the d=55 chord family.

In particular every odd multiple of P occurs in f(C(Q)): the image is a
coset containing P and with subgroup containing 2P. Direct addition gives

    2P=(16,-236),
    3P=(-91091/729,-10562552/19683),
    (r,y)=(91/135,10552/3645),   s=23684/3645.              (5b)

The new r lies just beyond the bounded census. It is an exact next
cross-family event, not an extrapolation. For any point of (5), one may
take a=(1+sqrt(-11))/2 and construct a contacting point at squared radius
s directly as

    b=(s+2+(y/4)sqrt(-11))/(2 bar(a)).

The equations |b|^2=s and |a-b|=1 follow by expansion; the same-radius
rotation of b still uses sqrt(-55). Thus all arithmetic points produce
actual field points, not merely abstract curve labels.

The next question is which of these contacts can coexist so as to restrict
a six-colour boundary relation, retaining actual phases and original ports.
A radius alone was already known to be insufficient in the type-M
experiments. Merely producing infinitely many contacts does not settle that
compatibility question.

**Hostile probe of the new construction (FINITE-EXACT).** Use both rational
points in (5b) and the seed, rotate their contacts by gamma^i for -3<=i<=3,
and take three consecutive points on each secondary circle (radial powers
-1,0,1). The complete induced graph has **49 vertices and 48 edges and is
a tree**. The receipt retains every coordinate, edge and a two-colouring.
Arithmetic is independently implemented in the degree-four basis
1,sqrt(-11),sqrt(-55),sqrt(-11)sqrt(-55). There are no hidden metric edges
at this depth. The points are a new test patch, not a context already glued
to the repository's eight fixed ports. Infinite supply is now proved;
additional compatible cross-edges are the missing coordinate.

## 3. A complete finite rotation-event reduction

**R3 (PROVED).** Let A and B be finite subsets of a conjugation-stable
number field K in C, both rotated about the marked centre 0. Keep A fixed
and replace B by gamma B, |gamma|=1. For each nonzero a in A, b in B put

    U=bar(a)b,       c=|a|^2+|b|^2-1.

A new unit cross-edge occurs exactly when

    U gamma^2 - c gamma + bar(U)=0.                         (6)

Hence there are at most 2mn unit-edge event angles, where m,n count the
nonzero points in the two pieces. Persistent edges involving the centre
are handled separately. Duplicate event angles are merged by exact equality.
Coincident vertices must also be excluded or explicitly identified; their
event equation is gamma=a/b.

**Proof.** Expand |a-gamma b|^2=1 and multiply by gamma. Since U!=0,
(6) has at most two roots. Its discriminant Delta=c^2-4|U|^2 is real.
For Delta<=0 the formula roots have norm one; Delta>0 has no norm-one
roots (equivalently |c|>2|U|). A nondegenerate allowed root is in K iff
sqrt(Delta) is in K. Tangencies and the zero-discriminant case are included.

For fixed pieces, every collision-free non-event rotation has the same
abstract graph: the two internal graphs and the persistent cross-edges.
Thus any stronger boundary-colouring restriction than this generic graph
must occur among the finite event list. This is complete for the chosen
pieces and hinge, with ports in the fixed piece; it is not a bound on the
number or size of pieces needed for six-colour forcing.

The half-angle variable gamma=(1+it)/(1-it) turns (6) into

    (2 Re(U)+c)t^2 + 4 Im(U)t + c-2 Re(U)=0.                (7)

Retain t=infinity (gamma=-1) separately. Over a real algebraic coefficient
field, field norms supply rational event polynomials; norm roots can be
spurious for the selected embedding, so authenticate the original equation.
The same method handles any finite rationally parameterized point family:
isolate zeros of distance, unit-edge and denominator polynomials, then test
one point per open interval and every algebraic event. In more than one
parameter, univariate Sturm alone is insufficient.

### Exact geometric positive and hostile control

Take the diamond vertices

    O=(0,0), A=(sqrt3,0), B=(sqrt3/2,1/2), C=(sqrt3/2,-1/2)

and a rotated copy sharing O. For 3/10<t<31/100 the program checks all 21
point pairs. There are ten persistent unit edges, no collision events, and
exactly one possible new edge, A to rotated A, at

    11t^2-1=0,       t=1/sqrt11.                            (8)

All other event norms have zero roots in this interval. Away from (8)
the graph has exactly 12 labelled three-colourings. At (8) it has none,
and exactly 384 labelled four-colourings: the classical Moser spindle.
This is an independently checked **event isolation** for a known gadget,
not a new spindle or a six-colour result. A rational-parameter-only search
would miss (8) entirely. In K=Q(zeta55), the angle survives but the diamonds
do not, as R1 already explains.

## 4. Rational Sturm chains as representation certificates

For a squarefree f in Q[X], the signed Euclidean chain is

    S0=f, S1=f', S_(j+1)=-rem(S_(j-1),S_j).

Its ordered transfer matrices satisfy

    (S_(j-1),S_j)^T = [[q_j,-1],[1,0]] (S_j,S_(j+1))^T.

The determinant is 1; the polynomial degree decreases at every nonterminal
division. This is a terminating algebraic recursion. It is not the
Fibonacci substitution recursion or the cutting-word meaning of *Sturmian*.
The useful connection is preservation of a signed boundary response, not
the shared surname or a common 2-by-2 matrix size.

For the eleven-squares paper's octic

    f=X^8-20X^7+178X^6-842X^5+1923X^4-496X^3
      -6754X^2+12420X-6865,

the exact chain has degrees 8,7,...,0 and variations 5,4,3 at -infinity,0,
+infinity. Thus it has one negative and one positive real root. The latter
is isolated in (3877/1000,1939/500). We also verify Appendix C's parameter
eliminant and the squarefree modular factor degrees (8), (7,1), (5,1,1,1)
at 29,7,131, using explicit Frobenius irreducibility certificates.

There is a second, representation-theoretic certificate. In E=Q[X]/(f),
let M_g be multiplication by g and form the rational symmetric matrix

    H_g(i,j)=Tr_E/Q(M_(g X^(i+j))),       0<=i,j<8.

Over R the algebra is R^r x C^c. A real factor contributes the sign of g
at that root; each nonzero complex factor contributes one positive and one
negative direction. Therefore

    signature(H_g)=sum_(real roots alpha) sign(g(alpha)).    (9)

This is the classical weighted Hermite signature identity, already used in
the repo's carry-corrected trace note. Here exact rational congruence
elimination, independent of Sturm division, yields

| multiplier g | positive / negative / zero inertia |
|---|---|
| 1 | 5 / 3 / 0 |
| X | 4 / 4 / 0 |
| X-3877/1000 | 4 / 4 / 0 |
| X-1939/500 | 3 / 5 / 0 |

The signature drop across the isolating interval is 2, hence exactly one
real root lies there. The rational multiplication representation certifies
the real root count without evaluating approximate eigenvalues.

It still does not select the physical root without an interval. For a
reducible defining polynomial, a nonzero quotient remainder can vanish at
the selected embedding; the inherited hostile
f=(X^2-1)(X^2-2), g=X^2-2, alpha in (5/4,3/2) is checked again.
The rational interval is the algebraic counterpart of the packing mark.

## 5. Why the golden and Paley fields meet here

On F11 let H={1,3,4,5,9}, and orient x->y exactly when y-x is in H.
There are no ties; multiplication by H and translation preserve this
actual binary predicate. This is the Paley tournament already present in
THM-4592. No arbitrary proof-obligation tournament is manufactured.

The group G=F11 semidirect H has order 55 and acts regularly on unordered
pairs. Consequently it has three ordered-pair orbits: diagonal, arcs,
reverse arcs. Every G-invariant matrix is in span(I,A,A^T), where A is
the adjacency matrix. The exact matrices give the standard Paley identities

    A+A^T=J-I,    AA^T=3I+2J,
    (A-A^T)^2=J-11I,
    char_A(X)=(X-5)(X^2+X+3)^5.                            (10)

Thus on the ten-dimensional mean-zero space, A generates Q(sqrt(-11));
over C the two five-dimensional channels have eigenvalues
(-1 +/- sqrt(-11))/2. The golden matrix Q=[[0,1],[1,1]] generates
Q(sqrt5). Their compositum contains sqrt(-55)=sqrt5 sqrt(-11), precisely
the second quadratic hinge family in R1. This is an explicit field map,
not an inference from the recurrence of 5,11 or55 in a table.

The same description also exposes a quotient loss. Symmetrizing A gives
J-I and loses the orientation completely. The three-dimensional invariant
operator algebra loses the chosen vertex/port response. Moreover a
G-invariant simple undirected graph on eleven vertices is empty or complete,
since G is transitive on unordered pairs. Nontrivial packing contact graphs
therefore cannot inherit this entire unmarked affine action.

The squared norms in (10) are controlled because of the **exact Paley
incidences**. They do not prove cancellation for a Collatz packet whose
source weights and terminal incidences are different. The field analogy
becomes a theorem only where the displayed maps preserve the target.

## 6. Boundary composition, recursion modes and the Lawvere boundary

For a finite graph G with named boundary B and q named colours, retain

    Z_G(beta)=#{proper q-colourings extending beta}.

If G1,G2 share exactly B, their interiors are disjoint and there are no
unrecorded cross-edges, then

    Z_(G1 union G2)(beta)=Z_G1(beta) Z_G2(beta).              (11)

Making some ports internal sums over their assignments. Gluing general
interfaces is tensor contraction; retaining only feasibility gives Boolean
relational composition. New metric cross-edges must be included as extra
constraints before contracting. Colour-permutation quotients need the
compatible assignments and multiplicities, not just a count of partitions.

This is the precise category-level construction: finite colour sets and
relations, or finite free modules and boundary tensors, with typed gluing.
It is not a claim that the chosen geometric category is cartesian closed.
For example on six labels equality and a cyclic-shift relation each have
both unary marginals equal to all six labels. Gluing equality with equality
has six solutions, while gluing equality with the shift has none. These
are abstract relation controls, not constructed six-colour geometric gadgets.

For tournaments, the analogous retained data depend on the operation.
Order-join permits scalar Hamiltonian-path multiplication; cyclic
substitution permits repeated visits to blocks and requires the inherited
path-cover profile. Fixed-Hamiltonian-path staircase presentations have
binom(n-1,2) independent bits, not that many isomorphism classes. Removing
both path endpoints keeps binom(n-3,2) bits and deletes 2n-5 bits. Boundary
labels must travel through this recursion.

The older three signed cell-count recurrences have characteristic factors
(X-1)^2, (X-1)^3, (X-1)^3(X+1). These are polynomial/parity modes (a root
1 of multiplicity m permits polynomial degree at most m-1). The Fibonacci
matrix instead has X^2-X-1 and exponential eigenvalues. Sturm has changing
quotients and decreasing degree. Their common algebraic language is an
ordered operator representation; their growth and termination mechanisms
are distinct. The historical claim that every recursion is cyclotomic is
not used.

Lawvere's [original theorem, section 1.1](https://conceptualmathematics.wordpress.com/wp-content/uploads/2025/03/diagonal_arguments.pdf)
requires an evaluation e:A x A->B representing every permitted map A->B
pointwise. Given h:B->B, represent a->h(e(a,a)); evaluation at its own
representing point yields a fixed point of h. For two-valued truth,
negation has no fixed point, so such total self-representation is impossible.
No such evaluation map is supplied by the tiling or tournament objects.
This preserves the correction MISTAKE-529 rather than reviving it.

The more useful fixed-point distinction here is **least versus greatest**.
For stopped Collatz T(1)=1, finite ROOT reachability is

    mu X. ({1} union T^(-1)(X)).                            (12)

A recursive parent-closure rule can certify indefinite compatibility; by
itself it need not produce a finite derivation in (12). The repo's
[chess trap analysis](chess_proof_system_20261012.md) already separates
greatest-fixed-point safety from ranked finite wins. For Chair44, the extra
weapon is that a translational period is an integer vector and unique
coarsening divides it by 2 at every level. Only zero survives all levels.
An arbitrary integer in a Collatz controller does not automatically divide
by 2 under its proposed refinement.

**R4 (PROVED exact-clock obstruction).** Suppose a finite quotient pi of
the positive integers intertwines stopped T with a deterministic map on
the quotient, and ROOT={1} is the pullback of a quotient predicate. Then
such a quotient cannot exist, even restricted to the powers of two.
If pi(2^j)=pi(2^k), j<k, applying the quotient transition j times would
identify ROOT 1 with non-ROOT 2^(k-j). Contradiction.

This is a future-language separation, not a claim that Collatz lacks a
finite proof or a finite instruction alphabet. A rule may carry an
unbounded exponent/clock, and an accelerated macro-edge may have variable
duration. The repair is exactly to retain that integer sidecar. A finite
diagram with discarded duration cannot be an exact step-for-step model.

## 7. Connection ledger and next decisive tests

| source -> target | map and preserved predicate | loss / repair | next decisive test |
|---|---|---|---|
| marked packing -> algebraic root | contact eliminant plus isolating interval; physical branch | polynomial alone loses branch | check every distance at the selected root |
| two geometric pieces -> boundary relation | enumerate (6), authenticate metric, contract (11) | independent lists lose joint compatibility | search all contact events, including algebraic isolated events |
| golden plus Paley representations -> hinge fields | sqrt5 times sqrt(-11)=sqrt(-55) | spectra lose ports and seed embedding | keep both the exact angle and the actual field graph |
| rational-radius gluing -> genus-one curve | discriminant square class gives (5) | radius loses cross-edge phase | compute rational points, then all marked cross-edges |
| recursive controller -> ROOT proof | sound steps plus source-preserving clock and well-founded rank | parent closure loses termination | exhibit a rank decrease or the independently positive bound |

For six colours, the current [component-palette theorem](six_colour_component_pairs_20261016a.md)
remains the strongest inherited rejection filter: all 210 boundary patterns
using at most four colours extend across any four-colourable exterior under
the stated incidence conditions. Event-generated candidates must first
survive this theorem, then supply a complete boundary relation.

For Collatz, [native Gram refinement debt](collatz_refinement_debt_20261016a.md)
is already paid uniformly in its stated scope. What remains is the uniform
paired-reset cancellation or an independent source-specific ROOT bound;
neither comes from Lawvere, a recurrence matrix, or an elliptic label.

For planar Jacobian, exact algebraic sign or pole checks can certify a
finite parameter family, but [inverse-clock compression](keller_inverse_clock_20261015f.md)
still owes its boundary return term. A trace representation and a faithful
polynomial inverse are different consumers.

Cards used: **search the statement before the method**, **retain local
incidence through counting**, **type every analogy**, and **correct the
object before sharpening the technique**. No new META-PATTERNS card is
needed: these computations strengthen existing triggers, especially the
need to test isolated contact events and retain a marked real embedding.

## 8. Reproduction and limits

Run:

```sh
python3 04-computation/experiments/recursive_interfaces_20261016b.py
python3 -O 04-computation/experiments/recursive_interfaces_20261016b.py
```

The two JSON outputs agree byte for byte. Dependencies: Python and SymPy;
the script explicitly selects exact Python ground arithmetic to avoid a
bundled SymPy/Flint modular-factor sorting incompatibility. It uses no
floating-point distance or spectral decisions, and checks are not disabled
by Python's optimized mode.

Frozen universes: the specified octic; all 21 hinge point-pairs on the
specified interval; all 3^7 and 4^7 labelled colouring assignments; all 55
affine maps and eleven vertices; the 33 powers 2^0,...,2^32; staircase
ranks for n=4,...,19; two 12-point cyclotomic orbits; and the stated reduced
rational mixed-shell parameters through denominator 100; good reductions
at 7 and 17 and the displayed elliptic additions; all 1176 pairs of the
49-point mixed patch. Positive controls
include the known spindle and the exact cross-edge (4); hostile controls
include reducible defining polynomials, a zero-diagonal inertia pivot,
identical unary marginals with incompatible gluing, and infinite legal
rotation orbits with only path graphs. Independent arithmetic paths are
Sturm division versus trace-form congruence, and symbolic geometry versus
exhaustive finite colourings. These do not replace an independent human
audit of the full note or any attached paper.
