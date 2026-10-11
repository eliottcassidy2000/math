# Forcing the pentagon boundary, and choosing an arithmetic geometry

2026-10-11. **PROVED, scoped and self-audited:** B1--B5 and the displayed
field-colouring deductions. **CITED:** the classical cycle-list fact and
Radchenko's Salem generator. **FINITE-EXACT:** list profiles and finite-field
colouring certificate. **OPEN:** a forcing six-colour context, a universal
six-colouring of the plane, and a seven-chromatic unit-distance graph.

## Inheritance and the interface

The [golden/taxicab boundary note](golden_taxicab_boundary_20261011.md)
already authenticates all25 algebraic vertices, all300 pair distances,
the induced pentagon with four private leaves at each vertex, and its
complete transfer-matrix formula. Repeating those coordinate checks is
unnecessary. The [quantized guards](unit_equality_quantized_guards_20261011e.md)
already give an equality-tester obstruction and a cited27-vertex necessary
floor for that different two-port problem. The [Moser-field audit](unit_equality_candidates_20261011e.md)
already excludes an entire four-colourable arithmetic field from producing
the desired six-colour equality tester.

The stronger inherited mechanism is **THM-4558, Lemma R**, in
[residue colourings of field planes](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md).
It credits the local colouring method to prior literature, proves its
denominator condition, and records entire four- and five-colourable fields.
The local constructions below are applications of that known mechanism.
In particular, the Polymath field obtained by adjoining the golden
quadratic to the Moser field is already known here to have chromatic5;
we do not replace that result with a weaker new conjecture.

Our new question is the *context* of the existing twenty ports. Let C be
a finite external graph whose intersection with the gadget is exactly those
ports. No additional edge may touch a hidden pentagon vertex. All actual
unit distances must be included when realizing C geometrically. The four
ports at hidden vertex i form B_i, and their used colours form S_i.
The available list is A_i=[6] minus S_i, so |A_i|>=2.

Closest mechanism: the authenticated list transfer. Hostile: five identical
two-colour lists. Corrected near miss: the five list cardinalities do not
record whether the omitted pairs agree. Least-used sidecar: a nonadjacent
pair of ports that can be recoloured with an unused colour.

## B1. Exact failure criterion, followed by a context barrier

A cycle with lists of size at least two is uncolourable precisely when
it is odd and all its lists are one identical two-element set.
This is a standard degree-choosability fact; see the
[primary lecture notes](https://www.math.tau.ac.il/~krivelev/teaching/co22/LectureNotes.pdf).

For completeness: if all lists are equal, two colours suffice exactly on
an even cycle, and three common colours suffice on any cycle. Otherwise
some adjacent lists differ; orient the cycle so L_0 has a colour absent
from L_last. Give vertex0 that colour and greedily colour the remaining
path. At the last vertex the initial colour is unavailable anyway, so
avoiding its predecessor is the only remaining restriction. This proves
the nonconstant-list case. If lists are equal with size>=3, the elementary
three-colouring settles the remaining case.

Consequently our pentagon boundary fails **iff** every B_i uses four
distinct colours and all five groups use the same four-colour set. It
does not fail merely because each hidden vertex has a two-element list.
One extra available colour at just one vertex repairs the common-pair
obstruction: the list profile ({0,1,2},{0,1},{0,1},{0,1},{0,1}) has exactly
two extensions. The common three-colour profile has30.

**B1, context theorem.** If C is five-colourable, then C joined to this
pentagon gadget is six-colourable.

**Proof.** Start with any proper five-colouring of C. Choose two ports in
one B_i that are nonadjacent in C and recolour both with the unused sixth
colour. They do not conflict with each other or any other vertex of C.
That group now uses at most three colours, making |A_i|>=3. All other
lists still have size>=2. The cycle-list criterion supplies an extension.

In the authenticated gadget each B_i is independent. More generally any
four distinct plane points contain a non-unit pair, since there is no
unit K4 in the plane, so the same argument applies to that broader
four-neighbour interface whenever its other hypotheses hold.

Thus a six-colourable C that forces the bad boundary relation must have
chi(C)=6; adjoining the hidden pentagon would give a seven-chromatic
unit-distance witness. This does not exclude other interfaces, additional
neighbours of the hidden vertices, or contexts already not six-colourable.
It shows exactly why repeatedly adjoining five-colourable gadgets under
this unchanged interface cannot accomplish the proposed forcing task.

## B2. The synchronization has a measurable information cost

Let a finite context be equipped with the uniform law on its proper
six-colourings, and suppose each B_i always uses four distinct colours.
Then its omitted pair A_i is uniform on the15 pairs, by permutation
symmetry of the six colour names. If the context forces all five pairs
equal, their joint law is the uniform diagonal law. Hence its total
correlation is exactly

    sum_i H(A_i)-H(A_0,...,A_4)=4 log15.                 (1)

Under the product of the five marginals, the diagonal has probability
15^-4. More generally, if its actual probability is p, data processing
gives

    D(P_joint || product_i P_i) >= d(p || 15^-4),        (2)

where d is binary relative entropy. This is an exact finite analogue of
the [supplied Kakeya paper's support comparison](kakeya_marked_context_intake_20261011.md).
It does not show that a geometric context can pay this information cost.
If one first conditions on all groups using four distinct colours, that
conditioning's probability must also be retained; it is not free forcing.

## B3. Golden Cartesian coordinates remain bipartite

Let K=Q(phi), phi^2=phi+1, and R=Z_(2)[phi], the ring with odd rational
denominators. Its residue ring modulo2 is F4, since T^2+T+1 is irreducible.

If a,b in K satisfy a^2+b^2=1, then a,b belong to R. To see this, suppose
2^s a=A, 2^s b=B are in R, at least one not in2R, with minimal s>=1.
Modulo2 the equation gives A=B, and both are nonzero. Write A=B+2C.
Then A^2+B^2=2(B^2+2BC+2C^2) is divisible by2 but not4, whereas it is
2^(2s). This is impossible. Now reduce the original unit equation:

    a_bar+b_bar=1 in F4.                               (3)

The F2-linear functional lambda(a+b phi)=a modulo2 has lambda(1)=1.
For a finite connected unit-distance graph with vertices in K^2, choose
an actual vertex origin. All vertex differences lie in R^2 because every
edge difference does. Colour a displacement(x,y) by lambda(x+y). Every
unit edge changes that colour by1. Therefore **every such graph is
bipartite**. The same local proof works for a real number field at an
unramified prime over2, with a chosen residue functional taking1 to1.

There is a sharper integral statement in the golden case. If both
coordinates lie in Z[phi], their two real embeddings satisfy the same
unit equation. Each coordinate a+b phi and its conjugate lie in[-1,1],
so |b|sqrt5<=2, forcing b=0. The only integral unit directions are the
four coordinate axes. Field coordinates do allow more directions, for
example (-phi/(phi+2), 2phi/(phi+2)); the denominator has norm5.
They still do not evade the bipartite colouring.

This does not make the golden pentagon bipartite. After fixing one edge
as the real unit, its **complex** vertices lie in Q(zeta_5), not in the
Cartesian plane Q(phi)^2. The sine coordinate and the metric embedding
are additional data. A golden norm unit such as phi is not automatically
a Euclidean unit vector.

## B4. A reciprocal family interpolates the golden and Salem cases

For an integer b consider

    P_b(z)=z^4-z^3+(2-b)z^2-z+1,
    P_b(z)/z^2=T^2-T-b,       T=z+z^-1.                 (4)

The discriminant of P_b is

    (1+4b)^2(2-b)(6-b).                                 (5)

| b | Trace discriminant | Geometry of the roots |
|---:|---:|---|
|1|5|P_1=Phi_10; the golden decagon, all four roots on the unit circle|
|2|9|(z-1)^2(z^2+z+1), a repeated-root boundary|
|3|13|two reciprocal positive real roots and two nonreal unit roots|
|4|17|the same mixed signature|
|5|21|the same mixed signature|
|6|25|(z+1)^2(z^2-3z+1), the second boundary|
|7|29|both trace roots outside[-2,2]; no nonreal unit-root pair|

For b=3,4,5 the quartics are irreducible over Q, and their unit-circle
root alpha is not a root of unity: a real conjugate has modulus>1.
Thus all alpha^j are distinct unit directions in the rank-four additive
group Z[alpha]. This is Radchenko's [Theorem1](https://arxiv.org/pdf/1807.03726)
at b=3; the adjacent b values and their trace interpretation are our
elementary comparison, not a new invention of the Salem construction.
The monic reciprocal polynomial supplies integral positive and negative
powers. This is a meaningful enlargement of the finite decagon generator.

The distinction from golden Cartesian integrality is exact: alpha is a
complex algebraic integer, but its real and imaginary Cartesian coordinates
involve division by2 and by2i. Integrality in Z[alpha] is not integrality
of both Cartesian coordinates in a totally real field.

## B5. Audit all local places: the four quartic fields have chi = 3,3,2,3

In the order b=1,3,4,5, the full complex field planes have chromatic
numbers3,3,2,3. The dyadic calculation first gives a useful but weaker
four-colour quotient for the three odd b values. Ramified primes below
then improve all three to3. This is an example of why one residue graph
must not be mistaken for the chromatic number of its source field.

For b=1,3,5, P_b modulo2 is the irreducible polynomial
T^4+T^3+T^2+T+1, giving F16. Complex conjugation alpha->alpha^-1 reduces
to z->z^4. If a field element z in Q(alpha) has Euclidean norm1, then it
has no factor2 in its denominator relative to Z_(2)[alpha]. Indeed take
2^s z integral at2 with nonzero residue; the product with its conjugate
has nonzero residue, contradicting 2^(2s) when s>0. Therefore its residue
satisfies

    z_bar^5=1.                                         (6)

Represent F16 elements by four-bit polynomials, with reduction polynomial
11111. The five allowed nonzero differences are

    {1,2,4,8,15}.

A proper four-colouring, listed on labels0,...,15, is

    [0,1,1,0,1,0,2,1,2,0,3,2,3,2,0,3].                (7)

All40 edges are checked in the replay. As in B3, choose one actual vertex
as origin in each finite connected field graph, reduce the displacements,
then apply(7). **Every unit-distance graph whose complex coordinates lie
in any of these three fields is four-colourable.** More directions, greater
degree, and larger Minkowski sums inside the same field cannot overcome
this entire-field upper bound. The finite quotient itself has independence
number5, checked exhaustively, hence is not three-colourable; that does
not assert every field graph realizes that quotient or has chromatic4.

For b=1 the field is Q(zeta_5). At its unique prime over5, conjugation
acts trivially on the residue F5. Every unit vector is a local unit,
and its residue is plus or minus1. The residue Cayley graph is C5,
so Lemma R gives a three-colouring.

For b=3 use the trace root T=1 modulo3 of T^2-T-3. Hensel's derivative
is nonzero, and (T+2)(T-3)=-3 gives v_3(T+2)=1. At alpha=-1+w,
the quadratic equation becomes

    w^2-(T+2)w+(T+2)=0.

It is Eisenstein at3, with residue F3 and a conjugation-fixed prime.
Norm-one residues are plus or minus1, so Lemma R again gives three
colours. For b=5, directly

    P_5(1+w)=w^4+3w^3-6w-3

is Eisenstein at3. The unique prime has residue F3 and is fixed by
conjugation, giving the same three-colouring. The induced odd cycles
below supply the matching lower bounds in all three cases. Thus these
fields are closed even for reaching chromatic4, despite their dyadic
quotient being four-chromatic.

At b=4 the residue polynomial is reducible, so that particular proof
does not apply. A different local place gives a stronger conclusion:
**the full Salem17 field is bipartite** in its unit-circle embedding.

Here is the denominator audit, rather than an inference from a ring
character alone. The trace field has T^2-T-4=0. At its even 2-adic root,
Hensel lifting is unique since2T-1 is odd, and T(T-1)=4 forces v_2(T)=2.
The quadratic polynomial of alpha is z^2-Tz+1. Setting z=1+w gives

    w^2+(2-T)w+(2-T),                                  (8)

which is Eisenstein at2. This defines a totally ramified quadratic local
extension with residue F2. The involution alpha->alpha^-1 fixes its
unique prime. For any field element z with z conjugate(z)=1, invariance
of the valuation gives2v(z)=0. Thus z is a local unit with residue1.
Reducing a connected graph's displacements at this place yields a proper
two-colouring. It handles arbitrary field denominators, not just the
integral module. Exact lifts through2^32 check the local arithmetic; the
Hensel/Eisenstein argument is the all-depth proof.

This strengthens the cheap initial test alpha->1 on Z[alpha]. It also
shows why “ramification breaks the earlier colouring proof” is not a
positive chromatic signal: another local place may colour more strongly.

### The same quartics provide authenticated five-, six- and seven-cycles

The coefficient relation P_b(alpha)=0 is a closed walk of unit vectors.
For b=3,4,5 take the successive steps

    1, -alpha, then (b-2) separate steps -alpha^2,
    then -alpha^3, alpha^4.

Their numbers are5,6,7. Partial sums give distinct vertices: before the
last step, any equality would be a nonzero relation among1,alpha,alpha^2,
alpha^3, contradicting degree4. The replay checks every pair and proves
that each of these three concrete graphs is an **induced** unit cycle.
At b=1 the analogous five steps are1,-alpha,+alpha^2,-alpha^3,alpha^4.
All coordinates and edges are stored in the output's power basis.

This is an exact recurring gadget connecting reciprocal trace arithmetic
to the5/6/7 cycle lengths. It is not an identification with spherical,
Euclidean or hyperbolic curvature. The b=3 and b=5 odd cycles witness a
lower chromatic bound3 in their fields; the b=4 even cycle is consistent
with its proved bipartite colouring.

For b=3,4,5, the infinite directions alpha^j also let the private-leaf
compiler stay inside the same field and integral module. At each step
only finitely many directions cause coincidences or an unwanted unit edge
with an existing point (two distinct unit circles meet in at most two
points). Choose another power and retain its exponent. This is a repeatable
exact gadget constructor, with its field-colouring upper bound attached.

Even combining these four quartic fields does not escape the four-colour
test. Their b=1,3,5 completions at2 are the same unramified degree-four
local field, with conjugation its second Frobenius power. Choose the
ramified degree-two place of b=4 from(8). The unramified and totally
ramified local fields are linearly disjoint, so their involutions combine
on the compositum, with residue action still the second Frobenius on F16.
Norm-one differences remain local units and reduce to the same five
differences. Thus the compositum is four-colourable too.

The [cyclotomic continuation](cyclotomic_palette_escape_20261011.md)
finds a genuine escape from this *particular local test* by adjoining
eleventh roots of unity to the golden field of fifth roots. It supplies
an explicit unit whose reduction collapses at one dyadic prime, and keeps
the distinction between breaking a colouring test and proving high
chromatic number.

## Useful next gadgets and declared stopping reasons

* Keep the twenty-port list relation as a test harness. A candidate context
  must retain its actual geometry, a proper six-colouring, and a complete
  boundary-relation certificate. Under this interface, a five-colouring
  immediately rejects it as a forcing candidate.
* Use the reciprocal quartic as a direction generator with an accompanying
  residue-colouring test. The b=1,3,5 full fields are three-chromatic,
  and b=4 is bipartite even with field denominators. To seek
  stronger forcing, combine geometries outside these closed fields and
  recompute their joint local colouring tests.
* An alternative interface can allow more than four neighbours at selected
  hidden vertices, or share internal constraint vertices. It must recompute
  the full list relation; the old cycle proof no longer settles that object.

[Exact checker](../../04-computation/experiments/marked_context_transfers_20261011.py),
[output](marked_context_transfers_20261011.json): all11^5=161051 list profiles
over a four-colour palette with list sizes>=2, independent literal
six-colour controls, the quartic identities, and all65536 subsets of F16
for the independence computation. The general six-colour list claim is
proved above, not inferred from the four-colour finite universe.
