# Same-clock closure and three different notions of selection

2026-10-31. **PROVED WITH INDEPENDENT PROOF AUDIT.** These are scoped
categorical/topological statements; no topos construction or unrestricted LRC
proof is asserted. Classical categorical structure is not claimed as novel.

## C1. A closed symmetric monoidal category that retains the clock

Fix a nondegenerate rational closed interval I. Let B_I consist of all finite
unions of intervals and points inside I with rational endpoints and explicit
open/closed endpoint flags. Regard B_I as a posetal category: there is one
arrow A→B exactly when A⊆B. It is a Boolean algebra under finite operations.

Set tensor A⊗B=A∩B and unit I. This is a symmetric monoidal category,
indeed a cartesian one. Its internal hom is

    [A,B]=(I\A)∪B,

because for every C,

    C∩A⊆B  iff  C⊆(I\A)∪B.

This proves the required adjunction rather than merely calling intersection
a tensor. The evaluation arrow [A,B]∩A→B preserves the actual clock t.
A witness t∈[A,B] and a different witness u∈A do not satisfy its domain;
the same-clock intersection requires t=u. Explicit endpoint flags make the
Boolean complement and hence the internal hom representable exactly.

This algebra has finite joins, not arbitrary completeness. Arbitrary unions
can require infinitely many components or nonrational boundary data. It is
not asserted to be a complete quantale. It is also not an ordinary topos:
Sub(I) contains at least two elements, but every hom-set of a posetal category
has at most one member, so no object Ω can classify all subobjects of I by
arrows I→Ω. In particular this is not, by itself, a realizability topos.

The closed-only fragment fails monoidal closure. With A=[a,(a+b)/2] and B=∅,
there is no largest closed C⊆I disjoint from A. The union of closed intervals
[(a+b)/2+ε,b] over positive rational ε is ((a+b)/2,b], whose closure hits A.
Any proposed greatest closed C would contain that union and its closure.
Thus retaining only closed interval unions is sufficient for conjunction
but not for the internal-hom adjunction.

## C2. Existential truth is not strong monoidal

Let E(A)=1 if A is nonempty and0 otherwise. This monotone map preserves
finite unions, but not tensor:

    E(A∩B) can be0 while E(A)∧E(B)=1.

For A={a},B={b} with a!=b this is immediate. In slice language, the product
of witness-bearing subobjects A→I and B→I is the fibre product A×_I B,
canonically A∩B, not the ordinary product A×B of independent times.
Discarding the clock before conjunction is the precise failed implication.

The same loss occurs even inside the exact bad-measure budget. On I=[0,1],
with radius1/4 and two unit speeds, phases(0,1/2) have common safe set
{1/4,3/4}; phases(1/4,3/4) have common safe set{0,1/2,1}. Both inputs have
identical speeds, interval width and individual bad measures(1/2,1/2), and
both satisfy the equality budget. Their witness sets are disjoint. Therefore
those numeric budget data alone cannot specify a correct common-time witness
for every phase lift. The constructive algorithm keeps the phase sidecar.

## C3. Nonempty phase fibres have no continuous global selector

Let T=R/Z, choose an integer m>=2 and0<r<=1/2, and define

    S_phi={t∈T: ||mt+phi||>=r},
    E={(phi,t)∈T²:t∈S_phi},  π:E→T.

Every fibre is nonempty: solve mt+phi=1/2 modulo1. Local continuous sections
exist by choosing a local real lift of phi and setting t=(1/2-phi)/m modulo1.
Nevertheless π has no continuous global section.

**Proof.** Suppose s:T→T were a section. The circle map

    g(phi)=m s(phi)+phi

would have image in the closed arc[r,1-r], which avoids0 and is contractible
(the point1/2 when r=1/2). Thus g has degree0. On the other hand

    deg(g)=m deg(s)+1,

which cannot vanish because deg(s) is an integer and m>=2. QED.

Equivalently, use a real lift s~ along phi∈[0,1]. Its endpoint increment is
an integer d. The lift m s~(phi)+phi has increment md+1, while remaining
inside a single translate of the proper safe arc forces increment0.
This gives the same contradiction without choosing a gauge across the loop.

At r=1/2 this is exactly the failure of a global section of a nontrivial
m-fold circle covering. For r<1/2 the fibres are unions of m closed arcs,
so the obstruction persists even with nonzero-width safe components.

Classical topology input: circle winding/degree and the degree-m power map;
see Allen Hatcher, *Algebraic Topology*, §1.1 and Example2.32 in §2.2,
https://pi.math.cornell.edu/~hatcher/AT/AT.pdf . The primary textbook was
checked for the circle power-map degree statement. The phase-family
application and elementary lift argument are given explicitly above.

## C4. What the constructive theorem does and does not supply

`CONSTRUCTIVE_CONDUCTOR.md` gives a polynomial-bit selector on finite
rational encodings satisfying its exact budget. It is permitted to branch
discontinuously on those encodings. C3 rules out a continuous section over
the entire circle of real phase data; it does not refute that rational-input
algorithm. Pointwise existence, a computable rational realizer, and a
continuous global section are three different predicates.

Source/target contract: full clock predicates -> truth values loses overlap;
full phase-coded interval data -> rational witness retains it; topological
phase families -> a putative continuous section adds a stronger regularity
requirement that fails by winding. None of these transfers identifies a
Boolean predicate category with an effective topos or proves a major
conjecture from a formal categorical vocabulary.
