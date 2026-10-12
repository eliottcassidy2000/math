# Unit-path coverage and marked colouring contexts in Q(zeta55)

2026-10-11. **PROVED, scoped and self-audited:** the connected-component
localization, integral three-colouring, explicit rational connectors,
field-preserving Cartesian construction, and preservation of its marked
boundary relation. **CITED:** the combinatorial high-girth theorem below.
**FINITE-EXACT:** the declared controls. **OPEN:** a six-chromatic forcing
context in this field, the field's full chromatic number, and the plane's
six/seven question under the separately stated H6 premise. No literature
priority is claimed for the algebraic lemmas.

## Inheritance and the construction target

The closest mechanisms are the nonintegral unit in
[cyclotomic palette escape](cyclotomic_palette_escape_20261011.md), the
complete boundary relation in
[boundary palette compression](boundary_palette_compression_20261011.md),
and the fixed-prime test in
[THM-4558, residue colourings](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md).
The CM/Kronecker argument already appears for conductor12 in
[THM-494, transverse resonance](../../01-canon/theorems/THM-494-transverse-resonance-is-rational-cosine-the-bisector-is-off-ladder.md).
The present application must keep the difference between integral and
nonintegral unit directions.

The hostile is the incoming
[four-chromatic field-escape family](unit_field_escape_recursion_20261011f.md):
escaping arbitrarily many small-prime screens does not force more colours.
The corrected near miss is transferring a residue graph's lower chromatic
bound to its source field. The least-used sidecar is now the **prime ideal
support of a displacement's denominator**.

The board is **marked source / denominator support / unit-path reachability /
full boundary relation / intentional cross-edge / clock divisor**. The anchor
is an independently paid source certificate, the niche is arithmetic plane
geometry, and the wildcard is the time-form comparison in the
[companion note](timeform_ramification_20261011.md).

Put L=Q(zeta55), A=Z[zeta55], and let c denote complex conjugation. Let
U={u in L: u*c(u)=1}. These are geometric unit directions; this does not mean
all units of the ring A.

## 1. The integral graph is exactly three-chromatic; the full graph has no triangles

If u is integral and u*c(u)=1, every Galois conjugate has complex absolute
value one, because c commutes with the Galois group. Kronecker's elementary
argument then makes u a root of unity: the monic integral polynomials whose
roots are the conjugates of u^n have uniformly bounded coefficients, so
there are only finitely many possible u^n. The roots of unity in L are mu110.

Set alpha=zeta55^11 and beta=zeta55^5. Since zeta55=alpha*beta^(-2), the
40 products alpha^i beta^j, 0<=i<4, 0<=j<10, form an integral basis of A.
Define an additive map ell:A->F3 by giving every such basis vector value1.
Its values on the five alpha powers are (1,1,1,1,-1); its values on the
eleven beta powers are (1,...,1,-1), with ten initial ones. Thus

    ell( +/- alpha^i beta^j ) = +/- ell(alpha^i)ell(beta^j) != 0.

Every integral unit difference is one of these110 roots. Hence ell gives a
proper three-colouring of the full unit-distance graph on A. The cumulative
alpha-edge pentagon is an odd cycle, so three is necessary.

The full field graph, despite allowing nonintegral directions, is triangle
free. An equilateral unit triangle would give w in L with w*c(w)=1 and
w+c(w)=1, hence w^2-w+1=0. This would put zeta3 in L. But adjoining zeta3
gives Q(zeta165), of degree80 rather than40. In particular the neighbours
of every point in L form an independent set. This is not a chromatic upper
bound: all additional colour forcing must come from the external context.

## 2. An exact universal description of all unit-reachable points

Let Gamma be the connected component of0 in the full field unit graph.
Then

    Gamma = sum_{finite} Z u, u in U.

It is a ring: multiplying two finite sums just multiplies their unit
directions, and U is multiplicatively closed. It contains A because roots
of unity generate A. For a prime ideal P of A fixed by c, every u in U has
v_P(u)=0, since 2v_P(u)=v_P(u*c(u))=0. Therefore Gamma is contained in A_P.

If c(P)!=P, choose a in c(P) but not P. The direction a/c(a) has a pole at
P. Thus Gamma is not contained in A_P for exactly those nonfixed primes.

Here is the needed overring fact, with its proof retained. For any ring
A subset R subset L, where A is Dedekind,

    R = intersection_{P: R subset A_P} A_P.                    (1)

If x is not in R, the nonzero proper ideal I={a in A: ax in R} lies in some
maximal ideal P. Then x is not in R localized at A\P. An overring of the
DVR A_P is either A_P or L: one element of negative valuation supplies the
inverse of a uniformizer. The latter alternative would contain x, so the
localized ring is A_P and R is contained in A_P. This proves (1), including
the direction that mere pole testing would otherwise leave missing.

Applying (1) gives the exact component theorem

    Gamma = intersection_{c(P)=P} A_P.                         (2)

Because L/Q is abelian, conjugation fixes either all or none of the prime
ideals over a given rational prime. Let S be the rational primes in the
second case. Then

    Gamma = A[S^(-1)],                                        (3)

and every component of the full graph is an additive coset of Gamma.
This is genuine universal **unit-path** coverage, with no colouring claim.
Translation identifies all these component graphs, so chi(L)=chi(Gamma).

For p not dividing55, the inherited decomposition calculation gives

    p not in S iff p = -1 mod5 and (p/11) = -1.                (4)

The primes5 and11 belong to S. For instance 2,3,5,7,11,13,17 are inverted,
whereas19 and29 are not. The two ramified cases follow from their unramified
projections: 5 has odd order5 modulo11, and11 is1 modulo5, so neither
decomposition group contains complex conjugation.

Consequently 1/19 is not unit-reachable from0 in this field, whereas 1/2,
1/3 and1/11 are. The description does not give a uniformly short path for
every allowed displacement. An exhaustive enumeration of finite unit words
does terminate for each point proved to belong to (3).

Any connected finite candidate can be translated into Gamma. Thus finite
searches may enumerate points of A with denominators supported on finite
subsets of S, authenticate every pairwise distance, and eventually include
every connected field witness. Fixing one denominator bank is a restriction,
not universal coverage by the entire field.

## 3. Small explicit connectors, and the price of leaving the integral lattice

The quadratic subfield Q(sqrt(-55)) lies in L. Put

    eta=(-3+sqrt(-55))/8,     eta*c(eta)=1.

The points

    0, 1, 1+eta, 1/4                                        (5)

form an induced three-edge unit path: its last edge is c(eta), since
1+eta+c(eta)=1/4. The nonadjacent squared distances are5/4,9/16 and1/16.
Powers of eta are distinct: if eta were a root of unity, its rational trace
-3/4 would be an algebraic integer, which it is not.

For alpha=zeta5, define

    u=(alpha-3)/(alpha^(-1)-3)
      =(-24alpha^3-32alpha^2-75alpha+72)/121.

All four conjugates u_i are unit directions, and exact trace gives

    sum_i u_i=419/121,     11 sum_i u_i - 38 = 1/11.          (6)

Equation (6) is an82-edge unit walk. Repeated points may be erased to obtain
a path; no claim that this particular ordering is induced is needed.
Together, (5)--(6) are actual constructive receipts for two denominator
changes, not just norm-one rotations.

There is also a clean obstruction to extending the integral additive
three-colouring. Since1/3 belongs to Gamma, every additive map Gamma->F3
vanishes: ell(x)=ell(3(x/3))=0. This excludes that colouring mechanism on
Gamma, not arbitrary three-colourings of Gamma. In particular it does not
establish four colours, much less six.

## 4. Field-preserving Cartesian refinement retains the marked boundary

Let G,H be finite induced unit graphs in L, with point sets A0,B0. Choose a
marked h0 in H and translate B0 so its coordinate is0. For a unit u, represent
(a,b) by a+u*b. With u chosen outside a finite exceptional set, the full
induced graph is exactly G square H.

Indeed one-coordinate changes reproduce precisely the original edges.
When both coordinate differences A1,B1 are nonzero, collision excludes at
most one u, and |A1+u*B1|=1 excludes at most two points on the unit circle.
The distinct powers eta^j therefore supply a deterministic search. The
first 3*binom(|G||H|,2)+1 powers suffice by this loose counting bound.
The entire original fibre G times {h0}, including its coordinates, is fixed.

Now fix k>=2 and suppose H has a proper k-colouring h. Every proper
k-colouring g of G extends to the product by

    C(a,b)=g(a)+h(b)-h(h0) mod k.                            (7)

Conversely restriction to the marked fibre gives a colouring of G. For any
ordered ports in that fibre, their complete boundary equality relation is
therefore **unchanged**. This is stronger than equality of chromatic numbers.
Deleting product vertices or edges while retaining the ports cannot remove
any of these extensions, so it cannot strengthen that boundary relation.

The exact control realizes C5 square P3 with15 vertices and25 edges in
Q(zeta5,sqrt(-55)), authenticates all105 pairs, and extends all30 labelled
three-colourings of the original pentagon. This remains a three-colourable
control, not a six-colour gadget.

**CITED connection.** Tomon's Theorem1.5 produces a same-chromatic high-girth
induced subgraph of a suitable Cartesian power of a finite graph. Combining
that preprint input with the field-preserving construction above retains
the coordinate field L as well. It does not guarantee retention of chosen
ports or their boundary relation. The geometric crowding price was already
proved in [the earlier girth packet](plane_girth_crowding_20261011b.md).
[Primary source](https://arxiv.org/html/2610.03517v1).

## 5. The deliberate cross-edge gate

To strengthen a boundary relation, leave the generic Cartesian regime and
authenticate at least one additional constraint. For nonzero A1,B1 in L,
set

    C=1-A1*c(A1)-B1*c(B1),   D=A1*c(A1)*B1*c(B1).

Writing w=c(A1)*B1*u turns |A1+u*B1|=1 and u*c(u)=1 into

    w+c(w)=C,     w*c(w)=D.

Thus the two possible orientations are the roots of w^2-Cw+D. Let F=L^c.
Since L=F(sqrt(-55)), an orientation belongs to L exactly when

    4D-C^2 = 55 h^2 for some h in F.                        (8)

In that case u=(C+h*sqrt(-55))/(2*c(A1)*B1), with the two signs when h!=0.
The zero-discriminant case is included. The sign at the actual real embedding
is automatic in (8). Every proposed orientation still needs all collision,
unit-edge and nonedge checks; (8) alone is not a graph certificate.

A cheap earlier test uses **every real embedding** of F: each must send
4D-C^2 to a nonnegative number. Equivalently the three lengths
|sigma(A1)|, |sigma(B1)|,1 must satisfy the triangle inequalities in every
conjugate plane. All Galois embeddings commute with complex conjugation,
so they preserve unit edges and nonedges. Passing these inequalities is
necessary but does not establish the square class in (8).

The cheapest hostile is two unit coordinate changes: D=1,C=-1, which would
require sqrt(-3) in L. Such a seam is impossible. One must use nonunit port
displacements or a different construction. This precisely identifies why
copying an equilateral/diamond trick into this field fails.

The next proposed search therefore has three independent obligations:
an allowed denominator bank, a square-class solution to (8) for a deliberate
seam, and a complete six-colour boundary relation that is genuinely smaller
after the seam. These are executable arithmetic and colouring questions;
their simultaneous success is still open.

### A complete small-family obstruction, and a surviving nonintegral seam

For the unit pentagon above, put phi=1+alpha+alpha^(-1). Its complete
nonzero ordered displacement set is

    D5=mu10 union phi*mu10.                                (9)

Consider all25 labelled points a+u*b from two such pentagons. If their map
is injective, their induced graph is **always** the Cartesian product,
for every unit u in L, not merely for generic u. To prove this, a possible
cross-edge has A1=r*epsilon, B1=s*delta with r,s in{1,phi} and
epsilon,delta in mu10. Write z=u*delta/epsilon. If r=s=1 then z+c(z)=-1,
which would introduce zeta3. In all other cases

    z+c(z)=(1-r^2-s^2)/(rs)=-phi,

so z=alpha^2 or alpha^3 and u lies in mu10. But (9) is invariant under
mu10. Hence D5 intersects -u*D5 and two distinct labelled product vertices
collide, contradicting injectivity. This rules out that entire25-point
Minkowski family as a source of strengthened pentagon constraints.

Replacing the second pentagon by the three-point unit path leaves one new
length pair (phi,2). Its square-class numerator in (8) is7phi-1, whose
quadratic conjugate is7(1-phi)-1 and whose norm is-55. It has opposite signs
at the two real embeddings, so it cannot equal55 times a square in F.
The pair (1,2) is degenerate and gives a root-of-unity collision; the other
cases were just treated. Therefore every injective C5/path product is also
exactly Cartesian. The finite control is not merely a lucky orientation.
The occurrence of55 here supplies an obstruction through its **negative
norm**, not permission to take a square root in the CM field.

A positive seam survives after changing the port distance. The points

    0, -eta, 3/4, -c(eta)                                 (10)

form an induced unit C4. Its diagonals have lengths3/4 and sqrt(55)/4.
This uses 4-(3/4)^2=55/16, so both conjugate circle intersections really
belong to L. It is a two-colourable connector, not a forcing gadget. It
provides an actual admissible nonintegral seam for searches beyond the
two closed product families, with every coordinate and nonedge retained.

## Reproduction and connections

[Program](../../04-computation/experiments/source_context_localization_20261011.py)
and [output](source_context_localization_20261011.json) give the finite
controls. The aggregate has 2,473 always-active exact checks with byte-identical
ordinary and optimized output. The infinite component theorem is proved in
(1)--(4), not inferred
from the tested primes. This note's proof audit is self-audit; no independent
agent audit is claimed.

The connection to ROOT is exact at the level of proof operations: retaining
a source fibre or integer through refinement preserves what has already been
proved, while a new forcing relation requires an additional authenticated
connection. The arithmetic meanings of unit paths and Collatz words are
different. Neither universal reachability (3) nor the Cartesian lift (7)
supplies a positive ROOT score.
