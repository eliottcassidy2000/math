# Order-ten prediction by canonical marked-root selection

24 October 2026. **PROVED:** canonical marked-root selection, completeness of
the reconstructed automorphism actions, and the exact prediction identity.
**FINITE-EXACT:** parent-atlas construction and order-ten upper-trace evaluation.
This is a canonical-construction method, without an algorithm-priority claim.

## 1. Result and source boundary

Starting with the certified6880 ordinary order-eight parent classes and their
complete vertex automorphism actions, the production pipeline obtains

    constructed order-nine parents: 191536;
    order-ten pointed attachment states: 97213472;
    order-ten stable-root upper trace U10: 97336068;
    independently proved structural correction E10: 5508;
    predicted a10 = (97336068 - 5508)/10 = 9733056.

Neither constructed count nor predicted target count occurs in the production
C++ sources. No known target value calibrates the computation. No production
canonicalization is performed on nine or ten vertices. Only2088 canonicalizations
of eight-vertex deleted parents are required while constructing the parent atlas.
The bulk order-ten scan uses no deletion, isomorphism, or target automorphism
calculation. It does retain actual parent adjacency and actual vertex actions;
it is not a calculation from scalar class counts or an unmarked metagraph.

## 2. Pointed extension states are exact

Fix a complete atlas of ordinary n-vertex tournaments P, with their actual
orientation-preserving permutation groups G_P. For each G_P-orbit on subsets,
choose its least binary mask S. Adjoin a root * with v->* exactly when v is in S.
Removing the root is inverse to this construction. Thus pairs(P,S) are in
bijection with ordinary vertex-pointed(n+1)-tournaments, not with individual
vertices and not directly with unpointed classes. The pointed states over one
child T are exactly its vertex orbits under Aut(T).

Stable directed refinement starts by degree and repeatedly records the old
color and outcounts into its current ordered colors. Lexicographic renaming
makes the numerical color order an isomorphism-invariant order. Therefore the
first stable singleton, when one exists, is a canonically chosen vertex.

## 3. Singleton selection and its complete action group

For every pointed extension T with a stable singleton, select the state if and
only if its root is the first stable singleton. Exactly one pointed state is
selected over each such unpointed child: the canonical vertex is unique and
fixed by every automorphism.

In a selected state(P,S), every child automorphism fixes the root and restricts
to a parent automorphism preserving S. Conversely every h in G_P with h(S)=S
extends uniquely over the fixed root. Thus the complete child group is exactly

    Aut(T) = { h in G_P : h(S)=S }, extended by h(*)=*.

No target isomorphism search is required. Notice that selecting an arbitrary
labelled singleton would not be valid; the invariant color order is essential.
A stable cell of size larger than one cannot in general substitute for a
singleton without supplying the missing root-orbit information.

## 4. Singleton-free residuals: parent-only marked keys

At order nine, stable cells have odd sizes. If none is a singleton, the only
possible partitions are9 and3+3+3. The production builder rejects any other
singleton-free pattern. The actual generated residual consists of232 pointed
states and56 unpointed classes, but these numbers are emitted results, not guards.

For any root v of a child T, delete v and canonically label the eight-vertex
parent. Let P_v be its supplied atlas index; let f_v map the original surviving
child labels to the canonical parent's labels; and let M_v be the image under
f_v of the vertices pointing to v. Define the pointed key

    K(T,v) = ( P_v, min_(h in G_Pv) h(M_v) ).

Canonical labellings f_v need not be unique. Two choices differ by an element
of G_Pv, so the minimized subset mask is independent of the chosen labelling.
Two pointed tournaments have the same key if and only if they are isomorphic:
a pointed isomorphism restricts to the corresponding parent isomorphism and
preserves the attachment; conversely that parent isomorphism extends over the
root. This is an exact equivalence, not a refinement signature.

Select a singleton-free state(P,S) precisely when its own key is minimum among
K(T,v) over all nine vertices. The minimum may be attained at several vertices,
but equality of keys says exactly that these vertices form one Aut(T)-orbit.
There is exactly one pointed class for that orbit and it appears once in the
parent/subset enumeration. Thus the rule selects exactly one representative
of every singleton-free unpointed child, with no duplication.

The singleton and singleton-free classes are intrinsically disjoint. Their
union proves the constructed atlas complete and pairwise nonisomorphic,
conditional only on the complete action-equipped input atlas.

## 5. Recovering the full residual action group

Let a selected state have source parent P and incoming-root subset S. For each
vertex v with K(T,v)=(P,S), enumerate exactly the parent transporters

    h in G_P satisfying h(S)=M_v.

Define a child permutation g by

    g(*)=v,       g(u)=f_v^(-1)(h(u)) for u in P.

The parent arcs are preserved because h and f_v are parent isomorphisms. Root
arcs are preserved exactly by h(S)=M_v. Thus every produced permutation is a
child automorphism. Conversely every child automorphism g has some image v of
the root and yields the unique transporter h=f_v*g on the deleted parent.
The construction therefore yields every automorphism exactly once.

For each output class the executable additionally verifies every supplied
permutation, every directed arc, the identity, distinctness, and full group
closure. These checks support the constructive proof; closure alone would not
certify group completeness without the preceding transporter argument.

The actual atlas has191480 singleton-selected classes and56 residual classes.
The sums of9!/|Aut(T)| equal2^36=68719476736, the entire labelled order-nine
universe. Its total refinement error is758, independently agreeing with the
Session50 classification, and its pointed trace identity reconstructs1724582.
The mass equality is an audit of a proved-distinct atlas, not a shortcut that
by itself proves distinctness or supplies a target count.

## 6. Refinement-only order-ten trace

For each constructed order-nine parent, enumerate its actual action orbits on
512 attachment subsets. For each resulting pointed order-ten state p, compute
the stable cell C(root_p) and accumulate

    U10 = sum_p |C(root_p)|.

The normal evaluator may stop once the root's color cell is a singleton,
because subsequent refinement cannot split a singleton. It never stops on a
nonsingleton heuristic. A separate FULL_REFINEMENT build disables this early
exit and recomputes the full partition; every per-parent output must agree.

For an unpointed child T with root orbits O, the pointed-state contribution is
sum_O |C(O)| = 10 + e(T). Thus U10=10*a10+E10. The independently proved correction
packet, byte-pinned under sources/, establishes E10=5508 by exhaustive structural
classification, not by comparing a known a10. Subtraction and division now give
a10=9733056.

The stable root-cell histogram is

    size1: 97153964; size3: 58118; size5: 1160; size7: 60; size9: 170.

The sum of multiplicities is97213472 and their weighted sum97336068.

## 7. Independent structural reconciliation by root-cell size

The total5508 correction separates into disjoint root-cell-size strata:

- Size3 contributes720. There are32 rigid strong order-ten triangle/singleton
  classes and eight order-ten classes obtained by padding the four order-nine rigid
  triple classes with a singleton on either side. These40 classes each have nine affected roots with gap2.
- Size5 contributes3360, exactly the84 rigid two-C5 classes, each with ten
  affected roots of gap4.
- Size7 contributes168. With a regular-seven cell and three singleton vertices,
  all quotient crossings are uniform. There are12 ordinary rooted four-vertex
  quotient classes, not only the four strong ones. Each bad-seven substitution
  contributes14. Among all three internal regular-seven types, there are
  1+1+3=5 root-orbit states, hence exactly12*5=60 size-seven pointed states.
- Size9 contributes1260. A nine-cell plus a singleton has uniform crossing and
  is one of two ordered SCC paddings of a regular-nine tournament. There are
  2*85=170 pointed nine-cell states; their total error is2*630=1260.

The root-orbit sums by sizes1,3,5,7,9 are therefore

    97153964; 173634; 2440; 252; 270.

They sum97330560, independently recovering10*a10. A zero error for a cell
means its cell is one orbit, not that its vertices are fixed.

For the size-seven quotient count, the unique strong four-tournament has four
rooted classes; the transitive four-tournament has four; the ordered joins
C3+1 and1+C3 each have two root orbits. Their total is12.

## 8. General weighted invariant-cell identity

For any invariantly chosen nonempty cell C of a child T, let r(p) be the actual
root-orbit size for a pointed state p whose root lies in C. Then

    sum_(such pointed states p) r(p)/|C| = 1.

This is simply the partition of C into actual vertex orbits. It supplies an
exact weighted count when root-orbit sizes are known, for example from a
proved cell-symmetry decoder. It is not by itself a method to output distinct
unpointed representatives or their complete permutation actions. The singleton
selection above has one summand and trivial root transport; residual selection
requires exact pointed keys and witnesses. Keeping weighted counting separate
from canonical representative construction prevents an illegitimate quotient.

## 9. Costs and limits

The normal atlas build took about2 seconds and the normal order-ten trace
about52 seconds in this environment. These are observed machine-dependent
runtimes, not an asymptotic or optimality claim. Input certification, artifact
replay, and independent audits have their own costs.

The input contains parent incidence and true vertex actions. The output atlas
is a full order-nine parent atlas obtained without nine-vertex canonicalization;
the target order-ten atlas is not enumerated. Future orders may have many more
singleton-free states, and structural E_n corrections remain unproved in
unrestricted higher orders. No all-order scalar recurrence follows.
