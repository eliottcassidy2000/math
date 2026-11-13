# Spatial filling, binary trees and Hamiltonian count fibres

**PROVED IN NOTES / CITED inheritance / FINITE-EXACT controls.**
Chan–Pak remains OPEN. The missing odd tournament counts tending to
exactly{7,21} is the conjecture itself, not an established input.
This session gives an exact spatial filling law inside existing fibres,
tests its possible extension to a construction tree, and recovers the
correct matrix and Collatz graph formulations.

## Three different meanings of perfect fill

| Space | Exact result | What remains missing |
|---|---|---|
| Vertex subsets of one tournament | Every cell carries a positive odd induced path count; the full array has an exact square-zero factorization | Realizing a prescribed numerical count |
| All neutral-growth representatives of a fixed strong-component word | An integer orthant, with exact uniform fill of each size simplex | Creating a previously absent core word or count |
| Binary addresses of odd target integers | A rooted tree of exactly the allowed labels after splicing out7 and21 | Actual tournament transitions that compose along every edge |

These maps preserve different information. A proof cannot replace a
missing physical transition with the fact that its scalar endpoints
are both present.

## Exact multiplicities and neutral spatial growth

The recovered complete inventories and independent enumeration give,
for ordinary tournament isomorphism classes on n>=1 vertices,

    N_n(5)=max(n−3,0),
    N_n(17)=2 max(n−5,0),
    N_n(7)=N_n(21)=0.

H5 has one strong atom of order4. H17 has two strong atoms of order6.
Every other component is a singleton. These formulas count actual
ordinary classes, not construction parameters or converse-merged types.

More generally fix an ordered word W of r nontrivial strong components.
The r+1 singleton gaps give coordinates g in N^(r+1). Every H-preserving
birth is one coordinate increment. There are g_i+1 literal attachment
cuts producing the quotient edge g→g+e_i. Consequently uniform choice
among literal neutral attachments gives

    Pr(g→g+e_i)=(g_i+1)/(d+r+1),  d=sum g_i.

After d births every cell of the size-d integer simplex has probability
1/binom(d+r,r), with exactly d! named construction histories per cell.
This is a perfect fill maintained at every stage, with the classical
Pólya weights implemented by actual tournament attachments.

The component count r is essential: its exact backward-difference law
recovers the number of primitive core words. Forgetting r makes a naive
birth counter negative already at H9. See
[the proof, controls and all-order coverage criterion](multiplicity/PROOF.md).

Let b_n be the least allowed odd value missing at order n. Padding makes
b_n nondecreasing. Chan–Pak is equivalent to b_n→infinity. The complete
small census gives5,9,17,35,63,611 at orders3 through8; those finitely many
values do not establish the limit. Neutral growth fills existing fibres
and cannot make b_n increase by creating a new count.

Even allowing ordinal products of a fixed finite core bank cannot close
this gap: all resulting counts have prime factors in one finite set,
so only O((log X)^s) distinct counts occur below X, where s is the number
of those primes. This excludes only fixed cores with products and
neutral padding. A finite alphabet of operations that creates new strong
cores is a different, potentially richer construction.

## Testing an arithmetic spanning tree

In H=2m+1 coordinates, binary children m→2m,2m+1 become2H−1,2H+1.
Delete forbidden7 and21 and splice their children to their grandparents:

    3→13,15;       11→41,43.

This is a proved arithmetic spanning tree on all allowed odd labels.
The [construction criterion](TREE_BLUEPRINT.md) specifies exactly what
physical states and transitions would turn it into a proof.

However its pure-birth lift is impossible. The complete induced-count
calculation at targets through73 finds precisely five failed edges:

    13→27, 25→51, 29→57, 31→63, 35→69.

The prime-source failures recover the earlier November2 theorem; the
full union calculation includes composite sources. These exclusions
hold at every order, without any unit-ready condition. For example,

    U(27)={1,3,5,9,15,17,27}

is the complete set of possible induced counts inside any H27 target.
No H13 source can grow into one. A cubic number of H27 representatives
does not overcome this absence of embeddings.

Four branches have physical multiplicative replacements:9→27,17→51,
19→57,23→69, by adjoining an ordered triangle component. The63 target
needs its separate strong carry atom because7 and21 are absent. This
supports a construction graph with several kinds of transitions:
ordinal products, binary growth where it is executable, and exact-target
core repairs. It does not prove that these operations cover all labels.
See [the complete reduction and independent audit](INDUCED_BOUNDARY.md).

## Matrices that remember Hamiltonian paths

Ordinary A^k counts walks of k edges, including repeated vertices.
Hamilton paths use n−1 edges; Hamilton cycles use n. A transitive
tournament has H1 and A^n=0, so the proposed unmarked power test fails.
The oddness theorem concerns the Hamilton count, not the sum or number
of nonzero entries of A^n.

Put X=diag(x_1,...,x_n), with commuting x_i and x_i^2=0. Then

    H(D)=[x_1...x_n] 1^T X(AX)^(n−1) 1.

Repeated visits vanish, while the endpoint matrix retains each starting
and ending count. A weighted matrix-tree cofactor gives the same endpoint
counts after squarefree extraction: branching trees vanish. Thus trees
and paths have a precise relationship, without identifying their totals.
The corresponding forest determinant recovers all unmarked path-cover
bin counts, including P2 and P3. Root-specific cover counts still need
their endpoint marks.

For tournaments the complete induced-count array is

    sum_S H(T[S])x_S = product_v(1+x_v) product_odd_cycles_C(1+2x_C).

It fills every subset cell modulo2. Yet an abstract positive odd array
need not come from a tournament: three triangles sharing an edge would
formally give H7, but completing the orientations forces a five-cycle
and raises the smallest completion to9. With m transitive pages this
same construction gives H=2^m+1, including5 and17.

General oriented graphs can attain7 and21: directed cycles of those
lengths do so. Even deleting one arc from a tournament yields small
examples with counts7 and21. The repair costs retain endpoint data;
restoring tournament completeness does not interpolate through every
integer. See [the matrix and tree identities](matrix/PROOF.md).

## Collatz as a rooted tree, and the literal relaxed rules

For ordinary Collatz, deleting the root cycle edge1→4 gives a rooted
tree on the integers known to reach1. Showing that this tree spans every
positive integer is exactly the remaining conjecture. Making every
original edge reversible is equivalent to ordinary Collatz connectivity.

The user's literal extra rules are different: even n may also double,
and odd n≡1 mod6 may go to(n−1)/3. The true reverse Collatz branch applies
instead at even n≡4 mod6. In the literal relaxation, added even doubling
changes neither reachability nor the minimum required peak: dyadic
detours can be replaced by halvings. The extra odd inverse is substantive.

There is an all-order geometric boundary. For even K>=2, every successful
relaxed path from2^K−1 must reach at least2·3^K−2. No fixed constant-factor
integer interval can therefore maintain root connectivity as size grows.
The finite policy check reaches1 from every start through100000, and
saved paths attain the exact lower peak for even K through40. Universal
reachability of this literal relaxation is not established here.

The requested four-number connection is exact under3n−1:

    5→7→5,
    17→25→37→55→41→61→91→17,
    21→31→23→17,

where arrows suppress even halvings. Thus5 and17 are cycle minima,
7 lies on the first cycle, and21 feeds the second. No known transfer
turns this basin information into tournament non-realizability.
For example25 lies on the second cycle and is attained by the ordinal
product of two H5 tournaments; cycle membership alone cannot explain7.
See [the guarded proofs and finite certificates](collatz/NOTE.md).

## Validation and provenance

[AUDIT.md](AUDIT.md) records the proof reviews and reproduction commands.
The numerical universes include all33,867 labelled tournaments through6,
all4,165 loopless labelled digraphs through4, exact neutral histories,
finite native/accelerated Collatz comparisons, and a complete inherited
strong bank with23,512 independently recounted subsets. The isomorphism
bank is inherited and pinned, not silently claimed as regenerated.

Relevant primary literature was checked: [Chan–Pak, Conjecture6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
[Busch's strong-order bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3),
[Irving–Omar's marked matrix treatment](https://arxiv.org/html/2412.10572v3#S2),
and [Lagarias's Collatz overview](https://arxiv.org/html/2111.02635v1).
The repository's prior proofs and corrections are linked beside the
claims that use them. No global metagraph, Collatz, Landau or Chan–Pak
classification is asserted.
