# The complete H=15 stratum and all its transitive cloning responses

Session 197. **All-order structural classification, exact unlabelled/labelled counting formulas, and complete one-vertex transitive-cloning responses.** The full {7,21} Hamilton-path spectrum remains open. This packet is independent of the bounded-slack and factorwise-unit proof routes, except for sharing elementary counting code.

## 1. A strong mixed extension strictly increases H

Let U be a nontrivial strong tournament, and add a vertex v mixed on U. Put O={x:v->x} and I={x:x->v}; both are nonempty.

For every Hamilton path P of U, insert v immediately before its first vertex in O. This is always legal because all preceding vertices belong to I. It gives a Hamilton path not ending at v, and deleting v recovers P. In addition, every Hamilton path of U ending at a vertex of I extends by appending v. Strongness supplies at least one such path. These two sets of outputs are disjoint. Hence

    H(U+v)>=H(U)+1, and therefore H(U+v)>=H(U)+2,     (1)

since both counts are odd. This proof does not assume that deleting v from every Hamilton path leaves a Hamilton path.

Use the [S195 universal strong-growth grammar](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/strong-structure/PROOF.md): from any induced triangle, a strong tournament grows by mixed singleton extensions or cyclic two-vertex ears. The latter have

    H(C(U))=3H(U)+F2(U), F2(U)>=2^|U|-2.             (2)

## 2. Exactly four strong H15 types

**Theorem 1.** Every strong H15 tournament has order five or six. There are exactly two isomorphism types at each order, with representatives

    A=[16,9,3,5,14],
    B=[24,17,3,6,12],
    C=[36,1,2,7,15,30],
    D=[44,1,2,6,15,30].                              (3)

Their automorphism-group orders are respectively 3,5,1,1.

**All-order reduction.** Start the universal strong-growth grammar from a triangle. A pair step at current order at least four has H(U)>=5, since a strong tournament has a positive Hamilton endpoint count at each vertex and its total H is odd. By (2) the child count is at least3*5+(2^4-2)=29, impossible below H15.

The only possible pair step is therefore directly from the triangle. It gives H=3*3+6=15 at order five. Any later mixed singleton step strictly increases H by (1), and any later pair increases it as well, so in this case the ambient tournament is exactly that five-vertex graph.

If there is no pair step, the tournament has a strong singleton chain from the triangle. Exhaustive labelled order-six enumeration gives minimum H=15 among strong tournaments. An order-seven member of the chain would have H>=17 by (1). Thus no H15 graph on this route has order at least seven.

It remains only to classify strong tournaments through order six. The verifier enumerates all labelled tournaments at these orders, independently computes H, and obtains exactly the four isomorphism classes (3). There are64 labelled strong H15 tournaments on five vertices and1,440 on six. The automorphism orders are checked by all permutations; the orbit totals5!/3+5!/5=64 and6!+6!=1,440 agree. This finite calculation is made globally complete by the preceding structural reduction, not by an assumed order cutoff.

## 3. Every H15 tournament, and its exact count by order

Strong-component Hamilton counts multiply, while transitive singleton components contribute one. The nontrivial positive factors of15 can therefore occur only in either of these configurations:

1. one strong H15 component of one of the four types (3);
2. one directed triangle of count3 and one strong four-vertex tournament of count5, in either order.

All other components are singletons. There are no further possibilities: every nontrivial strong component has odd H>=3. The H3 and H5 strong types are unique, as independently checked through their elementary order bound |C|<=H(C).

The SCC order is canonical. Singleton components before, between, and after the nontrivial components are just transitive padding, and their numbers in the separate gaps are complete isomorphism invariants once the active component types and order are fixed.

**Corollary 1.1.** For every n>=5 the exact number of unlabelled n-vertex tournaments with H=15 is

    a_15(n)=2(n-4)+2(n-5)+2*binom(n-5,2)
           =(n-3)(n-4).                            (4)

There are none for n<5. The two five-vertex types have n-4 padding placements each; the two six-vertex types have n-5 each. The two ordered (H3,H5) products have three padding gaps, giving binom(n-5,2) placements each when n>=7. The formula also evaluates correctly at n=5,6, where this product term is zero.

Every automorphism fixes each SCC setwise because its linear order is canonical. Each singleton is fixed; automorphisms multiply across the active components. The product type has automorphism order3, since the triangle has order3 and the H5 type is rigid.

**Corollary 1.2.** The exact number of labelled n-vertex H15 tournaments is

    n! * (5n^2-17n-32)/15, n>=5.                    (5)

Indeed sum 1/|Aut| over the classified isomorphism types is

    (n-4)(1/3+1/5)+2(n-5)+2*binom(n-5,2)/3.

Multiplying by n! gives (5). The first values of (a_15(n), labelled count) are

    n=5: (2,64),
    n=6: (6,2208),
    n=7: (12,31584),
    n=8: (20,408576),
    n=9: (30,5322240).

The verifier independently enumerates all labelled n=5,6 tournaments, including weak ones, to check64 and2208. These formulas concern one exact H stratum inside A000568, not the total number of tournaments or a spectrum-support theorem.

## 4. Every transitive cloning family from H15

Replace a marked vertex r by a transitive m-chain, m>=1, with all exterior incidences inherited. The [S196 clone-run theorem](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/spectral-moves/PROOF.md) gives

    H(T[r:m])=sum_k a_k k! S(m,k),                  (6)

and the largest k is1+nu, where nu is the matching number of the bipartite directed-triangle incidence graph through r.

For every root in the four strong classes (3), nu<=2. Hence only k<=3 occurs, and the entire function is determined by a1,a2,a3, with a1=15. These are literal clone-run counts; the verifier computes them by permutations and separated root-symbol gaps, and independently checks cloned adjacency matrices. Formula (6), not fitting three data points to an unproved recurrence, establishes the following complete list.

**Theorem 2.** Across every H15 tournament of every order and every marked vertex, the possible functions m -> H(T[r:m]) are exactly

    15,
    9*2^m-3,
    5*2^m+5,
    6*2^m+3,
    11*2^m-7,
    5*2^m+2*3^m-1,
    7*2^m+2*3^m-5,
    4*2^m+6*3^m-11,       m>=1.                    (7)

All are physically attained. The four binary functions occur at roots1,2,3,5 of C, respectively. The three mixed-base functions occur at root1 of A, root0 of B, and root0 of C. The constant function occurs by marking an added singleton SCC in any H15 tournament.

**Completeness beyond strong graphs.** Marking a singleton SCC only lengthens a transitive component, so H stays15. Otherwise cloning acts within one nontrivial strong component and all other strong-component counts multiply unchanged. For a strong H15 component, (3) and (6) give the listed functions. In the (H3,H5) product case, the triangle's cloning function is2^m+1, and the H5 type's possible functions are3*2^m-1 and2^(m+1)+1. These follow from the same literal run counts; their roots have matching number one. Multiplying by the untouched factor gives5*2^m+5,9*2^m-3, or6*2^m+3, already on the list. No other SCC configuration exists by Section3. QED.

**Corollary 2.1.** Duplicating one vertex transitively in an arbitrary H15 tournament gives exactly one of

    {15,25,27,33,37,41,59}.                          (8)

In particular it never gives29 or31. Neither target can be obtained by any positive transitive multiplicity: m=1 gives15, m=2 gives (8), and every nonconstant function at m>=3 is at least45.

Thus H15 cannot be normalized to a transitive-cloning root giving either literal binary child2H-1 or2H+1, even as a one-off step. This independently strengthens the no-returning-unit-root result from the factorwise bridge. It does not exclude other tournament constructions for29 and31, both of which are attainable.

## 5. Verification

Run `python3 verify_h15.py > receipt-h15.json` and `python3 -O verify_h15.py > receipt-h15.optimized.json`; compare stdout. The verifier imports only the frozen local `verify.py` for elementary DP and graph utilities. It enumerates all labelled tournaments through order6, computes isomorphism orbits and automorphism orders, checks the complete SCC count formulas and explicit representatives, independently counts clone-run words, checks their matching ranks, and builds the cloned graphs. All checks use explicit exceptions.
