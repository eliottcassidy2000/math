# Maximum-span majority selection and the mixed degree-three neighborhood

Date: 2026-10-19. Status: proved structural exclusion for a hypothetical poset with no balanced pair. This does not prove ALL7thin or assert that every rooted degree-three instance has a balanced pair inside its five-point neighborhood.

## Main conclusion

Let a finite non-chain poset have no pair whose comparison probability lies in [1/3,2/3]. Its >2/3 majority order is a linear extension. Choose an incomparable pair x,y of maximum index distance in this order, with x before y.

1. The ordered pair (x,y) is critical: D(x) is contained in D(y), and U(y) is contained in U(x).
2. Every common incomparability neighbor z lies strictly between x and y in the majority order.
3. It is impossible that |I(x)|=|I(y)|=3 and |I(x) intersect I(y)|=1.
4. More generally, if both selected root degrees are at most three, their common incomparability neighborhood must be empty.

The third assertion follows from a new elementary inequality for all seven possible mixed degree-three cores:

    p(x,z) + p(z,y) <= 4/3.

It holds with arbitrary finite exterior and no degree bound. The proof requires only independent prefix/suffix counts that are at least one. Rank-log-concavity and ideal-run moment bounds are not needed for this particular inequality.

## 1. Majority order and the safe selection lemma

For any three distinct points a,b,c, every linear extension satisfies

    1[a<b] + 1[b<c] + 1[c<a] <= 2.

Thus their three comparison probabilities cannot all exceed 2/3. If no pair is balanced, every pair has exactly one direction with probability >2/3. The resulting tournament has no directed triangle, hence is transitive. Denote its order by <_M. It extends P because a poset relation has probability one.

Choose an incomparable pair x<_M y maximizing index_M(y)-index_M(x).

If d<x but d is not below y, then d is incomparable with y: y<d<x would make y<x. But d<_M x<_M y, so (d,y) is an incomparable pair of larger span, a contradiction. Hence D(x) is contained in D(y). Dually U(y) is contained in U(x).

If z is incomparable with both roots and lies before x in M, the pair (z,y) has larger span. If z lies after y, the pair (x,z) has larger span. Therefore

    x <_M z <_M y                                         (1)

for every common neighbor z. In particular both p(x,z) and p(z,y) exceed 2/3.

This argument is structural: the critical-pair and common-neighbor conclusions hold for a maximum-span incomparable pair in any linear extension, even when that extension is not a majority order.

## 2. Why adjacency cannot replace maximum span without another argument

A strict transitive majority order need not contain an adjacent forward critical pair. Take the disjoint union of chains

    A0<A1<A2,       B0<B1.

There are ten linear extensions. The strict majority order is

    A0, B0, A1, B1, A2.

The ten forward comparison counts in this order are, row by row:

    A0: B0 6, A1 10, B1 9, A2 10
    B0: A1 7, B1 10, A2 9
    A1: B1 7, A2 10
    B1: A2 6.

None of the four adjacent pairs is forward critical. For example (B0,A1) fails U(A1) subset U(B0), and (A1,B1) fails D(A1) subset D(B1). Both middle adjacent comparisons have probability 7/10>2/3. The forward critical pairs are (A0,B1) and (B0,A2), neither adjacent.

This is not an example with no balanced pair: its two outer adjacent probabilities are 3/5. It disproves an adjacency-selection claim based only on transitivity of the actual strict majority order or on an individual adjacent >2/3 comparison. The stronger conditional assertion under the global no-balanced hypothesis is not settled by this example. Maximum span supplies the needed rigor without that assertion.

## 3. The seven rooted shapes

Suppose the maximum-span roots satisfy d(x)=d(y)=3 and have exactly one common neighbor z. Write

    I(x)={y,z,b},       I(y)={x,z,a}.

Criticality forces b<y and x<a. Indeed y<b would imply x<b by U(y) subset U(x), and a<x would imply a<y by D(x) subset D(y); both contradict the definitions of a,b.

Among b,z,a, the only possible relations are

    b<a,       b<z,       z<a.

The opposite directions respectively force x<y, z<y, or x<z. If b<z and z<a then b<a is compulsory. Consequently the exact list of relation-bit triples

    (b<a, b<z, z<a)

is

    000, 010, 001, 100, 110, 101, 111.                     (2)

The last is the previously settled Q5 core. The other six are not presumed eliminated by that theorem; the following direct inequality handles all seven together.

## 4. Direct exterior factorization

Assume only that I(x) union I(y) is the five-point set {b,x,z,y,a} and its induced core is one of (2). All exterior points are comparable with both roots. A mixed direction would compare x and y, so the exterior partitions into A below both roots and B above both roots. Moreover A<B pointwise, A<a, and b<B.

No core point is below an element of A: this would put it below both roots, inconsistent with membership in their incomparability union. Dually no core point is above an element of B. Thus the only core points that can interleave with A are b,z; the only ones that can interleave with B are z,a.

Fix a legal core word o. Let pre(o) be the prefix before its first root and suf(o) the suffix after its last root. Every A point precedes the first root; every B point follows the last root. Therefore the number of full extensions inducing o is exactly

    L(pre(o)) R(suf(o)),                                 (3)

where L is the extension count of A together with its ordered core prefix, and R is the extension count of its ordered core suffix together with B. The intervening core letters have a unique fixed order. Relations between these two factors and the intervening letters are automatically satisfied. Conversely their concatenation gives every extension inducing o exactly once.

Normalize L by e(A) and R by e(B). Empty prefix/suffix counts are one. Define the other normalized features by

    prefix: b -> c, z -> e, bz -> u, zb -> w;
    suffix: a -> d, z -> f, za -> v, az -> t.

Every used feature is at least one: concatenate any extension of A before the specified prefix, or the specified suffix before any extension of B. An order incompatible with the core never occurs and its unused feature can be set to one.

Let E0=e(P)/(e(A)e(B)), and let Nxz,Nzy be the similarly normalized comparison counts. They are exactly the sums of the weights in (3), with comparison indicators. There are respectively 30,20,20,25,18,18,11 legal core words for the shapes in (2).

## 5. Seven short nonnegative certificates

Set the nonnegative increments

    C=c-1, Z=e-1, U=u-1, W=w-1,
    D=d-1, F=f-1, V=v-1, T=t-1.

Define

    K = 2CD + 2CV + 2DU + 2C + 2D + 3U + 3V,
    X_r = CF + 2CT + C + rF + 3T,
    Y_r = DZ + 2DW + D + rZ + 3W.

Core-word summation gives the following exact identities for

    Delta = 4E0 - 3(Nxz+Nzy).

    shape 000: Delta = K + X_3 + Y_3
    shape 010: Delta = K + X_3 + 2
    shape 001: Delta = K + Y_3 + 2
    shape 100: Delta = K + X_2 + Y_2 + 4
    shape 110: Delta = K + X_2 + 3
    shape 101: Delta = K + Y_2 + 3
    shape 111: Delta = K + 2.

Every expression is nonnegative. This proves p(x,z)+p(z,y)<=4/3, independently of the size, degrees, or internal structure of A and B. The companion verifier derives these identities by standard-library polynomial arithmetic, enumerates legal core words, and checks the factorization against independent whole-poset extension counting on random exteriors.

Combining this inequality with (1) gives a strict contradiction. Thus no maximum-span pair in a poset with no balanced pair has parameters (3,3,1).

## 6. The one-sided-chain screen does give a balanced pair

The following elementary injection supplies the balance conclusion needed here.

Suppose x is incomparable with y, p(x,y)>=1/2, D(y) is contained in D(x), and U(x) minus U(y) is the chain s1<...<sr. Put s0=x. Then every si is incomparable with y. The last pair has no slaves, so the root-swap injection gives p(sr,y)<=1/2.

For i<r, in an extension with si<y<si+1, every point between si and y is incomparable with both: lower slaves do not exist, and any upper slave is at least si+1. Swapping si and y is injective, hence

    p(si,y)-p(si+1,y) <= 1-p(si,y),
    2p(si,y)-p(si+1,y) <= 1.                              (4)

If no pair (si,y) were balanced, the sequence starts above 2/3 and ends below 1/3; at its first such downward crossing (4) is violated. Thus one of these pairs is balanced. The dual proof handles lower-chain slaves. The empty and singleton slave cases follow immediately.

## 7. Small-degree consequence and remaining scope

For a critical pair, all exclusive neighbors belong to its forward slave set. Their two side sizes are

    |S^-(x,y)| = d(x)-c-1,
    |S^+(x,y)| = d(y)-c-1,

where c=|I(x) intersect I(y)|. If both degrees are at most three and c>=1, these sizes are each at most one. If their sum is at most one, Section 6 gives a balanced pair. If their sum is two, both degrees are three and c=1, which Sections 3-5 exclude for maximum-span selection. This proves the fourth main conclusion.

This is an actual structural reduction for any thinness bound, including seven: a maximum-span selected pair with both root degrees at most three can survive only with no common neighbor. It does not force all maximum-span roots to have small degree, does not dispose of the common-neighbor-free branch, and does not reconstruct Peczarski's six-thin configuration table.

From the repository root run `python 05-knowledge/experiments/verify_majority_maxspan_root3_20261019.py`. All certificate guards remain active under `python -O`; the program writes its receipt to standard output. The receipt is [majority_maxspan_root3_audit_20261019.json](majority_maxspan_root3_audit_20261019.json).
