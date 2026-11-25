# Changing the apex does not escape the flip obstruction; three 13-cycle ports

This follows [the common-apex midpoint calculation](E61-APEX-MIDPOINT-AND-CYCLE-PACKING.md) and [the order-13 Singer filters](SINGER-ORDER13-MODULES.md). It gives further obstructions in specified grammars, not an unrestricted Chan–Pak or degree-57 result.

## A different transitive deletion buys no nonzero jump

Suppose \(T,T'\) differ only on \(xy\), and \(T-a\) and \(T'-b\) are transitive. **Either they have a common transitive deletion, or both have exactly five Hamiltonian paths.**

Indeed, assume no common deletion. Then \(a\ne b\), and neither can be \(x\) or \(y\), since a deletion that removes the edit works on both sides. Delete both \(a,b\). The two transitive orders must differ only by an adjacent interchange \(x,y\): an intervening point would make the flipped triple a cycle. Put \(x<y\) in the order for \(T\). Inserting \(b\) in \(T-a\) must put it between \(x,y\), since otherwise the flipped order would still be transitive after deleting \(a\). Inserting \(a\) in \(T'-b\) similarly puts it between \(y,x\). All other vertices are therefore in the same transitive prefix or suffix, with uniform arrows to or from the four sites. On the four sites both tournaments contain \(x\to b\to y\to a\to x\). Either choice of each diagonal leaves exactly two directed triangles, and a four-site tournament with two directed triangles has five paths (direct enumeration, or the odd-cycle packing formula). Ordinal singleton padding preserves the count.

**Consequence.** If two tournaments straddle \(E_p=2^{p-1}(2^p-1)\) by a single flip, then, for \(p\ge5\), at least one of them has no transitive one-vertex deletion. The prior common-apex obstruction now excludes the whole class of pairs with individually one-apex endpoints, including the branch in which the apex was supposed to change.

## An individual one-apex endpoint has an exact active order

Write the transitive remainder as an ordered chain, with 1 where the apex wins. Strip leading zeros and trailing ones as ordinal padding. Let the resulting length be \(L\). The active word begins at 1 and ends at 0. Its count is
\[
H-1=\sum_{i<j,\ s_i=1,\ s_j=0} f(j-i),\qquad f(1)=2,\quad f(d)=2^{d-1}\quad(d\ge2).
\]
The endpoint inversion alone excludes \(L\ge2p+1\) when \(H=E_p\pm1\). Conversely, swapping adjacent 01 to 10 increases the sum by monotonicity of \(f\). The sorted split at \(k\) has sum \((2^k-1)(2^{L-k}-1)+1\); it is maximized by balanced halves. At \(L=2p-1\) the maximum \(2^{2p-1}-3\cdot 2^{p-1}+2\) is still below \(E_p-2\) for \(p\ge3\). The maximum is nondecreasing in \(L\) by padding. Thus \(L=2p\), and at exponent 61 the active tournament must have **123 vertices**.

Split the chain in half. Let \(u\) encode left wins and \(v\) right losses, place the nearest-center bits at the low end, and let \(G\) count the weighted inversions in a binary integer (higher 1, lower 0). The exact inverse equation is
\[
H-1=uv+G(u)+G(v)+(u\bmod2)(v\bmod2),\qquad
G(2w+b)=G(w)+(1-b)(w+(w\bmod 2)).
\]
Both half-words have exactly \(p\) bits. This is a coupled digit recurrence, not independent residue selection.

At \(p=61\), neither endpoint is attained if the apex wins or loses at most six active contacts (converse exchanges the cases). Write the losing positions \(A\subseteq\{2,\dots,122\}\), including 122. Then
\[
H-1=\sum_{a\in A}2^{a-1}-\sum_{b<a,\ a,b\in A}f(a-b).
\]
The replay checks 8,502,671 choices with up to four additional losses. For size six, fix the first four losses and the last loss at 122, and bisect the position of the remaining loss \(a\). The variable part is \(2^{a-1}-f(122-a)-\sum_{b<a} f(a-b)\), with the sum over the four fixed smaller losses. It increases strictly: the first term grows by \(2^{a-1}\), the second cannot shrink, and the last can decrease the total by at most \(2^{a-1}\sum_b 2^{-b}<2^{a-2}\). This uses exact integer comparisons. Dense individual one-apex endpoints are not excluded.

## Which odd closing blocks can an exact exponent-61 endpoint contain?

An isolated block obtained from a transitive \(\ell\)-site chain by reversing the long arc has \(2^{\ell-2}+1\) paths: those using the reversal partition the intervening vertices between the two increasing runs. In an ordered join its factor multiplies the other counts. For odd \(\ell\) put \(a=\ell-2\). Then:

- \(E_p+1\) for odd \(p\) is \(2\bmod3\), so it allows **no** such odd block.
- \(E_p-1=(2^p+1)(2^{p-1}-1)\). The exact odd exponents \(a\) allowed by divisibility are \(a=3\), the divisors of \(p\), and the odd divisors of \((p-1)/2\). For \(a\ge5\), use a [primitive prime divisor of \(2^{2a}-1\)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/770F320F7419E0E8494A7D109AC2D9C9/S1446788715000269a.pdf/primitive_prime_divisors_and_the_mathbfn_th_cyclotomic_polynomial.pdf). Its order \(2a\) forces \(p\equiv a\bmod 2a\) or \(2a\mid p-1\). These conditions are sufficient by divisibility of the corresponding powers. The exceptional \(a=3\) always works because both factors are divisible by 3.
- At \(p=61\) the admissible block lengths on the lower side are **3, 5, 7, 17, 63**. In particular the length-13 factor \(2049=3\cdot683\) cannot occur on either side. Moreover 7 divides \(E_{61}-1\) but no \(2^a+1\), since 2 has odd order modulo 7. An ordered join consisting only of odd closing blocks cannot realize even the lower endpoint.
- More generally the two odd endpoints at any jump of 2 are coprime, so **every common, unchanged ordered-join factor is transitive**. Disjoint cyclic supports alone do not produce a jump of 2. In an interacting grammar it is the signed count of marked cycles and their residual completions that must cancel to one.

## Hoffman–Singleton: disjointness beyond the translated family

Use the 50-site affine presentation \(A_{i,x},B_{j,y}\) with \(A\)-steps ±1, \(B\)-steps ±2 and \(y=x+ij\) across, over \(\mathbb F_5\). Index \(A_{i,x}=5i+x\) and \(B_{j,y}=25+5j+y\). In [Hafner's graph](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/JACO/Volume18_1/np3862514r41w675.fulltext.pdf), the earlier phase replay ruled out a triple of translates of the same 13-cycle under its specified order-5 action. Three different induced 13-cycles do fit:

~~~text
C: 0 1 2 3 28 8 7 6 44 15 16 17 30
D: 4 29 26 11 42 40 23 22 48 46 13 37 39
E: 5 9 36 38 35 19 32 10 41 24 45 47 49
~~~

Fix \(C\). Its complement has 968 induced 13-cycles, among which 207 unordered pairs are disjoint. From [transitivity on the 126,000 cycle supports](THIRTEEN-ODD-ALPHABET.md) there are 60,984,000 unordered disjoint pairs and 8,694,000 unordered triples. Four cannot fit in 50 sites. For a second cycle disjoint from \(C\), the numbers with \(0,1,2,3,4,5,8\) third-cycle choices are \(712,176,42,18,12,4,4\).

The exact full-group census is **245 pair classes and 35 triple classes**. Rooted multiplicities are \(2^6,4^{239}\) for pairs and \(3^1,6^{34}\) for triples, where the superscript counts classes. For each member of a pack, test all 26 polygon bijections to \(C\) using discrete color refinement and verify the resulting graph map on all edges; then minimize the sorted remaining supports. This covers every possible isomorphism of packs. The listed triple is the unique class with an order-two stabilizer. The reflection fixes precisely \(1,14,18,31,33,34\) and swaps \(D,E\).

The 35 triples give 29 distinct sorted cross-edge profiles. If \(s\) is the sum of the three cross counts, the complement on 11 vertices has \(s-59\) edges by the 7-regular degree sum. Its possible edge counts in the census are \(3,5,6,7,8,9,10,11,12,13,14\). The reflection example has cross counts \((24,24,23)\), leaving 12 edges. Thus even compatible cycle lengths plus the three cross counts are not a complete signature.

## A conditional degree-57 three-orbit stencil

Consider internally cyclic order-13 orbits in a hypothetical degree-57 Moore graph with an order-13 action. These are actual induced 13-cycles. For steps \(s_i\), let \(D_{ij}\subseteq\mathbb Z/13\) be the oriented shift list, \(D_{ji}=-D_{ij}\), and \(D_{ii}=\{\pm s_i\}\). The [two-orbit filter](SINGER-ORDER13-MODULES.md) says equal steps cannot touch and unequal internally cyclic steps have at most two shifts. In each ordered pair of orbits, adjacency phases and phases of two-step routes must be pairwise distinct, except for returns to the same vertex. This is precisely where a reciprocal row encoding stores more than the three cycle lengths.

For a pair joined by \(m\) shifts, the adjacency uses \(m\) phases and routes via either endpoint orbit use \(4m\). If a third orbit touches the two by \(r,t\) shifts, it uses \(rt\) more. Hence \(13-5m-rt\) missing phases must be supplied by routes outside the three orbits. In particular **a triangle of double contacts is impossible**: \(13-10-4<0\).

There is also a complete local classification of **2–2–1** triangles. After scaling phases, take the common double endpoint to have step 1. Up to swapping leaves and changing origins, the only class has leaf steps 4 and 6 and
\[
D_{01}=\{0,6\},\qquad D_{02}=\{0,3\},\qquad D_{12}=\{5\}.
\]
The check is small: the two contacts from the center can each contain zero by choice of origins, leaving \(12^2\cdot 13\) assignments per step triple. For each ordered pair, mark the direct phases and every sum \(D_{ik}+D_{kj}\), with multiplicity, omitting the zero returns on the diagonal. Reject any repetition. The four surviving origin choices at center step 1 are \((6,3,5),(6,10,2),(7,3,11),(7,10,8)\); all are the same pattern. In the displayed normalization the missing phases for pairs 01, 02, 12 are respectively **\(\{3\}\), \(\{8\}\), \(\{2,4,6,8\}\)**. Each requires exactly one outside two-step route. This gives a concrete constraint to propagate in a larger quotient search, not a claim that this stencil embeds in a degree-57 graph.

Replay all finite computations with **node apex-and-cycle-ports-check.js**. The exponent-61 test is the slow part.
