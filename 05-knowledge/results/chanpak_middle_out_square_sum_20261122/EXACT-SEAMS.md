# Middle-out perfect seams: exponent 11 and a global obstruction to the first induced step

Write \(E_p=2^{p-1}(2^p-1)\). A seam at exponent \(p\) is a pair
of tournaments with counts \(E_p-1,E_p+1\), differing on exactly
one unordered pair. An induced step means a common induced
embedding of the two orientations, retaining the flipped pair.
We allow an embedding to exchange the lower and upper orientations.

## Strong literal exponent-11 seam and its nested family

Take the nine-vertex lower deletion with outgoing bit masks
\[
 (256,265,267,337,295,143,55,95,224)
\]
and tag mask 84. Reverse only the arc between vertices 3 and 8
on the upper side. Attach an apex beating exactly the tagged
vertices; place an ordinal right suffix of \(m\) transitive vertices
after the deletion, all untagged with respect to the apex.
The core profiles \((h,c,z,k)\) are respectively
\[
 (1033,8190,-220,7457),\qquad(1079,8190,-244,7435).
\]
Here \(h\) is the core path count, \(c,z\) its unsigned and signed
weighted tagged-start reads, and \(k\) the count with the apex
attached but without a suffix. The actual suffix is the reverse
complement of the binary expansion of \(y=2^m-1\), namely \(m\)
zeroes. Hence the common binary read \(P(y)\) is zero, and \(y\)
is odd. In the exact transfer \(k+cy+hP(y)-z[y\text{ odd}]\),
the changed \(h\) is invisible while \(k-z\) changes by two.
The whole family therefore has exact counts
\[
 8190\cdot2^m-513,\qquad8190\cdot2^m-511 \quad (m\ge1).
\]
Both tournaments are strong, already at \(m=1\), and every
further suffix vertex reaches the apex and is reached from the
core. Keeping the core, apex and old suffix embeds each pair
jointly in its successor. This is actual induced growth of mixed
contacts, not just a rewiring or a congruence. The new vertex is
beaten by the deletion and beats the apex, so it is neither
a source nor a sink of the full tournament.

At \(m=8\) this gives the exact strong 18-vertex seam
\(2,096,127,2,096,129\) around \(E_{11}\). The replay counts each
full tournament independently by literal subset-endpoint DP.
The certificate does not claim order minimality. Together with
the earlier certificates this fills the consecutive exponent band
3 through 11.

The nested family is not an induced chain at consecutive Euclid
centers. For \(1\le m<8\), the midpoint has valuation \(m+1\),
so if it were \(E_p\) we would have \(p=m+2\) and then
\(4096=2^{m+2}+2^{8-m}\), which is impossible in that range.
For \(m>8\), the valuation is 9, forcing \(p=10\), but the
midpoint then exceeds \(E_{10}\). Thus \(m=8\) is its only
Euclid seam. The mechanism differs from equal-\(h\) positive
transmission: it preserves the gap on a distinguished address
rail despite changing the ordinary deletion count.

A further exact seam for the Euclid-form integer \(E_{21}\) uses
lower deletion \((240,141,153,177,194,150,174,0)\), tags 209,
flip \((0,4)\), no left flank, and right address \(1331008449\).
The profiles \((h,c,z,k)\) are \((159,1575,-159,567)\) and
\((159,1575,-159,569)\). The replay constructs strong order-40
pairs, checks that only the indicated arc differs and gives
exactly \(2,199,022,206,975\) and \(2,199,022,206,977\).
This exponent is not claimed to be a Mersenne prime.

## There is no common induced step from exponent 3 to 4, at any order

This is a global obstruction, not a negative search only up to
some externally imposed size. The [Busch lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
gives \(H(T)\ge5^{(n-1)/3}\) for strong tournaments. Thus a strong
27-count source has order at most seven, and a strong 119-count
target has order at most nine.

Recall the cross-component unit-seam restriction: if the lower
count is greater than three, both orientations have one common
nontrivial strong component containing the flipped arc, and all
other components are singletons. For completeness, the
[odd-cycle packing formula](https://arxiv.org/abs/2307.05569)
makes a reversal of a lower cross-component arc nondecreasing;
a gap of two forces precisely one new odd cycle with transitive
complement. The new triangles corresponding to \(a\to b\)
between \(A,B\) are counted by
\(\deg_A^+(a)+|M|+\deg_B^-(b)\), where \(M\) is the set of intermediate
vertices. If zero, no new cycle can exist. If one and its middle
vertex belongs to \(M\), the transitive complement and the sink
and source endpoints force the lower graph transitive. Otherwise
up to converse it comes from \(a\)'s unique outneighbor \(x\) in
\(A\), and \(b\) is a source in \(B\). The complement forces
\(R=A-\{a,x\}\) and \(B\) transitive. Any conquest by \(x\) of a
nonfinal vertex of \(R\) creates a second new odd cycle, the
five-cycle via that vertex and its successor. So \(x\) can beat
only the last vertex of \(R\), making the lower count at most
three. For an internal flip all other components contribute
odd factors dividing the difference two, hence are singletons.
The changed component cannot split, since the inverse would
then be a cross-component flip decreasing the count.

A strongly connected induced subgraph of a tournament lies in
one of its strongly connected components. Thus in a hypothetical
common induced step we may delete all singleton components of
the 27/29 source and of the 119/121 target. The same source
strong core, containing the flipped pair, lies in the target's
strong core, even if the two orientations are exchanged.
This proves the size cutoffs suffice globally.

The companion C++ replay independently enumerates all source
classes and possible target extensions. The source enumeration
uses Rédei's theorem: fix one spanning directed path and enumerate
all \(2^{(n-1)(n-2)/2}\) completions for each \(3\le n\le7\).
Filter exact counts 27 and 29 on both sides, and canonicalize the
marked lower tournament by relabelling along each of its directed
Hamiltonian paths and taking the least row/mark tuple. This is
complete because every isomorphism maps those paths bijectively.
The replay finds two marked classes at order six and two at order
seven, and none below six. For every class, in **both** assignments
of source orientations, it tests all common one- through
three-vertex extensions within the order-nine target bound
(with at most two new vertices for order-seven sources).
All 20 canonical source/size/orientation trials miss. No
semantics of induced lift that merely permits an orientation
swap can escape the conclusion.

## Joint weighted-cut search

For a fixed lower or upper substrate \(D\) and one new apex
beating \(S\subseteq V(D)\), define
\[
 w_{ij}=2\sum_{\substack{U\subseteq V(D)\\ |U|\ {\rm even}}}
 H(D-U)\,P_{ij}(D[U]),
\]
where \(P_{ij}\) counts spanning paths of the indicated graph
from \(i\) to \(j\) and \(H(\varnothing)=1\). The same cycle-packing
formula gives \(H(D+\mathrm{apex})=H(D)+
\sum_{i\in S,j\notin S}w_{ij}\). Any cycle using the new apex has
one even-length vertex support \(U\) left after removing it,
and factors uniquely from the disjoint cycles on \(D-U\).
Thus the two weighted cuts test the **same** contact mask on
both sides, not two separately attainable values.

The search uses ordered contact masks on the original core to
quotient by relabelling of newcomers. Every contact pattern to
earlier new vertices is still tested; the relative edge between
two newcomers is not fixed in advance. A partial branch is
pruned if either count has already exceeded its prescribed
target: every odd-cycle packing on an induced subgraph is also
a packing on its supergraph. At the final layer, the lower
weighted cut gives all counts at once, and the few lower
matches are checked by independent literal DP on the other side.
This is an exact bounded procedure, with neither stochastic
pruning nor a scalar approximation.

The stronger all-exponent existence of individual seams remains
open. A chain that requires a common induced step at *every*
adjacent exponent starting at 3 is impossible.

The exponent-11 family updates the literal center by \(M\mapsto2M+512\) per new suffix vertex and hence \(M\mapsto4M+1536\) per two. Its offset is an actual, reusable physical read, but it is constant rather than the changing \(2^p\) required by \(E_{p+1}=4E_p+2^p\). This distinguishes an exact gap-preserving growth rule from a rule that also tracks the prescribed exponent.

**Replay:** `python3 -B check_seams.py` and `c++ -std=c++17 -O3 -Wall -Wextra -Werror check_global.cpp -o /tmp/check_global; /tmp/check_global`. The latter takes no inputs and matches [the complete recorded output](global_expected.out).
