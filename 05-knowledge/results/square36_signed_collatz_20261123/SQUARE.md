# What a square label remembers on a Hamilton walk

Write \(Q_n\) for the simple graph on \(1,\ldots,n\) with an edge \(xy\) when \(x+y=k^2\). The earlier protected-path construction transports the edge \(19\!-!30\) under scaling. The order-36 test below asks when one can retain a prescribed edge in a **closed** walk visiting each label exactly once. The arithmetic observables are the square on an edge, the signs of its endpoint units at a prime, and the two alternating matchings of an even cycle. [Protected bridge and its lift](../chanpak_wheel_square_selector_20261117/square/PROOF.md) · [Square-path incidence and overlap](../session141_metric_aliquot_square_paths_20261102/square_paths/PROOF.md)

## The exact threshold for additive square-label observables

**PROVED.** Let \(A\) be any abelian group. For \(n\ge29\), suppose that a function \(g:[n]\to A\) and edge-label prices \(h(k)\in A\) satisfy
\[
                    g(x)+g(y)=h(k)\qquad(x<y,\ x+y=k^2).
\]
Then \(g(x)=a+xb\), and \(h(k)=2a+k^2b\) for every used root. These formulas also suffice. At \(n=28\) there is an extra independent parameter, already visible over the integers. Thus 29 is the least eventual threshold, and reducing the prices modulo any integer introduces no additional solutions above it.

**General sum sets.** The same problem on \([n]\) with an arbitrary permitted sum set \(S\) has a direct solution. Make a graph on the difference indices \(1,\ldots,n-1\) by joining \(i,j\) when \(|i-j|\ge2\) and \(i+j+1\in S\). Choose \(g(1)\in A\) and, independently, one element of \(A\) on each connected component. Use that element for every first difference \(g(i+1)-g(i)\) in the component. These are exactly all solutions.

Indeed the two edges \((i,j+1),(i+1,j)\) have the same sum, so their prices agree precisely when the corresponding first differences agree. Conversely, for a fixed sum \(K\) the smaller endpoints of all allowed pairs form the integer interval from \(\max(1,K-n)\) to \(\lfloor(K-1)/2\rfloor\). Comparing consecutive entries proves the common price throughout. This argument allows arbitrary torsion in \(A\).

For squares at \(n=28\) the difference components are
\[
 \{1,3,5,7,8,10,12,14,16,17,19,21,23,25,27\},\quad
 \{2,4,6,9,11,13,15,18,20,22,24,26\}.
\]
The new difference index 28 joins both 7 and 20, using 36 and 49. Each later index \(i\ge29\) joins an earlier index \(j\le i-2\). Take the square of \(\lceil\sqrt{i+2}\rceil\), which lies in \([i+2,2i-1]\), and put \(j=k^2-i-1\). The upper bound follows from \((\sqrt{i+2}+1)^2\le2i-1\) for \(i\ge29\). Connectedness now persists.

The non-affine choice through 28 has a short closed form:
\[
                 g(x)=\left\lfloor\frac{9x-3}{20}\right\rfloor,
                 \qquad h(k)=\binom{k}{2},\quad 2\le k\le7.
\]
Its differences are zero on the first displayed component and one on the second. A single representative edge per root gives the displayed prices. At vertex 29 the neighbor 7 would force \(g(29)=15-g(7)=12\), while neighbor 20 would force \(g(29)=21-g(20)=13\).

**Hamilton consequence.** For any spanning path with endpoints \(s,t\), the sum of the prices of its edges is \(2\sum_{v=1}^n g(v)-g(s)-g(t)\). Beginning at 29, all additive identities obtained by charging the two endpoints of each square edge depend only on the path length and \(s+t\). They cannot distinguish two endpoint pairs of the same sum. This quantifier concerns identities valid on **each edge**. Restrictions on the whole path and data attached to incidences can distinguish additional cases.

## Odd valuation layers distinguish split and inert primes

**PROVED.** Fix an odd prime \(p\) and odd positive \(j\). Partition the labels of valuation \(j\) into \(L_+,L_-\) according to the quadratic character of \(x/p^j\) modulo \(p\). Call an edge from this layer to another valuation an exit. If valuations differ on an edge with sum \(k^2\), the smaller valuation is \(2v_p(k)\), and the corresponding unit is a square modulo \(p\). Consequently every exit from this odd layer leads to a lower even valuation with square unit. If both valuations equal \(j\), then \(p^j\mid k^2\) implies \(p^{j+1}\mid k^2\). The two units sum to zero modulo \(p\). Their signs differ for \(p\equiv3\pmod4\) and agree for \(p\equiv1\pmod4\).

Write \(E_+,E_-\) for the exit counts in a selected subgraph and \(D_+,D_-\) for its total degrees over the respective classes. Counting internal incidences gives
\[
 \begin{array}{ll}
  E_+-E_-=D_+-D_- & (p\equiv3\pmod4),\\
  E_+\equiv D_+\pmod2,\quad E_-\equiv D_-\pmod2 &(p\equiv1\pmod4).
 \end{array}
\]
This works for any selected edges. A perfect matching gives \(D_\pm=|L_\pm|\). A spanning path gives \(D_\pm=2|L_\pm|-e_\pm\), where \(e_\pm\) counts its endpoints in the class. Thus the signed inert exit imbalance measures an endpoint charge. It is an incidence observable that the preceding rigidity theorem leaves available.

On \([36]\), for \(j=1\), the two class sizes at 3, 5, 7 are respectively \((4,4),(3,3),(3,2)\). Each perfect matching has equal exit counts at 3, odd exit counts on each side at 5, and exactly one excess positive exit at 7. At 7 the positive labels are 7, 14, 28 and the negative labels 21, 35. The sign law describes actual constraints even in cases where the same unsigned root color can be used by many matchings.

## The unique edge lost by alternating Hamilton closure at 36

**FINITE-EXACT.** Of the 57 edges of \(Q_{36}\), precisely \(6\!-!30\), with square label 36, belongs to no Hamilton cycle. It lies in 15 of the 295 perfect matchings, and in nine of the 70 spanning two-factors (subgraphs of degree two everywhere). In every one of those nine the component through this edge is odd. Its cycle, written starting at 6, is one of the following. Each occurs in three two-factors.

| Length | Cycle, with the final return to 6 understood |
|---:|---|
|3|6, 19, 30|
|13|6, 19, 17, 32, 4, 21, 15, 10, 26, 23, 2, 34, 30|
|13|6, 10, 26, 23, 2, 34, 15, 21, 4, 32, 17, 19, 30|

Exactly 59 two-factors consist solely of even cycles. They are the unions of the 95 unordered, edge-disjoint pairs of perfect matchings, with multiple decompositions when there are several components. Exactly 31 are single cycles. Thus a perfect matching can support the edge and a spanning two-factor can support it, yet the edge cannot appear in either member of a disjoint pair of perfect matchings. The adjacent edge \(19\!-!30\) occurs in all 31 Hamilton cycles.

The only edge between square vertices in \([36]\) is \(9\!-!16\), coming from \(3^2+4^2=5^2\). It serves 21 of the 31 Hamilton cycles. Vertex 27 has just two neighbors, 9 and 22. This forces \(9\!-!27\), so the remaining ten cycles all use \(7\!-!9\) in place of the Pythagorean edge. For an odd prime \(p\), the only possible square neighbor of \(p\) is \(((p-1)/2)^2\), by factoring \(k^2-t^2=p\). The prime–square ports on this horizon are therefore \(3\!-!1,5\!-!4,7\!-!9,11\!-!25,13\!-!36\). Inert primes here meet odd squares, split primes even squares.

A square edge acts on the current label by the reflection \(x\mapsto k^2-x\). For an odd closed walk with successive square labels \(K_1,\ldots,K_{2m+1}\), the return map is \(x\mapsto K_1-K_2+\cdots+K_{2m+1}-x\). For each of the three forbidden-edge components above the signed sum is 12, fixing the start 6. The triangle already gives \(25-49+36=12\). On an even cycle the alternating sum is zero. A zero translation tests closure of the labelled walk. The disjoint matching condition also enforces distinctness and evenness across all components.

The [certificate](square36.py) constructs the interval graphs, computes their components, enumerates regular factors by edge inclusion and exclusion with forced degrees, and checks all disjoint matching pairs against the independent two-factor enumeration. Run `python3 -B square36.py` in this directory with assertions enabled.
