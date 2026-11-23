# The 72-height cycle exists pairwise and fails at one tournament fiber

Let \(L=\{100v-1,100v+1:1\le v\le36\}\). Join distinct \(a,b\in L\) when \(a+b\) is a square and some two tournaments on the same vertex set with ordinary Hamiltonian-path counts \(a,b\) differ in exactly one arc.

**FINITE-EXACT, all tournament orders.** Exactly 101 of the 118 square-sum pairs are edges. This graph has a **unique undirected Hamiltonian cycle**. Order 13 is the least common tournament order at which the cycle is available edge by edge. There is **no coherent lift** through one representative per height, even allowing arbitrary isomorphisms, count-preserving flips and singleton padding at each height. The obstruction occurs at 199.

Here is the cycle, with the last entry joined to the first.

~~~text
  99  301  599 1001 2599 2301 1299 3601 2799 2101 1499 3401
2999 1901 1699 3201  399 1201 2399 2501 1099  501 1999 2901
3499 1401 2199 2701  899 1601 3299 3101 1799  701  199  201
 699 1801 3099 3301 1599  901 2699 2201 1399 3501 2899 2001
 499 1101 2499 2401 1199  401 3199 1701 1899 3001 3399 1501
2099 2801 3599 1301 2299 2601  999  601  299  101  799  801
~~~

The 17 excluded pairs, grouped by lower endpoint, are as follows. Every other arithmetic pair has an explicit outgoing-row witness. [101 one-flip certificates](physical_edges.rows)

| Lower count | Upper counts excluded at every order |
|---:|---|
|99|1,501, 2,401, 3,501|
|101|1,499, 2,399, 3,499|
|199|1,401, 2,301, 3,401|
|201|1,399, 2,299, 3,399|
|299|2,201|
|301|2,199, 3,299|
|401|2,099|
|499|3,101|

## An all-orders inverse for an edited arc

**PROVED.** For any positive odd \(a<b\), deciding the existence and minimum order of a one-arc reversal from \(a\) to \(b\) reduces to a bounded list of strong tournaments of height at most \(a\), ordered factorizations of \(a\), and a power-of-two divisibility test. The following construction makes the component masses needed for the inverse explicit. It uses the ordinal factorization of ordinary path counts and the preceding common-cut/strict-merge result. [Arc-fiber grammar](../chanpak_arc_fiber_grammar_20261122/README.md)

Write \(H(A)\) for the number of directed Hamiltonian paths of a tournament \(A\), with \(H(\varnothing)=1\), and set \(Q(A)=\sum_{S\subseteq A}H(S)H(A-S)\). Let \(s_A(u),e_A(u)\) be the full-path starting and ending counts at \(u\). Define \(S_A(u)\) by imposing that the first path in the two-slot count \(Q(A)\) starts at \(u\), and \(E_A(u)\) by requiring it to end there. Both \(H\) and \(Q\) multiply across ordinal sums.

Decompose a tournament of height \(a\) into its ordered strongly connected components (SCCs).

* If both arc endpoints lie in one component of height \(h\), the other components contribute the exterior factor \(g=a/h\). Test the internal change by reversing the arc there. The desired new component height is \(b/g\), requiring \(g\mid b\).
* Otherwise write the affected interval as \(A\to I\to B\) and reverse \(u\to v\), with \(u\in A,v\in B\). The exterior contributes factor \(g\). The change in height is
\[
 b-a=g\left(S_A(u)E_B(v)Q(I)-\mathbf1_{I=\varnothing}\,e_A(u)s_B(v)\right).
\]
Indeed the new paths through \(v\to u\) split into a path ending at \(v\) and another starting at \(u\). The endpoint components supply the two marked masses. Each intervening component can independently distribute its vertices between the two paths. In the original graph a path using \(u\to v\) has adjacent endpoint components and uses a full path ending at \(u\) and a full path starting at \(v\). Paths using neither arc persist.

The nonsingleton components have heights at least three. Their possible ordered height lists are therefore finite, and the strong minimum bounds their orders. Write the middle \(I\) as its nonsingleton components and \(r\) singletons. Then \(Q(I)=2^r\prod_D Q(D)\). If \(I\ne\varnothing\), the ratio \((b-a)/(g S_A(u)E_B(v)\prod_D Q(D))\) must equal \(2^r\). This fixes \(r\) exactly. If there are no middle nonsingletons, test the adjacent case \(r=0\) separately with its subtraction. Either endpoint component may itself be a singleton. All placements with the indicated counts give the same change. Exterior singleton padding has no effect and can be omitted for minimum order. These alternatives are both necessary and sufficient because an ordinal sum realizes each admitted list of components.

## Why the finite certificate covers every order

Each of the 17 exclusions has lower endpoint at most 499. The six deeper positive pairs in the following table have lower endpoint at most 699. Busch's minimum at order 14 is 1,125 and increases thereafter, so all nonsingleton components required by this calculation have order at most 13. This is a bound on individual components. The inverse above still accounts for arbitrarily many singleton components. [Strong minimum and spectral census](../path_spectrum_two_angles_metagraph_20261106/README.md)

| Pair | Minimum order among all one-flip realizations |
|---|---:|
|199 ↔ 701|12|
|299 ↔ 3,301|12|
|399 ↔ 3,201|13|
|501 ↔ 3,099|12|
|599 ↔ 3,001|12|
|699 ↔ 1,801|12|

To generate the finite catalogue, grow strong tournaments of height at most 699 from the directed triangle using every mixed incoming/outgoing neighborhood for a new vertex. Every strong tournament of order at least four has a directed cycle of length one less than its order by Moon's pancyclicity theorem. The induced tournament on that cycle is a strong predecessor. Its path count is at most that of the extension because each Hamiltonian path admits an insertion of the new vertex and distinct old paths have disjoint extensions. Thus pruning by height preserves all eligible cores. Canonical digraph labeling removes duplicates. The resulting class counts for orders 3 through 13 are \(1,1,6,35,353,6008,17748,12194,5369,1150,23\). [Moon pancyclicity as used in THM-4836](../../../01-canon/theorems/THM-4836-coherent-carrier-sets-of-tournaments-and-the-pair-or-triple-law.md)

The positive certificates directly count paths by their last vertex on subsets, once before the reversal and once after it. Exhaustively including or excluding each available edge of the 101-edge graph gives 37 spanning 2-factors. The counts by number of connected components 1 through 6 are \(1,6,8,14,4,4\). Only the connected factor is a Hamiltonian cycle. Its edge \(399\leftrightarrow3201\) requires order 13. All witnesses have order at most 13, and a common appended sink increases order while preserving the two counts and the physical reversal.

## Why 199 cannot carry both of its edges

The neighbors of 199 are exactly 201 and 701. Since 199 is prime, any tournament of that height consists of one nontrivial SCC of height 199 and two runs of exterior singletons. All its possible nontrivial components occur at orders 8, 9 and 10, with 24, 3 and 12 isomorphism classes respectively. The complete marked-component calculation finds that a move to 201 requires component order 8, while a move to 701 requires component order 10. These required orders hold with arbitrarily much exterior padding. [Stored strong-height catalogue](../path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt)

An isomorphism preserves component order. Any count-preserving reversal also preserves the orders of the nontrivial components. A reversal that genuinely merges components strictly increases the height by the prior strict-merge theorem, and the inverse of a strict merge strictly decreases it. Neutral flips between adjacent singleton components can only exchange their order. Inserting or removing singleton SCCs leaves the nontrivial order unchanged. Thus even a sequence of these operations at height 199 cannot connect the two required incidences. [Strict merge and split](../chanpak_arc_fiber_grammar_20261122/README.md)

Run `python3 -B physical.py` here for all 101 literal pairs, the 2-factor histogram and the unique cycle. Run `python3 -B physical.py --exhaustive` to regenerate the strong catalogue from the triangle, perform the all-orders exclusions and minima, and verify the incompatible component orders at 199. The exhaustive command uses a C++17 compiler and nauty's `labelg` for isomorph rejection. [Checker](physical.py) · [Strong-extension and port computation](strong_ports.cpp)
