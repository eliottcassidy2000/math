# The parity of the order of a vertex pair

**PROVED.** For a finite tournament \(T\) and distinct vertices \(i,j\), let \(N_{ij}(T)\) count directed Hamilton words in which \(i\) occurs before \(j\). Count directed simple cycles once up to cyclic rotation. Then
\[
 N_{ij}(T)\equiv [i\to j]+\#\{C\text{ a directed simple cycle in }T:i,j\in V(C)\}\pmod 2.
\]
The vertices may be arbitrarily far apart in the word or cycle. In particular, to form the tournament whose arcs record odd precedence counts, start from \(T\) and reverse **every pair** on the support of each simple directed cycle, taking the cycle list from the original tournament. Multiple reversals cancel.

## Marking two positions in a path

Use the commutative ring \(R=\mathbf F_2[x_v:v\in V]/(x_v^2:v\in V)\). A squarefree monomial records the vertex set of a simple path. With adjacency matrix \(A\) (rows are sources), \(X=\operatorname{diag}(x_v)\), and all-ones column \(e\), set
\[
 P=X(I-AX)^{-1},\quad r=Pe,\quad c=e^tP,\quad F=1+e^tPe.
\]
The series terminates. Its \((i,j)\) entry counts paths beginning at \(i\) and ending at \(j\), including singleton paths. Rédei's oddness theorem on each induced subtournament gives \(F=\prod_v(1+x_v)\).

Adding the all-ones matrix \(ee^t\) to \(A\) reverses all arcs in characteristic two and adds loops. Loops contribute zero because they repeat a variable. On the other hand, expanding the new path series by the number of added transitions gives
\[
 P^t=P+Pe\,(1+e^tPe)^{-1} e^tP.
\]
Indeed, a word with one or more added transitions is a list of two or more ordinary nonempty paths, separated by unrestricted transitions. Repeated vertices across those paths give zero.

Fix \(i\ne j\), put \(U=V\setminus\{i,j\}\), and strip the endpoint variables from \(P_{ij}=x_i x_j p_{ij}\). Write \(r'_j\) for the start-at-\(j\) path polynomial, omitting its mandatory \(x_j\) and forbidding \(i\); similarly \(c'_i\) ends at \(i\) and forbids \(j\). In the \((j,i)\) entry of the preceding identity, extract the coefficient of \(x_i x_j\) with coefficients in the squarefree algebra of \(U\). This gives
\[
 c'_i r'_j=\prod_{u\in U}(1+x_u)\,(p_{ij}+p_{ji}).
\]
A word with \(i\) before \(j\) splits uniquely into a prefix ending at \(i\), a middle path from \(i\) to \(j\), and a suffix starting at \(j\). With the two marked variables omitted, its enumerator is therefore
\[
 c'_i p_{ij}r'_j=\prod_{u\in U}(1+x_u)\bigl(p_{ij}^2+p_{ij}p_{ji}\bigr).
\]
Every nonconstant monomial squares to zero, so \(p_{ij}^2=[i\to j]\). Two internally disjoint directed paths from \(i\) to \(j\) and back specify exactly one simple directed cycle through both vertices. Their product counts these cycles. Extract the coefficient of \(\prod_{u\in U}x_u\). Every cycle receives a unique completion from the prefactor, proving the formula.

## The odd sign fiber and its arc derivative

Fix a total label order, let \(B(T)\) count arcs pointing backward in that order, and let \(O(T)\) count Hamilton words that are odd permutations of the labels. Summing pair inversions and using the theorem yields
\[
 O(T)\equiv B(T)+\sum_C\binom{|C|}{2}
       \equiv B(T)+\sum_{k\equiv 2,3\bmod4} C_k(T)\pmod2,
\]
where \(C_k(T)\) counts simple directed \(k\)-cycles. This computes the product of the signs of all Hamilton words as \((-1)^{O(T)}\). Every automorphism of a tournament is an even permutation, since an involution would reverse an arc between its exchanged vertices and the automorphism group has odd order. Thus each ordinary tiling orbit has odd size and constant sign. The same formula gives the parity of the **number of odd ordinary tilings** at the labeled viewpoint.

For a physical flip of the arc on the pair \(a,b\), delete that arc and let \(L_\ell(a,b)\) count simple directed paths of \(\ell\) edges from either endpoint to the other in the remaining digraph. If \(T'\) is the flipped tournament, then
\[
 O(T')-O(T)\equiv 1+
   \sum_{\substack{\ell\ge2\\ \ell\equiv1,2\bmod4}}
          L_\ell(a,b)\pmod2.
\]
The label term toggles. Every cycle gained or lost by the flip is a unique endpoint path closed by the old or new arc. Thus for a four-vertex tournament only triangles matter; starting at six vertices hexagons can affect the flip as well.

## Recovering cyclic viewpoints from linear ones

For any vertex set \(S\) of size at least three and any specified pair \(i,j\in S\),
\[
 \#\{\text{Hamilton cycles of }T[S]\}
  \equiv\sum_{\{i,j\}\subseteq U\subseteq S} N_{ij}(T[U])\pmod2.
\]
A cycle on support \(W\supseteq\{i,j\}\) contributes once for each \(U\) between \(W\) and \(S\), hence \(2^{|S|-|W|}\) times. Only spanning cycles survive. The fixed orientation contributes \(2^{|S|-2}\) times and vanishes. This Boolean inversion recovers cyclic parity using any single specified pair, with no closure condition imposed on the Hamilton words.

**Consequence for a physical cycle seam.** Suppose reversing the arc on \(a,b\) changes the Hamilton-cycle count by an odd number. For *every* pair \(i,j\), an odd number of subsets \(U\supseteq\{a,b,i,j\}\) have different odd-precedence directions for this pair across the flip. Subtract the two inversion formulas above. Subsets missing either flipped endpoint cancel since their induced tournaments agree. In particular each observation pair detects the seam somewhere in its induced views, even when the two full ordinary path counts are equal. The incoming seven-vertex seam \((C,H)=(6,81)\to(7,81)\) makes that last possibility literal. The same conclusion applies at every one of the 4096 extended binary Golay contact words in the 28-vertex construction, since the marked arc always changes the cycle count by exactly one. The Golay restriction preserves all marginals on at most seven contact bits; the witness subset in the inversion can have more vertices. [Unit cycle seams, exact path leakage and the Leech–Golay restriction](../../../session212/cycle-spectrum/PORTS-AND-GOLAY.md)

**FINITE-EXACT replay.** Run `python3 -B order_parity.py` from this directory. It enumerates directed simple paths in every ordinary class of orders two through seven, compares the pair theorem and (from order three) its Boolean inversion, then tests all 991 single-arc flips from the ordinary class representatives through order six. Each test uses literal vertex subsets, directed words and cycle closures.
