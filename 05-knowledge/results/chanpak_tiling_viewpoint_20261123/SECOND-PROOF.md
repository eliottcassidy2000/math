# When all Hamilton viewpoints have the same permutation sign

Fix a labelling and count the directed Hamilton vertex words according to even or odd permutation sign, giving \(H_0,H_1\). Changing labels can interchange the two counts. Their minimum and the predicate that one vanishes are intrinsic. The following proof uses strong-component structure at a cut vertex. For highly connected tournaments its external input is the classical arc-Hamilton theorem. It also gives a quantitative estimate for how many viewpoints of the minority sign survive.

## The defect of alternation pays for minority viewpoints

**PROVED.** Let \(T\) have order \(2m+1\ge3\). For every vertex \(v\) and every directed Hamilton cycle \(C\) of \(T-v\) up to rotation, color the cyclic positions alternately. Let \(p,q\) count the outneighbors of \(v\) in the two respective colors, and put \(\delta(v,C)=m-|p-q|\). This is the minimum Hamming distance from the cyclic word of contacts with \(v\) to either of the two alternating binary words. Then
\[
 2\min(H_0,H_1)\ \ge\ \sum_{v,C}\delta(v,C)
 \ \ge\ \sum_v \mathrm{HC}(T-v)\,|d^+(v)-m|.
\]

Cut \(C\) at each possible starting position, alternating permutation signs because \(|C|\) is even. Placing \(v\) first or last across the same ordered word gives the same sign since it passes an even number of letters. As the starting position ranges over the first color, we may prepend \(v\) at \(p\) starts, or append it after a second-color vertex beating \(v\) at \(m-q\) cuts. The other color gives \(q+m-p\) words of the opposite sign. Thus the minority count among these \(2m\) words is exactly \(m-|p-q|\). Each resulting word can be charged at most twice in the global sum: \(v\) is one of its two endpoints, and the remaining path determines its unique cycle closure when it closes. For the second bound use \(|p-q|\le\min(p+q,2m-p-q)\).

In particular, when all full paths have one sign, every \(v\) for which \(T-v\) is strong must have degree \(m\), and its contacts alternate along *every* Hamilton cycle of the deletion. We use the term **deletable** for such a vertex. The estimate also applies when other deletions are reducible.

## Strong tournaments of one sign

**PROVED using cited arc-Hamiltonicity and the finite base.** The only positive-order strong tournaments whose Hamilton paths all have the same permutation sign are \(K_1\), the directed triangle and \(Q_5=C_3[C_3,K_1,K_1]\). In the last example the outer blocks follow a directed triangle; write them as \(x\to C_3\to y\to x\). At even positive order a strong tournament has a Hamilton cycle, and cutting after two neighboring vertices gives paths of opposite signs. The [finite base](sign_balance.py) checks all ordinary classes of orders three, five and seven. At five the unique strong class of one sign has 15 paths. Its identification with \(Q_5\) also follows from the three cyclic internal orders of its triangle and the five sign-preserving cuts of each resulting cycle. [THM-4729, closed path counts](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md) · [Fully closed classes](../merged_metagraph_closure_reconciliation_20261020/CLOSURE_AND_NONSTRONG_ASYMPTOTIC.md)

Suppose by induction that an odd-order strong tournament \(T\) of order \(n\ge9\) has only one sign. There are two alternatives.

**No cut vertex.** Every vertex is deletable, so \(T\) is regular and \(T-v\) has even order at least eight with all outdegrees differing by at most one. Jakobsen's theorem says every arc of such an almost-regular tournament lies on a Hamilton cycle. Choose an arc between two vertices in \(N^+(v)\); some Hamilton cycle of \(T-v\) contains this arc. Its contact word fails alternation. [Jakobsen's exact theorem as stated in Theorem 1.8 of *On arc-traceable local tournaments*](https://doi.org/10.1016/j.disc.2007.12.042)

There is another derivation using only Thomassen's theorem that every arc of a 3-connected tournament is in a Hamilton cycle. In a regular \(n\)-tournament any cut of size \(k\) leaves nonempty parts \(D\to E\). The average outdegree in \(E\) is at most \((|E|-1)/2+k\), so \(|E|\ge n-2k\). The same result for \(D\) follows by converse. Hence \(k\ge n/3\). At \(n\ge11\), deleting a vertex leaves a 3-connected tournament, so the same contradiction applies. At \(n=9\), either the connectivity is at least four, or a minimum cut has size three. In the latter case equality holds in both average bounds: \(D,E\) each have order three and all cross edges with the cut \(K\) follow \(D\to E\to K\to D\). Regularity forces each of the three blocks to be a directed triangle. Delete one vertex from \(D\). Traversing the three remaining blocks as contiguous directed paths gives a Hamilton cycle on the eight remaining vertices, including three consecutive outneighbors in \(E\) of the deleted vertex. Again alternation fails. [Thomassen, *Hamiltonian-connected tournaments*](https://doi.org/10.1016/0095-8956(80)90061-1)

**A cut vertex \(v\).** Write the strong components of \(T-v\) in their forward order as \(B_1,\ldots,B_r\), with \(r\ge2\). Some first-component vertex is beaten by \(v\) and some last-component vertex beats \(v\). A strong component has a Hamilton path starting or ending at any specified vertex by cutting a Hamilton cycle. Fix a path of the last component ending at a vertex beating \(v\), and append \(v\) after the blocks. This extends all Hamilton words of every earlier component independently. Fixing a compatible path of the first component and prepending \(v\) does the same for the last. A component with paths of both signs would give paths of both signs in \(T\). Hence all components have a single sign and, by induction, odd size. Every Hamilton word of \(T-v\) has the same sign because all its words are concatenations of these fixed-order blocks.

Since some such word allows insertion of \(v\) in slot zero, in *every* such word all legal insertion slots have an even number of preceding letters. If \(v\) has mixed contacts with one nonsingleton component, choose a Hamilton cycle in that component. On the cycle some consecutive vertices \(b,a\) obey \(b\to v\to a\). Rotating its cut puts \(b,a\) at any prescribed consecutive internal positions of the block. There is an internal slot of odd absolute index because the block has at least three vertices. This contradiction shows that \(v\) is uniform toward each component.

If at least one component is nonsingleton, contract all the \(B_i\) and keep \(v\). The proper quotient is strong. It has one sign, since substituting a fixed path of each odd block into a quotient path preserves its permutation sign up to a common multiplier. Induction leaves the directed triangle as the only possible quotient. Indeed deleting \(v\) makes the quotient transitive, whereas no single vertex of \(Q_5\) meets all of its directed triangles. Thus \(T\) consists of \(v\to F\to L\to v\). If \(|F|>1\), every vertex of \(F\) is deletable because the three nonempty blocks still form a directed cycle of blocks after its deletion. Consequently each vertex of \(F\) has total outdegree \((|F|+|L|)/2\). Averaging its internal degrees gives
\[
 |L|+(|F|-1)/2=(|F|+|L|)/2,
\]
so \(|L|=1\) and \(F\) is regular. Its only possible nontrivial regular induction case is a directed triangle, which gives \(n=5\). Interchange the roles using converse if only \(L\) is nontrivial.

Finally, if all the components are singleton, write their transitive order as \(x_1,\ldots,x_{n-1}\). Strong connectivity forces \(v\to x_1\) and \(x_{n-1}\to v\). Every interior \(x_i\) is deletable by retaining the same rim through the remaining vertices. But its outdegree is \(n-1-i+[x_i\to v]\), which cannot be constant over the interior indices when \(n\ge7\). This finishes the induction.

## The ordinary isomorphism classes at every size

Every tournament has all viewpoints of one sign exactly when its strong-component word uses only \(K_1,C_3,Q_5\). In particular the joint series by size and the two nontrivial multiplicities is
\[
              (1-z-uz^3-wz^5)^{-1}.
\]
A word with \(a\) triangles and \(b\) five-blocks has ordinary Hamilton count \(3^a15^b\). This gives the equality case of the signed path triangle inequality \(|H_0-H_1|\le H(T)\) for all tournaments. Sign unanimity alone places no bound on the number of paths.

The tiling masses and converse merge are especially simple **on this complete family**. Both nontrivial strong blocks have automorphism group of order three. In \(Q_5\), the outer singleton degrees are one and three and distinguish them, and the only freedom is rotation of the inner triangle. The component order allows no interchanges between blocks. Thus a word with \(a,b\) as above has \(3^{a+b}\) automorphisms, giving exactly \(5^b\) ordinary tilings via the free action on Hamilton words. All three letters are self-converse. An ordinary component word therefore merges with its reversal. Palindromes have odd merged mass \(5^b\); a nonpalindrome and its reverse give even merged mass \(2\cdot5^b\).

More precisely let \(W(z,u,w)=(1-z-uz^3-wz^5)^{-1}\). The exact generating functions for the two sorts of merged nodes in this family are
\[
 P(z,u,w)=\frac{1+z+uz^3+wz^5}{1-z^2-u^2z^6-w^2z^{10}},
 \qquad D(z,u,w)=\frac{W(z,u,w)-P(z,u,w)}2.
\]
The coefficient of \(z^n u^a w^b\) in \(P\) counts the odd nodes and in \(D\) the doubled nodes at that size, content and prescribed mass. The palindrome consists of an arbitrary first half, its forced mirror and an optional central letter. For instance, among sign-unanimous presentations the first ordinary mass five is the single \(Q_5\) at size five; size six has its two placements beside a singleton, forming one doubled node of mass ten.

Run `python3 -B sign_balance.py` to regenerate the finite base from [the complete ordinary class bank](../../../session212/cycle-levels/CYCLE_LEVELS.md) and compute the cyclic-deletion bound on every class of orders three, five and seven.
