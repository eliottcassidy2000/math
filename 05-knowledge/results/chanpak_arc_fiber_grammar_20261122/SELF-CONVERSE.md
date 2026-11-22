# Converse symmetry detects parity of tournament birth and flip multiplicities

Write \(H(T)\) for the number of directed Hamiltonian paths. The converse reverses every arc and preserves this count. It acts both on the ordinary height fiber and on the one-missing-arc objects whose two completions have heights \(h,h+2\). In either fiber, the total number of isomorphism classes equals the number fixed by converse plus twice the number of distinct converse pairs. This applies at the smallest realizing order as well.

## Even order and the exponent-61 seam

**PROVED.** Every self-converse tournament of even order has \(H\equiv1\pmod4\). Hence every even-order scalar fiber of height \(3\bmod4\) and every even-order marked unit fiber splits entirely into distinct converse pairs.

A tournament automorphism group has odd order. An involutory automorphism would exchange two vertices and reverse the arc between them. For a self-converse tournament on at least two vertices, its group of automorphisms and anti-automorphisms has twice this odd order. Its involutions therefore reverse all arcs. Such a permutation \(\theta\) fixes at most one vertex, since an arc between two fixed vertices would reverse. At even order every orbit is a pair.

The weighted odd-cycle expansion reduces modulo four to
\[
 H(T)\equiv 1+2c_{\rm odd}(T)\pmod4,
\]
where \(c_{\rm odd}\) counts simple directed odd cycles of length at least three up to cyclic rotation. Reversing a cycle and applying \(\theta\) permutes these cycles. When every vertex orbit is a pair, an odd cycle cannot have invariant support, so \(c_{\rm odd}\) is even. For a self-converse missing-arc object the two distinct heights distinguish its completions, hence each completion is itself self-converse. Consecutive odd heights cannot both be one modulo four. [Grinberg–Stanley expansion in Irving–Omar, Corollary 20](https://arxiv.org/html/2412.10572v3#S3.SS4)

For \(E=2^{60}(2^{61}-1)\), every even-order coefficient of \(F_{E-1}\) is therefore even. A unique smallest lower witness must have odd order. At any even order a physical flip from \(E-1\) to \(E+1\) supplies a second inequivalent marked flip by converse. Any unpaired flip must occur at odd order.

## An odd-order quotient identifies the arcs that can change the residue

**PROVED.** Take a tournament on \(2m+1\) vertices reversed by an involution \(\theta\), with \(m\ge1\). Write \(f\) for its fixed vertex. Form a tournament \(Q\) on the \(m\) vertex pairs. Direct \(A\to B\) when a vertex of \(A\) has an odd number of outneighbors in \(B\). With the vertices of each pair listed in matched order, the incidence matrix from \(A\) to \(B\) is
\[
 \begin{pmatrix}a&b\\1-b&1-a\end{pmatrix}.
\]
Its two row sums agree modulo two, and either reverse row sum has the opposite parity. Thus the definition of \(Q\) is intrinsic.

For a pair \(J\), let \(b_J\) be the parity of the number of all directed simple paths in \(Q\) ending at \(J\), including the singleton. Let \(T^J\) reverse the single arc joining the two vertices of \(J\). Then
\[
             H(T^J)-H(T)\equiv 2b_J\pmod4.
\]
Only cycles that use this arc affect the difference. Nonfixed cycles occur in pairs. A fixed odd cycle through the arc consists of a directed path from \(f\) to one endpoint of \(J\), the crossing of \(J\), and the reflected reverse of that path back to \(f\). The outward path visits each pair at most once. Conversely each such outward path gives one fixed cycle in exactly one orientation of the crossing. There is no sign to retain modulo two. For any sequence of distinct pairs ending in \(J\), \(f\) has exactly one outneighbor in the first pair. At each subsequent step the parity of the number of continuations from either lift is the corresponding adjacency entry of \(Q\). The number of lifts is odd exactly for directed paths in \(Q\), giving the formula.

The same bit \(b_J\) equals the parity of the **Hamiltonian paths of \(Q\) starting at \(J\)**. Adjoin a vertex \(x\) with \(J\to x\) and with \(x\) dominating all other vertices. A Hamiltonian path using \(J\to x\) has a prefix ending at \(J\), followed by \(x\) and a Hamiltonian suffix on the unused vertices. Each possible prefix admits an odd number of suffixes by Rédei's theorem, including the unique empty suffix. Paths avoiding \(J\to x\) must start at \(x\) and continue along a Hamiltonian path of \(Q\) starting elsewhere. Both tournaments have odd total count, giving the asserted identity. Applied also to the converse, it recovers the equality modulo two of the counts of proper paths starting and ending at a vertex. [El Sahili and Abi Aad, path parity at a vertex](https://doi.org/10.1016/j.dam.2019.09.005)

Every Hamiltonian path starts in the first strong component of \(Q\), so \(b_J=0\) for every later component. For a transitive \(Q\) the source is the sole possible crossing for a physical unit flip. For a general quotient the test gives the full change modulo four.

**PROVED.** At every odd order \(2m+1\ge3\), exactly half of the unlabelled self-converse tournament classes have \(H\equiv1\pmod4\), and half have \(H\equiv3\pmod4\). Indeed \(\sum_J b_J\equiv 2^m-1\equiv1\pmod2\), since each nonempty subset of \(Q\) contributes its odd number of Hamiltonian paths. Reverse all \(m\) arcs joining paired vertices. The quotient and all \(b_J\) stay fixed, so this adds two modulo four. The operation is an involution on classes with a chosen reversing involution. Forgetting that choice is bijective on isomorphism classes. In fact the reversing involutions generate Sylow 2-subgroups of the extended automorphism group and are conjugate there. If a conjugator is reversing, composing it with the source involution produces an ordinary automorphism conjugating the same two involutions. Thus the residue-switching operation also gives the desired bijection on ordinary self-converse classes.

## An exact count with one fixed missing pair

**PROVED.** Put \(m=\lfloor n/2\rfloor\) for \(n\ge2\). Fix \(\theta=(0\ 1)(2\ 3)\cdots\) and omit \(e=\{0,1\}\). Let \(L_n(h)\) count the labelled partial orientations reversed by \(\theta\) whose two completions have heights \(h,h+2\) in either order. The number of self-converse marked classes is
\[
       \frac{L_n(h)}{2^m(m-1)!}.
\]
Any automorphism of the partial orientation fixes each completion since their heights differ. It acts freely on Hamiltonian paths of each, so its order divides \(\gcd(h,h+2)=1\). A self-converse unit object thus has exactly one anti-automorphism, whose square is the identity. This involution exchanges the endpoints of the omitted pair because it reverses the unique lower completion. All its other orbits have size two with at most one singleton. After labelling the involution and mark as prescribed, the possible relabellings are precisely the permutations commuting with \(\theta\) and preserving \(e\). They independently exchange the two vertices in each pair and freely permute the \(m-1\) unmarked pairs. Rigidity makes these labellings distinct.

Writing \(M=\binom n2\), the involution fixes \(m\) unordered vertex pairs and exchanges the other \(M-m\) pairs in twos. One orientation bit for each such orbit, omitting \(e\), gives \(2^{(M+m)/2-1}\) partial orientations. Each gives exactly two labelled scalar tournaments reversed by \(\theta\). In an ordinary scalar fiber known to be rigid, the number of its self-converse classes is the number of such scalar encodings of height \(h\) divided by \(2^m m!\). Prime heights at least five and \(E+1\) satisfy the inherited all-orders rigidity criterion. [Marked core grammar and rigidity proofs](README.md)

## Smallest classes fixed by converse

**FINITE-EXACT.** The table counts every self-converse marked unit class, summed over all heights. A central class has no exterior singleton run on either side. In the last column each listed lower height has one fixed central class.

| Order | Labelled partial orientations for the fixed involution and mark | Fixed marked classes | Fixed central classes | Central lower heights |
|---:|---:|---:|---:|---|
| 2 | 1 | 0 | 0 | |
| 3 | 2 | 1 | 1 | 1 |
| 4 | 8 | 0 | 0 | |
| 5 | 32 | 3 | 2 | 9, 13 |
| 6 | 256 | 0 | 0 | |
| 7 | 2,048 | 4 | 1 | 65 |
| 8 | 32,768 | 0 | 0 | |

For \(27\to29\), the full central polynomial is \(2z^6+2z^7\). Neither order has a self-converse central class, so the two central types at each order form a converse pair. At order seven the full fiber also contains four one-sided paddings from order six, giving six marked classes. Duality exchanges the lengths of the two exterior singleton runs. A fixed marked class would require both a fixed central class and equal runs. Consequently the self-converse \(27\to29\) series vanishes at **every** order. The same argument pairs the two \(11\to13\) central classes, even though the unique smallest tournaments at heights 11 and 13 are each self-converse. [All-order central numerators and exterior factorization](README.md)

Among the smallest Fermat-prime scalar witnesses, both rigid classes of height 17 at order six are individually self-converse. All 30 rigid classes of height 257 at order eight occur in 15 distinct converse pairs. Thus these two even birth multiplicities have different explanations. The comparison uses the full fixed-involution scalar enumerations of 512 and 65,536 labelled tournaments, respectively. [Fermat-prime birth counts](README.md)

Run `python3 -B self_converse.py --through 8` from this directory. The program ranges over every free orientation bit, computes both Hamilton path counts by subset and endpoint, tests centrality by the strong connectivity of the bidirected completion, and checks the even-order obstruction and the odd-order quotient formula on every enumerated partial orientation. Its default target column gives the counts for \(27\to29\). The two smaller exceptional unit fibers are included. [Enumeration](self_converse.py)
