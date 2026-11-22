# Transitive modules, matching poles, and an exact permanent reduction

**Status:** PROVED (combinatorial reductions using named prior work). This explains a structural feature of the Chan–Pak count and proves complexity results, not its full range. No literature-priority claim is made.

**What changes:** The earlier substitution kernel and balanced permanent coordinate become an exact one-call reduction from the 0–1 permanent. The known strong self-converse decoder then transfers hardness to the symmetric stratum. The intervening matching rank identifies the top live coordinate even in unbalanced contexts; we also give a unit-slope dyadic obstruction, cycle and palindromic variants, and a last-exit determinant for the general routing kernel. The familiar ingredients are included to make the proofs readable in isolation.

Let \(H(T)\) count directed spanning paths of a tournament, with \(H(\varnothing)=1\). Given a marked vertex \(v\), replace \(v\) by a \(k\)-vertex transitive tournament; every new vertex has the same edges to the exterior as \(v\) did. Write the resulting count as \(h_v(k)\), \(k\ge1\).

## 1. The divided-power coordinates of a transitive module (inherited)

This is the transitive-factor specialization of repository THM-3121 and Session 33’s boundary-observer law, reproduced for notation.

Collapse each maximal run of vertices from the new block to a symbol \(v\). The resulting words use each exterior vertex once and \(v\) some \(r\ge1\) times, never adjacently, and every transition is a base arc. Write \(a_r\) for the number of such words. For each word, assigning the \(k\) distinct block vertices onto its \(r\) ordered runs determines the path uniquely: inside each run they must follow the transitive order. Thus
\[
h_v(k)=\sum_{r\ge1}a_r\,r! \left\{\begin{matrix}k\\r\end{matrix}\right\}.
\]
The braces denote Stirling numbers of the second kind. In particular \(a_1=H(T)\), and the coefficients are recovered recursively from evaluations in increasing \(k\), dividing by the diagonal \(k!\) at step \(k\). Equivalently,
\[
h_v(k)=\sum_{j\ge1}c_j j^k,
\qquad c_j=\sum_{r\ge j}(-1)^{r-j}\binom rj a_r.
\]
This follows by inclusion–exclusion on missing runs. The first coordinates are nonnegative combinatorial integers; the exponential coordinates can be signed.

## 2. The largest pole is exactly a matching number

Let \(A=N^+(v)\), \(B=N^-(v)\), and \(G_v\) be the bipartite graph on \(A,B\) with \(ab\) an edge exactly if \(a\to b\). Let \(\mu=\nu(G_v)\). In a word with \(r\) occurrences of \(v\), each of the \(r-1\) interior gaps starts in \(A\) and ends in \(B\). Choose one forward \(A\)-to-\(B\) edge in each gap. The chosen edges are disjoint, so \(r\le\mu+1\).

Conversely, order the edges \(a_1b_1,\ldots,a_\mu b_\mu\) of a maximum matching in any order. Start with a Hamiltonian path of the unmatched vertices of \(B\), then use \(v,a_1,b_1,v,\ldots,a_\mu,b_\mu,v\), and finish with a Hamiltonian path of the unmatched vertices of \(A\). The prefix and suffix exist by Rédei's theorem and can be empty. The construction also works for any size-\(d\) matching with \(r=d+1\). Thus \(a_r>0\) exactly for \(1\le r\le\mu+1\). In particular
\[
\lim_{k\to\infty}\frac{h_v(k)}{(\mu+1)^k}=a_{\mu+1}>0.
\]
The formal ordinary generating function \(\sum_{k\ge1}h_v(k)z^k\) therefore has a simple pole at \(z=1/(\mu+1)\) and no pole closer to zero. This rate is not controlled by the number of incident directed triangles, but by the maximum number with pairwise disjoint exterior pairs.

More generally, if other vertices are replaced by nonempty blocks of sizes \(w_u\), independent of \(k\), then the same rate is \(1+\min_C\sum_{u\in C}w_u\), where \(C\) ranges over vertex covers of the *base* graph \(G_v\). To see this, apply König's theorem to the bipartite graph in which each vertex \(u\) is blown up to \(w_u\) copies with all bipartite edges dictated by \(G_v\). A cover must contain all copies of at least one endpoint of each base edge, and any base cover gives a cover of the stated total size. The internal orientations of the other blocks do not matter. In particular the weighted growth-rate data determine \(G_v\): give any tested pair of opposite-side vertices weight \(|A|+|B|+1\) each and all others weight 1. The minimum cover pays a large weight if and only if the tested pair is adjacent.

### The coefficient-one obstruction for dyadic tuning

Suppose the dominant base is \(2\), and the coefficient of \(2^k\) is \(1\). Then \(a_2=1\). Each edge \(a\to b\) of \(G_v\) gives a distinct two-run skeleton by using \(v,a,b,v\) between paths on the remaining \(B\) and \(A\). Hence there is exactly one crossing edge \(a\to b\). If some \(x\in A\setminus\{a\}\) beats \(a\), then \(v,x,a,b,v\) gives another skeleton with the same unmatched-vertex construction. If \(b\) beats some \(y\in B\setminus\{b\}\), use \(v,a,b,y,v\) instead. Therefore \(a\) is the source of \(T[A]\) and \(b\) the sink of \(T[B]\). The orientations now force the ordered sum
\[
T[B\setminus\{b\}]\ \Longrightarrow\ C_3(a,b,v)\ \Longrightarrow\ T[A\setminus\{a\}].
\]
The endpoint path product on the two outer blocks is at most \(a_2=1\), since all of those prefix/suffix choices give distinct skeletons. As each factor is positive, both equal 1. Consequently \(H(T)=3\), and
\[
h_v(k)=2^k+1.
\]
Thus allowing all marked hosts and cloning just one transitive module does not supply arbitrary constant corrections at unit dyadic slope. A proof of completeness by dyadic tuning would have to vary more than this single coordinate, or allow a different slope.

## 3. Exact counting hardness from the extremal run coefficient

The identity at the top coordinate is inherited from Session 33 (Boundary Observability §3a). The reduction below uses it in the inverse direction: fix a matrix encoding in the context, vary only the order of a transitive module, then extract the permanent from the count.

Let \(M\in\{0,1\}^{m\times m}\), \(m\ge1\), with \(|A|=|B|=m\). Put \(v\to A\), \(B\to v\), and \(a_i\to b_j\) if and only if \(M_{ij}=1\). Orient \(A\) and \(B\) transitively. The only way to have \(m+1\) occurrences of \(v\) is to fill all \(m\) internal gaps with exactly two exterior vertices each, with no prefix or suffix. Such a skeleton is an ordered perfect matching of the matrix; hence
\[
a_{m+1}=m!\operatorname{per}(M).
\]
This is valid when the permanent is zero, too. Thus counting directed Hamiltonian paths of tournaments is \(\#\mathrm P\)-complete under polynomial-time Turing reductions: evaluate the \(m+1\) tournaments of orders \(2m+k\), \(1\le k\le m+1\), then invert the triangular system. A permutation is a polynomial-size witness verifiable in polynomial time. Hardness is by Valiant's theorem for the zero–one permanent. The orientation needs only three internally transitive blocks, with the interblock edges \(B\to V\to A\) fixed and arbitrary edges between \(A\) and \(B\).

In fact a single count suffices by polynomial-time integer postprocessing. Put \(b=m+1\) and
\[
L=3m+4+2m\left\lceil\log_2(2m)\right\rceil,\qquad k=bL.
\]
Any skeleton is specified by an exterior permutation and a subset of the \(2m+1\) insertion slots; hence \(\sum_r a_r\le(2m)!2^{2m+1}\). The exponential coordinates satisfy \(\sum_j|c_j|\le 2^{b}\sum_r a_r\le(2m)!2^{3m+2}\). Using \(c_b=a_b\) and \((1-1/b)^b\le e^{-1}<1/2\) gives
\[
\left|h_v(k)/b^k-a_b\right|
 \le (2m)!2^{3m+2}(1-1/b)^k \le 1/4.
\]
Round to the nearest integer and divide by \(m!\). The tournament order and the bit lengths are polynomial in \(m\). This is a one-call metric reduction, not a parsimonious reduction: Rédei oddness precludes preserving all input values verbatim. By Zankó’s one-call completeness theorem for the 0–1 permanent, the tournament count is consequently complete under polynomial-time metric reductions.

## 4. The general module law: Euler flows paired against ordered path covers

The path-cover substitution law itself is inherited from THM-3121. We give its grouping by Euler flow and the outgoing-Laplacian coefficient.

For a tournament \(U\), let \(P_r(U)\) count ordered covers by \(r\) nonempty vertex-disjoint directed paths (singletons allowed). For a base tournament \(D\), replace each base vertex \(i\) by a nonempty module \(U_i\). Compress each maximal same-module run of a Hamiltonian path to \(i\). Conversely, for any directed walk \(W\) of \(D\) using every base vertex, its \(r_i(W)\) occurrences of \(i\) accept an ordered path cover. It follows that
\[
H(D[U_1,\dots,U_d])=\sum_W \prod_i P_{r_i(W)}(U_i).
\]
Only \(r_i(W)\le |U_i|\) can contribute, so the sum is finite. For transitive \(U_i\), these are the Stirling coordinates above.

There is an exact compressed formula: group the words by their start \(s\), end \(t\), and directed transition counts \(m_{ij}\). Only consider transition counts on base arcs; they must satisfy
\[
d_i^+-d_i^-=\mathbf1_{i=s}-\mathbf1_{i=t},\qquad
r_i=d_i^++\mathbf1_{i=t}=d_i^-+\mathbf1_{i=s}\ge1.
\]
Mark the parallel edges of each multiplicity temporarily and let \(\tau_t\) count the directed spanning trees into \(t\) in this multigraph. The number of compatible words is
\[
W_{s,t}(m)=\frac{\tau_t\,\prod_i(d_i^+-1+\mathbf1_{i=t})!}{\prod_{ij}m_{ij}!}.
\]
Proof: the last exits at nonterminal vertices of a labeled Euler trail form an in-tree rooted at \(t\); at each nonterminal vertex the tree edge is placed last in the local exit ordering, while all exits at the root are freely ordered. Conversely follow these local orderings from \(s\). By degree balance the walk cannot stop before \(t\); if it stops there with an unused edge, the unused edge balance and the remaining last exit tree edges force an unused exit from \(t\), a contradiction. More explicitly, the unused graph at that stop is balanced, so unused outgoing edges force unused incoming edges at the same vertex and vice versa. Following each unused last-exit edge toward the root therefore propagates unused outgoing edges all the way to \(t\). Each vertex-word then has exactly \(\prod_{ij}m_{ij}!\) labelings. The one-vertex, no-edge case uses the empty tree and \(0!=1\). By the directed matrix-tree theorem \(\tau_t\) is the \(t\)-minor of the outgoing Laplacian \(\operatorname{diag}(d^+)-m\). This separates the internal order symmetry of each module from the directed routing and last-exit tree on the quotient.

## 5. A strongly connected self-converse shell codes a square

The shell and its one-probe exact decoder are inherited from the repository’s November 18 strong blue probes. The complexity claim uses the permanent reduction in §3, not a new claim to the shell or decoder.

Given any nonempty tournament \(U\), begin with the ordered sum
\[
 a\ \Longrightarrow\ U\ \Longrightarrow\ C_k\ \Longrightarrow\ U^{\rm op}\ \Longrightarrow\ b,
\]
where \(C_k\) here means a transitive chain of length \(k\ge1\), not a cycle. Reverse just the arc \(a\to b\), so the new closing arc is \(b\to a\). Denote the result by \(S_k(U)\). It is strongly connected: every middle vertex can reach \(b\) and can be reached from \(a\), and \(b\to a\). It has an involutive anti-automorphism exchanging \(a,b\), matching the two copies of \(U\) vertexwise, and reversing the central chain. If \(k\) is odd the middle of that chain is the unique fixed vertex.

Write
\[
K(U)=\sum_{X\subseteq V(U)} H(U[X])H(U[V(U)\setminus X]).
\]
This counts an ordered cover by two paths where either may be empty. It is multiplicative for ordered sums: every split restricts independently to the blocks, and within each side the corresponding block paths are uniquely concatenated. It is also invariant under converse. For a transitive chain \(K(C_k)=2^k\).

A Hamiltonian path of \(S_k(U)\) avoiding \(b\to a\) is a Hamiltonian path of the original ordered sum: any such path uses at least one middle vertex, so never uses \(a\to b\). Those paths have count \(H(U)^2\). Every path using \(b\to a\) is uniquely \((P_X,b,a,P_Y)\), with \(X,Y\) any ordered partition of the middle vertex set and \(P_X,P_Y\) Hamiltonian paths on the parts. This is because the arcs from the middle to \(b\) and from \(a\) to the middle are universal. Thus
\[
H(S_k(U))=H(U)^2+2^kK(U)^2.
\]
This proves three consequences.

1. **Strong self-converse counting hardness.** Let \(n=|U|\), choose an odd \(k\) with \(2^k>(n!)^2\), for example \(2n\lceil\log_2(\max(2,n))\rceil+1\). Reduce \(H(S_k(U))\) modulo \(2^k\) and take the integer square root to recover \(H(U)\). This is a polynomial-size metric reduction. Combined with §3, exact tournament path counting remains \(\#\mathrm P\)-complete by one-call metric reductions even restricted to strongly connected self-converse tournaments of odd order with a displayed involutive anti-automorphism having one fixed vertex.
2. **Palindromic counting hardness, without a modular decoder.** For odd \(k\), no path fixed by reversal composed with the displayed involution can use the arc \(b\to a\). Indeed, that arc is mapped to itself, and no adjacency position in an odd-order vertex sequence is fixed by reversal. Using it would repeat the same arc twice. The unflipped paths are fixed exactly when their two outer paths correspond by reversed copy; the central chain already has the unique compatible order. The fixed-path count is exactly \(H(U)\). In particular counting the palindromic paths of a strongly connected odd-order tournament supplied with such an involution is \(\#\mathrm P\)-complete under the same one-call reductions.
3. **A norm-form boundary.** For even \(k\) the count is a sum of two integer squares; for odd \(k\) it is \(x^2+2y^2\), with \(x=H(U)\), \(y=2^{(k-1)/2}K(U)\). Thus inert primes of the relevant norm form can occur only to even valuation (for even \(k\), primes \(3\bmod 4\); for odd \(k\), odd primes \(5,7\bmod 8\)). For example, an inert odd prime dividing the norm forces both coordinates divisible by the prime, and one then removes its square recursively. This is a restriction on the shell, not on arbitrary self-converse tournaments. It makes explicit why symmetry in this particular construction supplies exact modular information while preventing any coprime inert-prime coverage.

For an OCF interpretation of the fork term, more generally let \(K_q(U)\) count ordered covers by \(q\) paths with empty paths allowed. If \(\mathcal C\) is a collection of disjoint directed odd cycles and \(u(\mathcal C)\) counts uncovered vertices, then
\[
K_q(U)=\sum_{\mathcal C}2^{|\mathcal C|}q^{|\mathcal C|+u(\mathcal C)}.
\]
Apply OCF on each part of an ordered split of \(V(U)\), then swap summations: each selected cycle must be monochromatic in the split, while uncovered vertices choose parts independently. In particular \(4\mid K(U)\) for \(|U|\ge2\): every term has power \(2|\mathcal C|+u(\mathcal C)\ge2\). This is a vertex-graded refinement of the ordinary odd-cycle packing enumerator, rather than evaluation at \(2q\) alone.

## 6. The same reduction counts directed Hamiltonian cycles

Use the matrix gadget of §3 and require a cyclic path, counted once up to cyclic rotation. Root it at the distinguished exterior vertex \(a_1\). Compress runs of the transitive block as before; exterior vertices remain singletons. A cyclic skeleton with \(r\) occurrences of \(v\) has \(r\) gaps from \(A\) to \(B\), so \(r\le m\). Each contributes \(r!\left\{\begin{smallmatrix}k\\r\end{smallmatrix}\right\}\) expansions. At the top \(r=m\), each gap is exactly a matched pair and the pair containing \(a_1\) is forced into the first position after the chosen cyclic cut. The coefficient is therefore \((m-1)!\operatorname{per}(M)\). Triangular inversion proves \(\#\mathrm P\)-completeness for cycles under Turing reductions; the same rounding argument gives a one-call metric reduction. In the latter argument the total number of skeletons is at most \((2m)!2^{2m+1}\) by choosing the exterior order and occupied cyclic gaps, so the same overestimate for \(L\) works with base \(m\) when \(m\ge2\); \(m=1\) is trivial. If desired, restrict the oracle to strongly connected tournaments: first test in polynomial time whether the bipartite graph has a perfect matching, return zero if not, and otherwise each pair of a perfect matching gives a cycle through the marked block. All vertices then communicate through that block, including its different copies.

## 7. A parity-preserving 2-adic residue proof using the same shell

This separates even- and odd-order witnesses, using the inherited strong blue probe with a padded central root. It is a companion approach to the earlier dyadic results (notably Session 197’s stronger constraints for exactly two Hamilton cycles on its admissible residue set, and the November 18 CRT construction for strong self-converse counts).

Replace the central chain in §5 by \(C_m\Rightarrow R\Rightarrow C_m\), where \(R\) is a fixed self-converse tournament with a displayed anti-involution. Exchange and reverse the two identical chains under the shell involution. The same proof yields a strongly connected self-converse tournament \(S\) with
\[
H(S)=H(R)H(U)^2+2^{2m}K(R)K(U)^2,
\qquad |S|\equiv |R|\pmod2.
\]
The second term can be switched off modulo any fixed \(2^s\) by increasing \(m\).

For \(s\ge3\), multiplication by \(25\) generates the subgroup of residues \(1\bmod8\) modulo \(2^s\). Indeed \(v_2(25^{2^j}-1)=3+j\): it holds at \(j=0\), and squaring adds exactly one because \(25^{2^j}+1\equiv2\pmod8\). Every power of \(25\) is \(H(U)^2\bmod 2^s\), by taking \(U\) to be an ordered sum of copies of a tournament of count 5; for exponent zero use a singleton. Thus any root of count \(e\) covers its entire residue class modulo 8, at every finer binary precision.

At even order take a transitive pair (\(e=1\)) and a transitive four-chain with its endpoint arc flipped (\(e=5\)); both are self-converse. This covers all residues \(1\bmod4\) with strong self-converse even-order hosts, meeting the known necessary even-order congruence (repository THM-791(a)). At odd order take a singleton and a directed triangle (\(e=1,3\)). Two five-vertex additional roots have \(e=13,15\), giving the other two odd classes. For explicit roots on \(0,1,2,3,4\), orient from larger to smaller except for the pairs below; the anti-involution is \(i\mapsto4-i\).

| Count | Forward pairs \((i,j)\) with \(i<j\) | # 3-cycles | # 5-cycles |
| --- | --- | ---: | ---: |
| 13 | \(03,14\) | 4 | 2 |
| 15 | \(02,03,14,24\) | 5 | 2 |

The directed 5-cycles in the two rows, each rooted at 0, are respectively \((0,3,1,4,2)\), \((0,3,2,1,4)\), and \((0,2,4,3,1)\), \((0,3,2,1,4)\); enumerate the possible first exits and continue along the available arcs to check the lists. The 3-cycle column follows from \(\binom53-\sum_i \binom{d_i^+}{2}\). The odd-cycle collection formula then gives \(1+2(\#3\text{-cycles}+\#5\text{-cycles})\) at order five, because two cycles cannot be disjoint there.

Reducing the constructions modulo \(2^s\) for \(s<3\) gives the compatible low-precision statements. This proves that no further finite binary congruence can obstruct an odd count from being realized by an odd-order strong self-converse tournament, or a \(1\bmod4\) count from being realized by an even-order one. This is density in the corresponding 2-adic open sets, not exact realization of arbitrary integers or a positive natural-density claim.

## References and scope

- Swee Hong Chan and Igor Pak, [Computational complexity of counting coincidences](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf), §6.5 (the tournament spectrum conjecture, attributed there to bof and Gordon Royle).
- Darij Grinberg and Richard P. Stanley, [The Rédei–Berge symmetric function of a directed graph](https://arxiv.org/abs/2307.05569), for the odd-cycle collection evaluation (OCF).
- Leslie G. Valiant, [The complexity of computing the permanent](https://doi.org/10.1016/0304-3975(79)90044-6), Theoretical Computer Science 8 (1979), 189–201.
- A [2025 public question](https://cstheory.stackexchange.com/questions/55659/is-counting-hamiltonian-paths-in-a-tournament-p-complete) formulates the exact tournament counting problem (priority not inferred from that question).
- Viktória Zankó, [#P-Completeness via Many-One Reductions](https://doi.org/10.1142/S0129054191000066), International Journal of Foundations of Computer Science 2(1) (1991), 77–82.
- Repository inheritance: THM-3121 (path-cover walk-content substitution), THM-3134 (endpoint jets), THM-3166 (ordered-join cover transform); [Session 33 boundary observers](../session33_continuation_states_20261023/tournament/BOUNDARY_OBSERVABILITY.md), especially §3a for the permanent identity; [November 18 strong blue probes](../chanpak_self_converse_core_20261118/encoder/PROOF.md), including the single-probe modular decoder; [November 18 CRT](../chanpak_self_converse_core_20261118/build/CRT.md); [Session 197 dyadic self-converse C2](../../../session197/dyadic-self-converse-c2/SELF_CONVERSE_C2.md).
- Independent smoke test: the accompanying Python subset-DP check verifies the path and cycle extraction on all 530 binary square matrices of sizes 1–3, and the reversal-shell identity on 225 cases. The infinite claims rely on the proofs, not the smoke test.
