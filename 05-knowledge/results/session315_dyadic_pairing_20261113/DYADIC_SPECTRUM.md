# Dyadic toothpick inflations: exact residue coverage and an absolute sparsity barrier

## 1. Fixed twin templates are too thin for the full spectrum

Let \(R\) be a fixed quotient tournament. Replace each vertex \(v\) by a transitive tournament \(TT_{L_v}\) of arbitrary positive order, directing every inter-module edge according to \(R\). Call a site *active* if its vertex belongs to a nonsingleton strong component of \(R\). At any inactive site the whole twin run is visited consecutively by every Hamiltonian path, because the strongly connected components are ordered and no path can return to a component after leaving it. Its length has no effect on the count.

At any active site \(v\), there is a directed triangle \(v\to a\to b\to v\). An exterior vertex in a strong tournament lies in a directed triangle: otherwise its out-neighbors have no edge to its in-neighbors, contradicting strongness. Replace \(v\) temporarily by two consecutive transitive twins \(v_1,v_2\). The sequence \(v_1,a,b,v_2\) is a directed path. Insert all remaining quotient vertices one by one into this path (any directed path in a tournament can be extended by insertion of one new vertex); the two appearances of \(v\) stay separated. In the final inflation, assign any nonempty proper subset of the \(L_v\) transitive vertices to the first \(v\) slot and its complement to the second. All other modules can be traversed consecutively. This gives at least \(2^{L_v}-2\) distinct Hamiltonian paths. In particular every active length is at most \(\log_2(X+2)\) in an output of count at most \(X\).

Hence any finite bank of fixed quotient templates with arbitrary independent transitive toothpick lengths realizes only \(O((1+\log X)^r)\) distinct counts at most \(X\), where \(r\) is the maximum number of active sites in a template. This does not restrict changing quotients, growing numbers of active sites, or nontransitive variable modules. If every active \(L_v\ge2\), one may split them simultaneously: traverse a Hamiltonian cycle of each nonsingleton quotient SCC twice (opening it to a path), then concatenate the SCC paths. This gives the optional stronger lower bound \(\prod_{\text{active }v}(2^{L_v}-2)\).

## 2. Exact one-root visit formula

Fix all sites but \(v\) at order one. For each \(j\), let \(a_j\) be the number of reduced path itineraries with exactly \(j\) nonempty runs of \(v\). The other vertices occur exactly once and \(v\)-slots are not adjacent. Then the count after replacing \(v\) by \(TT_L\) is
\[
H_L=\sum_{j\ge1}a_j\,j!\,\left\{\begin{matrix}L\\j\end{matrix}\right\}.
\]
The coefficient counts onto assignments of labelled twin vertices to the ordered slots; the traversal within each assigned set is forced. A \(j\)-slot itinerary forces a \((j-1)\)-matching from the out-neighbors to the in-neighbors of the root: each interior excursion has at least one such crossing and excursions are disjoint. Conversely a matching gives disjoint two-vertex excursions, and the unused in-neighbors and out-neighbors may be traversed as single prefix and suffix paths, as in [Session 303](../session303_ternary_normal_form_20261113/README.md). Hence maximal visits equal one plus the excursion matching number. Permuting the \(j-1\) internal excursions leaves the itinerary valid and acts freely; in particular \(a_3\) is even.

If the matching number is two, exactly the first three coefficients may occur. Put \(A=a_1-2a_2+3a_3\), \(B=a_2-3a_3\), \(C=a_3\). Then
\[
H_L=A+B\,2^L+C\,3^L .
\]
These are exact absolute counts at every length \(L\ge1\), not just formal expressions modulo powers of two.

## 3. Four strong templates cover all odd dyadic residues

Use the following outgoing bit rows, with vertex 0 the root. The excursion matching number is exactly two in each graph: one side has size two, and direct inspection gives a two-matching. Each quotient is strong, as certified by forward and backward reachability at 0. A uniform replacement of a vertex in a nonsingleton strong quotient remains strong.

| Rooted quotient outgoing rows | \((a_1,a_2,a_3)\) | \((A,B,C)\) | Odd residue classes mod 16 reached at large length |
|---|---|---|---|
| (20,1,2,7,14) | (9,10,2) | (−5,4,2) | 1, 13 |
| (44,17,2,6,13,30) | (25,28,6) | (−13,10,6) | 5, 9 |
| (36,1,2,7,15,30) | (15,22,6) | (−11,4,6) | 7, 11 |
| (34,48,3,7,13,28) | (41,41,10) | (−11,11,10) | 3, 15 |

One can recover the itinerary counts from direct path counts at \(L=1,2,3\): these are respectively \(a_1\), \(a_1+2a_2\), \(a_1+6a_2+6a_3\). The standalone replay checks the values and also the literal expansions at larger lengths.

For \(k\ge4\) and \(L\ge k\), the term \(B2^L\) vanishes modulo \(2^k\). Here \(C=2c\) with \(c\) odd in every row. The powers of 3 modulo \(2^{k-1}\) are exactly the units that are 1 or 3 modulo 8: their residues are in those classes and the order is \(2^{k-3}\), by the identity \(v_2(3^{2^t}-1)=t+2\) for \(t\ge1\) (successive squaring). Thus \(A+2c3^L\) reaches every lift modulo \(2^k\) of the two displayed classes \(A+C,A+3C\bmod16\). Choose the row by the target's last four bits, solve the finite discrete logarithm in the subgroup, and take a sufficiently large exponent in its class. Smaller powers follow by reduction. This gives every odd dyadic residue by literal strong tournaments, with only a choice among four templates and then one variable length.

Four templates are necessary within any bank of this specific single-variable, at-most-three-visit type if it covers all odd dyadic residues at arbitrarily high precision. Indeed \(C=a_3\) is even. For large \(L\) any such chart can supply at most two residues mod 16: with \(v_2(C)=1\) use the two classes of \(3^L\bmod8\), with \(v_2(C)=2\) use at most two classes mod 4, and larger valuation gives only one. Fewer than four fixed charts leave some odd class mod 16 absent at all \(L\ge4\); the finitely many smaller lengths cannot supply its unbounded number of dyadic lifts. The statement concerns this restricted grammar, not arbitrary toothpick banks.

These families give a dyadic analogue of the earlier 11- and 19-adic length bypasses: no site-by-site target selection is needed for modular reachability. The actual spectra remain sparse. In fact each row contributes at most logarithmically many absolute counts through \(X\), by the active-site bound. The exact guard for an absolute target is the original equation \(H_L=A+B2^L+C3^L\), not its residues at a chosen finite precision.

## 4. A stronger modular lift: every odd residue among strong self-converse tournaments

A single fixed-root self-converse chart with at most three visits has a genuine mod-4 restriction. Suppose its anti-involution fixes the root, and extend it by reversing the transitive order of the replacement run. The quotient has odd order because an anti-involution fixes at most one vertex. For even run length the expanded tournament has even order, so the extended anti-involution has no fixed vertex. It pairs all directed odd cycles: an invariant odd cycle would induce an order-two reflection of that directed cycle, with one fixed vertex. By the odd-cycle packing expansion \(H\equiv1+2c_{\rm odd}\pmod4\), where \(c_{\rm odd}\) counts individual directed odd cycles. Hence \(H\equiv1\pmod4\) for every even run length. In fact this congruence holds for every even-order self-converse tournament: its automorphism group has odd order, so if an anti-isomorphism has square of odd order $u$, its $u$-th power is an anti-involution, to which the same fixed-cycle argument applies. For at most three visits \(C=a_3\) is even, so \(A+B2^L+C3^L\) is constant modulo 4 for **all** \(L\ge2\). This grammar cannot supply arbitrary lifts of \(3\bmod4\).

A two-block mirror construction bypasses the restriction. Given a target odd residue \(c\bmod2^k\), with \(k\ge3\), choose a fixed self-converse tournament \(Q\) with \(H(Q)\equiv c\bmod8\). For the four possibilities use the singleton, directed triangle, unique strong order-four tournament, and regular cyclic order-five tournament, with counts \(1,3,5,15\). Every unit that is 1 modulo 8 is a square modulo \(2^k\): a root lifts from precision \(j\ge3\) by optionally adding \(2^{j-1}\), which flips precisely the next square bit. Choose \(T\) from the four-chart bank so that \(H(T)\) is a square root of \(c/H(Q)\) modulo \(2^k\). Form the ordinal sum
\[
S=T\ \to\ TT_r\ \to\ Q\ \to\ TT_r\ \to\ T^{\mathrm{op}}.
\]
It is self-converse by reflection of the block order. Write \(Q_U(t)\) for the ordered \(t\)-slot path-cover polynomial, including empty slots. In an ordinal join a path in each slot visits the components in their fixed order; hence both \(H\) and \(Q(2)\) multiply. Thus
\[
H(S)=H(Q)H(T)^2,\qquad
Q_S(2)=Q_Q(2)Q_T(2)^2\,2^{2r}.
\]

Finally take the cyclic join of \(S\) and two singletons. It is strong and self-converse (apply the anti-isomorphism on \(S\) and exchange the singleton blocks). The two singletons allow exactly three one-run itineraries in \(S\) and one two-run itinerary. The exact count is therefore \(3H(S)+[Q_S(2)-2H(S)]=H(S)+Q_S(2)\). Choosing \(r\ge\lceil k/2\rceil\) annihilates the latter term modulo \(2^k\), giving the prescribed odd residue. Lower precisions follow by reduction. The result concerns modular coverage, not arbitrary absolute counts; this wrap adds a real positive term that can be enormous even when zero modulo the target precision.

The original four-chart selector already needs long runs at high precision. Covering all \(2^{k-1}\) odd residues using at most four outputs at each length forces a maximal length of at least \(2^{k-3}\). The constructive upper bound is \(k+2^{k-3}-1\) using one full period per chart. Some residue therefore requires a count doubly exponential in the precision within this bank, by \(H\ge2^L-2\).

Within the displayed strong self-converse wrapper bank, the absolute spectrum is still only \(O((1+\log X)^2)\) through \(X\). The underlying template length is \(O(\log X)\) by the active-site lower bound on \(H(T)\), the pad length is \(O(\log X)\) by the positive \(4^r\) term, and all other choices belong to fixed finite palettes. Thus even imposing strongness and self-converseness does not make exact modular coverage approach absolute coverage.
