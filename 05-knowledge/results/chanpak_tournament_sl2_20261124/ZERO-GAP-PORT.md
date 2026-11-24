# A genuine unimodular interaction in an unweighted tournament

**PROVED.** A backward-key strip gives a rank-two conditioned transfer with determinant \(-1\) that counts ordinary Hamilton paths through the Rédei–Berge expansion. Concatenating its bare interfaces is not an \(\mathrm{SL}_2(\mathbb Z)\) action: physical spacer weights insert even determinants.

Place keys \(v_0,\ldots,v_m\) in a linear order, with \(g_i\) intervening spacers between consecutive keys. Orient all pairs forward except the \(m\) arcs \(v_{i+1}\to v_i\). The cut just after key \(v_i\) has one backward arc, from the next key. Every simple cycle crossing a cut must use that backward arc; crossing all cuts between its extreme keys forces one descending interval of keys. It then returns along an increasing subset of spacers in that interval, with no more keys available. If it crosses only one cut, the return must use at least one spacer (the direct forward arc was reversed).

The [Rédei–Berge expansion](https://arxiv.org/abs/2307.05569) sums over vertex-disjoint odd cycles, with each nontrivial odd cycle contributing a factor two. Choose a subset of the \(m\) indexed gaps. Each maximal selected run is one cycle, and distinct runs have disjoint key sets and spacers. For a run through \(r\) gaps and \(G\) spacers its weight is \(2^G\) if \(G>0\): exactly half the spacer subsets give odd length. If \(G=0\), its weight is 2 for positive even \(r\) and 0 for odd \(r\). This is an ordinary, all-path count, not a framed traversal restriction.

## Two ports, with a composition tax

For each positive-spacer gap let \(b=2^g\) and mark its state 0 or 1 according as it is unselected or selected. Between two positive markers a block of \(t\) zero-spacer gaps has conditioned interaction
\[
 K_t=\begin{pmatrix}f_t&g_t\\g_t&h_t\end{pmatrix},\quad
 g_t=1+\sum_{i<t}f_i,\quad h_t=1+\sum_{i<t}g_i.
\]
Here \(f_t\) is the free zero-only count (a selected run contributes 2 if its length is positive even and 0 otherwise). With one selected neighbor the adjoining run can have any length; split at the first zero, or take the all-selected word. With both neighbors selected, the first zero leaves a one-ended suffix. The resulting initial values are
\[
 K_1=\begin{pmatrix}1&2\\2&2\end{pmatrix},\quad
 K_2=\begin{pmatrix}3&3\\3&4\end{pmatrix},\quad
 K_3=\begin{pmatrix}5&6\\6&7\end{pmatrix}.
\]
The bare \(K_3\) is unimodular with determinant \(-1\). On the index-two lattice \(x\equiv y\bmod2\), a basis is \(v=(1,1)\) and \(K_3v=(11,13)\). The identity \(K_3^2=12K_3+I\) makes the matrix in this basis \(\left(\begin{smallmatrix}0&1\\1&12\end{smallmatrix}\right)\), the continued-fraction digit 12. It is not conjugate to that companion matrix on the full integer lattice: \(\det[v,K_3v]\) is always even for integer \(v\), so no cyclic basis of this form is unimodular.

**Composition qualification.** Each intermediate positive marker contributes its own physical \(\operatorname{diag}(1,b)\) with even \(b\). Thus even a unimodular bare port cannot generally be composed unimodularly on the native lattice. Positive matrices and formal inverses are different things.

## Exact full-strip dynamics

Let \(A\) mean no open run, \(O,E\) mean an open zero-spacer run of odd/even length, and \(P\) mean an open run containing a positive-spacer gap. Pending \(O,E\) get their factor only on closure. At a zero gap the update is
\[
 (A,O,E,P)\mapsto(A+2E+P,A+E,O,P).
\]
At a positive gap, of weight \(b=2^g\), it is
\[
 (A,O,E,P)\mapsto(A+2E+P,0,0,b(A+O+E+P)).
\]
Start with \((1,0,0,0)\) and finish with \(A+2E+P\). The zero update on the first three coordinates has determinant one; a positive gap resets to rank two. The replay compares the formulas with ordinary Hamilton-path subset DP over gap words of lengths up to six, spacings in \(\{0,1,2\}\), and tournament order at most ten.

Replay: run zero-gap-port-check.py with Python 3.
