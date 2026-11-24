# Literal strong tournament strips realize every finite clock

**PROVED.** For every modulus \(N\) and every residue class admitting an odd representative, there is a strong tournament from the one [backward-key strip architecture](ZERO-GAP-PORT.md) whose **ordinary** Hamilton-path count lies in that class. Modulo each odd \(m\), actual positive strip continuations implement both elementary \(\mathrm{SL}_2\) shears on the two conditioned count coordinates. For mixed moduli the same shears can be implemented with no drift in the already-set \(2\)-primary height. This does not give any prescribed absolute integer.

## Literal elementary shears at the odd clocks

Write \(D=\operatorname{diag}(1,2)\), and use the conditioned gap matrices
\[
 K_1=\begin{pmatrix}1&2\\2&2\end{pmatrix},\quad
 K_3=\begin{pmatrix}5&6\\6&7\end{pmatrix},\quad
 K_4=\begin{pmatrix}9&11\\11&13\end{pmatrix}.
\]
Their determinants, \(-2,-1,-4\), are units modulo every odd \(m\). An actual continuation at a positive-spacer marker appends \(t\) zero gaps followed by a positive gap with \(g\) spacers. It multiplies the conditioned row (state 0, state 1) on the right by \(K_t D^g\). The matrix \(D^g\) charges the new gap's weight \(2^g\) precisely when selected. Start with a positive gap giving row \((1,2^g)\); finish with column \((1,1)^{\mathsf T}\), so there is no lost condition at the ends.

Over the dyadic rationals,
\[
 U=\begin{pmatrix}1&1\\0&1\end{pmatrix}
    =K_4^{-1}K_3K_1^{-1}DK_1D,\qquad
 L=U^{\mathsf T}=DK_1DK_1^{-1}K_3K_4^{-1}.
\]
For example \(K_1^{-1}DK_1D=\left(\begin{smallmatrix}3&4\\-1&0\end{smallmatrix}\right)\); multiplying on the left by \(K_3\) gives \(K_4U\). Modulo odd \(m\), choose gap widths congruent to 1 or 2 modulo \(\operatorname{ord}_m(2)\), as needed; they can be arbitrarily large. The four physical steps with residues \(K_1D,K_1D^2,K_3D,K_4D\) are invertible in a finite group. Their inverses are positive powers of the **same physical steps** in that group. Thus \(D=(K_1D)^{-1}(K_1D^2)\) and \(K_t=(K_tD)D^{-1}\) are compiled by positive gap words, and so are \(U,L\). The result holds for arbitrary odd composite \(m\), not just fields. It does not implement their signed inverse matrices over the integers.

## Separating the dyadic clock by commutators

Write \(N=2^a m\) with \(m\) odd. For \(a=0\) skip the dyadic part; \(m=1\) needs only that part. Otherwise choose all appended widths \(g\ge\max(a,1)\) so the second coordinate is 0 modulo \(2^a\) after the step. On rows already of this form \(K_tD^g\) acts on the first coordinate by \((K_t)_{00}\), respectively \(1,5,9\) for \(t=1,3,4\).

To set any desired odd dyadic residue, first choose one of two starts. A single high-width positive gap gives \(x=1,y=0\bmod2^a\). Or begin with one spacer, row \((1,2)\), and append \(K_4D^g\) with high \(g\), giving \(x=31,y=0\bmod2^a\). These provide both signs modulo four. Repeating a high-width \(K_3D^g\) multiplies \(x\) by 5, which generates \(1+4\mathbb Z/2^a\) (by \(v_2(5^{2^j}-1)=j+2\)). This chooses the dyadic residue exactly, including the small exponents. All these steps remain invertible at odd \(m\): in particular the odd row is unimodular.

Track each high-width physical step in the finite product
\[
 (\mathbb Z/2^a)^\times\times\mathrm{GL}_2(\mathbb Z/m),
 \qquad K_tD^g\longmapsto((K_t)_{00},\,K_tD^g\bmod m).
\]
Each inverse is again a positive physical power in this product, so simultaneous equalities are respected by the compiler. Its odd projection includes \(D,U,L\) by the identities above. Take any lifts of these three elements. A commutator of lifts has first coordinate 1 because the first factor is abelian, while
\[
 D^{-1}UDU^{-1}=U,\qquad DLD^{-1}L^{-1}=L.
\]
Consequently both shears have **dyadically neutral positive physical words**, and generate \(\{1\}\times\mathrm{SL}_2(\mathbb Z/m)\). To see generation over the ring, use the unit pivot in the first column at each prime power, eliminate using row shears, and glue the shear parameters by CRT. Since this group is transitive on unimodular rows (complete the chosen row to determinant one by Bézout), advance the odd row to \((1,r-1)\bmod m\), whose sum is \(r\). This leaves the desired dyadic residue unchanged. The closing column therefore produces the requested count modulo \(N\).

The strip is strong: all keys reach backwards along their descending chain, and the first reaches the last by a forward shortcut when there are at least three keys. Each spacer lies on a forward two-edge path from first to last. A single positive gap is already a directed cycle closed through its key pair. Thus no separate ordinal padding or discarded path is needed. There is no claim of a small order bound, and the modular inverses used here do not give absolute integer descent: every added \(K_tD^g\) still has an even integer determinant.

The source of the finite classical symmetry is concrete here. The congruence action comes from composing actual conditioned reads of ordinary Hamilton paths and includes \(\mathrm{PSL}_2(\mathbb F_p)\) upon projecting for each odd prime \(p\ge5\); the Golay Mathieu action instead permutes a finite set of marked contacts. Neither is an action by automorphisms of one tournament.

Replay: run strip-crt-check.py with Python 3.
