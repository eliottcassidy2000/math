# An affine point-count lift of the backward-key tournament strip

**PROVED.** For the [physical backward-key strip](../chanpak_tournament_sl2_20261124/ZERO-GAP-PORT.md) there is a concrete affine-stratified \(\mathbb Z\)-scheme with count polynomial \(F_T(q)\) and \(F_T(2)=H(T)\). Its absolute specialization has an exact Fibonacci interface, but loses ordinary heights even at fixed tournament order. Conversely a fixed tournament and even its exact binary point set do not specify a unique absolute limit.

## Cells that recover the ordinary height at two

There are \(g_i\) transitive spacers between consecutive keys, with the one backward arc between the keys of each gap. In the odd-cycle expansion select a set of gaps and split it into its maximal consecutive runs. An \(r\)-gap run with \(G>0\) spacers has \(2^G\) realizations with decoration, by choosing an ascending spacer subset of parity \(r\bmod2\) and one binary decoration. With no spacers the weight is two for even \(r\) and zero for odd \(r\).

For \(G>0\), replace binary spacer choices by the integral affine hyperplane \(\sum_{i=1}^G z_i=[r\text{ odd}]\) and keep one free affine coordinate for the decoration. Over \(\mathbb F_2\) this is the parity rule verbatim; over any finite field with \(q\) elements the cell has \(q^G\) points. A zero-spacer even run uses an affine line, and an odd run is omitted. Take the product over runs and the disjoint union over admissible supports (including the empty point). The \(\mathbb F_2\) points are **decorated cycle collections**, whose cardinality is \(H(T)\) by the [Rédei–Berge expansion](https://arxiv.org/abs/2307.05569); no direct identification of individual Hamilton words is required.

For \(F_T(q)=\sum_jc_jq^j\), the usual finite-field zeta is \(\prod_j(1-q^j u)^{-c_j}\). The [Soulé limit](https://arxiv.org/abs/math/0304444) in this model is therefore
\[
 \lim_{q\to1^+}(q-1)^{F_T(1)}\prod_j(1-q^{j-s})^{-c_j}
   =\prod_j(s-j)^{-c_j}.
\]
First take real \(s>\deg F_T\), then continue the resulting rational function. This is a concrete count-polynomial specialization, not an assertion of an ordinary one-element field or a canonical arithmetic curve.

## Exact Fibonacci and a moving unimodular interface

In an interval of \(t\) zero-spacer gaps, write \(f_t(q)\) for the weight with neither neighboring positive gap selected, \(g_t(q)\) with one selected, and \(h_t(q)\) with both selected. A zero portion attached to a selected positive gap costs nothing extra. Splitting at the first unselected gap gives \(g_t=1+\sum_{i<t}f_i\) and \(h_t=1+\sum_{i<t}g_i\). The two-port matrix is \(K_t(q)=\left(\begin{smallmatrix}f_t&g_t\\g_t&h_t\end{smallmatrix}\right)\). For the free word, a maximal selected run is even and weighted \(q\), so
\[
 \sum_{t\ge0} f_t(q)z^t=
 \frac{1+(q-1)z^2}{1-z-z^2+(1-q)z^3}.
\]
Indeed begin with an unselected gap, or a positive even run followed by an unselected gap and a remainder, or a positive even run ending the word.

At \(q=1\) every allowed run has weight one; pair its selected gaps consecutively into dominoes. The other gaps are monomers. This gives \(f_t(1)=F_{t+1}\) with \(F_0=0,F_1=1\). Then
\[
 K_t(1)=\begin{pmatrix}F_{t+1}&F_{t+2}\\F_{t+2}&F_{t+3}-1\end{pmatrix},
 \qquad \det K_t(1)=(-1)^t-F_{t+1}.
\]
The determinant is Cassini minus the missing term. Exactly \(t=2\) is unimodular at the absolute specialization; positive-spacer multipliers \(\operatorname{diag}(1,q^g)\) specialize to the identity, so this port composes unimodularly across bridges. In fact \(K_t(1)+E_{22}=\left(\begin{smallmatrix}0&1\\1&1\end{smallmatrix}\right)^{t+2}\): one missing point in the both-selected read separates it from a pure Fibonacci power. Adding that point is a modification of the model, not a tournament edit.

In the ordinary \(q=2\) reading, the bridge has even determinant, and the bare **three**-gap port is instead unimodular. Over the whole count polynomial,
\[
 K_3(q)=\begin{pmatrix}1+2q&4+q\\4+q&7\end{pmatrix},
 \qquad \det K_3(q)=-(q-3)^2.
\]
Thus it is singular at field size three, unimodular at two and four, and not unimodular at the absolute specialization. A symmetry or inverse that exists in one reading need not exist in another. In particular the [physical mixed-modulus compiler](../chanpak_tournament_sl2_20261124/STRIP-CRT.md) at \(q=2\) is not an absolute inverse obtained by setting \(q=1\).

## Loss in both directions

1. **Same \(F_T(1)\), different ordinary heights, even at fixed order.** If all gaps have positive spacers then \(F_T(q)=\prod_i(1+q^{g_i})\). Fix total spacers \(S\ge2\) and take two gaps \((a,S-a)\) for \(1\le a\le\lfloor S/2\rfloor\). The hosts are strong, of common order \(S+3\), absolute count 4 and degree \(S\); the ordinary heights \(1+2^S+2^a+2^{S-a}\) are strictly decreasing in this range. Already \(S=4\) gives 27 and 25. The full polynomial retains the multiset of widths: the largest \(g\) contributes \(\Phi_{2g}\) in each copy of \(1+q^g\) and no smaller width contributes that factor. Strip those copies and repeat.
2. **Same tournament, different absolute models.** One gap with \(g\ge2\) spacers has the compressed lift \(1+q^g\). Alternatively decorate each of the \(2^{g-1}\) actual odd cycles by its own affine line, obtaining \(1+2^{g-1}q\). They agree at two and give different absolute zetas. More strongly, replace a positive-dimensional cell \(\mathbb A^j\) by two disjoint copies of \(\mathbb A^{j-1}\), taking the first binary coordinate as the label. On \(\mathbb F_2\) this preserves the finite set by explicit bijection; as schemes over the integers the replacement is not an isomorphism. Each split increases \(F_T(1)\) by one. Continuing until every cell is a point shows that the *same* tournament admits such models with all absolute counts from \(2\) to \(1+2^g\).

The robust interface is the full lift with its stated specialization, not an invariant assigned from the ordinary height alone. Compare [the prime-power Riesz interface](RH-F1-READING.md), where preserving a transform is likewise different from proving an upper bound, and [the exact finite spectral test](FAREY-SPECTRUM.md), where some graph zetas have a square-root bound and another provably does not.

Replay: run strip-f1-check.py with Python 3.
