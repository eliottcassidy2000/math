# A locally complete clock has a severe absolute-height cost

The [five-corridor clock](FULL_CONGRUENCE_CLOCK.md), with corridor lengths
\(4,2,5,3,7\) repeated and arbitrary positive widths, covers every
locally compatible odd residue at every modulus. But its number of
distinct ordinary outputs up to \(X\) is
\[
                            o(X^{203/500}).
\]
This is a stronger estimate than the bound for the whole terrace.
Consequently, uniformly over an increasing sequence of moduli \(m\),
almost every compatible residue needs a clock count bigger than \(m^2\).
The largest of their least heights is
\(\omega(m^{500/203})\) (the ratio tends to infinity).
This quantitatively rules out the hope that a locally supplied clock
witness could always be brought within a bounded multiple of its modulus.

## Two width types suffice for the bound

For a word with positive widths \(r_0,\ldots,r_k\), let
\(p_i=2^{r_i}/(1+2^{r_i})\), and select its positive-gap bits
independently with those probabilities. Write \(K_q(p,p')\)
for the expectation of the corresponding entry of \(C_q\):
\[
 K_q(p,p')=e_q(1-p)(1-p')+
 f_q\bigl(p(1-p')+(1-p)p'\bigr)+g_q pp'.
\]
The kernel \(e_q\le f_q\le g_q\) is increasing in either
selection bit. Positive association of increasing functions
of independent bits (prove by induction and the one-bit
covariance formula) gives
\[
 H\ge
      \prod_{i=0}^k (1+2^{r_i})
      \prod_{i=1}^k K_{q_i}(p_{i-1},p_i).
\]
This does not assert independence of neighboring interactions;
it is a one-sided inequality.

Separate \(r=1\) from \(r\ge2\), with probability floors
\(p^{(0)}=2/3\), \(p^{(1)}=4/5\).
For \(s=203/500\) let
\[
 w_0=3^{-s},\qquad w_1=\sum_{r\ge2}(1+2^r)^{-s}.
\]
Sum the inverse powers over all width words using the two-state
matrices
\[
 B_q(i,j)=K_q(p^{(i)},p^{(j)})^{-s}w_j .
\]
The contribution from \(a\) whole repetitions of the clock
is bounded by
\((w_0,w_1)(B_4B_2B_5B_3B_7)^a(1,1)^t\).
The integer replay certifies
\[
 B_4 B_2 B_5 B_3 B_7
 \begin{pmatrix}1000\\991\end{pmatrix}
 \le \frac{99}{100}\begin{pmatrix}1000\\991\end{pmatrix}.
\]
It bounds all inverse powers from above by rational ceilings and
the \(w_1\) tail by powers of \(2^{-s}\). Therefore the sum of
\(H^{-s}\) over all words converges; later numerical collisions
can only lower the sum over distinct outputs. The little-o bound
follows by splitting the summable series at a fixed height.

## Translate sparsity into an inverse obstruction

For every \(m\) denote by \(h_m(r)\) the least *actual count* of a
clock witness in the compatible residue \(r\). This is finite by
simultaneous local completeness. The number of residues with
\(h_m(r)\le Y\) is at most the number of distinct clock counts
at most \(Y\), because one count can cover at most one residue.
There are \(m\) compatible residues for odd \(m\) and \(m/2\)
for even \(m\). At \(Y=m^2\) this gives proportion
\(o(m^{406/500-1})=o(m^{-0.188})\).
At \(Y=Cm^{500/203}\) for any fixed \(C>0\), the number is
\(o(m)\), which proves the claimed worst-height bound.

For Chan–Pak, these are not conjectural "gaps at infinity":
they say precisely why a recipe that matches any fixed number of
digits is too weak to supply all exact targets even with a
generous height allowance. A genuine descent must use moves
whose absolute output can be chosen or reduced, not just
arbitrary modular inverses implemented by long positive words.
