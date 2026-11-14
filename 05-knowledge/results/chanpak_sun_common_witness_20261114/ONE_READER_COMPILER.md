# One exact scalar reader for all path-cover channels

Sun's [finite-chart compilation](https://arxiv.org/html/2611.14605v1#S8) preserves one common witness while eliminating dependent choices. Here the **entire** path-cover profile of an unknown fixed-order tournament is encoded in the ordinary Hamiltonian-path count of one fixed reader applied to it. This packages equality constraints rather than proving existence of their shared source.

For a nonempty tournament \(A\) let \(G_j(A)\) count ordered \(j\)-tuples of nonempty paths partitioning its vertices. Write \(G_0(A)=0\). Concatenating an ordered cover and recording its cut positions injects into permutations with \(j-1\) cuts. Thus for source order \(n\), \(0\le G_j(A)\le B=n!2^{n-1}\). Put
\[
K=1+\lceil\log_2(72n(n+1)B)\rceil,\qquad L=(n+1)K,\qquad S_j=G_j(TT_L),
\]
where \(TT_L\) is transitive. Thus \(S_j\) is the number of surjections from \(L\) vertices to \(j\) named slots; set \(S_0=0\).

Use the actual cyclic join \(A\to TT_L\to TT_L\to A\). Compress maximal runs of a directed Hamiltonian path to block labels. The labels follow the directed three-cycle. A run-vector \((k,k,k)\) gives three possible starting labels, while each unequal vector of entries in \(\{k,k+1\}\) gives exactly one. Consequently the ordinary Hamiltonian count is
\[
H=\sum_{j=1}^n w_jG_j(A),\qquad
w_j=S_{j+1}^2+2S_{j+1}S_j+3S_j^2+2S_jS_{j-1}+S_{j-1}^2.
\]
In the \(G_j(A)\) coefficient the first three terms count patterns with minimum \(j\), and the last two count those with minimum \(j-1\). This is a positive physical expression: each basis support expands to actual path covers.

For \(2\le r\le n+1\), \((1-1/r)^r<1/2\) by the binomial theorem, hence \(r(1-1/r)^L\le r2^{-K}<1/2\). Discarding nonsurjections by a union bound gives \(S_r\ge r^L/2\). Thus \(w_j\ge (j+1)^{2L}/4\), while \(w_i\le9j^{2L}\) for \(i<j\). For \(j\ge2\),
\[
\frac{B\sum_{i<j}w_i}{w_j}\le36nB\left(\frac{j}{j+1}\right)^{2L}
\le36nB\,2^{-2K}<1.
\]
The inequality is immediate at \(j=1\) with an empty lower sum. Descending division recovers \(G_n\) from \(\lfloor H/w_n\rfloor\), subtracts it, recovers \(G_{n-1}\), and so on. The physical reader has \(n+2L=O(n^2\log(n+2))\) vertices, and the decoder uses exact integers.

**Exact common-source reduction.** A proposed vector \(v\in[0,B]^n\cap\mathbb Z^n\) has one common source tournament of order \(n\) if and only if \(\sum_jw_jv_j\) is realized by a tournament **in this fixed cyclic-reader chart**. The reverse implication is the superincreasing decoding proof; it is invalid on arbitrary tournaments outside the chart. This is not a proof of Chan–Pak or even of local completeness: the remaining challenge is existence in the chart, not compression or ambiguity of simultaneous channels.

## Dependent coordinates from parity

The [odd-cycle packing expansion](https://arxiv.org/src/2307.05569v1/anc/redeiberge-long.pdf#page=14), specialized to \(q\) equal colors, gives
\[
F_q(A)=\sum_e a_e q^e
\]
with nonnegative \(a_e\), \(e\equiv n\pmod 2\), and \(a_n=1\). In each packing, a cycle of odd length \(\ell\) reduces the atom count by \(\ell-1\); all vertices not on a cycle give separate atoms. The largest exponent \(n\) is therefore the empty packing. Nonempty \(A\) has no constant term. The profile \(G_j\) and these coefficients determine one another, since
\[
F_q(A)=\sum_j\binom{q}{j}G_j(A),\qquad G_j(A)=j!\sum_e a_e \left\{\begin{matrix} e\\ j\end{matrix}\right\}.
\]
The second equality follows by collecting powers in the falling factorial basis. Thus the single scalar reader actually determines at most \(\lfloor (n-1)/2\rfloor\) nontrivial polynomial coefficients; it does not posit \(n\) independent physical degrees of freedom. This is dependent-coordinate elimination, not a converse saying that every coefficient vector of the prescribed parity is realized.

## Unbounded levels force higher channels

At a strong substitution layer choose one distinguished variable block; the other blocks are fixed nonempty tournaments. Nest \(d\) such layers over the unknown deepest \(A\). For each layer \(i\), let \(q_i\) be the ordinary Hamiltonian count of its strong quotient, assumed to have at least three vertices. Let \(U\) be the product of the ordinary Hamiltonian counts of every fixed auxiliary block used at all levels.

Collapse every maximal consecutive run of the deepest \(A\) to a symbol; outside \(A\) retain actual vertices. The exact ordinary count is a positive linear combination \(\sum_j c_j G_j(A)\), with coefficients depending only on the chart. To give a lower bound, on each layer either use any of the \(q_i\) quotient Hamilton paths and keep the old word contiguous, or use one fixed spanning directed cycle through the marked quotient vertex. In the latter case insert its traversal of the auxiliary blocks between two copies of a chosen deepest-\(A\) occurrence. There are \(j\) such positions for an old word with \(j\) occurrences. In both cases use any full directed path in each auxiliary block. The supports are distinct because the outer auxiliary vertices either leave the old block contiguous or split it. Removing them, then merging the split consecutive pair if necessary, recovers the prior word, split position and every auxiliary path.

The lower-bound enumerator evolves by \(P(z)\mapsto q_iP(z)+z^2P'(z)\) starting at \(P(z)=z\), with the auxiliary factor \(U\) put back at the end. Since the steps are polynomials in the same differential operator, for \(0\le k\le d\) this yields
\[
c_{k+1}\ge U\,k!\,e_{d-k}(q_1,\ldots,q_d).
\]
In particular \(c_{d+1}>0\). Testing on transitive \(A=TT_n\) shows that this count cannot be a fixed linear combination of \(G_1,\ldots,G_d\) valid for all sizes: the left grows at least as a positive constant times \((d+1)^n\), whereas the right is \(O(d^n)\). This is a channel-growth result, not a no-go for nonlinear compression such as the fixed-order scalar encoding above.

The top forced index and coefficient are sharp: the prior [singleton triangular wrap](../chanpak_absolute_reflection_20261114/README.md) has \(F_q(W)=q^2F_q+q(F_{q+1}-F_{q-1})\), where \(F_q=\sum_j \binom qj G_j\). Iterating this at \(q=1\) gives top index \(d+1\) with coefficient \(d!\), and no higher slot.
