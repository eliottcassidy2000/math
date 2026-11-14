# Density-one divisibility, subpolynomial exact spectrum

Use all the Fermat-indexed source components \(R_{2^j}\), \(j\ge1\), from
[the unbounded-rank bank](FERMAT_RANK.md). Glue a source with \(a_j\)
copies of these components and \(s\) singleton SCCs into the same cyclic
wrap. Write \(A(X)\) for the number of distinct resulting Hamiltonian path
counts at most \(X\). Despite density-one prime-divisor coverage, the exact
spectrum has the much smaller size
\[
\log A(X)=\frac{(\log\log X)^2}{2\log 2}
 +O\bigl(\log\log X\,\log\log\log X\bigr)
\]
as \(X\) tends to infinity (all logs in this display are natural).

For the upper bound write \(L=\lfloor\log_2 X\rfloor\). The formula
\(N=h+2^{s+1}\prod Q_{2^j}^{a_j}\) and \(Q_{2^j}\ge2^{2^j}\) imply
\(s+\sum_j a_j2^j\le L\). Thus \(s\) has at most \(L+1\) choices; the
coefficient vectors number at most
\(\prod_{j=1}^{\lfloor\log_2 L\rfloor}(1+\lfloor L/2^j\rfloor)\).
The logarithm of this product is at most
\((\log L)^2/(2\log2)+O(\log L)\), by summing
\(\log L-j\log2+\log 2\). Different orders of the SCCs do not change \(N\).

For the lower bound take an integer \(W\ge2\). Use any vector satisfying
\(\sum_j a_j2^j\le W\) and fix \(s=2W\). Since \(1+2^m<2^{3m/2}\)
and \(Q_m\le2^{3m}\) for \(m\ge2\), these counts obey \(N\le2^{5W+2}\)
and \(h<2^{2W}\). The residue \(N\bmod 2^{2W}\) recovers \(h\) exactly.
Distinct vectors have distinct \(h\), because \(h=\prod_j F_j^{a_j}\)
and the Fermat numbers \(F_j>1\) are pairwise coprime. This is a genuine
separation of absolute values, not of formal representations.

To count enough vectors, take \(J=\lfloor\log_2 W\rfloor\) and allow
\(0\le a_j\le\lfloor W/(J2^j)\rfloor\) independently for
\(1\le j\le\lfloor\log_2(W/J)\rfloor\), with the other \(a_j=0\).
This fits the budget and gives logarithmic cardinality at least
\[
\sum_{j=1}^{\lfloor\log_2(W/J)\rfloor}(\log(W/J)-j\log2)
=\frac{(\log W)^2}{2\log2}-O(\log W\,\log\log W).
\]
Take \(W=\lfloor(\log_2 X-2)/5\rfloor\) and combine the two estimates.

Consequently the bank misses nearly every exact odd target while giving
a physical divisible count for nearly every prime. Infinite quadratic
rank is necessary for this style of all-target proof, but it is nowhere
near sufficient without control of the absolute output.
