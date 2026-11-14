# Bounded source component size cannot yield a single-seam cutoff

Fix a grade limit \(b\), meaning every strongly connected component in the source has order at most \(b\), with no bound on the total order. For every \(j\ge1\), glue the ordered near-transitive block \(R_j\) at a marked root, not necessarily a Fermat exponent. The number of possible outputs through \(X\) is only a power of \(\log X\) whose exponent depends on \(b\). Thus even unbounded neutral growth, or stacking finite strong alphabets with arbitrary multiplicities, cannot yield a cofinite theorem for this seam.

Decompose the source at the root component \(C\) as \(P\to C\to D\). The inherited [low-edge calculation](../session328_fermat_gluing_and_phase_repair_20261114/FERMAT_LOW_EDGE.md) gives
\[
N=H(P)\bigl[H(C)H(D)+2^j L(C,v)Q_D(2)\bigr].
\]
All terms are positive. Counts for an ordinal join multiply in each ordered slot (including empty slots). Thus \(H(P),H(D),Q_D(2)\) are each at most \(X\) and have no prime factors above \((b+1)!\): the \(Q(2)\) of an order-\(m\) component is at most \((m+1)!\), because a two-slot cover is a permutation with one weak cut, and \(H\le Q(2)\). Each prime exponent is at most \(\lfloor\log_2 X\rfloor\). There are at most \(\sum_{m\le b}m2^{\binom m2}\le b^2 2^{\binom b2}\) choices of rooted \(C\) and at most \(\log_2 X\) choices of \(j\). An explicit upper bound on the image size is therefore
\[
b^2 2^{\binom b2} (\log_2 X)(1+\lfloor\log_2 X\rfloor)^{3\pi((b+1)!)}.
\]
This counts possible coordinates independently, which overcounts executable joined profiles. It is still sufficient. In particular for \(b\le (1-\epsilon)\log\log X/\log\log\log X\), fixed \(0<\epsilon<1\), Stirling's estimate gives \((b+1)!\le(\log X)^{1-\epsilon/2}\) eventually, and the displayed bound is \(X^{o(1)}\). Thus almost every integer up to \(X\) would need higher grade if this single-seam method supplied it at all.

This obstruction is about the maximum order of a single component, not the total order of the source. The two-port return theorem in the same folder is logically separate: iterated gluing may build one unbounded strong component, but those particular return ports still cannot cover primes.
