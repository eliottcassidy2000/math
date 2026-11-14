# Simultaneous divisibility for almost every fixed-length prime tuple

This extends the one-prime result for the Fermat-indexed cyclic-source
bank. Work with ordered tuples of distinct odd primes \(p_i\le X\),
with fixed tuple length \(t\). As \(X\to\infty\), the proportion for which
there is **one** physical Hamiltonian path count divisible by
\(\prod_{i=1}^t p_i\) tends to one. This is a prime-box density statement,
not a statement about all tuples or about the exact spectrum.

Choose \(r\ge t\) multiplicatively independent square classes of ratios
\(Q_j/h_j\) from the bank, and for this finite selection ignore the finitely
many primes dividing their numerator or denominator profiles. For a
tuple form a binary \(t\)-by-\(r\) matrix whose \(i,j\) coordinate
is 0 when \(Q_j/h_j\) is a square modulo \(p_i\) and 1 otherwise.

**Full row rank suffices.** In the product of the multiplicative
groups modulo the \(p_i\), let \(P\) be the product of their Sylow
2-subgroups. Project the generated ratio subgroup to \(P\), with image
\(B\) in additive notation. The binary matrix records the images of the
generators in \(P/2P\). Full row rank gives \(P=B+2P\), hence
\(P=B+2^aP\) for all \(a\) by iteration, so \(P=B\).
There is therefore a product of powers of the ratios whose 2-primary
projection is \((-1,\ldots,-1)\). Raise it to an odd common multiple
of all the odd parts of \(p_i-1\) to kill the odd-primary coordinates
without altering this projected involution. The resulting product is
exactly \((-1,\ldots,-1)\). Replace negative powers by nonnegative powers
using each generator's order modulo the whole tuple, and choose the
padding exponent \(s+1\) a positive multiple of the order of 2 modulo the
tuple, with \(s\ge1\). The actual wrapped count \(h+2Q\) is divisible by
all the \(p_i\) simultaneously. This also works for any separately chosen
exponents \(p_i^{e_i}\), since passing from the prime to the prime-power
unit group adds only odd-order kernel.

**Rank is asymptotically typical.** The independence of the \(r\)
square classes gives \(r\) independent quadratic characters modulo a
common conductor. Uniform prime distribution among its reduced classes
therefore makes all binary sign rows asymptotically equally frequent.
Ordered \(t\)-tuples of primes in the box give asymptotically independent
rows; discarding repeated primes and the finitely many bad primes has
vanishing relative effect. The proportion of \(t\)-by-\(r\) uniform
binary matrices with row rank \(t\) is
\[
\prod_{a=0}^{t-1}(1-2^{a-r}).
\]
For any \(r\) this is a lower bound on the limiting inferior success
proportion. Let \(r\to\infty\) at fixed \(t\).

A fixed palette still has an invisible arithmetic progression of primes,
and here the palette remains sparse in absolute output despite this
simultaneous divisibility result. There is no implication that every odd
modulus is hit, no uniformity as \(t\to\infty\) is asserted, and the
size of the output is not constrained by the prescribed product.
