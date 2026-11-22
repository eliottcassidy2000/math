# Two barriers to a full Liouville proof on the fixed +2 word law

**PROVED (method obstructions); cancellation of the full signs OPEN.**
These statements apply to the actual words `AB/BA` with the fixed finish in
[Fixed-gap sieve](../chanpak_prime_parity_inverse_20261122/FIXED-GAP-SIEVE.md).
Writing the choices as bits \(d_i\), the literal endpoints are
\[
N_\epsilon(\mathbf d)=B_n+\epsilon+D\sum_{i=0}^{n-1}a^i d_i,\quad
a=490032,\quad D=2424\cdot63618252914778017089536,\quad \epsilon=0,2.
\]
Here \(B_n\) is the lower count when all block increments are 360532.
The endpoint sets are at height \(X\asymp a^n\); each has precisely \(2^n\)
elements. All values at both endpoints are coprime to \(a\).

## 1. Quantifying the Fourier loss

Let \(e(t)=\exp(2\pi it)\). For the probability measure on the \(2^n\)
lower endpoints, the \(L^1(\mathbb R/\mathbb Z)\) norm of its transform is
\[
 I_n=\int_0^1\prod_{i=0}^{n-1} f(a^it)\,dt,\qquad f(t)=|\cos(\pi t)|.
\]
The translation contributes a unit phase, each choice contributes one cosine,
and multiplication by the nonzero integer \(D\) preserves Haar measure.
This is the same norm for the upper endpoints. Set \(m=2/\pi\).
For every \(t\in[0,1]\) the full inverse-branch average satisfies
\[
 \left|\frac1a\sum_{j=0}^{a-1}f((j+t)/a)-m\right|\le \frac{\pi}{2a}.
\]
Compare each sample to the integral over its own cell. Because \(f\) is
\(\pi\)-Lipschitz, the error after summing is bounded by
\(\pi (t^2+(1-t)^2)/(2a)\). Change variables on each branch of the
first multiplication by \(a\) in \(I_{n+1}\) to insert this inverse-branch
average into the integral defining \(I_n\). All integrands are
nonnegative. Iterating gives
\[
       (2/\pi-\pi/(2a))^n\le I_n\le (2/\pi+\pi/(2a))^n .
\]

A direct Fourier comparison using only
\(\sup_\alpha|\sum_{X\le k\le CX} g(k)e(k\alpha)|\ll X^{1-\delta}\)
bounds the word mean by \(O(X^{1-\delta}I_n)\). Here \(g\) would
be the full \(\lambda(k)\), \(\lambda(k+2)\), or
\(\lambda(k)\lambda(k+2)\). The method needs almost a full power of saving:
the norm has exponent \(-\log(\pi/2)/\log a+O(1/(a\log a))\),
about \(-0.03447\). Even a square-root bound does not make this
comparison tend to zero. In fact for any unit-modulus \(g\) on the interval,
Parseval forces its twist supremum to be at least the square root of the
interval length. This is an obstruction to the **supremum-times-norm**
method, not to a localized pairing of these two transforms.

For comparison, [Maynard's restricted-digit work](https://arxiv.org/abs/1604.01041)
uses arithmetic bilinear estimates, not just a single global bound.
The logarithmic two-form theorem of
[Tao](https://arxiv.org/abs/1509.05422) on ordinary consecutive input
indices likewise does not give a correlation on this sparse digital law.
Neither observation depends on the current status of an unweighted
ambient Chowla conjecture.

## 2. Exact multiplicative pairing cannot come from a signed bit matching

Take any two word cylinders, possibly at different depths, with the same
positive number of free bits; fix all other digits on each side and match
free bits by any permutation and any individual complementation. If the
two physical endpoint functions have a constant positive ratio \(r\)
on the cylinder, comparing the two values that differ only in a single
free bit gives \(r=\pm a^{i-j}\). Positivity excludes the minus sign; and
the fact that *all* physical counts are coprime to \(a\) forces \(i=j\)
(by valuation at any prime divisor of \(a\)). Thus \(r=1\).
The same output then has the same full sign. The two different endpoint
types in fact cannot coincide: the lower count is constant modulo
the fixed slope \(63618252914778017089536>2\) at every length,
whereas the upper has residue two more.

There is a companion obstruction for a **polynomial identity** asserting
constant ratio of the two-endpoint products. Regard all free bits as
indeterminates. The product in each cylinder factors over the rationals
into the distinct affine polynomials \(N_0,N_2\). Unique factorization
matches the factors up to scalar, and the two leading coefficient vectors
on a given side are equal, so the matched scalar is the same for both
factors. The one-factor argument then makes it one. This last conclusion
does not apply to an accidental equality of products asserted only on
Boolean vertices without a polynomial identity, nor does it rule out
nonlocal factor-dependent pairings.

**Next estimate.** Splitting the free bits into high and low groups
puts the old input in the exact form \(v_n=a^\ell u+v\), with \(u\)
and \(v\) independently distributed digital measures (each with its own
length). The missing input is cancellation of
\(\lambda(F(a^\ell u+v))\lambda(F(a^\ell u+v)+2)\) against these two
measures together, not an ambient scalar estimate, nor mixing at fixed
moduli. The two obstructions explain why a global Fourier bound or a
fixed-factor bit involution has not supplied it.
