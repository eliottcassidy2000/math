# A fixed-sum generating function and a Hardy–Littlewood test

This note concerns the physical [fixed two-tape carrier](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md)
with (B,Q,T,S,h)=(2855,5760,4034,3914,45) and
F(x,y)=B+Qxy+Tx+Sy+h[p(x)+p(y)]. Recall
p(2a+1)=p(a), p(2a)=p(a)+a+(a mod 2), and 0<=p(a)<=a.
The variables x,y range over nonnegative integers. No picture of
its local residues is assumed to describe its exact spectrum.

## An exact fixed-sum polynomial

For each s>=0 and e=0,1 set

    W_s^e(u,z) = sum_{0<=x<=s, x=e mod 2} u^x z^[p(x)+p(s-x)].

Use W_0^0=1, W_0^1=0, and W_-1^e=0. Write W_s=W_s^0+W_s^1,
and b_t(e)=(t-e) mod 2 in {0,1}. The complete recurrences, for s>0, are

    W_(2t)^0(u,z)   = z^t sum_e z^[e+b_t(e)] W_t^e(u^2,z)
    W_(2t)^1(u,z)   = u W_(t-1)(u^2,z)

    W_(2t+1)^0(u,z) = sum_e z^e W_t^e(u^2 z,z)
    W_(2t+1)^1(u,z) = u z^t sum_e z^b_t(e) W_t^e(u^2/z,z).

These follow by writing x=2a+e and s-x=2b+e'; s and e fix
e' and a+b. The odd-sum first line pays a+(a mod 2) when x is
even; the second pays b+(b mod 2) when y is even. The even-sum
first line pays the sum of these surcharges; when both are odd the
two high prefixes sum to t-1 and no surcharge is paid.
The Laurent substitution in the odd last line has nonnegative
final exponents, since a<=t.

For a literal address interval [L,U], first discard monomials
of u-exponent outside the interval. Then use the linear monomial
map

    u^x z^j --> Z^[B+S s+Q x(s-x)+(T-S)x+h j].

The resulting polynomial V_(s,L,U)(Z) enumerates physical path-count
presentations in that window. Its coefficient [Z^N] is exactly the
number of witnesses (not isomorphism classes); it is positive exactly
on the local exact spectrum. Hitting an integer tolerance band
[N-D,N+D] is tested by summing the corresponding coefficients.
For modular counts, reduce exponents modulo m after the map, or
equivalently apply a roots-of-unity filter before coefficient extraction.
This *is* a positivity test for a finite window, not an assertion
that a congruence is an exact witness.

## Congruent in every marginal does not justify an Euler product

In the exact-sum window s=9, 0<=x<=7 the literal outputs are

    38261, 84641, 119051, 142301, 154031, 154151, 142661, 119651.

They occupy both odd classes mod 4 and every class mod 7.
Modulo 28 they miss 9,15,17,19,21,27. This remains an arbitrarily
high physical window: put (X,Y)=(16H+x,16H+9-x) with H a multiple
of 28. By the exact binary split, at fixed x+y the quadratic and
linear parts are shifted by constants, and the p-sums change by

    2p(H) + H[30-rev_4(x)-rev_4(9-x)].

The boundary-bit term vanishes as H is even. All classes are
therefore translated by 90p(H) modulo 28. Also the lanes x=0,1
have equal reverse-address sum, so these two counts remain
exactly 46380 apart. Pick a target within their interval in any
shifted missing class. It passes mod 4 and mod 7 separately,
has a witness below and above, but fails the joint mod 28 test in
the actual same window. The corresponding x and y differ by at
most nine and the windows attain unbounded height.

This is a failure of independence of *window-conditioned*
marginals, not an obstruction to a theorem using all windows or
a statement about simultaneous p-adic admissibility at every
precision. It tells us to use the joint polynomial V when building
local factors for a chosen window; multiplying marginals loses
which pair supplied each residue.
