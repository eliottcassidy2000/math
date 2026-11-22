# The entire cube-root controller window is empty at exponent 61

Put \(E=2^{60}(2^{61}-1)\), \(N=E-1\). On canonical Farey spines
4 and 8 the lower chart is
\[
 H(x,y)=k+tx+cy+qxy+r[x\text{ odd}]+e[y\text{ odd}]
       +h(f(x)+f(y)),
\]
and reversing the designated arc always gives \(H+2\).
The [inherited complementary construction](../chanpak_decimal_kernel_20261122/CUBE-ROOT-PHASE.md)
places both addresses at a common sum \(S=x+y\). Set
\(z=\lfloor\sqrt{N/q}\rfloor\); its absolute placement selects
a binary scale \(B\), a divisor \(g\mid h\) from the odd part
of the modulus, and a sum in
\[
 S\equiv-1\bmod B,\qquad
 2z+2B\le S<2z+(g+2)B.
\]
The same \(z\) is obtained by using the midpoint \(E\) instead of
\(N\), or perturbing either by \(1/2\). The formula follows by
spacing the allowed sums by \(gB\) and taking the first at least
 \(2z+2B\).

**Finite-exact.** No pair in this entire placement can give
\(H=N\), at either spine. The test allows all powers of two \(B\),
not just the cube-root size; it replaces \(g\) by \(h\), drops
the remaining controller congruences, and allows **all** word
pairs with the resulting sum. Thus this exclusion applies to
shorter complementary controllers and to choosing the modulus
adaptively, as long as the same placement window is used.

| spine | \(z\) | possible powers \(B=2^b\) | sum presentations | sums with a nonempty integer envelope | candidates after the digital divisibility check | exact hits |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 1076566223924186 | \(0\le b\le99\) | 227300 | 62987 | 8485 | 0 |
| 8 | 135481075480154 | \(0\le b\le93\) | 3056602 | 526568 | 185 | 0 |

A sum can occur in more than one scale, so the table counts
presentations, not unique sums. The closest residuals **among the
phase-screened envelope candidates of these windows** were respectively
\(-320791632154190\) and \(-177097323900534090\). These are not a
record for all flank pairs. The unrestricted 4 and 8 search is
not exhausted by this result.

## Exact inverse and all-aspect long-terminal exclusions

Any pair at the target satisfies \(S\le (N-k)/\min(c,t)\) by
positivity. At a prescribed sum \(S\), put \(A=qS+t-c\) and
\(D=A^2-4q(N-k-cS)\). Since \(0\le f(u)\le u\),
\[
 D\le (2qx-A)^2\le D+4q(hS+r+e).
\]
Take exact integer square roots, intersect both root intervals
with \(0\le x\le S\), and then require the residual after the
bilinear, linear and parity terms to be divisible by \(2h\)
because \(f\) is even. For large intervals the implementation
solves the quadratic congruence separately at each of the
three prime factors of \(h\), in each parity, and applies CRT.
For short intervals it checks the same divisibility directly.
Every survivor is evaluated in the original integer chart
with \(y=S-x\), with no rounding or bounds on word complexity.

Independently, without the placement window or any aspect ratio
restriction, the same inverse excludes **all** target pairs
with \(S\) ending in at least **79 ones** at spine 4 or **73
ones** at spine 8. This covers respectively 6,168,852 and
4,886,502 sums, with 23 and zero final phase candidates,
none exact. It includes the axes and arbitrary imbalance in
these particular terminal sectors.

Reproduce with Python from this directory. The script imports
the existing integer chart and receiver generation from the
sibling packet. Run `python cube_controller_exact.py` without arguments for planted controls;
run with `--spine 4` or `--spine 8` for the full placement
window. Add `--terminal 79` or `--terminal 73`, respectively,
for the all-aspect scans.

## Why making the inherited bank's modulus larger is not enough

For the [existing bank implementation](../chanpak_decimal_kernel_20261122/complementary_controller.py)
with modulus \(2^aM\) and common suffix width \(\ell\), any
output at or below \(N\) must have \(v_2(S+1)\ge\ell+a+1\).
The cap on \(S\) above bounds this valuation by 101 at
spine 4 and 95 at spine 8. The bank formula alone gives
\(\ell+a+1\ge2(a+1)+3(M-1)2^a(9+a)\). Thus no mixed modulus
with \(a\ge1,M\ge3\), and no odd modulus \(M\ge5\), can even
fit into the absolute budget as this bank is written.
For \(M=1\) the exact power-of-two bank uses one switch per
precision bit, one final loader, and its leading zero block,
so \(\ell=2a^2+5a+2\). The surviving literal bank moduli by
this length test are only \(1,2,3,4,8,16,32\); at 32
the required terminal-one length is 83, already excluded
globally above. For \(M=3\) the bank needs 71 terminal ones.
These are limitations of that implementation, not lower
bounds on the length of every modular algorithm. In
contrast the empty placement window has no algorithmic
length assumption.
