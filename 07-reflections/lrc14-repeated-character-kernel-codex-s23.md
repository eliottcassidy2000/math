# LRC14 Repeated Character Kernel

**Source:** codex-2026-06-19-S23, HYP-2632 / T880.

The useful progress is that the HYP-2630 repeated tail is no longer just a
list of coimage masses.  It has a small signed integer kernel.

The finite transform identity

```text
S_d(a) = sum_{a.r=0} C_d(r) = (1/7) sum_t C_hat(t a)
```

checked over all `159` projective classes, turns the problem into a character
kernel before any reciprocal hyperplane estimates start.

For `4+2`, the QR/NQR split is exact:

```text
2*S/U = -43 - 7*chi_7(a).
```

For `4+1+1`, the important correction is that the old multiplicative
signature is not complete.  The zero lane is affine:

```text
a+b = 2 mod 7.
```

That line kills `(0,2)`, `(3,6)`, and `(4,5)`.  This is the new coordinate the
next proof should keep.  The finite signed ledger is already smaller than the
absolute ledger:

```text
4+2:   -108 U
4+1+1:  +54 U
net:    -54 U
abs:    162 U
```

So the next analytic theorem should not spend effort bounding the repeated
packets by absolute mass.  It should attach this signed `chi_7`/affine table
to the two-large reciprocal hyperplane sum and then prove a cotangent or
Dedekind-style cancellation bound.
