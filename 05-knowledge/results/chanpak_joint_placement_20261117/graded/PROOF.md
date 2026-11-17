# The graded image of the three-contact repair: exact shifts, octal holes

[Zehavi's paper](https://arxiv.org/pdf/2611.18218), §§1.2 and 3,
suggests a method rather than a transferable theorem: filter the parameters,
bound what survives the kernel constraint uniformly, and only then compare
scales. Here the constraint is an elementary affine repair identity, not
an Iwasawa module. Keeping the starting count and shift in the *same*
fiber exposes a base-8 obstruction invisible in either projection.

Let G=2^t>=4, r odd in [1,G), s=G-1-r, and use the literal
[three-contact move](../../chanpak_free_seam_second_move_20261117/three_contacts/PROOF.md).
Its positive prefix integers have x even, z=x-y, and

    u=G(16x+5)+r, w=G(16y+4)+s,
    Delta=15y-(16384G²-5)z-B_r,
    B_r=4096G²+1004G+30+2048G r+40 rev_t(r).

The B formula follows from the previous reverse-complement expression
using s=G-1-r. Here rev_t is bit reversal of exactly t digits.
The ordinary lower count is H=81+142(u+w)+256uw+5(F(u)+F(w));
the reflected partner is always H+2.

## Leading joint image, with a true integral right inverse in each band

For sequences with bounded z, r/G approaching theta, and Delta/G²
approaching lambda in [-1,1], positivity forces z>=0 at large G, and

    15 sqrt(H)/(524288 G³) = 8z+2+r/G+Delta/(2048G²)+O_Z(1/G).   (1)

Indeed y/G²=(16384z+4096+2048r/G+Delta/G²)/15+O_Z(1/G);
the positive prefixes differ by bounded z, so u,w=16Gy+O_Z(G).
The 142 and F terms are lower order (1<=F(v)<=v+1), and
\sqrt{256uw}=16\sqrt{uw}. The error is uniform in r and the
allowed shifts for fixed z bound Z. Thus the normalized source and
specified fine correction are not independent: for fixed z and lambda
the possible source radii converge to the closed interval
[8z+2+lambda/2048,8z+3+lambda/2048].

Conversely there are real, integral witnesses throughout each limit
interval. This does not require choosing z from Delta modulo 15.
Fix any nonnegative integer z and any even Delta with |Delta|<=G².
The condition for integral y is Delta+(16384G²-5)z+B_r=0 mod15.
For t>=5 freeze bit 0 of r to 1 and freeze bits 5 and above
arbitrarily. After removing the common unit G modulo 15, the
changes to B_r from turning on bits 1,2,3,4 are respectively
11,7,14,13. Indeed bit i changes B_r by
2048G 2^i+40*2^(t-1-i); divided by G this is 8*2^i+5*2^(-i)
mod15. The subset sums of those four coefficients cover all
15 phases. So any 32-wide block of suffix values contains
a valid odd r, y obtained by exact division, and x=y+z.
Even Delta and B_r, and odd coefficient of z, force x even.
For fixed z>=0 these prefixes are positive at large G.

In particular for any prescribed theta in [0,1] one can select
a suffix with r/G within O(1/G) while preserving the **exact**
requested even Delta. Formula (1) gives the convergent graded
image with both directions, not just an upper envelope.

## Across sizes: a nested octal obstruction

Compare source counts at ambient scale G0=2^k, writing

    xi = 15 sqrt(H)/(524288 G0³),    8 <= xi < 64.

For a physical scale G=2^(k-j), relation (1) puts xi near

    8^(-j) [8z+2+lambda/2048, 8z+3+lambda/2048].

For a fixed even shift (so lambda tends to zero), each j>=0
therefore reads one octal digit: ignoring interval endpoints,
the band condition is floor(8^j xi)=2 mod8. At j=-1 the
additional strip in [8,64) is [16,24], from z=0; j<=-2 cannot
hit this fundamental annulus. For unrestricted z in the
*iterated graded model* (each z fixed before taking its limit),
the points invisible to all positive levels are those avoiding
digit 2 in the checked octal positions, with [16,24] removed.
This is a missing-digit set of Hausdorff dimension log(7)/log(8):
fix the integer prefix on a surviving open interval and count
7^n equal cylinders of length 8^-n for the upper bound; giving
these cylinders uniform mass 7^-n gives the matching lower bound
by the interval mass estimate. Endpoints have no effect on
dimension. This describes the model, not the closure of actual
unrestricted-scale witnesses; the union of the finite-level
bands is dense.

More robustly, **every fixed near-balance bound leaves recurring
open intervals of absolute targets missing from this template**,
even allowing all even corrections |Delta|<=G². Fix a constant C
and impose |u-w|<=C G. Then z is bounded independently of size.
The source ratio sqrt(H)/G³ is bounded above and below by positive
constants for large G by (1), so for xi in [8,64) only a finite
range of j can occur; as k grows, all admissible G grow too.
The point xi*=27+3/7 has octal expansion 33.333... . For every
j>=0 it misses the digit-2 bands even after expanding either
endpoint by 1/2048 in the 8^j coordinate: the fractional part is
3/7, and the units digit is 3. It also lies outside the possible
j<0 bands, including the 1/2048 expansion. The finite family
has a positive minimum margin there. Hence an open interval
around xi* depending on C, scaled back by G0^6, contains no
such source count for every sufficiently large k. Adding the
reflected +2 and repaired count changes this normalized radius
by o(1), so the same conclusion holds on a slightly smaller
interval for all four endpoints.

This explains what a joint construction has to defeat:
programming every local shift in every suffix cylinder still
does not let this near-balanced grammar choose its own coarse
octal band. Neither this obstruction nor the paper implies an
obstruction to arbitrary tournament counts or to a different seam.
