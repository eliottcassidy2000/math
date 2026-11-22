# Collatz parity cylinders and Mahler rounding shells at E61

Throughout, the target is the canonical lower height \(N=E_{61}-1\) on
spine 4 or 8, with the [inherited chart](../chanpak_decimal_kernel_20261122/FAREY-HYPERBOLA.md). We check positive
flanks in both physical orientations. Neither dynamics is assumed to
govern the unrestricted two-word chart.

## Absolute inversion of a Collatz cylinder

The standard shortened map is \(T_+(x)=x/2\) at an even input and
\((3x+1)/2\) at an odd input. Write \(T_-\) for the same map with \(-1\)
instead. Use the [Bernstein–Lagarias parity coding](https://doi.org/10.4153/CJM-1996-060-x)
and the repository's [signed 9/8 mirrors](../chanpak_mirimanoff_unit_lifts_20261122/MIRROR-UNIT.md), but don't
freeze the odd-step count or the number of halvings in advance. On
each residue \(x=2^d Z+B\) the exact \(d\)-step value is \(C Z+D\)
with \(C=3^m\), where \(m\) counts odd steps. To extend the cylinder
write \(Z=2Z'+b\). The next odd/even choice depends only on \(Cb+D\)
modulo two; the new constant is \(T_\pm(Cb+D)\) and the new slope
is \(C\) or \(3C\). This is a binary branch tree with no input-size
bound.

Put \(P(x,y)=k+tx+cy+qxy\). The literal height \(H\) obeys
\[
 P(x,y)\le H(x,y)\le P(x,y)+h(x+y)+r+e<P(x+1,y+1).
\]
The last strict inequality uses \(q>h\) and \(q+t+c>r+e\). In any
one cylinder both coordinates strictly increase with \(Z\). So
one first solves the positive affine-part quadratic inequality
\(P(2^d Z+B,CZ+D)\le N\), and tests only its last solution \(Z\).
One more test modulo \(2h\) is valid because \(f\) is even. The
floating value in the replay is only a root seed: integer
comparisons and a fallback binary search certify the last \(Z\).
This proves a uniform exact inverse per parity branch, not a
bounded search over starting integers.

**FINITE-EXACT.** In all 268,435,440 depth/sign/spine/orientation
presentations through depth 24 there were eight parity-corrected
envelope hits and zero residues modulo \(2h\). Thus neither signed
orbit can link the two addresses up to this depth. The eight
presentations need not represent distinct pairs. This does not rule
out longer connections or different signs chosen within one script.

~~~sh
g++ -O3 -std=c++17 -Wall -Wextra -Werror farey-collatz-shell.cpp -o farey-collatz-shell
./farey-collatz-shell 1 check
./farey-collatz-shell 24
~~~

## Exactly three fixed rows capture every rounding script

Let \(L(x)=\lfloor3x/2\rfloor\) and \(U(x)=\lceil3x/2\rceil\).
All the endpoints from arbitrary interleavings of \(L,U\), after
exactly \(d\) steps from \(x\), are **all** integers
\([L^d(x),U^d(x)]\). Indeed there is no hole
between \([L(s),U(s)]\) and the next interval at \(s+1\):
\(L(s+1)=U(s)+1\). Induct on \(d\).

Set \(\lambda=(3/2)^d\). The accumulated rounding error is at most
\(\frac12\sum_{j<d}(3/2)^j=\lambda-1\). Write \(z\) for the largest
integer with \(P(z,\lambda z)\le N\), using the physical chart or
its transpose. If \(x\ge z+2\), then \(y\ge\lambda(x-1)+1\) and
\(P(x,y)>N\). If \(x\le z-2\), then \(y+1\le\lambda(x+1)\), so
\(P(x+1,y+1)\le N\) contradicts \(H=N\). Thus only
\(z-1,z,z+1\) are possible, however many rounding scripts there
are. One computes the exact reachable interval in each of these
rows and inverts it against the original \(f\), first fixing
parity and then using the [low-word split](../chanpak_e61_fibonacci_target_20261122/SUFFIX-INVERSE.md).
There is no uniform cap on the length of an address in a tested
row, except the one implied by its absolute height.

**FINITE-EXACT.** Both spines miss at every depth through 96,
allowing all mixes of up and down, both orientations. If only
one rounding direction is used, the maps are strictly increasing
and the one-chain inverse works at each depth. For the ceil map
the value from 1 goes to infinity; for the floor map this is true
from 2, while 1 is fixed. The cutoff is the first depth where the
smallest growing start already has \(P>N\). This checks **all**
depths. The tested numbers of depths in (floor, ceil) order are
(168,169), (168,169) in the two orientations at spine four, and
(157,159), (158,159) at spine eight; all are empty.

In [Mahler's original problem](https://doi.org/10.1017/S1446788700005371)
a Z-number \(\alpha>0\) has fractional parts
\(\{(3/2)^j\alpha\}<1/2\) for every \(j\). If
\(g_j=\lfloor(3/2)^j\alpha\rfloor\), then \(g_{j+1}=U(g_j)\).
For even \(g_j\) a fractional remainder \(<1/2\) cannot change the
next integer part; for odd \(g_j\) the next unrounded fraction
starts at \(1/2\), so the next Z-condition forces a carry. Also
\(g_0>0\): any positive \(\alpha<1/2\) first crosses \(1/2\) below
\(3/4\), and \(1/2\le\alpha<1\) already fails. The calculation
therefore excludes every pair of integer parts from one such orbit
as the two E61 flank addresses. It says nothing about whether any
Z-number exists, and it does not justify demanding either relation
on independently chosen Farey words.

~~~sh
python3 farey-mahler-roundings.py --check
python3 farey-mahler-roundings.py --depth 96
python3 farey-mahler-roundings.py --deterministic
~~~

## Where the fraction 3/2 is a physical multiplier

**PROVED.** There is a different tournament \(3/2\) clock that does not
require the two Farey flanks to lie on a special curve. For an arbitrary
nonempty tournament \(A\), let \(q_j(A)\) count assignments of all
vertices into \(j\) ordered directed-path slots, allowing empty slots,
and write \((h,Q,R)=(q_1,q_2,q_3)(A)\). Ordinally prefix \(m\) singletons
to obtain \(P^m A\). A prefixed singleton can go first in any of
the \(j\) slots, so the reads are \((h,X,Y)=(h,2^m Q,3^m R)\).

Let the closed wrapper be the triangle \(a\to A\to b\to a\)
with singleton anchors \(a,b\). Opening means reversing \(b\to a\);
the resulting tournament is the ordinal join \(a\to A\to b\).
Use two nested wrappers with independent anchor reversals.
This is a literal four-corner square on fixed vertices. By the
[ordered-slot transfer](../chanpak_farey_complementary_exact_20261122/TRIANGLE-SELECTOR.md),
the closed wrapper gives \(H'=h+Q,\ Q'=4Q+2R-2h\), whereas the
open wrapper multiplies \(q_j\) by \(j^2\). Write inner bit first
and 1 for closed. The four heights are
\[
 H_{00}=h,\qquad H_{10}=h+X,\qquad H_{01}=h+4X,\qquad
 H_{11}=-h+5X+2Y.
\]
Thus the first flip reads \(X=H_{10}-H_{00}\) and the mixed
difference reads \(Y=h+(H_{11}-H_{10}-H_{01}+H_{00})/2\).
Adding one ordinal pad multiplies the ratio \(Y/X\) by exactly
\(3/2\). No individual flip distance is asserted to stay fixed.
This is an exact 2/3 separation on an arbitrary interior, unlike
the numerical \(9/8\) endpoint ratio on the [square-sum cut rail](../chanpak_middle_out_square_sum_20261122/SUPPORTED-SELECTOR.md).

More generally \(j-1\) nested wrappers give a cube with \(j-1\)
literal flips. Only its all-closed corner can take an upward shift
at every step of the tridiagonal transfer, so its top Boolean
difference has coefficient \((j-1)!\) at \(q_j\), and depends
otherwise on lower reads. Those are recoverable from the lower-rank
cubes. The resulting padded \(j\)-slot registers need not be
freely selectable, so the clock is a diagnostic, not a proof of
Chan–Pak or a scalar filling rule.
