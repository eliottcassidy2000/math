# Selective cut by absolute factorization, with an independent-switch law

Let $h=H(B)$ and $Q=P_B(2)$ for a nonempty tournament $B$. For a
mark $U\subseteq V(B)$ define:

* $(E_s,O_s)$: weighted counts of nonempty even/odd directed paths
  **starting in $U$**;
* $(E_e,O_e)$: weighted counts of nonempty even/odd directed paths
  **ending outside $U$**;
* $c$: weighted count of nonempty even directed paths starting in $U$
  and ending outside $U$.

Each path $P$ has weight $H(B-P)$ (empty complement has count one),
and set $W_s=E_s+O_s$, $W_e=E_e+O_e$. The unpadded apex insertion,
beating exactly $U$, gives $h+2c$ by the odd-cycle packing identity.

## Two tapes, two independent heights, an exact inverse

Build $TT_\ell\to B\to TT_m$ by ordinal joins, with one apex that
beats the whole left tape, beats exactly $U$ inside $B$, and loses
to the whole right tape. For $X=2^\ell$, $Y=2^m$, $\ell,m\ge1$,
the resulting count is
$$
  QXY+(W_e-Q)X+(W_s-Q)Y+h+Q+2(c-E_s-E_e).
$$
Indeed, classify the odd cycle containing the apex by which block
supplies its first and last vertex after deleting the apex.
Packings avoiding the apex contribute $h$. The internal case gives
$2c$. A nonempty subset of a transitive tape has a unique spanning
directed path. Its even/odd subset counts are $(X/2-1,X/2)$ on the
left and $(Y/2-1,Y/2)$ on the right. Consequently left-to-core
contributes $XW_e-2E_e$, and core-to-right contributes $YW_s-2E_s$.
The left-to-right case can use an arbitrary, possibly empty subset
of $B$; the signed two-cover convolution on $B$ vanishes. To see
the vanishing, use odd-cycle packings and color every singleton or
cycle block independently by which side of the convolution it
belongs to. Each block is odd, so its two signs cancel. Thus the
left-to-right contribution is $Q(X-1)(Y-1)$. This proves the formula.

Put $b=Q-W_e$, $d=Q-W_s$ (both positive), and
$z=h+Q+2(c-E_s-E_e)$. For a prescribed target $N$ in this chart the
exact inverse is to factor the integer
$$
     Q(N-z)+bd=(QX-d)(QY-b).
$$
Try positive factorizations until $(f+d)/Q$, $(g+b)/Q$ are both powers
of two of exponent at least one. The mask comes from the offset
$z$ of the chosen physical source. This is an integer test, with
no approximations or choices of phase after constructing the graph.

### Three odd counts at every pair of heights

Two actual six-vertex cores and three marked realizations have a
common profile $(h,Q)=(45,432)$, $(E_s,O_s)=(114,144)$,
$(E_e,O_e)=(57,72)$. Outgoing row masks are zero-based.

| outgoing rows | $U$ bit mask | $c$ | $z$ |
|---|---:|---:|---:|
| 40,37,17,6,11,28 | 23 | 45 | 225 |
| 48,49,3,7,12,28 | 15 | 46 | 227 |
| 48,49,3,7,12,28 | 27 | 47 | 229 |

So the window is $432XY-303X-174Y+\{225,227,229\}$, beginning
at $999,1001,1003$. To choose the mask directly from $N$, test
$z$ from the three displayed offsets, factor $432(N-z)+52722$,
and use factors $432X-174$ and $432Y-303$ as above. Changing
either tape length does not erase the three-way local choice.

A smaller example uses the circulant $C_5$ whose outgoing
differences modulo five are $3,4$. All per-vertex profiles are
$(9,12)$ for starting and also for ending paths, by rotation and
the literal subset DP. Its $H=15$, $Q=120$. Mark $U=\{0,4\}$ or
$\{0,3\}$, yielding respectively $c=13,14$. Then the chart is
$120XY-57X-78Y+\{71,73\}$. With only a right tape its counts are
$42Y+\{5,7\}$ ($m\ge1$). Both core switches use the same source
as an unpointed tournament.

## An exact independent-switch law in arbitrary surrounding tapes

Take an ordinal word of modules, with arbitrary apex masks in every
module. Record each pointed module by
$(h,Q,E_s,O_s,E_e,O_e,c)$. If some module changes only its
last coordinate $c$, *all crossings over that module are fixed*:
if the module is strictly between the first and last block of
the apex cycle, it supplies the unrestricted $Q$ and signed
convolution zero. If it contains one endpoint, it supplies one
specified endpoint aggregate, even or odd. If it lies outside
the cycle's block interval, its whole factor is $h$. Only the
cycle wholly inside the switched module changes. The exact total
change for independently switching any collection of such modules is
$$
  2\sum_i (c_i'-c_i)\prod_{j\ne i} h_j.
$$
In particular the three-way order-six switch above remains valid
with arbitrary apex contacts in arbitrary transitive tape before,
after, or both; those tape factors have $h=1$. Independent copies
of the same core give the entire additive digit interval in the
displayed step, rather than only individual choices. They do not
give unit spacing if there are multiple copies: each other core's
factor $h>1$ scales every local contribution.

For mutually coprime $h_i$, one can immediately filter the digit
choices by CRT: for $K=\prod h_i$, the target half-gap fixes
$(K/h_i)(c_i'-c_i)$ modulo $h_i$ individually. Even a compatible
residue vector is not enough: the weighted integer sum must equal
the half-gap, not just agree modulo $K$. This is another
obstruction to reading off universal coverage from independence
alone.

The positive family and its factor inverse are all-order, but
still sparse: finitely many fixed two-tape charts of this form
produce $O((\log N)^2)$ counts at most $N$ (both parameters are
individually logarithmically bounded because the apex meets a
nontrivial cut and each tape has contacts). No claim of universal
target coverage is made. The next mathematical step is to choose
a physical three-coordinate profile as a function of $N$, not
merely choose among finitely many fixed offsets once the profile
has been chosen.
