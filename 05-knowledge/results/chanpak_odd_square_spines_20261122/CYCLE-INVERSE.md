# Two flank histories meeting at one integer

Let \(H(x,y)\) be the canonical lower tournament count on spine 4 or 8 and let
\(N=2^{60}(2^{61}-1)-1\). A positive solution gives a literal single-arc flip
from \(E_{61}-1\) to \(E_{61}+1\). The earlier signed Collatz search tested
whether one flank was an iterate of the other. The construction here lets each
flank take its own steps. Changing the sign by a two-bit rule provides an
unconditional, finite, unrestricted search on positive flanks.
The actual computation below stops short of its worst-case depth. [Canonical
chart](../chanpak_decimal_kernel_20261122/FAREY-HYPERBOLA.md) · [Earlier
one-sided inverse](../chanpak_farey_two_three_clocks_20261122/COLLATZ-MAHLER.md)

## Why an even input modulus and an odd output modulus fit together

**PROVED.** Suppose a history admits exactly the inputs \(x=B+A u\) with
outputs \(D+C u\), for integers \(u\ge0\), where \(A\) is a power of two,
\(C\) a power of three, \(0\le B<A\) and \(0\le D<C\). Take two histories
and prescribe the right output minus the left output, an integer \(g\).
Their exact intersection is empty unless
\[
                  D_R-D_L\equiv g\pmod{\min(C_L,C_R)}.
\]
If compatible, put \(M=\max(C_L,C_R)\). Choose the smallest
\(w\ge\max(D_L,D_R-g)\) congruent to \(D_L\bmod C_L\) and
\(D_R-g\bmod C_R\). Both addresses are now on the one increasing line
\[
 x=B_L+A_L(w-D_L)/C_L + (A_LM/C_L)z,\qquad
 y=B_R+A_R(w+g-D_R)/C_R + (A_RM/C_R)z,\qquad z\ge0.
\]
Indeed the moduli divide one another, and all common output pairs are
\((w+Mz,w+g+Mz)\). The congruence records precisely the information
discarded by removing the low input bits.

At any positive integral target there is at most one candidate on this line.
Write \(P=k+tx+cy+qxy\) for the polynomial part of the physical chart.
The [digital envelope](../chanpak_farey_two_three_clocks_20261122/COLLATZ-MAHLER.md)
gives \(P(x,y)\le H(x,y)<P(x+1,y+1)\). Thus only the last \(z\) with
\(P(x(z),y(z))\le N\) can work. An integer square root finds it. The final
test uses the literal binary read \(f\), including the parity of both
addresses. The root can be zero, and an address can be zero on a boundary.
The reported computations count positive address pairs.

## Ordinary negative Collatz and the optimal two offsets

Use the shortened map \(T(n)=n/2\) on even inputs and \((3n-1)/2\) on odd
inputs. Negation conjugates it to the usual \(3n+1\) map on negative
integers. Write \(C_1=(1)\), \(C_5=(5,7,10)\) and
\(C_{17}=(17,25,37,55,82,41,61,91,136,68,34)\) for its three known
positive cycles. Completeness of this list as basins of all positive
integers is conjectural. [Bell and Lagarias, §1](https://doi.org/10.4064/aa170-2-1)

**PROVED.** For any two positive inputs that eventually enter these cycles,
there exist independent iteration times at which their absolute difference
lies in \(\{0,4,24\}\). Among sets containing zero and covering all three
basins by this rule, two positive offsets are necessary, and the least
possible maximum with two positive offsets is 24.

Inputs in the same basin can meet. Between distinct cycles, the possible
absolute differences for \(C_1,C_5\) are exactly \(\{4,6,9\}\), and none of
these occurs between either of the other two pairs of cycles. The two
remaining sets of differences intersect in \(\{24,36,54,81\}\), by direct
subtraction from the displayed cycles. Starting the flanks on the cycles
themselves forces these bounds. The choice \(4,24\) is realized by the
phase pairs \((1,5),(1,25),(10,34)\). Here 24 is an extremal **integer
gap**. No tournament congruence is required.

All of \(2,3,11\) reach the fixed point through \(2\to1\), \(3\to4\to2\)
and \(11\to16\to8\to4\). The exceptional tournament counts have different
destinations under this map. The number \(7=2\cdot3+1\) is already on
\(C_5\), and \(21=2\cdot11-1\to31\to46\to23\to34\) enters \(C_{17}\).
This proves a three-basin observation about these five numbers. The argument
does not imply that a tournament omits either count.

**PROVED.** The fifteen displayed phases are all distinct modulo \(3^5\),
and \(1\equiv82\pmod{3^4}\), making five sharp. A word copying \(d\) steps
of a cycle starting at phase \(a\) admits \(x\equiv a\bmod 2^d\) and has
output congruent to \(T^d(a)\bmod 3^m\), where \(m\) counts odd steps.
This follows by iterating the two affine branches. If two such words each
have at least five odd steps and their outputs coincide, their terminal
cycle phases agree. The shorter word is the terminal segment of the
longer one, since a marked cycle has a unique preceding cycle phase.
On one fixed parity word the iterate is injective. Cancelling the
common segment shows that one input is an iterate of the other.
The two-sided inverse is genuinely wider when these conditions fail.
At four odd steps the claimed *phase separation* already fails.

## A contracting sign choice gives an unconditional cover

**PROVED.** Define
\[
 S(2z)=z,\qquad S(4z+1)=3z+1,\qquad S(4z+3)=3z+2.
\]
At odd \(n\), choose the sign making \(3n\pm1\) divisible by four and take
both halvings. Every positive integer reaches 1, by
\(4(S(n)-1)\le3(n-1)\) and \(S(1)=1\). Moreover each prefix of this
deterministic history has the stated cylinder form. From an old branch,
split \(u\bmod2\) at the unique residue making \(D+Cu\) even. Its
successor has strides \(2A,C\). On each of the other two residues mod 4,
the odd rule gives strides \(4A,3C\). The canonical output constants remain
in \(0\le D'<C'\), since the least old output in these classes is below \(2C\)
or \(4C\) respectively. The three new cylinders partition the old one.

The fiber at a positive \(v\) always has its even predecessor \(2v\).
There is precisely one additional odd predecessor when \(3\nmid v\).
it is \((4v-1)/3\) when \(v\equiv1\bmod3\) and \((4v+1)/3\) when
\(v\equiv2\bmod3\). Hence every genuine split in the inverse tree
takes place at a unit modulo three. All preimages at arbitrary depth
of a multiple of three form one doubling ray. Removing the self-loop
at 1 makes this a rooted tree on *every* positive integer. The merging
vertex and the branch words retain the information needed to recover
the actual two flanks.

For a positive-coordinate bound \(U\ge1\), let \(d(U)\) be the first
nonnegative \(d\) with \(3^d(U-1)<4^d\). After \(d(U)\) steps every
\(n\le U\) is at 1. For a prescribed target \(N\ge k\), use
\(U_L=\lfloor(N-k)/t\rfloor\) and \(U_R=\lfloor(N-k)/c\rfloor\),
discarding a nonpositive bound. It suffices to match exact depths
\((d(U_L),d(U_R))\) at zero offset and run the one-candidate inverse.
At \(E_{61}-1\) the sufficient depth pairs are **\((245,243)\)** and
**\((230,229)\)** on spines 4 and 8. This is complete on positive flanks
regardless of the ordinary Collatz dynamics. If testing a coordinate
axis, the corresponding zero stays at zero, so offsets \(+1,-1\)
between its terminal value and the positive flank suffice. The pair
\((0,0)\) has height \(k\). The depth bound is logarithmic, but the
number of words grows as three to the depth. It gives a stopping proof,
not a polynomial-time complexity claim.

## Completed E61 searches

**FINITE-EXACT.** The following tests found no exact pair. “Compatible”
counts matched word presentations before the absolute test. Distinct
presentations can describe the same pair. The band includes the two
parity corrections and \(0\le f(x)+f(y)\le x+y\). Every band survivor
failed the following necessary congruence modulo \(2h\).

| Families at each spine | Signed offsets | Compatible at spine 4 | In band | Compatible at spine 8 | In band |
|---|---|---:|---:|---:|---:|
| Arbitrary \(T\)-histories, each depth 0–12 | \(-24,-4,0,4,24\) | 10,115,548 | 7 | 10,115,548 | 4 |
| \(T\)-histories copying any of the 15 cycle phases, **all depths** | same five | 452,936 | 0 | 411,122 | 2 |
| Arbitrary contracting \(S\)-histories, each depth 0–8 | \(0\) | 3,304,509 | 1 | 3,304,509 | 0 |

For the all-depth calculation, if \(2^d>\max(U,a)\), the congruence
\(x\equiv a\bmod 2^d\) and \(1\le x\le U\) force \(x=a\) or an empty
branch. From that depth one more cycle period exposes every possible
terminal phase. The right and left sides may start on different phases and
run for different times. In the bounded zero-offset checks a pair already
meeting earlier can advance together until one clock reaches the depth cap.
Nonzero offsets require every pair of depths. The full \(S\) cover
outlined above has not been exhausted.

Reproduce the exact searches and planted unequal-depth, sign and offset
controls from this directory. The script uses Python integer arithmetic.

~~~sh
python3 -B cycle_inverse.py --check --cycles --depth 12 --offsets -24 -4 0 4 24
python3 -B cycle_inverse.py --controlled 8
~~~

## Compressing all controlled words into \(O(2^D)\) unrestricted rows

**PROVED.** Fix any positive lower-chart target \(N\ge k\) and a depth cap
\(D\). All pairs of positive addresses that meet after at most \(D\)
controlled steps on each side can be tested with fewer than
\(68\cdot2^D\) whole one-variable inversions in each physical orientation.
Each inversion leaves the second flank completely free. Only the count of
rows has this bound. This gives a reduction independent of the magnitude
of the prescribed target.

One controlled step is \(\lambda x+\epsilon\), with
\((\lambda,\epsilon)=(1/2,0),(3/4,1/4),(3/4,-1/4)\).
Since \(|\epsilon|\le1-\lambda\), the errors telescope under composition.
After \(d\) steps containing \(m\) odd events,
\[
 \left|S^d(x)-l x\right|\le1-l,\qquad l=2^{-d}(3/2)^m .
\]
For a pair of histories let their multipliers be \(l,r\) and interchange
physical sides to make \(s=l/r\ge1\). If the histories meet then
\(\lvert sx-y\rvert\le (2-l-r)/r\). Define
\[
 z=\max\{j\ge0:P(j,sj)\le N\},\qquad
 W=\left\lceil(2-l-r)/l\right\rceil.
\]
Every admissible first address belongs to \([z-W,z+W]\). If
\(x\ge z+1+W\), then \(y\ge s(z+1)\) and \(P(x,y)>N\).
If \(x\le z-1-W\), then \(y+1\le sz\) since \(s\ge1\), giving
\(H(x,y)<P(x+1,y+1)\le P(z,sz)\le N\). Thus the target supplies
one integer center per slope, and one need only keep the largest
radius for that slope. A rectangle of whole rows replaces all of its
parity histories.

To bound their total, put \(c_0=\log_2(3/2)\). Every slope has the form
\(s=2^{-i}(3/2)^j\), with \(i=d_L-d_R\), \(j=m_L-m_R\) and
\(i\le c_0j\). The cap on \(d_R\) and nonnegativity of \(m_R\) give
\[
 l^{-1}\le 2^D 2^{\min(i,0)}(2/3)^{\max(j,0)}.
\]
For \(j\ge0\) the sum of the last two factors over all permitted \(i\)
is at most \((2+c_0j)(2/3)^j\), totaling at most \(6+6c_0\). For
\(j<0\) it is at most \(2(3/2)^j\), totaling at most 4. Since
\(2W+1\le5/l\), the combined number of rows is at most
\(5(10+6c_0)2^D<68\cdot2^D\) per orientation, even if overlapping
ranges are counted separately. In particular a long string of even steps
can enlarge a radius, but the slopes that can have this large radius
become correspondingly fewer.

**FINITE-EXACT.** At depth 18 the compressed search used 514 slopes
and **11,906,304 whole positive rows** in each of the two orientations on
both spines. Every one of those rows is empty for \(E_{61}-1\), regardless
of how the other flank is chosen. At spine 4 the left-fixed and right-fixed
scans gave respectively 55,107 and 50,137 parity-corrected envelope pairs,
with 101 and 89 passing the \(2h\) congruence. At spine 8 the figures were
18,665 and 11,954, with 1 and 2 passing. The 193 surviving presentations
all failed the literal binary height. A raw address may participate in both
orientations. This clears any controlled coalescence through depth 18,
as well as every other pair represented on the indicated rows. The
worst-case depths above remain substantially larger.

~~~sh
python3 -B cycle_strips.py --check
python3 -B cycle_strips.py --depth 18
~~~
