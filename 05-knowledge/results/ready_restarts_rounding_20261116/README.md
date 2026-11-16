# The restarted rows: persistent physical windows and exact integer selection

**Proved here:** two explicit rooted tournament windows \(13/15\)
and \(41/43\) renew forever, and each member of the first window
is an induced one-vertex parent of its corresponding member in
the second, at every depth. A separate fixed edit realizes the
positive \(+2\) correction from a triangle-scaled 13 ray. Exact
rounding formulas identify the necessary scalar carry and the
binary host of a classified mirror patch. None asserts a repair
for every target.

## 1. What fails at 11, 29 and 35

[THM-5272](../../../01-canon/theorems/THM-5272-cut-form-and-bank-rigidity-for-one-head-over-a-base-the-commutator-carry-of-the-dyadic-and-triadic-steps-and-27-63-1473-1729-in-the-descent-tree.md)
uses \(D(h)=2h-1\). A port in its finite census means some literal
\(h\)-tournament has an induced one-vertex child at \(2h-1\).
It is not a claim that an arbitrary selected host can do this,
or that the same port can continue indefinitely.

In half-count coordinates \(n=(h-1)/2\), a fork is \(n\equiv2\bmod3\)
and a fork whose odd inverse child is sterile has \(n\equiv5\bmod9\).
These translate to \(h\equiv5\bmod6\) and \(h\equiv11\bmod18\).
The odd inverse formula is \((2n-1)/3\), so sterility here means
that this child is divisible by 3. The forks through 45 are
\(5,11,17,23,29,35,41\); intersecting with the no-\(-1\) list gives
exactly \(11,29,35\). Only 11 and 29 are sterile among these three.
Over the whole archived census through 104 the failed forks are
11, 29, 35, 47, 53, 71, 101; the sterile list is 11, 29, 47,
65, 83, 101, and 65 and 83 have the doubling port. Thus a modular
descendant type and a physical missing-port diagnosis are distinct.

Write \(T(h)=3h\), implemented by ordinal join with a directed
triangle. The arithmetic identities
\[
D(7)=13,\qquad D(T(7))=D(21)=41,\qquad T(D(7))=39
\]
are the first restart of the commutator. The other dyadic child
\(2h+1\) gives 15 at 7 and 43 at 21; after tripling the first
two-child window, its endpoints 39 and 45 move inward to 41 and
43 at the same midpoint 42. None of these identities starts from a
physical tournament at 7 or 21.

## 2. The literal, indefinitely extensible restart diagram

Use outgoing-mask notation, where bit \(j\) of row \(i\) records
\(i\to j\). Take the common six-vertex base
\[
(34,36,41,51,39,0).
\]
Add root 6 with out-mask 3 or 6. These are respectively
\[
(34,36,105,115,103,64,3),\qquad (98,36,41,115,103,64,6),
\]
of counts 41 and 43. They differ in exactly the old arcs on
\((0,6),(2,6)\). Delete vertex 3 in the lower graph and vertex 1
in the upper, and read surviving vertices in increasing order.
The resulting rooted six-vertex rows are
\[
(18,20,49,55,32,3),\qquad (48,21,57,51,32,2),
\]
of counts 13 and 15. Their roots are both index 5: each is the
restriction of vertex 6. In these displayed labelings they differ
on exactly three old pairs \((0,1),(0,5),(2,3)\). This is not a
claim that a different choice of witnesses cannot use fewer edits.

For a root \(x\), record \((h,g,r,s)\), where \(g=H(T-x)\) and
\(r,s\) are the sums of \(H(T-V(P))\) over even-order or odd-order
directed paths \(P\) starting at \(x\) (include the singleton in
\(s\)). The four profiles are
\[
(13,3,6,9),\ (15,3,6,9),\ (41,9,20,29),\ (43,9,20,29).
\]
The retained-root append from the [ready-port archive](../inverse_fences_landing_20261116/INVERSE.md)
adds a new vertex beaten by all old vertices except \(x\), which
it beats. Its exact profile action is
\[
(h,g,r,s)\longmapsto (h+2r,\ g,\ r+s-g,\ r+s).
\]
For the first coordinate, all new odd cycles pass through the
new vertex and use an even-order old path starting at \(x\).
For the last two coordinates, any old rooted path can end with
the new vertex except the singleton at \(x\); that extension
reverses its parity. Paths omitting the new vertex have neutral
sink padding in their complement. Removing \(x\) leaves the new
vertex as a sink, so \(g\) is unchanged.

When \(s-r=g\), this appends the rule \(r\mapsto2r\) indefinitely
and holds the odd intercept \(h-2r\) fixed. Apply the same append
\(k\) times at all four roots. The complete diagram has rows
\[
W_{13}(k)=(1+12\cdot2^k,\ 3+12\cdot2^k),\qquad
W_{41}(k)=(1+40\cdot2^k,\ 3+40\cdot2^k).
\]
The induced deletions commute with every append; all new vertices
keep the same arcs in the two members of their own window.
Consequently the horizontal edit budgets stay 3 and 2, and the
vertical links stay induced one-vertex growth. The lower members
use the \(-1\) doubling digit; the upper use \(-3\), not \(-1\).
For example 43 goes to 83 at this root, not 85. This is consistent
with the no-\(-1\) at 43 in THM-5272.

A particularly sharp width/transport distinction is that **no**
cross-only edit of the ordinal join of a directed triangle with
an intact \(H13\) core can give 41. Nor can one of a triangle with
an intact \(H15\) core give 43. Both target counts are prime, so
all nontrivial source components and the triangle would have to
lie in the unique nontrivial strong component of the output.
The intact 13 or 15 host has at least five nontrivial vertices
(the order-four maximum is 5; decomposition of 15 into 3 and 5
uses even more). The target would require at least eight strong
vertices, but [Busch's sharp floor](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)
is at least 45 there. The literal diagram grows the unscaled
inputs instead of trying to trim an intact triangle product.

## 3. A scale repair that is actually physical on one entire ray

For the lower 13 ray, \(3W_{13,\mathrm{lower}}(k)=3+36\cdot2^k\).
Take \(H13\) rows \((12,21,24,18,1)\), root 0, and put an ordinal
triangle before them. The resulting order-eight graph has root
index 3, rows
\[
(250,252,249,96,168,192,144,8),\quad (h,g,r,s)=(39,9,18,27).
\]
Change exactly the five old arcs on
\((0,1),(0,5),(2,6),(3,6),(4,5)\) to obtain
\[
(216,253,185,32,136,209,156,8),\quad (41,5,18,23)
\]
at the same root. Both profiles have \(s-r=g\) and the same \(r\).
Thus the five edits change \(3+36\cdot2^k\) into \(5+36\cdot2^k\)
on all lengths, even after the tail has been appended. This is an
exact realization of the **positive** scaled carry \(+2\) on
every member of the 13 ray. Its target 41 witness is not the
\(41\) witness of the matched 41/43 window; its intercept is 5
rather than 1. For the negative scaled carry no analogous result
is claimed.

There is a useful alternate construction with a marked predecessor.
For an old vertex \(y\) let \(R_y\) count all weighted splits
ending the first run at \(y\); precisely
\(R_y=\sum_{A\ni y} E_A(y)H(T-A)\).
It is at least \(H(T)\) by splitting a full path at \(y\).
To add a near-source beating all old vertices except \(y\), start
from an ordinal source and reverse just the edge to \(y\).
The old paths lost are the ones starting at \(y\) after the new
vertex, and the gains are exactly \(R_y\). If there are \(t\ge1\)
ordinal source pads before \(T\), the loss term is zero and
\(R_y\) is multiplied by \(2^t\): assign each pad to either run
while retaining the unique linear order within each assignment.
The resulting count is
\[
(1+2^t)h+2^t(R_y-h).
\]
If \(y\) has just one predecessor \(z\) and \(T-y\) has exactly
one full path starting at \(z\), then \(R_y-h=1\): the only
unconcatenable splits put \(z\) first in the second run and have
the singleton \(y\) as their first run. This gives physical
\(3h+2\) at \(t=1\) and the maximal positive even correction for
each multiplier \(1+2^t\) on those hosts. The verifier has an
explicit 13→41 instance. This does not preserve an intact
triangle or certify a reusable marked predecessor at every count.

In particular the two earliest sterile forks do have this different
kind of growth port: the 11 witness (18,4,9,19,6) marked at 3 gives 35,
and the 29 witness (58,4,9,18,38,14) marked at 0 gives 89. Both satisfy
R_y=h+1 and use one neutral pad followed by the near-source extension.
Thus even a missing doubling port at a sterile fork is not a missing
ternary carry port. A separate rooted lower lift reaches 125 from
41 on two additional vertices: masks (256,5,9,323,271,287,439,319,6)
at root 8 have profile (125,27,62,89), and their restriction to
(0,1,2,3,4,6,8) has (41,9,20,29). The same profile update
(g,r,s) → (3g,3r+2,3s+2) recurs here on the lower side only;
this does not supply the paired 127 or indefinite ternary renewal.

## 4. Three independent ways to select the correct integer

### Nearest odd parent and a finite carry

At depth \(k\ge1\), write \(L=2^k\) and use the balanced dyadic
row \(C_k(h)=\{Lh+d:d\text{ odd},\ |d|<L\}\). Here \(h\) is
odd. For any odd multiplier \(c\) and any \(P\in C_k(ch)\),
**there is a unique** odd \(H\in C_k(h)\) and even \(e\) with
\[
P=cH+e,\qquad |e|\le c-1.
\]
Take the nearest odd to \(P/c\). It stays inside the walls
\(L(h-1),L(h+1)\) of the original row. A tie would put \(P/c\)
at an even integer, which is impossible. Multiplying the
strict distance \(<1\) by \(c\) gives the asserted error bound.
For \(c=3\) the reselected source needs only \(-2,0,+2\),
independent of depth. Holding the same signed address instead
would incur the growing commutator error \((c-1)d\).

The positive carry already succeeds at two large prescribed primes:
30011 selects the unique odd parent 10003, and 30089 selects 10029.
In the [previous binary inverse](../prescribed_prime_30000_20261116/README.md)
their even addresses are respectively 6246 (also 6558) and 5272.
To see why it preserves the rounded parent, start with the
one-apex word RC(bin v), where the apex z has one full start
and R_z=1+v. Add y with z as its only predecessor, making
H(B)=f(v)+v and R_y=H(B)+1. The one-pad carry above now gives
exactly the prescribed prime 3H(B)+2. This is the crossed-triangle
apex inverse viewed from a different marked vertex, not a new
independent positive family. It identifies a literal rounded
ancestor at these targets. The neighboring primes 30013 and
30091 instead select 10005 and 10031 with negative carry;
the positive port does not supply those lifts.

Equivalently multiplying \(Lh+d\) by 3 and re-normalizing gives
parent \(3h+2j\) and digit \(3d-2Lj\) with \(j=-1,0,1\).
Those are the three nearest-parent stages. The thresholds in
\(d/L\) are \(\pm1/3\), each avoided by at least \(1/(3L)\):
\(3d\) cannot equal \(\pm L\). An approximation with certified
error strictly below \(1/(3L)\) therefore selects the exact
carry. The inner lane has \((L+2(-1)^k)/3\) addresses and each
outer lane has \((L-(-1)^k)/3\), by counting the odd integers
subtended by the thresholds. This is scalar selection; the
physical graph must still implement the selected correction.

### Root profiles and rounding the stable slope

For a general near-sink root, from \(k=1\) onward the exact
pencil is
\[
H_k=2^k(r+s-g)+(h-2s+2g).
\]
The new state's \(s-r=g\) after just one append, so the
formula follows by summing the ensuing geometric sequence.
A uniform bound on the difference of two pencils along
arbitrarily large common \(k\) requires equality of their
integer slopes. If the error from the intercept is below
\(1/2\) after division by \(2^k\), rounding selects the slope
exactly; the intercept remains a separate target coordinate.
This explains why a one-vertex first hit is not enough for
an exact persistent window.

### Square-root address recovery

The [classified two-adjacent-swap bank](../prescribed_prime_30000_20261116/TWO_SWAPS.md)
has an exact root decoder, not a search over outside hosts.
For a nonempty permitted middle word of length \(m\), set
\(L=2^{m+4}\). Both window endpoints satisfy
\[
L x^2\le H<L(x+1)^2.
\]
Indeed the exteriors match at \(X=Z=x\); the leading term
is \(Lx^2\), the central profile sum is at most \(2L-3\),
\(f(x)\le x+1\), and the central count plus two is less
than \(L\) for all listed centers. The concatenation formula
gives the bracket. So \(x=\lfloor\sqrt{H/L}\rfloor\), and
the exact residual of the formula tests the two endpoints.
For touching pairs \(x\) is odd and the lower count is
\(16x^2+31x+f(x)+f(x+1)+10\); the width is 2 and both endpoints
are in \([16x^2,16(x+2)^2)\). Choose the largest odd \(x\)
not exceeding \(\lfloor\sqrt{H/16}\rfloor\).
The number of center tags is logarithmic in the target.
This accepts 10007 and 10009 at \(x=12\), rejects 30011 and
30013 in this restricted family, and does not claim every
prescribed prime has some tag.

## 5. Two-corner selection and why scaling can miss the center

More generally a \(+2\) pair is specified by **two** integers:
its midpoint and its half-width. Tripling 10007/10009 gives
30021/30027, midpoint 30024 and half-width 3. Trimming inward
without changing midpoint gives 30023/30025, respectively
\(7\cdot4289\) and \(25\cdot1201\). The independently
constructed 30011/30013 has midpoint **30012** and half-width
1. It uses the width-reducing correction and a different
center, unlike 39/45→41/43 where the center agrees.

There is an exact common local law for the physical 41/43
switch and the 30011/30013 switch. Over a fixed base with
a new head and two optional contacts, write the four counts as
\[
F_{00}=a,\quad F_{10}=a+2u,\quad
F_{01}=a+2v,\quad F_{11}=a+2(u+v-C).
\]
By THM-5272's cut form \(C=W(e,f)+W(f,e)\ge H(B-e-f)\).
The opposite diagonal is a prescribed \(p,p+2\) exactly if
\(v-u=1\) and \(a+u+v=p+1\). The quadratic interaction cancels
from those two equations, even though the intermediate counts
can be far below them. In particular their mean is at most
\(p+1-H(B-e-f)\); an exact narrow diagonal is not a narrow
four-corner bank. Exact or certified sub-half-unit approximations
to the two integer selectors suffice, but not a bound on just
one selector.

| Pair | \(F_{00}\) | \(u,v\) | \(C\) | \(F_{11}\) |
|---|---:|---:|---:|---:|
| 41/43, masks 3 and 6 over common root | 27 | 7, 8 | 2 | 53 |
| 30011/30013, switches at vertex 4 | 29895 | 58, 59 | 1246 | 27637 |

The large-prime square-root and contact-switch facts complement
each other: rounding can find a unique address **inside a chosen
family**, and an exact context change can give a small carry
at that address. It is still open how to choose a family and
context that closes the two integer equations for arbitrary
prime targets.

## Replay

Run python3 verify.py, python3 extra.py, and python3 square.py. The root-profile and
row checks use direct subset path DP; the dyadic lane and square
root checks use exact integers.

Further exact construction: [a two-head template reaches 13, 41, 125 and several of the prescribed primes](TWO_HEAD.md), with a proved limitation on matching its tail slopes.
