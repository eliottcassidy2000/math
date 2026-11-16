# Selecting center 126: a new physical interface for 125/127

**Exact result.** There is a two-arc \(125\leftrightarrow127\) switch on eight
vertices that keeps its +2 difference under every subsequent append at
the same root. Its midpoint is precisely three times the earlier
41/43 midpoint at every depth. It does not have the previously blocked
full profile: the correct continuation coordinate is smaller.
The natural \(5^3\) and \(10^2+3^3\) recipes also lead to explicit
physical constructions and explicit obstructions, depending on what
is kept fixed.

## 1. Retain the response, not the whole deletion profile

For a root \(x\), use the [earlier notation](../ready_restarts_rounding_20261116/README.md)
\((h,g,r,s)\): count, deletion count, even and odd weighted root-path
sums. Put \(b=s-g\). The append of a new vertex beaten by every old
vertex except \(x\), which it beats, acts on the **observable** triple by
\[
(h,r,b)\longmapsto(h+2r,\ r+b,\ r+b).
\]
Thus \(g\) and \(s\) are not separately needed to predict the entire
future count sequence. Conversely the initial count and the first
two increments recover \(r\) and \(b\). Two rooted graphs retain a
constant count difference (starting from depth zero) exactly when
their \(r,b\) agree. On the stable locus \(b=r\) the count is
\(h-2r+2r\cdot2^k\).

There is also an exact interpretation of stability. Let \(e\) be the
number of Hamiltonian paths ending at \(x\) and \(R=r+s\) the total
weighted split response. Count a near-sink append in two ways:
it adds \(2r\) by the odd-cycle expansion, and it adds \(R-e\)
by splitting paths on the single backward arc. Hence \(s-r=e\).
Thus stability is \(e=g\): **every** Hamiltonian path in the deletion
ends at a vertex that beats \(x\). In particular a stable root lies in
the last strong component. The rooted response of any preceding
ordinal block is a common multiplier of \(g,r,s\), since every
root path leaves that block entirely in its complement.

### Literal two-edit window

In outgoing-mask notation, take
\[
\begin{array}{c|c|l}
h & (g,r,s)\text{ at root }2 & \text{rows}\\ \hline
125 &(9,60,69)&(4,5,72,35,15,87,155,63)\\
127 &(11,60,71)&(4,5,72,99,47,71,147,63).
\end{array}
\]
Reverse the arcs on pairs \((3,6),(4,5)\). Neither touches the root;
all later near-sink appends at root 2 are compatible with the edits.
The two counts after \(k\) appends are \(5+120\cdot2^k\) and
\(7+120\cdot2^k\). The old proposal required \(g=27,r=62,s=89\) on
both sides. This new interface keeps \(r=60,b=60\) on both sides but
does **not** require either deletion count to stay 27.

The single-edit corners make the two cancellations visible. Their
profiles (in lower, first-only, second-only, upper order) are
(125,9,60,69), (141,13,66,79), (105,9,50,59), (127,11,60,71).
All four are stable. The count changes +16 and -20 receive mixed
correction +6; the slope changes +12 and -20 receive mixed
correction +8. Thus the same physical two edits settle both
the +2 in the intercept and the zero net change in slope.

### The center comes from an old row, not a trial approximation

Tripling the old 41/43 row gave \(3+120\cdot2^k,9+120\cdot2^k\).
Our new row is the **inward trim** by \(+2,-2\), with the same
midpoint \(6+120\cdot2^k=3(2+40\cdot2^k)\) for all \(k\).
This is also literal in a common labeling: put an ordinal
triangle before each of the earlier rooted 41/43 hosts and keep
the root at index 9. The four order-ten rows are:

| Count | Outgoing masks |
|---:|---|
| 123 | (1018,1020,1017,272,288,840,920,824,512,24) |
| 129 | (1018,1020,1017,784,288,328,920,824,512,48) |
| 125 | (1022,1020,952,768,296,968,412,792,512,80) |
| 127 | (1022,1020,952,768,360,840,396,824,512,80) |

The horizontal moves are both two old arcs. The vertical edits
123→125 are \((0,2),(2,6),(3,4),(3,9),(5,7),(6,9)\);
129→127 uses \((0,2),(2,6),(3,4),(4,6),(5,9),(6,9)\).
All four are stable with \(r=60\); the four deletion counts are
27, 27, 9, 11, respectively. Hence this entire edit diagram
persists at every depth. No minimality for the six vertical
edits is claimed. In particular this is **not** an induced
extension of the old 41/43 graphs or an edit keeping the
triangle and core intact.

An exact inverse for the center of this family is elementary:
\(C-6\) must be positive, with odd part 15 and dyadic valuation
\(k+3\) for \(k\ge0\). Then the two linked graphs give \(C\pm1\).
This selects 126 and the subsequent centers of the very same
pencil; it does not select an arbitrary center, e.g. 30012.
It is different from a modulus-only test.

## 2. The slope 100 recipe works on the cube side only

Ordinally join three copies of the rooted \(H5\) graph
\((2,4,9,3)\). Mark vertex 1 of the **last** copy (global index 9).
The profile of a single root is \((5,1,2,3)\); the earlier
two copies multiply the profile by 25. Thus the \(5^3\) graph
has \((125,25,50,75)\) and a literal slope-100 pencil
\[
5^2+100\cdot2^k.
\]
The partner suggested by \(127=3^3+100\) would require a stable
\(H127\) root with \(r=50\), **allowing any** deletion count.

There is no such root at any order. Here is the reduction to a
finite exact check. A putative root is in the unique nontrivial
strong component after removing leading singleton padding (its
positive \(e=g\) pins it to the last component, and \(H=127\)
is prime). This removal leaves the profile unchanged. The next
near-sink append makes the result strong and has \(H=227\).
Since \(227^3<5^{11}\), the
[Busch bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
restricts the output to order at most 11, hence the input to
order at most 10. In a ready root the odd-path sum identity is
\[
s=(h+g)/2+\eta,\qquad \eta\ge0,
\]
where the singleton contributes \(g\), the closed nonsingleton
odd paths contribute \((h-g)/2\) by the odd-cycle identity, and
\(\eta\) is the remaining open odd-path sum. Thus \(g=h-2r+2\eta\ge27\).
Pad to order 10 with leading singletons. Deleting the mark gives
an order-nine base of count in [27,127]. The complete
[order-eight atlas](../actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv)
and all augmentations with count at most 127 cover those bases
by induced monotonicity.

For each base, the even weight \(r\) and the nonsingleton part
\(s-g\) are two subset sums of the outgoing contacts of the
mark: select the first vertex of a rooted path. This allows exact
screening before doing full DP. The replay finds 17,090 base
presentations, 84 masks with \(r=s-g=50\), and no \(H127\).
This is stronger than testing only \(g=27\). It does not exclude
the slope-120 pair above. Nor can any cross-only edit of the
cube of three intact \(H5\) components produce 127: the prime
count would force all twelve vertices of those components
into one strong component, beyond its Busch bound.

## 3. The same arithmetic is physical in other interfaces

### A square seam, an asymmetric seam, and path pools

For an ordinal pair with one reverse arc between marked vertices,
the exact count is
\[
h_1h_2-e_1e_2+R_1R_2
\]
when the second marked profile is taken in the converse. For the
one-apex word 1001, the marked profile \((h,e,R)=(5,0,10)\) gives
\[
125=10^2+5^2.
\]
In fact for \(a=2^j+1,\ j\ge1\), the word of binary address \(2a-1\)
has profile \((a,0,2a)\), giving \(5a^2\) at the square seam.
The ordinal cube gives \(a^3\), meeting this diagonal exactly at
\(a=5\). That is the exact reason these two graph recipes can both
present 125. It is also [Cayley's](../sumner_t5_petersen_collatz_20261005.md) \(K_5\) spanning-tree count
\(5^{5-2}\); that is a distinct enumeration, not a bijection here.

A symmetric diagonal with the same rooted graph on both sides
always gives \(1\bmod4\), since \(R-e=2r\), so \(R^2-e^2\equiv0\bmod4\).
Even allowing **unequal hosts with equal responses** cannot give 127.
Indeed \(R_i\ge h_i\ge e_i\) and the join with common response \(R\)
lies in \([R^2,2R^2]\). Only \(R=8,9,10,11\) need be considered.
Prime 127 forces the product of components preceding either mark
to be 1, since that factor divides \(h,e,R\) and hence the join
and \(R\). Strip their singleton pads. A near-sink append on each
trimmed host makes it strong at count at most \(2R\le22\), so the
host has at most five vertices. Pad before the mark to five
and check all \(2^{10}\) orientations. The output sets at common
response 8, 9, 10 and 11 are respectively
\(\{65\}\), \(\{105,123,125,153,159,161\}\),
\(\{125,145,165\}\), \(\{193,199,201\}\).
Thus the obstruction is not just a failure to try another root
of a 5-path witness.

The **off-diagonal** one-arc seam with words 11 and 11010 has
profiles \((1,0,4),(19,1,27)\) and does yield \(127=4\cdot27+19\).
The 27 survives as a physical response, not as the leftover
after a squared response. For this asymmetric slot there is a target-first inverse at 127.
With the left word 11 fixed, an even right address v must satisfy
P=4(v+1)+f(v). The bounds 1+v/2 ≤ f(v) ≤ 1+v give
(P-5)/5 ≤ v ≤ (P-5)/(9/2). At P=127 this leaves exactly one even
integer, v=26, and its exact residual verifies f(26)=19. So the
new scalar is selected from the target without guessing a host.
This interval can admit several candidates at larger P.
There is even a literal \(100+27\) in a different interface.
For the 127-path rows \((2,64,11,19,39,79,157,63)\), deleting
vertex 4 leaves exactly 27 paths. The one-head cut formula adds
exactly \(2\cdot50=100\) at that vertex. In the very same graph,
there is also a different path partition: exactly 27 full paths
use arc \(5\to3\), and 100 avoid it. The two interpretations count
different sets, and we do not identify them.
That does not give an intact ordinal
triangle cube. In fact no H127 witness contains an induced
\(C_3\to C_3\to C_3\): its strong core has at most 10 vertices
by Busch, so there is at most one extra vertex. Each path can
touch at most two triangles by edges incident to that extra
vertex. Cyclically rotating the least untouched triangle gives
free path triples, forcing divisibility by 3. Here it is the
[triangle contact mechanism](../berggren_target_return_20261111/arithmetic/TRIANGLE_ALL_K_PROOF.md),
not the numerical value 27, that decides whether the block fits.

There is a second direct path from the earlier H27 work to
the chosen center. Use the actual [27→53 source](../inverse_fences_landing_20261116/INVERSE.md)
with rows (62,124,120,112,32,64,17). Add a new head with
out-mask 54 or 100. Subset DP gives 125 or 127 respectively
on these two extensions of the **same** base. The one-head
cut sums are 49 and 50. Reversing the three head contacts
at pairs (1,7), (4,7), (6,7) exchanges the two counts.
Thus here 27+98 and 27+100 give center 27+99=126 literally,
but the H27 base is not the forbidden ordinal triangle cube.
This is independent of the persistent two-edit interface in §1.

An exact joint selector comes from holding contacts 2 and 5
outgoing and assigning \(u,v,w\) to contacts \(1,4,6\).
The eight corner checks give
\[
(H-27)/2=46+46u-18v+4w-25uv-6uw-5vw.
\]
Set \(u=v=t\), a synchronized two-arc switch. Then the selector is
\(46+3t+4w-11tw\). Its opposite corners give 49 and 50
with the same base and hence the requested midpoint 126.
The other corners have counts 119 and 111; a narrow diagonal
does not imply a narrow cube.

### A physical \(5^2\)–\(3^3\) swap at a different scale

For a nonempty intermediate spacer \(M\), orient \(A\to M\to B\)
and reverse one arc from the mark in \(A\) to the mark in \(B\).
The ordinary SCC paths lose nothing (the exposed arc skips \(M\)).
A path using the reversed arc splits each of \(A,M,B\) into
the portions before and after that crossing. Therefore
\[
H=F H(M)+G Q(M),\quad F=h_A h_B,\quad G=R_A R_B,\quad
Q(M)=\sum_I H(M[I])H(M-I).
\]
In particular \(G\ge F\); for a transitive \(k\)-spacer,
\(H(M)=1,Q(M)=2^k\). With strong one-apex outer modules marked at their apices \(G\)
is odd, so the proposed raw 127 split with this transitive
spacer would require \(F=27,G=25,k=2\), impossible.

More generally \(Q(M)\ge2^{|M|}+2(H(M)-1)\): expand the odd-cycle
packing formula in the two factors and assign each chosen cycle
as a whole to either part. A nonempty \(j\)-cycle pack with
\(\ell\) total vertices contributes \(2^{|M|-\ell+2j}\)
instead of its weight \(2^j\) in \(H(M)\). Each such term grows
by a factor at least 2. For \(H(M)=3,9,27\), use sizes at
least \(3,5,6\), respectively (the small-order maxima are
1, 5, 15 at order 2, 4, 5). Together with \(G\ge F\) this
already excludes \(GQ=100\) when \(F H(M)=27\) and \(H(M)>1\).
This is an obstruction for the proposed positive spacer
split, not a claim that 127 is missing.

But change 25 to the **common response product 35**, using
words 100 and 110 for \(F=5^2,G=5\cdot7\), or word 100010
and a singleton for \(F=27,G=35\cdot1\).
Both are literal ordinal seam constructions. This gives
\(25+35\cdot2^k\) and \(27+35\cdot2^k\) for every \(k\ge1\):
independently of the rooted slope-120 window, it gives a
physical way to swap the two intercepts at a common scale.
Its centers are \(26+35\cdot2^k\), beginning with 96 and 166,
not 126. The exact inverse is to test that \(C-26\) has odd
part 35 and positive dyadic valuation.

## 4. What the recurring numbers do and do not imply

The archived [Eisenstein viewpoint](../eleven_thirds_boundary_20261012.md) gives the identity \(127=7^2+7\cdot6+6^2=1+3\cdot42\)
contains the old midpoint 42 directly. Also
\(19^2+19+1=3\cdot127\). The equality
\((z^2+z+1)/3=z+4\cdot3^{2m-1}\) at \(z=1+2\cdot3^m\)
gives \(19+4\cdot27\) at \(m=2\), exactly the off-diagonal
port values above. At \(m=1\) the tentative core is the
forbidden 7. This is an arithmetic derivation of the candidate
tags, not a graph existence theorem at other \(m\).

A separate exact arithmetic catalog comes from
\(x^2-3y^2=-2\). Multiplying \(x+y\sqrt3\) by \(2+\sqrt3\)
gives the sequence \((1,1),(5,3),(19,11),(71,41),\ldots\).
These are all positive solutions: for \(y\ge2\),
\(3y/2<x<2y\) and multiplication by the inverse unit gives
a positive smaller pair, reaching \((1,1)\).
For any \(z\), counts \(z^2+x^2\) and \(z^2+3y^2\) would
share the exact center \(z^2+x^2+1\) with width 2. At
\((z,x,y)=(10,5,3)\) this is the user's identity at 126.
The next integers 19, 11, 41 recur in earlier lanes, but the
Pell catalogue only selects candidate integers, not physical
ports. Even freely varying \(z\) gives only
\(O(\sqrt B\log B)\) such centers below \(B\), by exponential
growth of the Pell sequence.

At 30000 scale, the [prescribed-prime construction](../prescribed_prime_30000_20261116/README.md)
showed that finding a width-two edge move and finding its
midpoint are different jobs. Here we can choose that midpoint
**exactly** on a real, infinitely renewable row, by changing
the deletion profile rather than the retained response.
The slope-100 no-go also shows why an excellent numerical
decomposition can fail to implement the requested physical
coordinates. No arbitrary-prime selector is claimed.

As a warning one row farther out, there is a two-edit 377/379 pair
with the same retained response r=186 and stable profiles
(g,r,s)=(23,186,209),(29,186,215). The two old arcs are (3,7)
and (4,6), at root 2 in the attached replay. This selects the same
depth-zero midpoint 378 as tripling 126, but its subsequent
center is 6+372·2^k, not 18+360·2^k. Their difference is
12(2^k−1). The new slope 372 is exactly three times the first
proposed 125 slope 124, not three times the actual centered
slope 120. Selecting the correct first integer can still choose
the wrong future slope. We do not claim a centered trim at
every ternary stage.

### Reduce the size of the restart seeds

In both literal two-edit pairs, the first two numbered vertices
are actually the last two near-sink appends at the common root.
Removing them identifies a smaller matched seed and preserves
the same two edit pairs on the retained vertices:

| Seed pair | Same \(r\) | Same \(b=s-g\) | After one append | After two |
|---|---:|---:|---:|---:|
| 31/33 | 17 | 13 | 65/67 | 125/127 |
| 89/91 | 51 | 42 | 191/193 | 377/379 |

These small seeds are not yet stable (\(b\ne r\)). They need only
match the observable pair \((r,b)\) to keep their width; the first
append puts each on the stable locus. If selecting a center
after two steps, the exact scalar to match is
\(C_0+4r+2b\), and the future slope is \(2(r+b)\).
This separates the two targets in a smaller source graph rather
than asking for the full blocked \((g,r,s)\) at the large count.

## Replay

    python3 verify.py
    python3 norm.py
    c++ -O3 -std=c++17 slope100.cpp -o /tmp/slope100
    /tmp/slope100 ../actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv

The Python checks use literal subset path counts, including
the two-arc repairs and centered trims. The C++ check is the
finite part of the all-order slope-100 exclusion.
