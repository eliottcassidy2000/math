# Figurate coordinates give overlapping physical center ranges

**Status:** exact two-center overlap at all sufficiently large depths
and at all balanced two-tape addresses described below. Both centers
have the same continuation interface. Not universal Chan–Pak coverage.

## The two geometric pictures

Put S=(X+Y)/2 and D=(X−Y)/2, recording their common integer or
proper-half-integer layer. The sum is 2S and the two doubles
are 2S±2D. The product is S²−D²; the mean of the squares
is S²+D², and their signed half-difference is 2SD. At fixed
sum, separation therefore loses exactly a square in the product.
Dilation acts linearly on the sum and quadratically on the
product. For a^X,a^Y the product becomes a^(2S) and their mean
is a^S cosh(D log a): the same coordinates but with the
multiplicative center and ratio distinguished.

The [Vandermonde truncation](../../../07-reflections/the-vandermonde-truncation-law-polygonal-vs-polyhedral-opus-S317.md)
shows that the polygonal entry i+(j−1)binom(i,2) is precisely the
first two layers of the simplex entry binom(i+j−1,j). Its discarded
layers are explicit positive binomial modes; summing rows or
shallow diagonals gives a controlled difference, not interchangeable
counts. The [figurate-stagger atlas](../../../01-canon/theorems/THM-1360-figurate-stagger-atlas-pisot-ladder-proth-mersenne-relational.md)
gets doubling, Fibonacci and supergolden growth from the same
Pascal array by different shears, and uses triangular numbers
as the exponent counting free orientations. These give three
distinct tools here: centering in sum/difference coordinates,
retaining the quadratic mode under dilation, and isolating the
boundary term whose variation was lost in the projection.
The [additive-basis note](../../../07-reflections/goldbach-polygonal-zeckendorf-additive-bases-s501.md)
warns that a polygonal basis theorem needs independent summands
or an exact carry; sparse atoms alone do not give a covering.

## A physical sum-product center with one signed boundary term

Use the [arm transfer](HALF_GRID_REPAIR.md), with signature
(h,g,B,U,V,Ubar,Vbar). Place the converse of a right transitive
tape of address X on the left, a core in the middle, and a right
tape of address Y on the right. The literal count is
```
h+BXY+Ubar X+U Y+g(C(X)+C(Y))−Vbar (X mod2)−V (Y mod2),
```
where C is the binary carry checksum of
[the joint decoder](JOINT_CENTER_WINDOW.md). This follows by two
applications of the arm law: the left tape contributes a start
read X, and the start read carried into the right tape is BX+U.
The signed term sees only the neighboring boundary. When U=Ubar,
the unsigned center is B(S²−D²)+2US. Swapping the addresses
keeps that part and the two checksum terms fixed; if X is odd and
Y even it exchanges the signed contribution. Actual outside
continuation can change under address swapping. Instead we will
reverse the signed contribution by reversing two caps at the
**same** addresses, preserving the final interface.

Use cap C on masks (0,5,17,7,43,79,159,63), root 3. Its profile is
```
(29,17,448,91,-17,91,-15).
```
Thus C and its converse have identical h,g,B,U,Ubar, but
different signed orientation. Let A be the old neutral
65/69-count pair with deletion count 17. Use a new 201/203-count
pair D with common signature (65,976,206,−28,570,−52), given by:

| H | outgoing masks | root |
|---:|---|---:|
| 201 | (8,33,35,70,15,153,55,95) | 1 |
| 203 | (8,33,19,38,75,149,47,95) | 1 |

Let W be an m-cell all-low word (A D)^m. Fold the core as
C W W^op C, where converse reverses the arm order and its arcs.
The second version uses C^op instead of C at *both* ends.
The all-low folded middle is self-converse, so the two all-low
cores are converse and have the same h. Since the cap's unsigned
reads equal each other, they have identical unsigned profiles;
in particular both have U=Ubar. To verify the last assertion,
the unsigned-profile transfer depends on only the unsigned data,
so replacing both caps by converse leaves that data fixed, but
also exchanges starts and ends by full core converse.

Put G=17·65=1105. In the folded core
g=17²G^(2m), B=448²(224·976)^(2m),
and the signed reads of the first cap orientation are
```
(V,Vbar)=17 G^(2m) · (−17,−15).
```
Thus at any X odd and Y positive even, reversing both caps
changes the center by exactly 34G^(2m). The final glued
signature (all six reads apart from h) is identical in the
two versions: the unsigned reads agree, and after the left
and right tapes are in place the signed reads depend only on
the outer tapes and the unchanged products of deletion counts.
If needed, append one ordinal sink to set the ordinary root end
count to zero. The retained-root continuation then agrees too.

## A literal overlapping pair, not merely overlapping bounds

Allow independent protected upgrades in all 2m A slots
and all 2m D slots of the folded word. Upgrading j D arms
and k A arms adds
```
P·(17j+130k),     P=2·17²·G^(2m−1),    0≤j,k≤2m.
```
The two center choices differ by 65P. Since 17 and 130
are coprime, the normalized repair palette contains every integer
in [2193,294m−2193] for m≥65. Choose j0=23r mod130 (in 0..129),
then k0=(r−17j0)/130. The lower endpoint 2193=17·129 makes
k0≥0. If k0>2m advance j by 130 until k (falling by 17)
first reaches ≤2m. The largest permitted j in that residue is
at least 2m−129 and gives k≤2m by the upper bound. The first
such k is nonnegative, since 2m≥130. This is an actual
bounded digit test, not just a congruence.

The two shifted palettes therefore share the full interval
[2258,294m−2193], relative to the lower center in units of P.
Their union certifies [2193,294m−2193+65]. Any target in the
intersection has two physical constructors with the same outside
continuation. At depth 65 the intersection is [2258,16917].
The assertion holds for **every** positive odd X and positive
even Y, not just a numerical coincidence or a prescribed ratio.
At m=0 the small shared-interface check (without the palette)
directly recounts the two centers 538419 and 538453, including
the retained-root state after a sink reset.

## Reading the absolute height

On the diagonal ray X=u,Y=2u for positive odd u, the lower
center simplifies, since C(2u)=C(u)+u+1:
```
H_0(u)=h+2B u²+3U u+g(2C(u)+u+1)−Vbar.
```
All additional shifts are nonnegative and less than B
for m≥1. Indeed their upper limit is P(294m+65), and its
ratio to B is
```
[2·17²/(448²·G)] · (G/(224·976))^(2m) · (294m+65) < 1/2.
```
This follows already at m=1 and the expression decreases.
Since U≤B−g, C(u)≤u and |Vbar|<g, the height above
h+2Bu² is less than 3Bu+3B, which is at most the difference
to h+2B(u+1)² for u≥1. Thus at a fixed depth the exact
integer-square-root address is u=floor(sqrt((N−h)/(2B))).
Reject an even address, then test the two signed centers and
their bounded repair digits. For u=2^a(z+1)−1 with a≥1,
C(u)=C(z) from the trailing-one rule. At fixed z this gives
an exact quadratic in 2^a, with modes 4^a,2^a,1; the cap
shift acts only on the constant mode. That is the physical
analogue of the figurate projection, not a guessed asymptotic.

The lock Y=2X is optional. Allow X odd and X≤Y≤2X with Y even.
At a fixed depth and trial X define
A_X=h+UX+gC(X)−Vbar and B_X=BX+U. Read the only possible
Y as floor((N−A_X)/B_X), then check parity, balance, and the
exact debit for each cap and repair choice. This division is
valid because the remaining debit is at most 2gX+P(294m+65)
< BX/4+B/2 < B_X, using g≤B/8. Every eligible X is at most
sqrt((N−h)/B), since BXY is one positive term. Thus there is an
exact decoder for a two-address family which need not be on the
doubling ray. It does not imply the existence of an address
for every absolute N.

## Another consequence: signed choices are a nearest-neighbor layer

Expanding the signed part of the arm transfer over a word gives
exactly
```
sum_i (product of g_j for j outside {i,i+1}) · V_i Vbar_(i+1).
```
The unsigned terms can see distant positions through products
of B, but the signed term always reads only the last old
arm and the incoming arm. For identical caps with either
orientation, put p=−16, q=−1, and write (V_i,Vbar_i)
=(p+q s_i,p−q s_i) with s_i=±1. The mixed term
pq(s_i−s_(i+1)) telescopes when the ends are fixed.
A switch in orientation adds 2q²17^(n−2) relative to a
same-orientation bond. With both ends fixed to the base
orientation there can be 0,2,4,...,2 floor((n−1)/2) switches,
all realized. Consequently an n-cap word, n≥2, gives the
protected interval H_base+4·17^(n−2)·[0,floor((n−1)/2)]
on its lattice. The endpoint ports and all other signature
reads are unchanged. This is an exact interaction-depth
analogue of retaining and discarding layers in the figurate
archive, and uses a smaller exponent in the repair step than
independent replacements of one ordinary count at a time.

**Scope.** The common interfaces and arithmetic overlaps above are
proved. The doubling ray alone stays sparse (O(sqrt(X)) counts
through X even summing depths, since B grows exponentially and
each address gives O(m) choices). The broader two-address family
still needs a theorem selecting a balanced absolute address;
a reflection always produces a neighboring center on the correct
local progression, not a universal inverse. The next useful
question is whether many *different products* can bring their
certified ranges into the same absolute band without paying an
incompatible residue or changing the required continuation.

Replay: `python3 folded_centers.py`; it uses the earlier
`joint_decode.py`, `reflection_replay.py`, and `fixtures.py`.
The large compiled graphs use the all-order transfer for exact counts;
the literal small graphs are recounted by subset DP.
