# An exact scaled-swap decoder for prescribed pairs of counts

**Proved.** For every positive multiple G of 30, the simple scaled swap
below realizes two windows of three counts at exact common-coordinate
gap G. For a prescribed G and absolute target N, there is an exact
inverse within this construction: one integer square root, and at most
one candidate address for the specified endpoint and center. The
possible endpoint density for each fixed G is zero. This is a scoped
exact pair decoder, not global Chan–Pak.

## A physical involution with a linear gap

Use the [literal six-site center bank](../chanpak_target_address_20261115/INVERSE.md).
Write its three counts as

    H_j(x,y) = 432xy + 129x + 258y + 15 e(x) + 30 e(y)
               + 45[p(x)+p(y)] + 135 + 2j,

with e the parity bit, p(2u)=p(u)+u+e(u), p(2u+1)=p(u),
0≤p(u)≤u. Restrict to an even left address and put
x=2u, y=v. The involution exchanges u and v; both sides are literal
tournaments, not signed formal combinations. The leading product
and the linear baseline remain identical, and the checksum gives

    H_j(2u,v) − H_j(2v,u) = 45(u−v)+15[e(u)−e(v)].

In particular set d=u−v≥1. For d even the gap is 45d; for d odd the
two parities of v give 45d−15 and 45d+15. Each positive multiple
G of 30 is obtained at exactly one d:

    G=90t     → d=2t,     either parity of v (t≥1);
    G=90t+30  → d=2t+1,   v odd;
    G=90t+60  → d=2t+1,   v even.

All three center choices survive the swap. One can also match different
center choices on the two sides and add an offset −4,−2,0,2,4 to the gap;
this is not a claim about the other congruence classes.

## Solve the absolute endpoint by a single square root

For either specified endpoint at center j, let

    B_j(v)=864v(v+d)+258(2v+d)+135+2j.

The exact endpoint lies in the band

    I_v=[B_j(v), B_j(v)+135v+90d+75].

For the upper endpoint its remainder is
45[p(u)+p(v)+u+e(u)]+30e(v), and for the lower swap u,v in that
expression. Both are nonnegative and satisfy the displayed upper bound.
The bands are pairwise disjoint and increasing: the baseline step is
864(2v+d+1)+516, which exceeds the entire allowance for v.
Consequently pick the greatest integer v for which B_j(v)≤N.
If v is nonpositive or its parity is wrong, there is no witness in this
lane. Otherwise test the precise checksum, not an interval surrogate.

Here is the exact square-root formula. Put A=864, B=864d+516,
C=258d+135+2j. If N<C, reject. Else the unique candidate is

    v = floor((isqrt(B²+4A(N−C))−B)/(2A)).

The floor is unaffected by replacing the square root by its integer
floor, since the divisor and subtracted offset are integers. For a
target odd N, j is determined by N modulo 3. This is an all-size
decider with physical witnesses; no unrestricted numerical census
is required.

There are Θ(sqrt(X)) endpoint values at most X for any fixed gap G,
since the corresponding endpoints have leading term 864v² and each
allowed parity class has positive density in v. Their spacing grows
on the order of v: the difference between consecutive allowed
endpoints is bounded above and below by constants times v+d+1,
using the disjoint bands and step of the baseline. Thus even exact
gap control at all scales does not place its endpoints densely.
The missing variable is a reader or source that moves B_j itself
while preserving the exact linear difference.

## What the old structures contributed

The earlier [square-sum absorber/binary shell work](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/square_sum_berggren_absorber_bijection_20261019.md)
keeps the host edge and scale as coordinates, not just the existence
of an insertable path. In the scaled swap, the gap likewise fixes the
difference coordinate, but the baseline depends on a second positive
integer; that is why the discriminant, parity and checksum tests are
all needed. This is an exact analogue of an inverse branch with an
explicit domain, not an isomorphism of the unrestricted odd spectrum
with the Berggren tree. Unlike the old target-first inverse for all
two-address graphs, this paired subproblem has no address enumeration.
