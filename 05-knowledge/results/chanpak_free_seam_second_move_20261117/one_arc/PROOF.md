# A genuine one-arc second switch, and where it runs out

Use the [free reflected seam](../../chanpak_wheel_square_selector_20261117/FREE_SEAM.md)
with (h,P,L,R)=(5,256,142,142). Write the lower orientation with
u positive odd and w positive even. Its count is

    H(u,w)=91+5(F(u)+F(w)-2)+142(u+w)+256uw.

Swapping u,w gives the other orientation and adds 2. We will
toggle **one actual apex arc** in the even flank of the lower
orientation and its corresponding arc in the reflected upper
orientation. Both switches survive for arbitrary positive addresses,
because the toggle is not in their least or most significant bits.

## An infinite one-arc switch of any prescribed even size

Take w=4a, with a positive, and toggle its second-lowest bit to
make w'=4a+2. The same single apex adjacency is reversed in that
binary flank; on the upper side it is reversed in the reflected
flank. The common horizontal change is

    delta = 512u+294-5(a+[a odd]).

Indeed F(4a)=F(a)+3a+[a odd] and F(4a+2)=F(a)+2a+2.
For any prescribed even delta, choose a positive odd u with
512u+294-delta divisible by 5. There are infinitely many such u.
Set a=(512u+294-delta)/5. This is a positive even integer once u
is large, and gives exactly delta with no parity adjustment.
The resulting four counts are

    H, H+2, H+delta, H+delta+2.

This is a physical commuting square, not a comparison of
unrelated existence witnesses. All four tournaments have the
same order and remain strong. In particular delta=-4, u=1,
a=162 gives exact four consecutive odd counts
261163,261165,261167,261169 on 19 vertices. A subset-DP replay
checks the actual four adjacency matrices.

The operation requires a relation between the two high addresses.
It is not a free delta applied after an arbitrary initial placement.

## Count of potentially useful one-bit moves

Consider every one-bit flip that preserves the two flank parities
and leading bits of this fixed core, on either flank. At bit j>=1
let Y=2^j, and write the changing address as
2Y a+r, or 2Y a+Y+r, with 0<=r<Y and a>=1. Write w for
the other positive flank address (the roles can be reversed).
The signed difference in the 0-to-1 direction is

    delta = PYw + 142Y + 5(RC_j(r) + g_j(r) - a - [a odd]).

Here RC is reverse complement in exactly j bits and g is 1
if the leading bit of that j-bit suffix is zero. This follows
by first flipping the low bit of 2a, which changes F by
-a-[a odd], and then applying the j-digit suffix law. In
particular a near-perfect cancellation requires a very skew
split, a about (P/5) Y w even when j=1.

**Bound.** Fix X>=2, 0<=E<=sqrt(X). Among these literal one-bit
edits, the number of possible endpoint count values at most X
that belong to moves with |delta|<=E is O((E+1)sqrt(X)log X).
The absolute constants are fixed by the core; multiple representations
can only lower the number of endpoint values.

Proof: a >= (P/5)Yw - E/5 - 1. If the last two terms cost
at most half the first, a >= (P/10)Yw. On either endpoint the
physical count is at least P(2Ya)w, so Yw=O(sqrt(X)).
In the other case Yw=O(E+1), again O(sqrt(X)). Thus for each j
there are O(sqrt(X)/Y) possible positive w, and only O(log X)
possible j. There are at most Y r's. For fixed (j,w,r),
the delta constraint puts a+[a odd] in an interval of
length 2E/5, with at most two preimages per integer, hence
only O(E+1) possibilities for a. Sum the product, and multiply
by a fixed factor for direction, flank and endpoint.

For example if a separate selector first reaches an absolute
error O(N^(1/3)), the endpoints accessible by one such fine
edit still have density zero: use E=O(X^(1/3)), obtaining
O(X^(5/6)log X) endpoint values below X. The same argument
also excludes an O(X^(2/5)) residual strategy relying on only
one such arc. This result does not forbid several coordinated
edits, a change of core, parity-changing switches with a new
reader, or unrelated Chan–Pak realizations.
