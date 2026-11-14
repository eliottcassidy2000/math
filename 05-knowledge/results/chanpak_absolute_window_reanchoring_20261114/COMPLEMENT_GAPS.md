# Every gap in the complementary collar, and a dry band at every dyadic square

**PROVED:** the full fixed-slice gap formula and an infinite family of
square-adjacent bands poorly covered even by the union of all widths
of the complementary construction. This is a limitation of that chart,
not tournament nonexistence.

Use the physical complementary chart from
[TIGHT_COMPLEMENT](../chanpak_tight_mod_six_20261114/TIGHT_COMPLEMENT.md).
Let M=2^s, s>=4, u>=0, 0<=q<M/2. Its two addresses have a common high
prefix u and complementary length-s tails a=M/2-1-q, b=M-1-a.
In the actual contact word, the central 2s bits are the palindrome
a followed by rev_s(a). Inherited notation gives

    H(q) = B(u,s) + 2g(q)-q²-2q+p(q)-2 top(q),
    B(u,s)=M²u(u+1)+M²/4-u+2g(u)+p(u)-3,

where top(q) is the leading bit of the declared width s-1, and p(q)
is its last bit. All q except the u=0, a=0 boundary yield a strong
tournament. No restriction on u beyond nonnegativity.

Let k be the number of trailing ones of q and v=floor(q/2^(k+1)).
For 0<=q<M/2-1, the exact positive gap is

    H(q)-H(q+1) = 6v + 2 + 2p(v) + 2(top(q+1)-top(q)).       (1)

For k=0, g(2v)-g(2v+1)=v+p(v) gives the formula directly. For k>=1,
q=2^(k+1)v+2^k-1 and q+1=2^k(2v+1). Removing the k terminal ones
or zeroes gives g(q)=g(v)+v+p(v) and
g(q+1)=g(v)+(2^k-1)(2v+1)+1. Substitution gives (1).
This proves sharp monotonicity, not only the estimate that the
gap is at least two.

Gap two occurs exactly at q=2^j-1 for 0<=j<=s-3:
v must vanish, and the next candidate j=s-2 crosses the leading-bit
boundary and has gap four. Therefore the only three-point step-two
window in a fixed (u,s) slice is q=0,1,2, and there is no four-point
window in any slice. The other two-point windows are exactly the
dyadic-carry pairs listed in the older construction. The existence
formulas are inherited; exhaustion of the list is the point here.

## The union across all widths still leaves an explicit dry region

Let A=2^r, r>=4. Consider the A-1 odd numbers

    I_r = {A²+3, A²+5, ..., A²+2A-1}.

Across **all** values of u and M, at most 2(r-1) members of I_r
are hit by the complementary chart. In particular, I_r contains
at least ceil((A-2r+1)/(2r-1)) consecutive missed odds
(in step-two order).

For N in this band, floor(sqrt(N))=A. By the inherited square-band
inequality M²u²<H<M²(u+1)², any representation has u=floor(A/M).
If M<=A, then u=A/M and the inherited lower bound
H>=M²u(u+1)-u-1=A²+AM-A/M-1 exceeds A²+2A because M>=16.
If M=2A, u=0 and its maximum is A²+1.

All remaining widths have M>=4A and u=0. At a>=M/4,
the product a(M-1-a) alone exceeds A²+2A; it is increasing
up to M/2-1. Thus any hit has q>=M/4. Every second gap of this
part of the chart, at even q, is at least 3M/4+2 by (1).
An interval of length less than 2A cannot therefore contain three
outputs from one width. Finally a>=1 is necessary for H>1, and
a(M-1-a)>=M-2, so M<=A²+2A+2<2A². The widths 4A through A²
are exactly r-1 powers of two. This proves the count bound; the
missed odds split into at most 2r-1 runs. These are pinned bands
near explicit squares, not just gaps obtained by an average-density
argument. Other tournaments may supply every one of the misses.

Replay: python3 complement_gaps.py. It verifies the exact gap formula
and the equality classification on all offsets at widths 4..13 and
performs the independently decoded dry-band test through r=14.

## Adding the free-height four-point bank still does not fill these bands

The [linear bank](LINEAR_WINDOWS.md) has four offsets at every positive
integer u, but even it cannot be added to a finite set of dyadic
charts to fill these particular bands. Here is a bound on how much it
can cluster locally, not just a global count bound.

Write F(u)=908u+9g(u). For A=2^r, r>=8, choose the least power L
of two with 908L>=5A. Then L<10A/908. Split each u=Lv+b with
0<=b<L. Any preimage of a point in I_r using any listed offset
has v<=(A²+2A)/(908L)<=(A+2)/5, since g and the offsets are positive.
The exact suffix identity says

    F(Lv+b) = (908L+9(L-1-rev(b)))v + 9g(v) + 9 e(b)p(v) + c_b,

where e(b) is zero or one. For v'>v in a given tail b, their
centers differ by at least

    908L - 9v - 9 >= 5A - 9(A+2)/5 - 9 > 2A+6,

using 1<=g(v)<=v+1. The four offsets span 6, so among targets
in I_r at most one v can work per b, and contributes at most four
points. Thus the free-height bank contributes less than 40A/908
points, about 4.41% of the A-1 odd numbers in this band.
The two-axis dyadic bank adds O(r²) points, and the previous
fixed affine/quadratic dyadic banks add O(r). These specific
families, even combined, miss a positive proportion of every
sufficiently large I_r. This does not bound the full spectrum.
The repair must change the physical family, not only increase
the number of addresses allowed in this four-point bank.
