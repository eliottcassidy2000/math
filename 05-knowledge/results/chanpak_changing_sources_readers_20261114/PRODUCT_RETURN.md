# A four-odd window with two heights and an exact neutral return

Let A=TT_k and B=TT_m, and take the cyclic source S=C3[r,A,B]
in that order. Form S→R, reversing just the edge from r to the
entry w of the reader R. Put a=H(R), let q_e,q_o be the weighted
sums of even/odd paths *ending* at w (weighted by the count of the
deleted complement), and let f=a+q_e+q_o. The [source projection](
SOURCE_SWITCH.md) gives
[
 H=(2^k-1)(f2^m-a)+2a+2q_e.
]
The following four literal readers, discovered in the pointed
reader search, all have entry 6, a=73 and f=146. The numbers in
the two right columns are available from subset path DP.

| outgoing row masks, zero-based | q_e | q_o |
|---|---:|---:|
| 8,5,9,82,39,15,183,63 | 27 | 46 |
| 2,8,3,53,39,71,159,63 | 28 | 45 |
| 2,4,9,51,39,71,159,63 | 29 | 44 |
| 2,12,9,81,7,31,183,63 | 30 | 43 |

Consequently at **every pair of positive heights** k,m there are
four consecutive odd counts
[
 73(2^k-1)(2^{m+1}-1)+\{200,202,204,206\}.
]
The first window is 419,421,423,425. An exact target inverse
chooses one of the four offsets, requires that the positive
remainder divided by 73 factors as two Mersenne numbers with
exponents k>=1 and m+1>=2, and then uses the literal corresponding
reader. This applies only in this window family.

## Exchange the exposed height, keep the absolute count

When k>=2 and m>=1 the map
[
 (k,m)\longmapsto(m+1,k-1)
]
is a literal involutive rebuild of the source with both slabs still
nonempty. By commutativity of the two Mersenne factors it preserves
H for each of the four readers separately. The total order k+m+9
is identical. More than H is preserved: the number of full
Hamiltonian paths starting at the retained root r is a, and the
number ending there is the number of full reader paths ending at
w. These claims hold on both arms. A full path starting at r
cannot use the backward seam w→r and must exhaust S before R;
its full source part is unique. A full path ending at r must use
w→r at the end, and exhaust S−r before R since it cannot return to
earlier source vertices.

Despite this agreement, the rooted **weighted** path budget at r
changes. Let W_R(w) be the sum of Hamilton counts after deletion
of paths *starting* at w in R, and put Q_R=P_R(2). If W(T,r) is
the sum of complementary Hamilton counts over all paths starting
at r in the glued tournament, then
[
 W(T,r)=[1+(2^k-1)2^m]Q_R-W_R(w).
]
Indeed a path starting at r has a rooted source part and possibly
a nonempty reader suffix. Its complement has no backward seam
because r is gone, so complement counts factor. Summing the
optional reader suffix over all subsets gives Q_R. The only
excluded crossing is directly from the one-vertex source path
(r) to w. Its removed weight is W_R(w) times H(S−r)=1. The
unseamed rooted source weight is 1+(2^k−1)2^m by the cyclic
source formula. Subtracting after the height exchange yields
exactly
[
 W(T,r)-W(T',r)=Q_R(2^{k-1}-2^m).
]
The difference vanishes only on the fixed diagonal k=m+1.
Thus a window can have two literal sources of the **same size and
absolute count, with the same two full endpoint counts**, yet
with different continuation information at the retained root.

The difference can be read by physical continuations without
assuming any converse-endpoint identity. Append an order-two
reader TT2 by one reversed seam at r. Reversing into the sink
has (a,q_o,q_e)=(1,1,1); its output is H+2E_r+2O_r=H+2W(T,r).
It distinguishes the two rebuilds at every nonfixed pair of
heights. The replay directly checks this third stage on literal
tournaments at (k,m)=(2,2) versus (3,1).

This is a target-preserving return **within** a physical
Mersenne-product family, not an all-target selector. The geometry
is the same compression lesson as in the older
[Green-conductance sidecar](
../../../07-reflections/lrc14-green-conductance-connectivity-codex-20260628.md):
two observations can agree exactly, while a different observation
retains exchange data. Here the changing coordinate is a concrete
rooted path weight, not an analogy in place of the calculation.
