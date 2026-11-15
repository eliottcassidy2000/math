# Absolute steering by paired binary tapes

These are **physical** approximations from above and below, not a proof of
Chan–Pak, and not a claim that a finite target passed merely a residue test.

## Fix the carrier; change both independent tape addresses

Use the height-zero [unit-slide carrier](../chanpak_reservoir_slide_20261115/UNIT_SLIDE.md).
Its painted eight-vertex core has outgoing rows
`(4,1,130,7,143,159,63,75)`, mark mask 35, a leading 1 guard
and a terminal 0 guard, all joined in ordinal order. Add one apex
beating the marked vertices. A left binary address x is the ordinary
binary mark string; a right address y is the little-endian binary
string with 1 meaning unmarked. Empty tapes are permitted. The first
and last actual spine marks are 1 and 0, so the full tournament is strong.

The inherited two-tape formula here reads
```
F(x,y) = B + Qxy + Tx + Sy + h[p(x)+p(y)],
(B,Q,T,S,h) = (2855,5760,4034,3914,45).
p(0)=0, p(2u+1)=p(u), p(2u)=p(u)+u+(u mod 2).
```
In particular 0 <= p(u) <= u. These are counts of different physical
tournaments, not values of a signed surrogate. A fixed one-sided address
has the other tape as an essentially linear address, with bit-reversed
variation of its slope. Simply rounding a one-sided address does not
control the checksum. Two tapes give a direction where the product is
stationary to first order.

## All-order target-first bracket: error O(N^(3/8))

For every sufficiently large N there are literal counts L < N <= U
from this *fixed* carrier with N-L and U-N both O(N^(3/8)). There is
an explicit decoder using only O(log N) bisection steps on integers
of O(log N) bits. It does not enumerate graphs or factor N.

1. Put a = ceil(sqrt((N-B)/Q)), choose a power q=2^r of order a^(1/4),
   and write a=qu+v, 0<=v<q.
2. Let x=a+qi, y=a-qi, i=0,...,u. The two addresses have the same
   low bits, so the reversed low-word phase is *common*.
   Write c=q-1-rev_r(v) and delta=1-bit_(r-1)(v), using r>=1.
   The [binary split identity](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md)
   gives exactly
```
p(x)+p(y) = 2p(v)+2uc + p(u+i)+p(u-i)
            + delta[((u+i) mod 2)+((u-i) mod 2)].
```
3. Thus F(x,y)=G(i)+E_i, where
```
G(i)=B+Qa^2+(T+S)a+h[2p(v)+2uc] - Q(qi)^2+(T-S)qi;
0 <= E_i <= D=h(2u+2).
```
   Use the first i with G(i)<N-D for the lower witness and the
   last i with G(i)>=N for the upper one. Each cutoff is found by
   binary search: for either threshold the nonnegative i in the
   upper level set form an initial interval, by concavity and
   G(0)>=N. At i=u the product term is Qv(2a-v)=O(aq),
   and the other terms are O(a), so G(u)<N-D for large N.
4. At either threshold, solving the concave quadratic shows qi=O(sqrt a),
   also at the first index beyond the threshold since q=O(a^(1/4)).
   Indeed G(0)-N=O(a), D=O(a/q+1), and the linear coefficient in qi
   is a fixed constant. A single index step changes G by O(q sqrt a+q^2).
   The physical witnesses are therefore within O(q sqrt a+q^2+a/q+1)
   of N. Choosing q of order a^(1/4) balances the two dominant terms
   and yields O(a^(3/4))=O(N^(3/8)) on both sides.

The small residual is not automatically a multiple of a repair
step, nor within the O(log N) width of the known slide bank. One
cannot substitute this bracket for exact placement. It is an all-scale
way to tune the center without enlarging or changing the cyclic core.

## Second redesign: allow all low bits to vary but freeze their carry

For any x+y=s, write c_j for the binary carry into position j. The
sum of the two bits at position j equals
`bit_j(s)+2c_(j+1)-c_j`. Thus a carry word, together with s, fixes
the sum of the bit-reversed fractions. It also specifies a full cube
of physical choices: wherever this value is one, allocate that bit
to either address. Where it is zero or two, allocation is forced.
Conversely those choices all give s and the same carry word
(column addition proves this from c_0=0). This is a literal
two-strand analogue of fixed overlap in the tournament frame work:
information about the *choice* of strand can be discarded without
discarding information about what the two addresses jointly read.

For a first application take a power q=2^r, r>=1, any u>=0, and
```
x=qu+v, y=qu+(q-1-v), 0<=v<q.
```
The low bits split the all-ones word without carries. The large
low-phase terms cancel exactly:
```
p(x)+p(y)= u(q-1)+2p(u)+p(v)+p(q-1-v)+(u mod 2).
```
The first two terms are constant, and the remaining sum lies in [0,q].
With z=v-(q-1)/2 and a=qu+(q-1)/2, the smooth part of the resulting
count is B+Q(a^2-z^2)+(S+T)a+(T-S)z+h[u(q-1)+2p(u)].
The exact count differs from this by at most hq. Let L be the actual
value at v=q-1 and U at v=q/2. Then
`U-L = Q(q^2-2q)/4 + O(q)`. Consecutive v in this walk change the
actual count by at most (Q+|T-S|+h)q: use the bound on the smooth
difference and the interval [0,hq] for the error. Consequently every
target in [L,U] is bracketed within that step bound (no monotonicity
assumed, just first crossing). For u in [q,2q), these are infinitely
many bands of width of order sqrt(N) and mesh O(N^(1/4)), at heights
N of order q^4. Their dyadic centers are not universal. The result
explains how to lower the scale further when a target is already
near a suitable shell: freeze the **carry**, not the individual
low bits. Seeking an absolute selector for general carry cells is
a different problem from extending the original unit slide.

## Reproduction

`python replay_steering.py` (also with `-O`). Eight examples use
the literal subset-endpoint path counter, independently of the tape
formulas. The large cases check only the exact integer decoder and
unknown-target comparison, not any broader claim about the spectrum.
