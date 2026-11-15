# All primitive global transports form a free repair language

This extends the [two-letter orbit law](WINDOW_ORBIT_SERIES.md)
using every primitive +2 chart, at all widths. A chart denotes the
physical identity F(x,A z+V)=F(A x+U,z)+delta for every x,z>=0.
The two boundary codes and the loss compose: consuming successive
V suffixes appends their corresponding U suffixes in the opposite
tape's concatenation order, and widths and losses add.

## Exact cuts and primitive suffix codes

Use the [unique ternary word theorem](ALL_WIDTH_TRANSPORTS.md).
At width 12k the free ownership word, read in the V low-to-high
direction, is represented by e=P^(k−1)R with
P=110222011212 and R=11022. At a putative cut after j periods,
1<=j<k, the shorter consumed V block must have top seven bits
89. In the last P period of that block, these are positions
5,...,11 read from low to high. On these positions e is
(2,0,1,1,2,1,2); 89 has low-first bits (1,0,0,1,1,0,1).
Thus the required V ownership at the three relevant free
positions (7,8,10) is exactly (0,1,0).

Conversely, at this marker the low V block has its top seven
bits 89 and the following data really split. On the retained
lower field of the first factor one sees P^(j−1)R, because
the first five positions of the cut P are R. The mirrored
high U block has top seven bits 87 inherited from the whole
chart. In the second factor the high V block still has top
seven bits 89, and the top seven positions of the low U block
are the mirrored complements e−V on positions 5,...,11:
read from high to low they are (1,0,1,0,1,1,1), the binary
digits of 87. The retained data are P^(k−j−1)R. Thus both
parts are legitimate charts. Requiring only the shorter U's
top value forces the very same ownership marker.

Therefore cut exactly at all marked periods. This is a
unique factorization into nonempty primitive transports,
and the primitive U codes and V codes are each **suffix-free**
by the corresponding shorter-block marker test. There is at
most one longer inverse or forward arrow at a given pair,
even when **all** primitive +2 moves are included.

## Primitive shift generating function and the actual length law

Let t mark units of 12 bits in a chart and Z its signed loss.
Define
```
A(Z) = Z^(-16) (1+Z^(-32))(1+Z^50),
C(Z) = Z^16 (1+Z^(-32))²(1+Z^50)²(1+Z^(-68)).
```
The full nonempty transport series from the width theorem is
`T(t,Z)=t A(Z)/(1−t C(Z))`. By free factorization its primitive
series is exactly
```
P(t,Z)=T/(1+T)=t A(Z)/(1−t(C(Z)−A(Z))).
```
It has 4*28^(k−1) terms with multiplicity at width 12k.
The first eight +2 primitive multiplicities, from [the replay](replay_primitive_language.py),
are 1,0,178,0,76140,0,40625112,0. The selected alphabet of
all primitive +2 arrows is infinite. For example, C−A has both −16 and +16 as period losses,
with strictly positive coefficients. Appending these two
allowed primitive periods repeatedly to any base term of
loss +2 preserves that loss, yielding infinitely many widths.

Let rho be the total digit-cylinder density of the selected
alphabet:
```
rho = [Z^2] P(1/4096,Z) = 0.0002441432153069537... .
```
An exact way to bound the decimal is truncation at eight
widths: the uncounted tail is at most
`4*28^8/(4096^9*(1−28/4096))<4.7e−21`.
The total primitive code density (all losses) is 1/1017.
The no-suffix address sets have actual natural density 1−rho:
for a cutoff Y, a code with radix A>2Y has no occurrence
since its top bit is 1. Each remaining suffix contributes at
most Y/A+1<=3Y/A to the count. This bounds long-code errors
by the tail of a convergent series, and finite code cylinders
have exactly their expected natural densities.

Maximal consecutive chains defined by the full primitive +2
alphabet are still finite, canonical and have all the same
endpoint recurrences as the [two-letter law](WINDOW_ORBIT_SERIES.md).
The endpoint/window GF holds verbatim after replacing the two
letters by the entire primitive unit alphabet. Their length law is
`(1−rho)^2/(1−rho t)^2`, where t now counts the number of
edges, not bit-width. To justify the real factor for infinitely
many letters, use that a set of natural density 1−rho has
harmonic density 1−rho by partial summation. For each fixed
word of letters the hyperbola proof gives the same leading term
as before, with its own inverse product of radices. Uniform
large-letter tails are controlled by the sum of those inverses:
the leading suffix block is always at least half its radix,
so H>=c A_w (x+1)(z+1) uniformly, without any uniformity in a
fixed-word error term. This proves the exact limiting PGF.

This is an infinite menu of **literal adjacent-count corrections**
and a generating function for its maximal windows. It still
does not produce positive coefficients for every target. The
[width tariff](SHIFT_MULTIPLICITY_AND_TARIFF.md) explains why
adding longer corrections while retaining the same reader
cannot replace the remaining absolute placement theorem.
