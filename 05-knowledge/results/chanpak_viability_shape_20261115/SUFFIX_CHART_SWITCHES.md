# An exact unit slide by transporting a binary suffix across the core

Use the fixed height-zero two-tape carrier from
[absolute steering](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md):
`F(x,y)=2855+5760xy+4034x+3914y+45[p(x)+p(y)]`, for all nonnegative
addresses. Both tape addresses remain physical when their lengths
change. In particular this construction does **not** add a second cyclic
core or assume an arbitrary new interface can read only F.

For any binary suffix b of length d, write L=2^d and r(b) for its d-bit
reversal. The exact surcharge split is

```
p(Lw+b)=p(w)+p(b)+(L−1−r(b))w+[b<L/2](w mod 2).
```

For both a,b in the upper half of the address block, the parity term
vanishes. The count with y constrained to a suffix is

```
F(x,Lz+b) = QLxz + (Qb+T)x + (SL+45(L−1−r(b)))z
            + B+Sb+45[p(x)+p(z)+p(b)],
Q=5760, T=4034, S=3914, B=2855.
```

Thus if Qb+T=SL+45(L−1−r(a)) and Qa+T=SL+45(L−1−r(b)),
swapping the upper addresses x,z and using suffix a gives the exact loss
`S(b−a)+45(p(b)−p(a))`. This proves, for instance, that the full cylinder
generating functions `G_b(w)=sum_(x,z>=0) w^F(x,Lz+b)` and G_a are
monomial translates, not independent term banks.

The two exact choices needed here are:

| Side carrying suffix | d | a | b | r(a) | r(b) | p(a) | p(b) | loss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Right | 28 | 182684216 | 184270779 | 29695093 | 232775157 | 162475172 | 24479448 | 2 |
| Left | 8 | 179 | 179 | 205 | 205 | 36 | 36 | 0 |

For the left lane simply transpose T and S before applying the test.
The two right matched coefficients are 1061399691074 and
1052261088194; the left matched coefficient is 1034954. The coefficient
and translation tests are integral equalities, not checked only modulo
some M. The checksum values follow from its binary recurrence.

## Transfer and exact windows at any depth

Set

```
A = 2^36,   U = 256*182684216+179 = 46767159475,
V = 2^28*179+184270779 = 48234217403.
```

First use the right move on suffix 184270779, then the left zero move
on suffix 179. The physical result is the exact identity

```
F(x,Az+V) = F(Ax+U,z) + 2      (x,z>=0).
```

Write U_k=U(1+A+...+A^(k−1)), V_k likewise (empty sum 0). Starting
from any x,z, take the initial pair (x,A^k z+V_k) and move the bottom
block of the right tape to the left k times. At step j the addresses
are `(A^j x+U_j, A^(k−j)z+V_(k−j))`, and their counts are **every
odd integer** from `F(x,A^k z+V_k)−2k` up to the initial value.
Both outer addresses remain free at every height. For example,
x=z=0, k=1 gives the adjacent counts 189077072195365 and
189077072195367. Thus this is not just an adjustable residue or
a k=1 coincidence.

The construction gives an exact window-viability generating function
for these moved-suffix windows. If v_k(N) counts ordered outer
pairs whose depth-k window contains N, then

```
sum_N v_k(N)w^N =
  (1+w^(−2)+...+w^(−2k)) sum_(x,z>=0) w^F(x,A^k z+V_k).
```

It attaches one finite geometric kernel to a full two-address bank,
not to a single dyadic center. This does not assert that every
coefficient is positive, or that two ordered pairs cannot encode the
same tournament or count.

## A precision tariff on this positive construction

At fixed k, solving in the quotient pair has aspect
`F(x,A^k z+V_k)=Q A^k xz+O_k(x+z+1)` with positive linear
coefficients, so the total ordered center presentations through X are
`(X/(QA^k))log X+O_k(X)`: compare with the divisor region with
coordinates translated by sufficiently large fixed nonnegative
integers, and handle the two axes separately. A periodic sublattice
version of the same comparison gives this main term directly.
For any fixed precision M the center counts are locally uniform:
varying the unrestricted x yields uniform residues in each growing
dyadic x-interval by the prefix-state mixing; the coefficient of x
after dividing (F−B) by 2 is a unit modulo gcd(45,M), for *every* z.
The geometric kernel also preserves that fixed-modulus leading
uniformity.

There is nevertheless an all-length capacity restriction. For k>=1,
V_k/A^k stays above V/A>0, and F(x,y)>=Qxy+Sy. Hence
`F(x,A^k z+V_k)>=c A^k(x+1)(z+1)` for an absolute c>0.
A window containing N<=X has center at most X+2k, so k=O(log X).
For k>=K>=1 the number of targets covered by these windows, even
counting every presentation separately, is bounded by

```
O( X log X sum_(k>=K) (k+1) A^(−k) )
 = O( X log X (K+1) A^(−K) ).
```

If K is the ceiling of (1+epsilon) log_A log X, this is o(X).
So this literal unbounded unit return does **not** give a uniformly
wide repair at all targets merely by increasing k. This is a
statement about this moved-suffix grammar, not a no-go for moving
cores, selective contacts, or combining these windows with a
different coarse selector.

## Replay

`python replay_unit_suffix.py` verifies the literal matched
coefficients, the two-tape identity over specified pairs, all the
count differences over repeated transports, and the complete
width-16 shift census. The proof, not the sampled addresses, gives
the all-order identities.

## Relation to the inherited 12-bit slide; a genuinely different 36-bit return

The general unit-return mechanism itself is inherited from the
[guarded unit slide](../chanpak_reservoir_slide_20261115/UNIT_SLIDE.md).
With its words w=010111000111 and z=000100011001, retain the central
guards and absorb the other guard positions into the tapes. The left
suffix is 1 followed by w without its last bit, hence U12=2787.
The right suffix is the reversal and complement of z without its first
bit followed by 0, hence V12=2875. It satisfies precisely
`F(x,4096z+V12)=F(4096x+U12,z)+2`, for arbitrary high x,z.
Thus this earlier slide already has the same two free outer addresses
and better bit cost. The window generating function and depth tariff
apply to it with A=4096 as well.

What is different about the new 36-bit move is its **composition from
two unequal transpositions** and its nonfactorization into shorter
global transports of this same kind. To make that statement precise,
call a suffix pair U,V at width d a global transport if
`F(x,2^d z+V)−F(2^d x+U,z)` is constant for all x,z>=0.
The split shows that the two linear coefficients must cancel
(test x and z in steps of 2), then that both top bits must be 1
(test the parities). The remaining necessary and sufficient tests are

```
128V + rev_d(U) = (2^d−1)(T/45+1),
128U + rev_d(V) = (2^d−1)(S/45+1).
```

Since T is coprime to 45, d is divisible by ord_45(2)=12. At d=12
the complete list (U,V,loss) is
`(2787,2875,2), (2795,2873,−48), (2803,2874,34), (2811,2872,−16)`.
This list and width divisibility are also related to the earlier
marked-core period calculation: no new minimal unit period is claimed.

Any factorization of the 36-bit address move into shorter global
transports must either begin with width 12 (and consume the bottom
12 bits of V) or end with width 12 (and write the bottom 12 bits of U).
Those bottom bits are 4027 and 2227, respectively, neither on the
required lists. Thus the new 36-bit move cannot be assembled from the
earlier 12- or 24-bit global transports, even though the alternating
transposition calculation above proves it in two steps. This is a
local grammar obstruction for the exact tape-coordinate transport,
not a statement about all ways to obtain the same integers.
