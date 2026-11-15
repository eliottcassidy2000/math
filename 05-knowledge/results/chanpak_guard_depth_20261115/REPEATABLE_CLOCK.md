# Repeatable absolute digits from an intentionally different reader

The [mixed-radix physical interface](../chanpak_outlet_blocks_20261115/MIXED_RADIX.md)
gives closed state ((D,C,N,z)), where $N=H+2E$. One result of this
session is a physical adjacent pair that locks the outgoing continuation
at one initialization, but not at any positive incoming address. The
following different pair gives an indefinitely repeatable branch. It
uses a deliberately *different* response to the already accumulated
address. Equal response, though natural, would erase the positional
information of repeated choices.

## The two literal patches

Use $K$ with outgoing rows ((8,1,3,6)), $H(K)=5$,
$P_K(2)=32$, and five singletons. Place them in either arrangement:
```
A: singleton, singleton, singleton, singleton, singleton, K;   marks 110110
B: K, singleton, singleton, singleton, singleton, singleton;   marks 100000
```
As in the mixed-radix model, keep the spine in ordinal order and let the
same apex beat the entire marked blocks. Place a final outlet beating
only that apex when the construction is finished. Both patches start
with mark 1 and end with mark 0. In particular, the initial state and
all inter-patch boundaries have $z=0$.

Composition of the literal block matrices gives
$$
(D',C',z')=(5D,1024C+864D,0)
$$
for both. Their two exact count readers are
$$
 N'_A=5N+1903C+1618D,\qquad N'_B=5N+2011C+1706D.
$$
Thus the construction can switch freely, with **identical continuation
coordinates at every depth**. The switch makes precisely the difference
$108C+88D$, and reading the switch in a new position gives it a new
weight. At one patch the outputs are 1623,1711. This is not the
adjacent initialized pair of the companion note.

## Positional weights, exact inverse, and scale budget

Start $(D,C,N,z)=(1,0,1,0)$. After $j$ patches
$$
 D_j=5^j,\qquad C_j=864(1024^j-5^j)/1019.
$$
Write $N_{A^k}$ for the all-A baseline at depth $k$.
The entire spectrum at that depth in this family is exactly
$$
 N_{A^k}+\sum_{j=1}^k \epsilon_j W_{j,k},\qquad
 W_{j,k}=5^{k-j}(108C_{j-1}+88D_{j-1}),\quad \epsilon_j\in\{0,1\}.
$$
Every word gives a strong tournament: the apex goes to all of its
first block, all subsequent blocks are accessible in ordinal order,
and the unmarked last block goes back to the apex; the final outlet
is also in this component. Thus these are genuinely physical digits,
not integer vectors assumed realizable.

The weights are superincreasing. In fact, $A_j=C_j/D_j$ has update
$A_{j+1}=(1024/5)A_j+864/5$, so
$$
 (88+108A_{j+1})-\tfrac{1024}5(88+108A_j)=728>0.
$$
Hence every weight is more than 204 times its predecessor.
The $2^k$ outputs are distinct. For an absolute target at a
prescribed depth, subtract the all-A baseline and greedily decode
from the largest weight; accept only if the residual is exactly
zero. At depth at least two the gcd of the weights is exactly 8:
the first two, after dividing by $5^{k-2}$, are (440) and
(93752) with gcd 8, and the last is not divisible by 5.
But clearing that congruence plainly does not imply membership.
The $2^k$ counts have size $O(1024^k)$, giving an unbounded
polynomial-size, but far from full, set of exact baselines. Conversely
the all-A minimum is bounded below by a positive constant times
$1024^k$ once $k\ge2$, so the union up to $X$ has
$\Theta(X^{1/10})$ distinct counts. The upper bound follows by
summing $2^k$ over possible depths and the lower by choosing the
largest depth for which the whole shell fits below $X$.

## Why matching one more coordinate blindly is counterproductive

In any two patches with common transformations of (D,C,z) and
common direct multiplier $h$ on $N$, write
$\Delta N=\Delta \beta\,C+\Delta\alpha\,D$.
If $\Delta\beta=0$, a choice at position $j$ contributes
$\Delta\alpha\,h^{j-1}D_0$ immediately, and therefore contributes
$\Delta\alpha\,h^{k-1}D_0$ at final depth $k$, independent
of position. Even a perpetually legal normalized unit switch
of that type gives at most $k+1$ different counts.
Here $\Delta\beta=108$ makes the choices positional.
The remaining design problem is to keep that leverage but fill a
much larger digit range: the growth base (1024/5) dwarfs the
current two choices. The local gcd already improved to 8 without
closing that absolute gap.

### An alphabet-size ceiling for one movable core

Fix the sequence of ratios $r_i=P_i/h_i\ge2$.
The normalized outgoing address at incoming $C=0$ is
$$
\sum_i w_i(r_i-1)\prod_{j>i}r_j .
$$
It is injective in the marks: the contribution at the first
differing position exceeds the total of all later contributions,
since that later total is $prod_{j>i}r_j-1$.
Thus an alphabet sharing an outgoing address needs distinct
ratio sequences, not just different mark words on one sequence.
For one fixed strong core and $l$ otherwise singleton blocks,
there are at most $l+1$ such choices at any common address,
whereas the normalized scale multiplier is
$(P_K/H(K))2^l\ge 2^{l+1}>l+1$.
Even exhausting all marked labels in this moving-one-core family
cannot build a filled positional interval by repeating one matched
address. It must either admit several independently movable
strong profiles, change contact geometry, or overlap with a
different normalized clock and absolute baseline.
