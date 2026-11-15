# A sharp matched repair of three odd half-steps

By [endpoint reflection](ENDPOINT_REFLECTION.md), for the painted
right-block interface with ordinary count h, marked unsigned
starts s and **unmarked** unsigned ends t, one has
`t≡sigma(s+h)−h mod 3`, where sigma=(-1)^m. Indeed the total
weighted ends are Q(2)−h, Q(2)≡sigma h, and the weighted
end sum over the marks is −sigma s. Thus equality of h,
order parity and marked starts mod 3 fixes the right end
residue even if Q(2) itself is not matched. Holding its signed
end channel fixed as an integer then gives the additional
factor two. In this matched end-only channel, the actual
repair moves by three places on the odd-count lattice.

More generally, if two rooted start polynomials of equal order
parity agree in their first k Taylor/Hasse coefficients at a=1,
the resulting terminal unsigned-end change is a multiple of
3^k, and its next digit is the signed difference of the k-th
coefficients. The signed-end lock adds the factor two.
The coefficients count actual rooted odd-block deletions, so
this is a stronger interface requirement, not freedom to choose
the new digit independently.
## Sharp physical step on an exact common continuation state

The two order-6 cores below mark only their designated root.
Use h,Q,S,T for H,Q(2),weighted starts,weighted ends,
and d and l for the full end and full start at that root
(the two opposite signed readouts).

| Outgoing bit rows | root | h | Q | S | T | d | l |
|---|---:|---:|---:|---:|---:|---:|---:|
| (52,17,2,7,12,30) | 4 | 29 | 320 | 43 | 50 | 5 | 4 |
| (42,20,1,6,13,30) | 3 | 29 | 320 | 43 | 56 | 5 | 4 |

Prepend a single marked vertex to the rooted right core,
in ordinal order, and add one apex beating the marks. The
local apex count of the right core is h+S−d=67. The entire
count by the [exact seam transfer](../chanpak_selective_target_20261114/CUT_STATE.md)
is
```
67 + (Q−h−T) + (h−l) = 333 or 327.
```
The two assembled tournaments are strong by the mixed
root contact and the strong cores. Their common outgoing
state is (h,S,d)=(29,363,5), so any common transitive
continuation keeps the difference exactly 6, allowing both
heights and orders to grow uniformly. This is a literal
sharp first-digit lane, not a claim about the width of
a consecutive repair window or coverage of new counts.

## Two triadic chains and what they do not imply

With v in one factor A of an ordinal sum A→B, its
rooted S_v(a) gains factor Q_B(a+1) and T_v(a) gains
Q_B(a); with v in B the roles swap. This follows by
splitting a path at the cut and summing a distinguished
optional extension slot. In particular Q_(C3)(a)=a³+2a.
Each directed triangle SCC multiplies both ordinary endpoints
of another factor by exactly one factor of 3 in valuation:
Q_(C3)(1)=3, Q_(C3)(2)=12. The triangle's own S and T are
3 at every root. Thus the product family of a count-5 core and
r triangles has H=5·3^r and total endpoint content exactly
3^r in valuation (the 5-core has a unit endpoint).
But the same scalar need not prescribe that content.

The independent DP controls give the following:

| Shape | order | H | Q(2) | Minimum 3-adic valuation of weighted endpoints |
|---|---:|---:|---:|---:|
| strong 5-core | 4 | 5 | 32 | 0 |
| same core plus C3 | 7 | 15 | 384 | 1 |
| same core plus two C3 | 10 | 45 | 4608 | 2 |
| regular strong cyclic core | 5 | 15 | 120 | 1 |
| strong core rows (48,49,3,7,12,28) | 6 | 45 | 432 | 0 |
| strong core of the unit slide | 8 | 45 | 1440 | 0 |
| coupled strong core rows (238,220,216,16,225,206,8,72) | 8 | 63 | 1872 | 0 |

In particular zero endpoints mod3 do not establish triangle
ancestry (the strong 15 core), and a scalar 3² factor need not
erase any low endpoint signal (strong 45 and 63). More generally
if two SCCs have Hamilton count divisible by 3, every root sees
at least one such external factor and all roots' endpoints vanish
mod3; if exactly one such SCC exists, unit endpoints can occur
only inside it. These tests say nothing against 63. The latter
is a coupled core, not the unattainable 7 or 21 with triangles
attached. See [its two-path certificate](../chanpak_collatz_exact_addresses_20261114/READOUT_63.md)
and the [small-exclusion certificate](../chanpak_window_descent_20261115/SMALL_EXCLUSIONS.md).
The shared lesson for a constructive proof is to distinguish
a strong seam that creates the scalar factor from an ordinal
product that simply propagates an existing count.

## Replay

Run the standard and optimized Python replay. Literal subset DP
checks the paired graph counts, the 5/15/45 and 63 controls, and
both polynomial identities on every induced subtournament of the
sample graphs. The identities and continuation claims follow
from the proofs, not from extrapolating the numerical checks.
