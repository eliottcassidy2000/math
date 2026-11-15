# Let the target choose the complementary modulus

The [complementary phase construction](COMPLEMENT.md) gives a uniform
error bound O(N^(1/3)). Keeping the exact upward alignment cost, rather
than bounding it by the modulus, makes the bound on certain targets
smaller.

Write L=sqrt(N/432) and s0=ceil(2L). For each power of two
1<=g<=sqrt(s0), put

    delta_g=(-1-s0) mod g, S_g=s0+delta_g.

Use the same signs and the same exact three-center values as in the
complementary construction; now record the cost of its alignment for
this target. For d=x-S_g/2≥0, the deterministic upper envelope is

    H_2 <= 108 S_g² + (387/2)S_g - 432 d² - 129 d
                       + 45 S_g + 45 + 139.

At d=0 the excess of this upper envelope over N is
O(L(delta_g+1)), because S_g−2L≤delta_g+1, S_g=O(L), and all remaining
terms are O(L). Thus the first grid point certified below N has
d=O(sqrt(L(delta_g+1))) and still lies with positive addresses for
all sufficiently large N, uniformly over permitted g. At the near
equal split, the whole window is above N as before. The exact
complementary decomposition controls the checksum fluctuation along
this shorter interval by O(d+g+L/g+1). Binary sign crossing therefore
gives actual triples on both sides, each within

    O(sqrt(L(delta_g+1)) + g + L/g).

One may minimize the displayed bound over the dyadic choices using
integer arithmetic and then perform the two sign bisections. There are
O(log N) candidate powers, so the same O(log² N) arithmetic-operation
bound applies. Taking g of order L^(1/3) recovers the uniform result.

A strict target-specific power saving is guaranteed when the upward
alignment is O(1) for g comparable to sqrt(L): then the displayed
bound is O(N^(1/4)). In particular, take S=2^(2r)−1 and g=2^r. For all
integer targets

    108(S−1)² < N <= 108S²

we have s0=S and delta_g=0, and g is a permitted choice up to a
constant factor. Strictly under g²<=s0 the neighboring power g/2
is permitted and equally aligned. These are infinitely many entire
intervals of length Θ(sqrt(N)) where the exact three-count windows
are closer than O(N^(1/3)). There is no corresponding claim of
density one, and the improved residual is not a zero residual.

Replay: `python adapt.py` evaluates 771 large aligned targets, with
exact sign and physical address checks, plus representative examples.
