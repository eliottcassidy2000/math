# Polynomial-precision phase placement with the physical +2 bank

Combine the [linear three-packet decoder](THREE_PACKETS.md) with the
[complementary bracketing geometry](../chanpak_qrh_residual_20261115/COMPLEMENT.md).
The linear decoder is crucial; sending the previously isolated digits
through the same proof does not give polynomial precision.

Fix an odd N, put L=sqrt(N), take the same center index j=2N mod 3, and
let B be the total pre-address guard length (including resets and the
terminal choice for mod 4). Set q=2^B. In each length-q block of a
controlled pre-address U, the other pre-address V is determined by the
fixed literal sum S, and the suffixes are complementary. Choose a
dyadic number g comparable to L^(1/3) q^(−2/3) with q·2^k dividing g.
Choose S≡−1 mod g just large enough that 2S+2≥sqrt(N/108). In particular
S−sqrt(N/432)=O(g+1). The quotient prefixes in each block have
sum ≡−1 mod 2^k because q·2^k divides S+1. The linear guard applies to
*every* pair of prefixes and realizes the target modulo 3·2^k in one
member of the physical three-center window, after the both-one closure
x=2U+1, y=2V+1. The other members still differ by exactly +2.

Set d=U−S/2≥0 and R=2S+2. The main count term is

    108 R² + (387/2) R − 1728 d² − 258 d.

The remaining terms are nonnegative and O(L). In the first full
address block, d≤2q and the linear term dominates the loss when
q²=o(L). A block at d≈C sqrt(Lg) is entirely below N for an absolute
large C, using 0≤p(U)+p(V)≤S. Uniformly through this band, the
complementary checksum has a common main part plus
O(|d|+g+L/g+1), since the low reversed coefficients cancel at
fixed sum. Two adjacent block representatives differ in address by at
most 2q, so their output jump is

    O(q sqrt(Lg) + q² + g + L/g + 1).

The selected g balances the first and fourth terms, giving
O(L^(2/3) q^(2/3)), provided q≤N^(1/16); the other terms are smaller.
Use sign bisection on the lowest and highest count separately
(using strict inequality on the lower side) to find whole windows
on opposite sides. This does not use checksum monotonicity. All
constants are independent of k when the size conditions hold.

In particular, use the linear decoder with B≤10k−7. Set
k=floor((log₂ N)/160)−C₀ with a large absolute constant C₀. For large
N this gives q≤N^(1/16), q²=o(L), and a large power-of-two margin
for q·2^k dividing g. Thus for every sufficiently large odd target:

- two physical +2,+2 windows are placed on opposite sides at absolute
  distance O(N^(3/8));
- in both, the same designated count agrees with N modulo 3·2^k,
  where 2^k is at least a positive constant times N^(1/160);
- the resulting signed integer quotients have absolute size
  O(N^(3/8−1/160)).

This improves the previous 2^(c sqrt(log N)) modulus to a power of N,
not the absolute O(N^(3/8)) bound itself. It does not yet force the
integer quotient to vanish. Unlike a near-hit sieve, the proof
constructs the very same physical windows on which the congruence and
the absolute bracket hold. A decoder can enumerate each packet's
support and store representatives; a pairwise sumset and a lookup
find a three-packet witness without searching tournaments.
