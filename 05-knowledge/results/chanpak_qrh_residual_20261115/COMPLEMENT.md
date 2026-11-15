# Complementary low digits give exact triplets at distance O(N^(1/3))

For every sufficiently large real target N (the code takes integers), the previous independent-address window has one whole triplet on either side of N within O(N^(1/3)). All outputs are literal Hamilton-path counts; both triplets still have gaps exactly +2. This is not an exact selector for unrestricted targets.

The prior [six-vertex center](https://github.com/eliottcassidy2000/math/blob/2aa5cd874debced736c3461a65bd544f41d17f95/05-knowledge/results/chanpak_target_address_20261115/INVERSE.md) gives, for each positive integer pair x,y,

    H_j(x,y) = 432xy +129x +258y +45[p(x)+p(y)]
               +15(x mod 2)+30(y mod 2) +135+2j,  j=0,1,2.

This is the exact graph construction, not an interpolation of the three counts. Here p(0)=0, p(2u)=p(u)+u+(u mod 2), p(2u+1)=p(u), and 0≤p(u)≤u. The latter follows by induction using p(u)≤u−1 for odd u.

## Replace a common low-bit phase by complementary low bits

For g=2^s, b in [0,g), the inherited low-bit identity is

    p(gu+b)=p(u)+A_b u+p(b)+z_b(u mod 2),

where A_b is the s-bit reversal of g−1−b and z_b indicates that the leading digit of b is zero (A=z=0 when s=0). It follows by appending the digits one at a time; every appended zero charges the current prefix and its parity. In particular

    A_b + A_(g−1−b) = g−1.

Instead of requiring x and y to share the same low block, require their low blocks to be complementary. This is automatic under the single condition S=x+y ≡ −1 mod g; it puts **no restriction on the step of x**. Write

    x=gu+b, y=gv+(g−1−b), u+v=(S−g+1)/g.

With d=x−S/2≥0 the product and linear terms give

    P(x,y)=108 S²+(387/2) S−432 d²−129 d.

Moreover

    p(x)+p(y) = (g−1)(u+v)/2 + R,
    R=(A_b−A_(g−1−b))(u−v)/2
      +p(u)+p(v)+p(b)+p(g−1−b)
      +z_b parity(u)+z_(g−1−b) parity(v).

Since |A_b−A_(g−1−b)|≤g−1 and g|u−v|≤2d+g−1, we have

    |R| = O(d+g+S/g+1).

The large linear checksum term depends only on S and g and is therefore exactly fixed across the unit grid, even as the low block changes. This is the new compensation.

## Place the grid from the target, then bisect signs (not monotonicity)

Write L=sqrt(N/432). Choose g as the largest power of two at most L^(1/3), and S as the least integer congruent to −1 modulo g with S≥2L. Then S=2L+O(g), and the nonnegative imbalances d=x−S/2 start at 0 or 1/2 and increase in **unit** steps.

At the starting imbalance, H_0>N: 108S²≥N, the positive linear term dominates the bounded half-unit loss, and all checksum terms are nonnegative. At an imbalance d=O(sqrt(Lg)) the entire triplet must lie below N: indeed P(ceil(S/2),floor(S/2))−N=O(Lg), the whole additive checksum is at most 45S+45, and the quadratic loss 432d² can dominate both by a fixed factor. Such d=o(L), so both addresses remain positive. The code selects the first grid point where the universal upper bound P+45S+45+139≤N by an ordinary monotone binary search on P.

Throughout this interval, consecutive values of each H_j differ in absolute value by

    O(sqrt(Lg)+g+L/g) = O(L^(2/3)) = O(N^(1/3)).

For the smooth part this follows by differencing its quadratic at unit steps; for the residual the preceding uniform bound is enough. We do not assume H_j is monotone. If two indices have opposite signs relative to a threshold, bisecting by sign still finds two adjacent indices with opposite signs. Bisect on H_0 vs N and keep the strictly above endpoint, then bisect on H_2 vs N and keep the at-or-below endpoint. At each endpoint all three counts are on the same side and differ from N by at most the successive-value bound plus 4. Both windows retain all residues modulo 3.

The method uses O(log N) bisection evaluations, each O(log N) arithmetic operations on O(log N)-bit integers, and each final graph has O(log N) vertices. It reduces the prior uniform O(N^(3/8)) error while keeping the exact three-choice center. Reaching equality for all N would still require a growing repair or a different argument: neither complementing the low digits nor sign bisection turns an O(N^(1/3)) bound into an exact interval.

Replay: `python complement.py` checks the identity, reconstructs brackets across 40,000 consecutive targets and additional targets through 512 bits, and prints exact examples. These checks supplement the uniform mathematical proof.
