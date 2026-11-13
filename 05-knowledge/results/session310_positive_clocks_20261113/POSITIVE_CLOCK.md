# The odd-cycle expansion reads hidden bits of the determinant clock

Let Q_T(t) count ordered path covers in t labeled slots, allowing empty slots, with H(T)=Q_T(1). For a packing F of vertex-disjoint directed odd cycles of length at least three, put k(F)=|F| and u(F)= the number of vertices not covered by its cycles. The identity

    Q_T(t) = sum_F 2^k(F) t^(u(F)+k(F))

is the all-ones specialization of [Grinberg–Stanley, Theorem 1.39](https://arxiv.org/abs/2307.05569). It also follows directly from the repository’s [odd-cycle formula](../../../01-canon/theorems/THM-002-ocf.md): color each vertex by its path slot, apply the odd-cycle formula within each color class (an empty slot has weight 1), and reverse the order of summation. Every packed cycle is monochromatic and every uncovered vertex is free. Each packed cycle contributes two and admits t choices of slot, as does each uncovered vertex. This gives an elementary proof from the inherited one-variable count. The expansion retains the number of uncovered vertices as well as the number of cycles. Two packings with the same weight at t=1 can therefore supply different clock exponents.

## Lowest exponents and the two colors

At t=2 a packing contributes exactly 2^(u+2k). If n≥2 this is divisible by 4 (the empty packing gives 2^n). For even n≥4 it is divisible by 8: for k=1 the odd cycle leaves at least one uncovered vertex, and for k≥2 the exponent is already at least four. More precisely:

- For odd n≥3, Q_T(2)/4 modulo 2 equals the parity of the number of directed Hamilton cycles, since u+2k=2 forces one spanning odd cycle.
- For even n≥4, Q_T(2)/8 modulo 2 equals the parity of the number of directed cycles of length n−1, since u+2k=3 forces one cycle and one uncovered vertex.

For t=3, the exponent u+k is congruent to n mod 2, because every cycle has odd length. It is positive for nonempty T. Thus Q_T(3) is divisible by 9 for even n. More precisely Q_T(3)/3 modulo 9 for odd n≥3 is twice the Hamilton-cycle count modulo 9: after exponent 1 the next possible exponent is 3. For even n≥4, Q_T(3)/9 modulo 9 is twice the number of (n−1)-cycles plus four times the number of unordered spanning pairs of disjoint odd cycles modulo 9. The terms of exponent two are exactly those two types. This calculation explains why one versus two visits is not just a port distinction: it changes which near-spanning cycles can be observed in a small-prime residue.

The all-order repository [half-grid polynomial](../../../01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md) already gives order parity; the positive packing expansion additionally identifies the precise combinatorial coefficients that control these residues.

## A clock-dependent determinant modulus

Use the inherited [polarized repair](../target_collatz_transitions_20261113/POLARIZED.md) with h=H(M), j=H(N), p=Q_M(2), q=Q_N(2), determinant e=j(6p+h)−6hq, and output Y=36pq+6hq+18pj−hj. Exact cancellation gives

    Y+e = 12 p(3q+2j).

If |N|≥2, q is divisible by 4 and j is odd, so 3q+2j has valuation exactly one at 2. If |N|=1, q=2,j=1 and it has valuation exactly three. Consequently

    v2(Y+e)=3+v2(p)    if |N|≥2,
    v2(Y+e)=5+v2(p)    if |N|=1.

If |M|≥2 and e=±1 then Y≡−e mod 96; when |M| is even and at least four the modulus strengthens to 192. The singleton boundaries cannot evade the modulus. For M=TT1, e=13j−6q; since q≥(|N|+1)j by cutting a Hamilton path into two slots, a nonsingleton N forces e≤−5j. For N=TT1 the equation is e=6p−11h, and a nonsingleton M forces e≥7h. If both modules are singletons, e=1 and Y=191. Hence **every** determinant ±1 output is ±1 modulo 96, including the boundary. In particular prime residue classes outside these two are necessarily missing from this specific last-step rule; by Dirichlet this is an infinite obstruction, not an overall spectral gap.

The known negative-det repair with p=3424 has v2(p)=5, and 435779329−1 indeed has valuation 8. This is a precise least significant bit obstruction, not just a finite filter in a response bank. With |N|≥2, a determinant output for which v2(Y+e)=5 must have |M|=2 or odd order with an odd number of Hamilton cycles. In the inherited clock Π, the stronger intercept modulus from the same expansion is 2^4·3^2·5=720 when |N|≥2, and 2160 when |M| is also even; it follows from 60Q_M(3)Q_N(2).

The proved constraints suggest that searches should pair output valuation with the physical cover data of the corresponding module, rather than fill gaps in a residue list and infer a continuation. They do not preclude a different chart from reaching omitted residue classes or prove a supply of modules realizing each permitted combination.
