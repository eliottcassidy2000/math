# The first two source coordinates can agree while the next one differs

Two strong order-eight sources have exactly the same (H,Q2)=(345,5640), but distinct Q3:

    A=(128,129,99,71,15,27,51,124), Q3=44505
    B=(36,97,50,7,11,216,29,95),   Q3=45225.

The replay counts by literal subset DP. For either [selective quotient](../session324_busch_inverse_and_gap_structure_20261113/COMMUTATIVE_INVERSE.md), replacing A by B with any fixed nonempty reservoir V raises the final Hamilton count by exactly 720·γ·Q_V(2), with γ=2 or 6. Thus equal H and slot budget Q2 really can conceal the next response. These graphs are not in the depth-one carry case.

There is a small-order reason their difference has just one direction. At fixed order n, two such slot polynomials have equal leading terms (the empty cycle packing), identical parity of exponents, and positive minimum exponent δ=1 if n is odd, 2 if even. Their difference vanishes at z=1,2 and hence at z=−1,−2, and its degree is at most n−2. It is therefore divisible by z^δ(z²−1)(z²−4). For n≤8 the quotient is necessarily constant (or zero), since only degrees of the same parity can appear. At order eight the displayed pair gives exactly

    Q_B(z)−Q_A(z)=2z²(z²−1)(z²−4).

At order nine the arithmetic direction space first allows more than one dimension, but this polynomial observation does not certify physical examples with multiple directions at that order. The factorization also implies the full ordered-cover dominance in this small-order regime, by the same falling-factorial argument as in the [first-carry theorem](SLOT_CARRY_INVERSE.md).
