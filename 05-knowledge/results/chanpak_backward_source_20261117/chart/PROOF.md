# Recovering a one-arc source backwards from the contacts it must have

The exact target inverse of the [signed-height packet](../chanpak_signed_height_20261117/README.md) produces candidate H and P₂ for the source. This note starts by assuming a witness with this profile exists, deletes the two endpoints of its switched arc, and asks which equations recover the missing contacts. It identifies two sharply different channels (interior and fringe) even for a simple transitive remainder. Its modular construction does not solve the exact integral lift.

## 1. An exact inverse chart on a fixed remainder

Write the original arc x→y, and retain the tournament A on the other m vertices. Let X_i,Y_i be the indicators of x→i and y→i. For q=1 or 2, let K_q(i,j) be the weighted number of directed paths on an odd nonempty vertex set S in A whose endpoints are i,j, with weight P_q(A−S). Singletons are included. The [Grinberg–Stanley odd-cycle expansion](https://arxiv.org/abs/2307.05569) gives, with the difference reversed-minus-original,

    ΔP_q = 2q D_q,
    D_q = X^T K_q(1−Y) − Y^T K_q(1−X).

Indeed every affected odd-cycle packing has exactly one odd cycle using the switched arc. The rest of that cycle is an odd-vertex path in A, reached from x and returning to y after the reversal, or vice versa before it. All other pieces lie wholly in A−S. This keeps the **actual** orientations and does not approximate P₂ by H.

Put d=X−Y in {−1,0,1}^m. At d=+1 the value of Y is forced to 0, and at d=−1 it is forced to 1. On the common-contact set d=0 the bits Y are free. Then

    D_q = d^T K_q 1 + Y^T(K_q−K_q^T)d.

Thus a fixed signed disagreement pattern converts the exact repair condition into a **two-dimensional 0–1 subset-sum** over the shared contacts. The native (+10,−8) switch asks for D_1=5 and D_2=−2. In particular X,Y must be incomparable under inclusion: if X≥Y, the pair weight is d_i(1−Y_j)+Y_i d_j≥0 everywhere; if X≤Y both differences reverse sign. This excludes a port with only one direction of disagreement, however freely its shared contacts are changed.

The absolute profile is a different condition. For fixed A,d and a fixed orientation of x→y, both H and P₂ of the whole source are multilinear polynomials of degree at most four in the free Y bits. In any ordered path cover at most four arcs incident with x,y enter its indicator product, and Boolean powers reduce to their first power. If the prescribed reader is aH+bP₂, the same holds for N. The exact polynomial is obtained by evaluating the source on contact supports of size ≤4 and applying Boolean subset inversion. Therefore, on a fixed remainder and disagreement pattern, the backward problem is: find a 0–1 point solving two linear repair equations and a quartic absolute equation (or both quartics if H,P₂ are specified separately). This is an exact compiler, not an efficiency or universal existence assertion; the common-contact choices that solve the linear constraints must not be deduplicated before checking the quartic.

A necessary check before looking for A: at even source order n and canonical ΔP₂=−8, the switched arc toggles the parity of the number of odd cycles of length n−1. If C counts such cycles and E counts unordered packs consisting of two nontrivial odd cycles covering all vertices, expanding P₂ modulo 32 gives

    ΔP₂ ≡ 8 ΔC + 16 ΔE (mod 32).

There is no other contribution below 32: the number of pack pieces is even, and four pieces with a nontrivial cycle already give a factor at least 32. The common all-singleton pack cancels. Likewise ΔH modulo 4 sees only individual nontrivial odd cycles. For (+10,−8), the shorter odd cycles therefore change by an **even** number in total, while length n−1 changes by an odd number. This points the inverse toward a change in a near-spanning cycle, not merely a local triangle edit.

## 2. A transitive remainder separates low bits in two different ways

For A=TT_m in increasing order, the diagonal is K_q(i,i)=q^(m−1); below the diagonal it is 0; above, with g=j−i−1,

    K_q(i,j)=q^[i+(m−1−j)] * [(q+1)^g−(q−1)^g]/2.

There is just one path on a chosen subset. The g interior vertices are placed either inside the selected path or in the complement, retaining odd cardinality of the inside subset. At q=1 the off-diagonal coefficients are 0 at distance 1 and 2^(j−i−2) at greater distances. Thus D_1 modulo 2^L reads only the diagonal and distances ≤L+1. At q=2 the explicit exterior factor 2^[i+(m−1−j)] means D_2 modulo 2^L reads only edges crossing from the first L positions to the last L positions, provided m−1≥L (so the diagonal vanishes). The latter is a **fringe** read, not a short-distance read. The two conditions can be programmed separately.

**All-order two-adic lift, with a bounded spine.** For any L≥2 and any requested residues d_1 and *even* d_2 modulo 2^L, there are actual contact words on any sufficiently long even spine, of length O(L²), that give exactly these two residues simultaneously. D₂ even is necessary because the total source order is even, so P₂ at fixed order varies by multiples of 8.

First set both contact bits equal to 1 throughout the last L vertices except the final two, whose contact bits are both 0. Other contacts begin at (0,0). Only the final two vertices now contribute to the fringe read. For each row i in the first L with signed disagreement d_i, its coefficient v_i is K_2(i,m−2)+K_2(i,m−1). For i even the first summand has valuation i+1 (its g is odd), while the second has valuation at least i+2 (its g is even). For i odd the second summand has valuation i, while the first has valuation at least i+3. Here one uses 3^g−1 equal to 2 mod 4 for odd g and divisible by 8 for even g. Thus each pair of rows 2r,2r+1 has two coefficients equal to 2^(2r+1) times odd integers. Disagreements −1,0,+1 independently realize all four residues modulo 4 after division by 2^(2r+1): either odd coefficient gives the two odd residues with sign, and the sum or difference of both gives 2. Choose the disagreements greedily in base 4, starting from the lowest two available bits after the forced zero bit. This programs the full even d_2.

Then work well away from both fringes. For each r=0,…,L−1, put one common (1,1) contact at distance r+2 to the right of an optional switch position i_r, with all other nearby Y bits zero. Separate the switch positions by 2L+4 and keep each more than L+1 from the fringes. The optional switch changes X_i from 0 to 1 with Y fixed. Its contribution to D_1 modulo 2^L is exactly −2^r: the diagonal and all allowable forward powers would otherwise sum to 1+∑_(d=2)^(L+1)2^(d−2)=2^L, but the common (1,1) excludes the distance r+2 term. No earlier Y bit in local range contributes. The toggles are independent because D_1 is linear in X for fixed Y. They realize any correction to d_1 and do not alter the fringe read D₂. The actual arcs x→y can now be included and reversed on the same labelled vertices.

In particular the pair (5,−2) passes every dyadic test in this small structural class; a bounded no-witness search cannot become an all-order obstruction just by increasing that modulus. Conversely the high bits produced by this procedure are uncontrolled. An exact backward solver has to use the quartic baseline and still control the signed high digits of both reads.

## 3. Exact finite control and boundary

The four-symbol transitive-spine recursion in `spine.cpp` examines all contact pairs through m=12 and finds no exact (5,−2), even though the modular result is unconditional for sufficiently long m. The two-opposite-disagreement subfamily was separately tested by the exact two-dimensional meet-in-middle using `scan.cpp` on even m≤30; it finds no witness in that range. Neither finite result rules out arbitrary spine lengths or more disagreements; cores with nontransitive remainder do exist elsewhere in the packet. Transitive-spine roles are measured here against literal subset DP in `replay.py`, which also replays the simultaneous modular construction.

The next target-driven operation is **not** to pick a modulus and add arbitrary length. Keep the near-spanning-cycle parity and the two linear constraints while changing the fixed remainder A enough to put the quartic baseline at the prescribed count. The contact chart tells precisely what must be preserved when the remainder changes and what data the modular construction still lacks.
