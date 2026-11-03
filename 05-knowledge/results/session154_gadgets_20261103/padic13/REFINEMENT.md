# Exact higher 13-adic image and 19-adic fibres

The notation and literal strong-tournament construction are as in TRANSPORT T3 and GENERAL.md. In particular z_k=(25,34)A^k=(Y_k,X_k), A=[[11,18],[4,7]], h_k=Y_k+6; the two companion counts are 3h_k+2 and 3h_k+4. All claims concern congruence images, not existence at a specified integer or primality.

## Complete image at every 13-power

**Theorem.** For m≥2, h_k mod 13^m hits exactly those residues whose reduction mod 169 is not one of

    6, 19, 58, 71, 84, 97, 110, 123, 162.

Thus the image has 160·13^(m−2) elements. For the first companion, the nine forbidden residues mod169 are the images 3a+2; for the second, they are 3a+4. Each lifted residue that is attained is attained for infinitely many k. At m=1 the image is all 13 residues.

**Proof.** Integer computation gives T=ord(A mod13)=56 and

    M=A^56 ≡ 144 I + 169 D (mod 2197),
    D ≡ [[7,12],[7,0]] (mod 13).

The 56 base phases visit every h mod13. All phases with Y nonzero have derivative (z(M−I)/13)_1≡11Y nonzero mod13 and so supply every higher lift of their particular h mod13 (the general first-order lift lemma). The four phases with Y=0 mod13 are

    k:                 13  27  41  55
    h_k mod 169:       45  32 136 149
    (h_(k+56)−h_k)/169 mod13: 8  1  5 12.

Here the displayed divisions are integer divisions, since each difference is divisible by 169. The only candidates mod169 above h=6 mod13 are therefore these four, and the nine other lifts are precisely those listed in the theorem.

It remains to show that each of those four does lift without a further obstruction. Let z be any state in one of these phases. Its Y is a multiple of 13 and its Y mod169 and z mod13 are unchanged by any M step. For s≥2 and n=13^(s−2), the scalar decomposition yields

    (z M^n−z)_1 / 13^s ≡ 11(Y/13)+(zD)_1 (mod 13).

Indeed, 144^n−1 has 13-adic valuation s−1 and quotient 11 mod13, while the first-order D term contributes n·169 D; all higher D terms vanish modulo 13^(s+1). The displayed unit therefore remains the same along the phase. The four values in the table are units. Varying the number of such increments j=0,...,12 makes the next digit arbitrary, and induction fills every residue above the phase's mod169 value. Conversely, no sequence can leave its established mod169 image. These are actual nonnegative indices: start at the listed phase and add nonnegative increments. Periodicity modulo each finite modulus, together with growth of h, gives infinitely many k and infinitely many distinct counts.

Multiplication by 3 is a unit mod13^m, giving the two companion image statements. The two companion choices at a given k remain correlated: they differ by 2. No claim is made about mixed-prefix words, other seeds, or all strong tournaments.

## Exact 19-adic fibre size

For 19, T=171 and C=(A^171−I)/19 mod19 = [[2,18],[4,17]]. The finite certificate shows that every one of the 171 base phases is regular, and each of the 19 residues of h mod19 appears exactly 9 times. The first-order lift lemma then gives **exactly 9 distinct index residues modulo 171·19^(m−1) for every h residue mod 19^m**, for every m≥1. The same exact nine-index count holds separately for each companion because 3 is invertible. This is a statement about index residues modulo a valid period, not a claim that distinct indices produce distinct integer values modulo the same modulus.

## Provenance and checks

The p=7 transport, p=19 order, seed and positive-offset theorem are inherited from TRANSPORT.md at f6d5304849e50b29df6e368239e92fffa92e57e3. The mod169 obstruction was independently audited in this session. The second-order lift and exact 19-adic fibre refinement are deductions here, not a claim of novelty for p-adic lifting. Run `python check.py` for integer checks of cycle, singular derivatives, exact image sizes through 13^4 and the 19-base-phase distribution.
