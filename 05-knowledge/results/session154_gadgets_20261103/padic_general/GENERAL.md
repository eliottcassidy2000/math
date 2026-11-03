# Prime-power coverage and a lift obstruction in the q2 clock

Use the supplied literal-seed and even-triangle transport of TRANSPORT, section T3: z_k=(25,34)A^k, A=[[11,18],[4,7]], h_k=Y_k+6. The actual strong-tournament companions have counts 3h_k+2 and 3h_k+4. The following are congruence statements only, not spectrum completeness or primality assertions.

## An exact finite test and a sufficient lift criterion

Let p be an odd prime other than 5. Put T=ord_{GL_2(F_p)}(A), so T exists because det(A)=5, and define the integral matrix C=(A^T−I)/p. Write z_k=(Y_k,X_k), d_k=(z_k C)_1 mod p. Then:

1. For every m≥1, A^(T p^(m−1))≡I mod p^m, so enumerating 0≤k<T p^(m−1) is an exact finite decision procedure for the residues reached by h_k mod p^m. This can be expensive and is not a uniform classification over primes.
2. Mod p², the complete reached set is the union over 0≤k<T of {h_k+p j d_k : 0≤j<p} mod p². This is an exact finite test, including singular branches.
3. If for each a in F_p there is k<T such that h_k≡a and d_k≠0 mod p, then every residue mod p^m is attained infinitely often for every m. In fact each regular base phase gives exactly one k mod T p^(m−1) for every compatible target residue. This condition is sufficient, not claimed necessary.

Proof. With M=A^T=I+pC, binomial expansion gives M^(p^(s−1))≡I+p^s C mod p^(s+1) for odd p. Hence adding j T p^(s−1) changes h by j p^s d_k modulo p^(s+1); d_k is unchanged mod p on the fixed phase k mod T. A nonzero d supplies exactly the desired next digit. The m=2 union is immediate from (I+pC)^j≡I+jpC mod p². Unboundedness of positive A gives infinitely many distinct counts. For p≠3, multiplication by 3 is invertible modulo p^m, so the same criterion gives both companion families all residues. For p=3, the companions are fixed at 2 and 1 respectively mod 3, so neither is surjective even if h were.

## A positive case at 19

At p=19, T=171 (as in TRANSPORT T3) and

    A^171 mod 361 = [[39,342],[76,324]],
    C mod 19 = [[2,18],[4,17]].

The following table gives, for each a=0,...,18, one (k,d_k) with h_k≡a and d_k≠0:

    (20,18) (5,12) (2,18) (19,2) (18,15) (15,13) (4,14)
    (3,10) (9,15) (39,10) (16,3) (6,13) (0,15) (14,12)
    (58,14) (22,14) (34,3) (12,10) (1,2).

The finite criterion proves that h_k and each companion separately hit every residue mod 19^m for every m, infinitely often. Thus both attain every exact 19-adic valuation. The supplied program verifies order, entries and the witnesses using modular integer arithmetic.

## p=13: mod-prime surjectivity does not lift

Here ord(A mod 13)=56; h_k visits all 13 residues. Nevertheless

    A^56 mod 169 = 144 I,   C mod 13 = 11 I.

Exactly k=13,27,41,55 mod 56 have h_k≡6 mod 13. Their respective values mod 169 are 45,32,136,149. For those phases Y_k≡0, hence d_k=11Y_k=0 mod 13; adding any multiple of 56 leaves their Y, and hence h, unchanged modulo 169. Therefore those four are the only values among the 13 possible lifts of residue 6. In particular h≡6 mod 169 is impossible; the corresponding missing companion residues are 20 and 22. All claims follow by the displayed finite cycle. The example is a failure of this specific pure-q2 seed family, not of the full mixed prefix language or of all strong tournaments.

## The excluded determinant prime and a 3-divisibility warning

At p=5, A is singular, so GL order cannot be used. Direct iteration mod 5 gives h_0≡1, and thereafter h_k runs periodically through 2,4,0,3. Thus h≡1 mod 5 occurs only at k=0, where h_0=31≡6 mod 25. In particular h≡1 mod 25 is impossible and the whole mod-5 image arises partly from a one-time transient. For p=3, each companion has the fixed residue already stated. These examples are explicit barriers to extrapolating a mod-p clock blindly.

## Scope and provenance

The seed, transfer, positive-offset condition, p=7 field computation, and p=19 clock order are inherited from [TRANSPORT.md at f6d530](https://github.com/eliottcassidy2000/math/blob/f6d5304849e50b29df6e368239e92fffa92e57e3/05-knowledge/results/golden_inverse_defects_20261103/TRANSPORT.md). The derivative criterion, p=13 lift obstruction, p=19 coverage certificates and p=5 transient are deductions here, not claims of novelty in p-adic dynamics. Do not infer pairwise independent choices: at fixed k the two companion counts differ by 2. Do not infer closure under a further ternary replacement.
