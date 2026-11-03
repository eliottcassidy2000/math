# Simultaneous 7- and 19-adic coverage by one q2 index

Use the actual source seeds and transfer in TRANSPORT T3/T8. Let A=[[11,18],[4,7]], z_k=(25,34)A^k=(Y_k,X_k), h_k=Y_k+6. The transported strong-tournament counts are 3h_k+2 and 3h_k+4, with the pure q2 prefix automatically satisfying even triangle parity.

**Theorem.** For every a,b≥1 and every residue r mod 7^a·19^b, infinitely many nonnegative k satisfy h_k≡r modulo that product. The same is separately true for each companion. More precisely, each target has at least 18 distinct index classes modulo

    L=2736·7^(a−1)·19^(b−1).

At a=b=1 the exact number of index classes per target is 18 if r≡6 mod7, and 21 otherwise. Thus all 133 residues occur. These are congruence statements, not global integer spectrum or primality claims.

**Proof.** The period orders modulo 7 and 19 are 48 and 171, and the associated derivative matrices are

    C7=(A^48−I)/7 mod7 = [[6,6],[6,0]],
    C19=(A^171−I)/19 mod19 = [[2,18],[4,17]].

The exact finite certificate verifies:

1. For each target h mod7, exactly six of the 48 phases have that observation and nonzero (z_k C7)_1 mod7. The entire 48-cycle has six phases observing 6 mod7 and seven phases for any other observation.
2. For each target h mod19 and each k mod3, exactly three of the 171 phases have that observation and nonzero (z_k C19)_1 mod19. In fact every phase is regular.

For a regular p phase k0 with period T, binomial expansion of A^T=I+pC shows that adding j T p^(s−1) changes h by j p^s(z_k0 C)_1 modulo p^(s+1). It yields exactly one index residue mod T p^(a−1) lifting k0 to any prescribed h modulo p^a. In particular the index stays congruent to k0 modulo 3 because both base periods are divisible by 3.

The two resulting index moduli have gcd exactly 3:

    gcd(48·7^(a−1),171·19^(b−1))=3.

For any one of the six regular 7 phases, precisely three of the nine 19 phases for the desired observation have the same index mod3. Their unique lifts remain compatible, and the generalized CRT gives one class modulo their lcm L. Distinct phase pairs give distinct classes, so there are at least 6·3=18. For a=b=1, all phases are visible: the 19 distribution is three in each mod3 class, and the total 7 count is six or seven. Consequently the exact counts are 18 or 21. Periodicity and monotone growth of the positive A orbit provide infinitely many distinct integer counts for each residue.

As 3 is invertible modulo 7^a19^b, an affine change of target proves the separate companion statements. The pair is not arbitrary: for each k the two companion values differ by 2. This is a synchronized single-index result, not the invalid inference that individual modular coverage alone implies simultaneous coverage. No ternary state closure or prime-quotient inference is made.

## Verification and provenance

The seed, transport, and two base clocks are inherited from TRANSPORT.md at f6d5304849e50b29df6e368239e92fffa92e57e3. The regularity checks and compatible-phase CRT argument are deductions here. Run `python check.py` to verify the finite phase assertions, all 133 base fibres, and the larger mixed modulus 7²·19². The latter computational check illustrates, but does not replace, the all-exponent proof.
