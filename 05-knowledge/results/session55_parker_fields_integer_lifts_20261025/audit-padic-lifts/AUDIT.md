# Independent audit of prime-power lifts

Status: PASS. The all-level proofs, thresholds, and scope restrictions are correct. No correction requested. Normal and optimized independent controls agree; both producer verification modes pass.

For the 3-power threshold, the opposite-pair constraint modulo 3 forces all entries into the center's residue class. A zero center forces root divisibility by 3 and entry divisibility by 9, permitting descent. A unit center forces all entries into the class 1 modulo 3, which has exactly 3^(k-1) unit-square residues. Thus nine distinct entries require k>=3. The base cases have insufficient square residues; division at k=2 may terminate at the trivial ring, which also cannot support nine distinct entries.

For the 2-power threshold, opposite-pair sums modulo 4 force all entries 0 modulo 4 or all 1 modulo 4. The former descends by 4; the latter consists of odd squares, hence entries 1 modulo 8. Nine distinct entries then require 2^(k-3)>=9, equivalently k>=7. Descent preserves distinctness because multiplication by p^2 is an additive bijection from the lower quotient to the indicated ideal. It need not be a unital ring isomorphism and the proof does not require one: divided roots supply the square predicate.

The two affine Lo Shu templates have exactly the listed entries and line sums. Every A entry is 1 modulo 3 and every B entry is 1 modulo 8. They remain distinct from moduli 27 and 128 respectively; strict no-wrap starts at 81 and 128. Entries such as 22 and 57 are not ordinary integer squares. There is no integer-square witness here.

The 2-adic induction shifts the old root by 2^(j-1), not 2^j. Its square changes by 2^j times an odd factor modulo 2^(j+1), with the correction-square term zero for j>=3. It therefore establishes levelwise existence but does not establish that an arbitrary chosen root at level j lifts while preserving that root's entire residue. This caveat is explicitly present in the proof and is essential. The odd-prime unit-root induction does preserve the previous residue and is unique.

For Theorem 2 each fixed root modulo p determines exactly one root for each lifted unit entry. The magic entry lattice has three free pivots over every quotient; hence p^(3(k-1)) root vectors and p^3 children per vertex. This works even without distinctness, and distinctness, if present modulo p, persists. The Jacobian rank explanation is valid since multiplying all columns by units preserves rank. The result is a local dimension law, not a global integral existence theorem or an all-prime local-global counterexample.

Characteristic-two fields and 2-power quotient rings are different: in a characteristic-two field opposite entries must agree, immediately precluding nine distinct entries. The quotient rings have characteristic 2^k, not 2. No field-cardinality substitution is licensed.

## Reproduction

Run python3 check.py and python3 -O check.py; compare output to normal.json and optimized.json. The independent checker enumerates every square-entry pivot triple at moduli 2,3,4,8,9,16,27,32,64,128, checks both template families at five levels, and directly enumerates all 3^9 lifts of the fixed root vector (1,...,1) modulo 3, finding precisely 27 magic children modulo 9. This final test does not use the pivot parametrization to generate candidates. All checks raise exceptions rather than using assert. Output certificates pin the three producer files by SHA256. These finite controls supplement, not replace, the audited all-level proofs.
