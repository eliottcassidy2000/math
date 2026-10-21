# Independent audit: all-prime-power signed Paley markers

Verdict: APPROVED for the exact producer files hashed in receipt.json. No mathematical correction required. Review performed 2026-10-21, with mathematical repository read-only and no Git commands.

## Approved scope

For every prime power q=p^e congruent to 3 modulo 4, the full signed Paley group in field coordinates is AΓL(1,q). Its regular subsets have minimum size three except q=3,7, where no regular subset exists. For e>1 (necessarily odd), the full-degree normalized parameter criterion, both Frobenius–Möbius obstructions when 3 divides e, existence bound, and exact Möbius formula for free signed three-set orbits are valid. The e=1 case is separately proved in the prime-field predecessor; the extension-field formula is not asserted at e=1.

Status separation: full group identification uses a CITED classical theorem; elementary consequences are PROVED; explicit replay and independently implemented controls are FINITE-EXACT. This audit makes no novelty, arbitrary-tournament, Hamiltonian-gap, or LRC claim.

## Proof obligations checked

1. Primary source independently retrieved: Bruen–Levinger, Canadian Journal of Mathematics 25 (1973), Theorem 1 on page 1060, DOI https://doi.org/10.4153/CJM-1973-113-4. Publisher PDF text explicitly states a proper divisor d of q−1, difference ratios in its order-d subgroup, and exactly the claimed affine-Frobenius maps with (q−1)|d(p^j−1). For d=(q−1)/2 and odd p the condition is automatic. No excluded small field is present. Ordinary and anti-automorphisms are correctly distinguished by square versus nonsquare slope.
2. Affine kernel: checking all five nonidentity induced permutations gives t∈{−1,2,1/2} or t²−t+1=0. All have degree at most two. Characteristic three merely merges exceptional conditions; no division by three is used.
3. Two different injections are needed and supplied. The projection to the cyclic Galois group is injective because the affine kernel is trivial. The action on the three-set is faithful because fixing 0,1 forces b=0,a=1, and fixing a full-degree t kills the Frobenius power. Consequently the stabilizer order divides both e and six, and is one or three. No order-nine or hidden pointwise stabilizer is omitted.
4. If order three occurs, one nonidentity element projects to p^(e/3)-Frobenius. The two cycle orientations yield precisely the two displayed equations; inversion covers the other nonidentity Galois power.
5. The union-bound existence proof is strict already at r=3: r³−7r/2−2>0. Summing p^d for all integer d≤e/3 is an overestimate of proper-subfield sizes, not an erroneous assertion about subfield degrees.
6. For e=3u and d|e, d/gcd(u,d) divides three. In the order-one case the fixed quadratic has D roots in each odd-degree subfield, including the single repeated root at p=3. In the order-three case both polynomial equations are separable, have degree p^g+1, and every algebraic-closure root lies in F_(p^(3g)); excluded poles 0,1 cannot occur. The inverse-Frobenius conversion is valid because rho has coefficients in the prime field. Möbius inversion therefore counts full-degree roots correctly. The two obstruction root sets intersect only in degree at most two, so subtraction of 2L_e is exact.
7. Affine normalization preserves regularity under the full group. A normalized parameter in a proper subfield has a nonidentity pointwise Frobenius fixer. Each labelled unordered triple has six ordered normalizations, and division by e*q*(q−1) applies exactly to free orbits. Thus N_e/(6e), not the affine-only quotient, is correct.
8. Small subsets fail by reflections; complementation plus the explicit five parameters proves q=3,7 have no regular subset. The prime predecessor's {0,1,3} certificate is valid for the remaining prime orders and is not wrongly transferred to extensions.

## Exact controls

Producer verify.py was rerun normally and under -O. Both outputs byte-match each other and the producer's normal.json and optimized.json. Results:

- q=27: full-degree stabilizers {1:18,3:6}; one free signed triple orbit.
- q=243: full-degree stabilizers {1:240}; eight free signed triple orbits.
- q=343: full-degree stabilizers {1:324,3:12}; eighteen free signed triple orbits.

Independent arithmetic in independent_q27.py imports none of the producer implementation. It enumerates all 2,106 distinct semilinear maps and all 2,925 three-subsets. Exact stabilizer histogram: {1:2106,3:702,18:117}. In F3[z]/(z³+2z²+1), the triple {0,1,2+z} has trivial affine stabilizer but order-three semilinear stabilizer. With base-three coefficient encoding its nonidentity maps (a,b,j) are (7,5,1) and (4,1,2), inducing images (5,0,1) and (1,5,0). This directly refutes transferring affine regularity alone to the full signed group. Normal and -O independent outputs match.

A second independent algorithm, independent_rootcounts.py, counts roots via polynomial gcd of (1−X)X^(p^u)−1 with X^(p^d)−X, without constructing extension fields. Every M_d matches for (p,e)=(3,9),(7,9),(3,15), covering repeated factors in e, characteristic three, nontrivial proper subfields, and inverse-Frobenius restriction (e=15,d=3).

Reproduce from workspace:

    python3 session24-paley-powers/verify.py
    python3 -O session24-paley-powers/verify.py
    python3 session24-paley-powers-review/independent_q27.py
    python3 -O session24-paley-powers-review/independent_q27.py
    python3 session24-paley-powers-review/independent_rootcounts.py

## Boundaries and approval

No missing assumption was found within the stated q≡3 mod4 and e>1 extension scope. General odd or even e outside that scope is not approved. The full group identification is sourced, not independently discovered by the finite computations. The existing exception classification is inherited, not advertised as new. The exact approved allowlist and source/replay hashes are in this review's receipt.json. Any subsequent source edit requires rechecking its affected claims.
