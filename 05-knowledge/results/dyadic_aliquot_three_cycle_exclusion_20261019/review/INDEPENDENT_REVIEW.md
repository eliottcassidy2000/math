# Independent review: dyadic-times-prime aliquot three-cycle exclusion

19 October 2026. **Verdict: PASS for the stated all-exponent arithmetic theorem and its scoped aliquot corollary.**

The argument proves that for three pairwise distinct powers of two at least 2, the unique rational solution of the cyclic affine system has a nonintegral variable at the largest core. The odd-common-denominator propagation lemma then excludes integrality of every other variable individually. Together with the separately checked repeated-core classification, this proves the proposed hypotheses A3 and D3 in their original restricted meanings.

No finite search is being substituted for the universal proof. This review makes no historical-priority claim and no assertion about general aliquot three-cycles or trajectory termination.

## Reviewed sources

The proof reviewed is `DYADIC_THREE_CYCLE_EXCLUSION.md`, SHA256:

`3af6d6b9cd0e501f102d1cf009731f9f791b5a2854191d8d9915d2c7a333f129`

The inherited repeated-core classification is Section 3 of `aliquot_dyadic_frontier_20261014e.md`, SHA256:

`07d9d18c360f639922b8c9ed9ecfb14887d88e29ff3ab729f987a1fad4697348`

## 1. Elimination, uniqueness, and orientation coverage

For an ordered core word `(x,y,z)`, substitution through the three edges independently gives the numerator of the variable at `z` as

    B(x,y,z) = (y−1)(x−1)(2z−1)
               + x(y−1)(2x−1) + xy(2y−1).

The coefficient matrix has determinant

    D = xyz − (x−1)(y−1)(z−1) > 0.

Thus there is exactly one rational cyclic solution, with largest-core variable `B(x,y,z)/D`. Expanding `B(c,d,e)` gives the stated forward numerator. Swapping `c,d` gives exactly

    B(d,c,e) − B(c,d,e) = (c−d)(2(c+d)−1).

Sorting the three distinct cores loses no directed case: the six permutations form exactly two classes under cyclic rotation, represented by `(c,d,e)` and `(d,c,e)`. Both formulas retain the variable at `e` in the last position.

Because all cores are even, `D` is odd. Every variable has a denominator dividing `D`. If any one is integral, its outgoing edge makes a power of two times the next variable integral; an odd denominator cannot be canceled this way. Integrality therefore propagates around the whole cycle. This correctly connects the largest-core obstruction to the predecessor's original first-variable nondivisibility formulation of A3.

## 2. The strict size trap is valid for both orientations

The four slope and boundary identities in Section 2 were independently expanded and checked. On `c≥2`, `d≥2c`, `e≥2d`, all required signs are strict. In particular, the lower slope is positive even at the smallest case `c=2,d=4`, and neither endpoint of the trap can survive.

Since the reverse numerator is smaller than the forward numerator, proving the lower bound for the reverse orientation and the upper bound for the forward orientation proves

    c−1 < r < 3c−1

simultaneously for both. The divisibilities `c | d` and `c | e` yield `D≡1 (mod c)` and both numerators `≡−1 (mod c)`. Hence an integral `r` satisfies `r≡−1 (mod c)`. The only permitted value in the strict interval is exactly `r=2c−1`.

This step does not use primality or any assumption that the cyclic variables are initially positive. The exact rational system and strict bounds apply to every putative integer solution.

## 3. Both valuation contradictions survive the edge cases

The forward substitution gives exactly

    e(2c²−c+d−1) = c(2d²−1).

Both parenthesized factors are odd. Therefore `v₂(e)=v₂(c)`, which contradicts `e>c` because both are powers of two.

The reverse substitution gives exactly

    e(2c²−c+d−1) = 2c(c−1) + 2d²(c−1) + d.

For `c=2^a,d=2^b,b>a+1`, the three valuations on the right are `a+1,2b+1,b`. The minimum is unique, so `v₂(e)=a+1<b`, contradicting the required `e>d`.

For the remaining case `d=2c`, the equation is exactly

    e(2c²+c−1) = 2c²(4c−3).

The odd factors force `e=2c²`; cancellation would require `2c²−3c+2=0`. Its discriminant is `−7` and its leading coefficient is positive. Thus this case, including `c=2,d=4`, is impossible. No boundary exception is omitted.

## 4. The inherited repeated-core classification is sufficient and correct

Every three-node word with a repeated core is a cyclic rotation of `(C,C,T)`. Independent elimination confirms the predecessor's expressions for `k=(p+1)/C` and for the variable at `T`, including all interval-difference identities.

- If `T<C`, then `0<k<2`, incompatible with the positive even integer forced by positive odd cyclic variables.
- If `T=C`, all variables are `2C−1`; the source nodes coincide, so this is a fixed point rather than a primitive three-cycle.
- If `T>C`, the strict interval and congruence force the third variable to be `C−1`. The resulting equation is `(C−1)T=C²(3C−2)`. Reduction modulo `C−1` forces `C−1 | 1`, hence `C=2,T=16`.

The exceptional formal variables are exactly `(23,13,1)`, up to cyclic rotation. They produce the formal sources `(46,26,16)`, but `s(16)=15`, not `46`. The prime guard fails at `1`, so this is not a native aliquot cycle.

For actual nodes `2^a p` with `a≥1` and `p` an odd prime, multiplicativity gives precisely the displayed affine system. The distinct-core theorem and the repeated-core classification cover all patterns. This establishes the stated native corollary with no exponent or height cutoff.

## 5. Independent executable controls

The companion `independent_audit.py` imports no author checker. It uses exact rational Gaussian elimination as an independent numerical oracle and SymPy only for exact symbolic identities. It raises explicit exceptions on failure, so Python optimization cannot disable its checks.

The normal and optimized runs produce byte-identical receipts:

- 25 symbolic identities, including the matrix determinant, all three rational formulas, both orientations, the strict-bound identities, the valuation equations, and the repeated-core elimination identities.
- 8,000 complete ordered exponent triples in `1..20`, checked by rational Gaussian elimination. The sole integral controls are the 20 equal-core fixed templates and the three rotations of the formal ghost.
- 83,328 oriented distinct-core triples with sorted exponents in `1..64`, checking the strict trap, retained congruence, and largest-variable nonintegrality.
- 729 Gaussian-elimination triples drawn from exponents `{1,2,3,17,128,129,512,1024,4096}`, including large gaps and values beyond the predecessor's finite frontier.
- 8,128 valuation-pair controls with `1≤a<b≤128`, explicitly separating the adjacent-exponent case.
- Two hostile controls retaining the native-prime and exact-dyadic guards.

The nondyadic hostile control is worth preserving: the same formal affine system with even cores `(2,44,158)` has integral variables `(393,9,3)`. Thus an extension from powers of two to arbitrary even cores would be false, even though the displayed size order and divisibility by the smallest core still hold. The final valuation step uses the exact powers-of-two hypothesis essentially.

### Reproduction

Python 3 and SymPy are required. From the directory containing the checker:

    python -B independent_audit.py > independent_audit.json
    python -B -O independent_audit.py > independent_audit_optimized.json
    cmp independent_audit.json independent_audit_optimized.json

Checker SHA256:

`5bef0a1465e3fac392704980aa9468fefaddedd0f08591213906061059e1f74a`

Both receipt files have SHA256:

`a588f8ac9dd1f02beff61ef132c91fa2ed1c34b8e05e3151425c267ff4861d16`

These calculations are independently reproducible regression controls. The unrestricted quantifiers are supplied by the algebraic proof audited above.
