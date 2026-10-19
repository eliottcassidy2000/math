# Independent review: perfect-core twin edges

19 October 2026. Reviewed `../PERFECT_CORE_TWIN_EDGES.md`, its verifier, and the relevant inherited guarded identity and dyadic-classification statements. The review below refers to the pre-clarification source; both wording recommendations were incorporated in the published note, as recorded in `../SOURCE_PROVENANCE.json`.

## Verdict

The three central mathematical assertions are correct on the inherited positive-integer core domain. There is no prime-infinitude or historical-novelty claim. Two small wording clarifications are recommended before publication:

1. State **positive integer C** explicitly in the opening theorem. The inherited divisor-sum and aliquot functions have positive-integer domains.
2. In section 4, name the globally classified dyadic dynamics as the **auxiliary cofactor-map dynamics**. The inherited global classification is for F, which agrees with aliquot dynamics only while its prime guard holds. The current phrase “its exact changing-core dynamics were globally classified” is unnecessarily ambiguous.

## Mathematical audit

- For prime p not dividing C, multiplicativity gives s(Cp)=(sigma(C)-C)p+sigma(C). Equal shifts h at two distinct legal primes imply (sigma(C)-2C)(p1-p2)=0; hence sigma(C)=2C and h=2. Perfectness conversely gives the stated translation. In fact, if the shift is already specified to be 2, one legal prime suffices.
- For C=6 and p>=5 prime, the guard is automatic. If the actual target is 6q with q>=5 prime, cancellation gives q=p+2. Thus both directions and injectivity are exact. The sole excluded ordinary twin pair with lower prime below 5 is (3,5), and its attempted source 18 does not obey the guarded formula: s(18)=21.
- For other fixed perfect cores, retaining both p not dividing C and q not dividing C removes only finitely many twin pairs. The target guard is important if the target is to remain in the same admissible class: 28*5 -> 28*7 is a genuine edge, but its target marking violates coprimality. Finite removal preserves the equivalence of infinitude assertions.
- The source-cutoff identity counts sources, with the exact condition 6p<=X. It is a counting identity, not an asymptotic estimate.
- The fixed-6 exit bound is correct and sharp. If p=1 modulo 3, the first target cofactor is composite. If p=2 modulo 3, either it exits on the first step, or the legal second step has cofactor p+4 divisible by 3 and greater than 3. Example: 30 -> 42 -> 54. Each run contains at most one internal prime-to-prime edge, although the union over starts can have arbitrarily many such edges if twin primes are infinite.
- The two-step argument never applies the affine formula to a composite source cofactor. The accidental equality s(54)=66=6*(9+2) does not restore the guard; s(150)=222 != 6*(25+2) is a useful complementary control. No subsequent aliquot behavior is inferred.

## Independent controls

The original verifier passed both normally and under `python -O`, reporting 901 guarded edge checks and 60 twin edges with lower prime from 5 through 2000.

`independent_controls.py` imports no production code. It uses factor-product divisor sums and an Eratosthenes prime sieve, with direct proper-divisor enumeration checking the divisor-sum engine on 1 through 1000. It also checks two-prime shift rigidity on every positive core through 1000, perfect cores 6, 28, 496, 8128, actual fixed-6 edge sets against twin-pair sets through lower prime 20000, cutoff counts, actual two-step exits, and source/target/composite hostiles.

Normal and optimized independent outputs agree: 341 fixed-6 twin edges; 1919 starts exit after one step and 341 after two, among all 2260 prime starts from 5 through 20000. All controls pass. These are bounded regression checks supporting the elementary proofs, not evidence promoted to an unbounded prime-distribution theorem.

Reproduce from the package directory:

    python -B review/independent_controls.py
    python -B -O review/independent_controls.py
