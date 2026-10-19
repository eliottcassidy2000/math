# Independent review of Wronskian and scalar Frobenius results

2026-10-19. Status: PASS for the stated scope.

The root reviewer read both `WRONSKIAN_FROBENIUS_TRANSFER.md` and
`SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md` and checked the proof mechanisms,
including the weighted coefficient-space echelon bound, discarded-node jet
budget, ordinary/Hasse differential losses, sorted assignment inequality,
commutative power lifting, and exact minimal off-grid coefficient.

The main all-precision lift iff is sound. Necessity iterates the assumed
Frobenius congruence in the same fixed modulus and uses the independently
proved exact grid valuation at exponent p^(m-1). Sufficiency constructs
A=X^d B^(p^(m-1)); modulo p, the coefficient field is F_p, so this really
lifts the prescribed coefficient polynomial. The theorem concerns existence
of an integral Laurent-polynomial lift at each fixed precision. It neither
makes every fixed lift compatible nor supplies a compatible nonmonomial
lift at all precisions.

The minimal-exponent argument is exact because after removal of monomial
shift every exponent is nonnegative: to contribute to the first exponent
not divisible by p, only one copy of that term and otherwise constant terms
can occur. This also explains why discarding the integral lift before the
stated time threshold would be invalid.

Independent computation, using no imports from the producer's scripts,
exhaustively enumerates every coefficient lift in four declared polynomial
boxes and compares existence with the support-coset criterion:

- p=2, m=3, degree<=4: 32,768 lifts, 7 admissible reductions.
- p=3, m=2, degree<=3: 6,561 lifts, 13 admissible reductions.
- p=5, m=2, degree<=2: 15,625 lifts, 13 admissible reductions.
- p=2, m=4, degree<=3: 65,536 lifts, 5 admissible reductions.

Total: 120,490 exact complete-lift checks. Normal and optimized runs agree
byte-for-byte. The finite boxes are controls for the universal proof; they
do not independently exclude lifts of larger degree. The proof does.

Reproduce:

    python3 independent_lift_check.py > independent_lift_check.json
    python3 -O independent_lift_check.py > independent_lift_check_optimized.json

Scope boundaries retained: scalar and commuting ring, integer coefficient
lifts, literal polynomial support before any seed convolution or torus
quotient, separate integral versus pointwise modular lattices. The matrix
and Rule30 hostiles prevent expanding the result beyond its hypotheses.
No conclusion about a prime race, a Jacobian conjecture case, or a major
conjecture's resolution is inferred from these instruments.
