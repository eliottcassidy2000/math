# Independent unit-coefficient firewall audit

Status: ACCEPTED counterexample, not a new lower-bound theorem.

The checker uses full antisymmetrization over all permutations for the
Pfaffian and Leibniz determinant sums, independently of the producer's
recursive Pfaffian and Gaussian-elimination methods. Four rational skew
edge values verify Pf(K)=3x-2 and det(K)=Pf(K)^2. The complementary
Pfaffian is 3. At x=2/3 the full Pfaffian vanishes but x²=4/9. Exact
weighted coordinates are u=32/9, v=-5/3, u-|v|=17/9. Eleven active checks
pass identically normally and with Python -O. No floating arithmetic is used.

The failed implication is odd coefficient implies unit coefficient. The
scoped unit-coefficient survivor is valid. The total-margin refinement is
not disproved by this example and remains open. The accepted independence
subset proof uses entirely binary principal sets; the multiplicative
robustness proof uses the full Gram sum and Cauchy-Schwarz, with the actual
coefficient energy. Neither accepted argument assumes unit coefficients.

Reproduce: python3 -B check.py and python3 -B -O check.py.
