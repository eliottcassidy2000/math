# Independent review: ternary Frobenius constructions and cyclic lift obstructions

2026-10-19. **PASS in the stated scopes.** This review independently checks the arguments, reconstructs both displayed universal polynomials, and reproduces the finite rank claims. It makes no external-priority claim. The initially suggested all-prime minimum formula is refuted by the independently verified characteristic-29 exception below.

## Frozen targets

- `session6-periods/TERNARY_FROBENIUS_RECONSTRUCTION.md`, SHA-256 `27ab222ffc25bf84d535f76f0ad2dfff5b0bb5d7f30345b4de537b826e9bc60f`.
- `session6-periods/CYCLIC_NORM_LIFT_OBSTRUCTIONS.md`, SHA-256 `81de0a3f14c55e90a3910cf25d359c0dc9aae7fead5ddff128aebe8875bdbc02`.
- `session6-periods/EXCEPTIONAL_PRIME_29.md`, SHA-256 `93555994c2f5dbb23f6123950c6f054d57a2b3529f54fd554ea40915bbceabe8`.
- `next_thresholds.json`, SHA-256 `cf74dd17c3d82bd486a6be2ee6e3f822211b4ddd1736a1cc18c0958155bdfa08`.
- Symmetric numerator certificate `rational_kernel_d12_k2_eps1.json`, SHA-256 `cc432fd1d4f7c7a654ec6225ad10554a0ef53249335a61a81b7e8ac2d5f4e234`.
- Alternating numerator certificate `rational_kernel_d3_k3_eps-1.json`, SHA-256 `74bbc4145000b840aecaec75340a65e09fb9364a4ab570cadd81cbe965ec0c29`.

The alternating saved certificate is the **negative** of the displayed numerator. This was independently discovered, reported, and is now explicitly documented in the target. It does not affect the kernel identity, but direct formula comparison must retain this sign. No mathematical correction remains outstanding.

## 1. Odd supported classification

For primes p >= 5 and odd n, the claimed dimension floor(log_p n) is correct. The support and cyclicity reduction is coefficientwise: after removing pure powers using the invertibility of three, every supported cyclic polynomial is uniquely Phi(q), with q zero on the coordinate axes. Splitting the transposition action is valid because two is invertible. The transposition conjugates A to A squared, so it preserves the norm kernel.

The symmetric odd sector vanishes. In the alternating sector, the residual ternary norm identity forces the odd-in-the-first-variable part of q to be additive. Over F_p, additive exponents are exactly powers of p. Antisymmetry therefore gives the independent generators

    Q_a = X^(n-a) Y^a - X^a Y^(n-a),  a = p^r < n.

Independence holds because a is odd whereas n-a is even. For every available a, the binary norm is the same polynomial F_n = X^n - Y^n - (X-Y)^n. If n is a p-power, F_n vanishes and all coefficients are free; otherwise their sum must vanish. These two counts both give floor(log_p n). The exclusions p = 2,3 are substantive; at p = 3 even the degree-one cyclic sum supplies an exception.

The six-alternant construction from two additive polynomials is also valid. The rational determinant identity cancels the norm before specialization, and row divisibility supplies a polynomial after division by xyz. Taking distinct positive Frobenius powers gives nonzero interior-support solutions. The characteristic-two identity Alt(1,3,7) = xyz times the square of the Dickson determinant was independently expanded.

## 2. Alternating family of degree 3p+3

The displayed seven-alternant formula was reconstructed directly, without reading coefficients from the saved certificate. Its six-variable numerator P has 42 monomials. Exact integer expansion verifies

    (x-y)P - x A(P) + y A^2(P) = 0,

where A acts simultaneously on the original and auxiliary pairs. This is exactly the denominator-cleared rational norm identity: the signs follow from A(xyz) = -y(x-y)z and A^2(xyz) = x(x-y)z. Cyclicity, alternating symmetry, and invariance under (X,Y,Z) -> (X+c x,Y+c y,Z+c z) also hold as exact polynomial identities.

Substituting auxiliary variables by Frobenius powers is valid because the substitution intertwines A in characteristic p. Every resulting exponent is positive for p >= 5, and the seven unordered exponent triples are distinct. The coefficient of x y^3 z^(3p-1) in 2G is one. Thus this family is nonzero for every p >= 5, with no exceptional large prime hidden in a coefficient denominator.

## 3. Symmetric family of degree 2p+12

The compact displayed H formula and cyclic expression were independently expanded into a numerator with 363 monomials. The saved symmetric coefficient certificate matches exactly. The same integer norm identity, cyclicity, full symmetry, and gauge invariance all verify exactly.

After division by xyz, the original-variable exponents lie in [-1,12], and auxiliary exponents in [0,2]. For p >= 17 the map (a,b) -> a+p b is injective coordinatewise on this range, since its width is 13 < p. The coefficient 11 on the term with original exponents (1,12,-1) and auxiliary exponents (0,0,2) cannot disappear modulo such a prime. This proves nonvanishing uniformly.

The remaining primes need separate checks because this injectivity bound does not cover them. Direct specialization verifies coefficient 10 at (1,6,27) for p = 11 and coefficient 11 at (1,12,25) for p = 13. The respective specialized polynomials have 336 and 315 nonzero monomials. At p = 3,5,7 the entire symmetric specialization is zero. Thus the claimed p >= 11 range, including both exceptional small cases, is justified.

## 4. Exact finite scopes and minima

A separately written norm-matrix constructor expands integral binomials on distinct cyclic monomial orbits. It imports no production or author verification code. Independent forward modular elimination gives first nonzero boundary-restriction kernels at:

    p:       2   3   5   7  11  13  17  19
    degree: 11  11  18  24  34  38  46  50

Every lower nonnegative degree was checked, and the first kernel has dimension one in each case. Orbit sums contain each distinct monomial once; this matters for size-one orbits in characteristic three. Independent supported odd dimensions were also checked through degree 39 and at degrees 49,51,121,125,127 for p = 5,7,11,13.

The rational ansatz bounds were independently reproduced with a different orbit enumerator and forward elimination modulo 100003:

- Auxiliary degree one: no solutions in either symmetry for original rational degrees 0 through 30.
- Auxiliary degree two: none through degree 11; degree 12 has a one-dimensional symmetric kernel and no alternating kernel.
- Auxiliary degree three: none through degree two; degree three has a one-dimensional alternating kernel and no symmetric kernel.

These are 96 exact rank cases. Full column rank modulo a prime proves full column rank over Q. For either displayed nonzero rational certificate, the modular rank one below full column rank bounds the rational dimension by one, while the exact certificate supplies the matching lower bound. The denominator and interior-support constraints used in the enumerator agree with the stated ansatz.

These bounded ansatz results do not exclude other rational carriers, and the initial eight finite minima do not prove an all-prime minimum. The statement delta_p = min(3p+3,2p+12) is **false**, as the following separately verified extension shows.

### 4.1 The characteristic-29 hostile and extended minima

The independent forward-elimination implementation checked every degree from zero through the new first kernels:

    p = 23: first degree 58, nullity one;
    p = 29: first degree 42, nullity one;
    p = 31: first degree 74, nullity one.

Every lower degree has zero interior kernel. All three supplied witnesses were also expanded literally to verify cyclicity, strict interior support, homogeneity, and norm zero. These are exact coefficient identities, not evaluations at field elements.

The characteristic-29 note's 23 signed alternants were separately reconstructed directly from its displayed table with SymPy. They produce 138 terms, agree modulo 29 with the census witness, and have coefficient one at x^2 y^12 z^28. The integral norm has 334 nonzero coefficients whose greatest common divisor is exactly 29. Thus degree 42 gives a rigorous counterexample to the proposed minimum 70, independently of the lower-degree rank census. The exception does not affect either universal upper bound, the odd supported classification, or its lift theorem.

The stated integer factor is valid: write u=x^2, v=y^2, w=z^2. The polynomial E is alternating in u,v,w and divisible by uvw. Its remaining alternating polynomial is divisible over Z by the Vandermonde (u-v)(v-w)(w-u), because the three distinct monic linear factors divide it. The quotient is symmetric. Its degree is 15 in u,v,w, or 30 in x,y,z. An independent exact division additionally finds the integral symmetric quotient has 127 terms.

Finally, for the cyclic all-interior degree-42 integral norm matrix, independent construction gives shape 946 by 274 and rank 273 modulo 29. Adjoining the exact divided norm -NE/29 raises the rank to 274. Hence this class has no cyclic all-interior lift modulo 29 squared. Rational injectivity follows from the characteristic-zero restriction theorem. Together with mod-29 nullity one, the failed lift of this nonzero kernel generator proves that the sole 29-divisible Smith invariant has valuation exactly one. The resulting 29-primary torsion is Z/29. This is a finite-prime even-degree interior result, not an extension of the all-degree odd supported theorem.

## 5. All-precision supported lift theorem

The proof over R = Z/p^m, p >= 5, is correct, including at degrees n = p^r. The additive-polynomial lemma follows by inspecting every middle binomial coefficient. For an exponent e > 1 that is not a p-power, Lucas supplies a coefficient prime to p. At e = p^r, all middle coefficients are divisible by p, while binom(p^r,p^(r-1)) has valuation exactly one. Therefore all nonlinear additive coefficients are multiples of p^(m-1), but the linear coefficient may initially be arbitrary.

Repeating the supported reduction over R is legitimate because two and three remain units. It gives

    q = c Q_1 + p^(m-1) sum_(a=p^r<n, r>=1) c_a Q_a.

The exceptional linear coefficient is then controlled by the norm: N_2 Q_1 = F_n integrally, and N_2 Q_a = F_n modulo p. If n is not a p-power, F_n has a unit coefficient, forcing c to be divisible by p^(m-1). If n is a p-power, the minimum coefficient valuation of F_n is exactly one, which again forces c to be divisible by p^(m-1). Degree one is zero by axis normalization. Consequently the full supported kernel over R is exactly p^(m-1) times the mod-p supported kernel, and is killed by p.

The Smith corollary follows without extra assumptions. The integral supported cyclic lattice has the distinct-orbit-sum basis, and its rational norm map is injective in odd degree. The mod-p dimension counts the positive p-valuations of its nonzero Smith factors. The mod-p^2 kernel theorem excludes valuation two or greater. Thus the p-primary torsion of its cokernel is

    (Z/p)^(floor(log_p n)),  p >= 5.

Its prime-to-six torsion order is exactly the prime-to-six part of lcm(1,...,n), and its prime-to-six exponent is the product of primes from five through n. Neither statement controls the 2- or 3-primary torsion.

Independent complete integer Smith forms were computed at every odd degree through 31 and additionally at degree 49. Their p-valuations agree for p = 5,7,11,13, including the repeated p-torsion at the prime-square degrees 25 and 49.

## 6. Full cyclic nonliftability: finite only

For the specified degree-p supported class, lifting cyclically modulo p squared is equivalent to solving the augmented norm equation with right side -F_p/p. Independent literal integer columns and forward elimination reproduce rank increase one at every prime 5 <= p <= 61. At p = 59 the full mod-p norm kernel has dimension two; at the other tested primes it has dimension one. This does not invalidate the obstruction to the specified class.

The support-preserving theorem does not rule out arbitrary full-support cyclic corrections at untested primes. The all-prime full-cyclic nonlift question remains open. Conversely, dropping cyclicity allows lifting through the projector I-N/3; this is consistent with the obstruction rather than a contradiction.

The carry identity and finite-logarithm substitution in the lift note follow directly by expanding (1+T)^p. A single vanishing evaluation, including the Wieferich evaluation at T=1, does not imply that the polynomial obstruction vanishes or that it lies in the cyclic norm image.

## 7. Primary-source attribution

The primary paper Ilin–Ono, *Ternary polynomial relations and periods of cusp forms*, arXiv:2610.19309v1 (13 October 2026), was inspected directly: https://arxiv.org/pdf/2610.19309 . The target correctly attributes the characteristic-zero supported classification and odd-degree vanishing to Theorem 1.1, the normalized binary restriction isomorphism to Corollary 1.3, and the positive-characteristic invariant-product example to Section 7.4. The characteristic-two rank-eight Dickson-module statement is Lemma 5.3, and the integral alternating divisibility statement is Theorem 5.2. The new certificates here are not attributed to that paper. This review verifies the attribution and the elementary constructions used here, not every proof in the source paper.

The finite-logarithm definition was additionally checked directly in the introduction of Elbaz-Vincent–Gangl, *On Poly(ana)logs I*, https://arxiv.org/pdf/math/0008089 . Its attribution is correct; the preprint includes an appendix by Kontsevich.

## 8. Reproduction and implementation separation

Run from this review directory (Python with SymPy and NumPy installed):

    python3 verify_period_certificates_sympy.py > period_certificates_sympy_receipt.json
    python3 verify_period_ranks.py
    python3 verify_period_ansatz.py
    python3 verify_period_lifts.py
    python3 verify_period_next_minima.py
    python3 verify_exception29_sympy.py

The first verifier reconstructs formulas using SymPy exact polynomial operations, independently of the producer's sparse implementation. It records 4,083 explicit checks. The second constructs integer monomial columns and uses forward modular elimination. The third uses an independent six-variable orbit enumeration. The fourth uses the reviewer's own integer norm columns, augmented ranks, and SymPy integer Smith normal forms. No verifier imports any producer code. The formula, extended-minimum, and exception verifiers accept `--base PATH` for their input directory, and otherwise discover either the publication parent directory or the original sibling packet; the other verifiers have no input-path dependency. Normal and optimized reruns of all six pass with byte-identical receipts. Checks use explicit failures, not Python assertions.

Receipts: `period_certificates_sympy_receipt.json`, `period_ranks_receipt.json`, `period_ansatz_receipt.json`, `period_lifts_receipt.json`, `period_next_minima_receipt.json`, and `exception29_sympy_receipt.json`. The corresponding normal/optimized controls are retained beside them. Exact finite computation supports the stated finite claims; the all-degree supported classification, both infinite nonzero families, and all-precision supported theorem rest on the proofs above.
