# Independent audit: cubic excellent-Frobenius lift and OPEN next-digit probe

2026-10-24. Verdict: ACCEPTED within the explicit boundaries below. No repository or frozen producer packet was modified.

## 1. Cubic theorem

The packet proves F(t) = (1+4t+48t²) F(t^sigma) modulo 27, with sigma the excellent lift t(q³), not naive t³. This closes only Open Problem 2 of the inherited note, given its modular-polynomial framework. The independently reviewed modular-input audit is preserved in sources/modular-input-audit.md. No external priority claim is certified.

The inherited pole and Fricke bounds force T3 to have coefficients 27,0,A,Bs,C+Dt. The constant product (1−x)(3+x)³ gives A=−18, B=−8, C=−1. The first q coefficient of the norm quotient gives D=24. This is a finite-dimensional modular identity argument, not an inference from 80 matching terms. The exact root R=G/Gsigma is a unit. Dividing the root equation by R² is legitimate modulo 27. Substitution R=sX with s²=(1+t)(1−27t) yields

Q[(1−24t)X²+8X]+18 = 0 modulo 27.

Modulo 3, X(X−1)=0 and X is a unit, so X=1+3a. Reduction modulo 9 gives 3(1+t)(a+t)=0; hence X=1+6t+9b. Expanding gives the residual 9[(1+t)b−t−t³]. Therefore b=t(1+t²)/(1+t) modulo 3. Squaring yields R²=1+4t+21t² modulo 27, which equals the claimed truncation modulo 27. All divisions are by units or by known coefficientwise powers of 3. No unproved algebraic ansatz is used.

The t³ coefficient of F/Fsigma is 760−4=756=27·28, nonzero modulo 81. This establishes sharpness. Neither the next-digit Cartier correction nor the higher-depth tower follows from this theorem. An independent binomial/recurrence negative control explicitly rejects replacing the excellent lift by t³: s7(5)=273504, while the t⁵ coefficient of (1+4t+48t²)F(t³) is 48·4=192. Their difference 273312 is 18 modulo 27. This is distinct from the t³ obstruction modulo 81 for the correct excellent lift.

Independent check.py uses direct lattice enumeration for theta, not the producer's divisor sums. It computes the three residue classes A(q)+zB(q)+z²C(q), z³=q, and the cubic norm A³+qB³+q²C³−3qABC, not the producer's Newton trace method. It obtains t from the Euler eta product E/theta³, and s from q dt/dq = t theta² s, independently of the producer's h expression. All T3 coefficients agree through q119. The excellent mod27 congruence is separately checked through q120. Exact recurrence composition G(t(q)) equals the direct lattice theta through q120. These finite checks support, but do not replace, the infinite proof above.

## 2. Next-digit finite packet

The separate claim C_p²G = G+p²a_p theta_t²G modulo p³ remains OPEN. The frozen design calibrates a_p only at n=1 and uses n² as its fixed shape, with no fitting to the holdouts. This audit verifies that code follows that design; it does not independently establish the chronological creation time of the preregistration.

A separate exact-integer recurrence generation of c(0)..c(40000), checking every division, reproduces all 33 odd-prime calibrations for 3..139 and all 7,995 holdouts n≥2, p²n≤40000. The overlapping strata have 1,983 with p dividing n; 7,162 with p dividing c(n); and 1,864 with both. No stratum is excluded or canceled. Every calibration, per-prime subgroup count, and first linear-model failure agrees with the frozen fixture. The linear model fails for every tested prime; at p=3,n=2 its residual is 18 modulo 27. Normal and optimized independent audit results are byte-identical.

This does not prove the formula for all primes or indices. It does not prove an additional digit at arbitrarily deep Cartier iteration. The p=2 case is excluded correctly because the stated calibration divides by 2. The proof text explicitly identifies the missing second-order lift-error cancellation rather than claiming extrapolation as proof.

## 3. Frozen packages and reproduction

Audited cubic manifest SHA256: 73e83450a128e083bcaa790bec57fa9827c1711d919d4a9ba3f42e4d97b24537.

Audited next-digit manifest SHA256: 44bf84733e7bef8b40acf2f0912d27ef19842a24d827309d4b3eb3f6fbd9338c.

Both producer replay scripts passed recursive exact membership, payload hash verification, normal and optimized recomputation. Both producer hostile harnesses passed under normal and optimized Python, rejecting corrupted proof, missing fixture, and undeclared nested file. These producer harnesses test those three conditions, not every possible adversarial manifest encoding. Frozen manifest hashes were rechecked after audit.

The standalone independent checker and copied evidence are in this audit package. Run python verify.py, or python -O verify.py, from any directory to check manifest closure and rerun both checker modes. Run python hostile.py for copied-package missing, modified, extra, duplicate-manifest-key, and unsafe-path rejection in both interpreter modes. No assert statement controls acceptance.
