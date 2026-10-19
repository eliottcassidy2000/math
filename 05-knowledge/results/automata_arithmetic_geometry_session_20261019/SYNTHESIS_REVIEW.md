# Independent synthesis audit

2026-10-19. Reviewer: tournament integral-lattice producer, independently checking the synthesis against its proved integral decomposition and elementary linear algebra.

**Verdict: PASS on the mathematical content reviewed below, with the scope/wording refinements below incorporated.** This is a proof audit, not an independent reproof of the supported ternary norm torsion theorem or scalar Frobenius lifting theorem. Section 6 is checked conditional on the stated supported odd-degree norm classification.

## 1. Matrix lifting lemma

The map beta(v) = [M vtilde / p] into coker(M mod p) is well-defined and F_p-linear. Changing the integer lift by p w changes the quotient by M w. Its vanishing is exactly the solvability of Mbar w = -M vtilde/p modulo p, hence exactly the existence of a lift of v in ker(M mod p^2). No additional rank or saturation hypothesis is needed.

The Smith-coordinate interpretation is correct. For precision, if a nonzero diagonal entry has valuation e >= 1, a nonzero mod-p coordinate on that entry admits a compatible kernel lift modulo p^k **if and only if k <= e**. It fails at p^(e+1). A zero Smith entry permits compatible lifts to every precision. For a mixed mod-p class, the maximum permitted k is the minimum e among its nonzero coordinates on nonzero Smith entries; zero-entry coordinates impose no bound. This is about existence of lifts of a fixed mod-p class, not every chosen intermediate representative being extendible.

Thus beta alone detects exactly e=1, and cannot distinguish e=2 from e=3 or an actual integral-kernel direction. The full integral diagonal data does distinguish them. The synthesis now includes the explicit modulus statement above, avoiding an off-by-one or arbitrary-lift interpretation.

## 2. Tournament application and integral normalization

The entries binom(n-r+m,m), with multiplicity t_r-t_(r-1), agree with equations (24)–(26) of PADDING_FROBENIUS_PROOF.md. They are a unimodular diagonal presentation, not necessarily ordered Smith invariant factors. This distinction causes no problem: each coordinate already gives a direct cyclic cokernel summand, and prime-by-prime valuation multisets determine the actual Smith factors.

Consequently the first-lift and complete valuation claims hold for precisely the displayed map D_m:C_n(Z)->C_(n+m)(Z). Kummer's carry count is the valuation of that entry, so one carry means no lift of a nonzero mod-p coordinate to p^2; e carries permit exactly p^1 through p^e. Repeated ordinary deletion contributes the additional scalar m!, exactly as stated.

The adjoint normalization is literal: identify a primitive string D_j kappa with x^[j], pair x^[j] with X^k by delta_(j,k), and note that multiplication by x^[m] is adjoint to the Hasse derivative of order m. Its m!-multiple is adjoint to the ordinary derivative. This identifies an integral operator and its normalization, not the additional evaluation/ideal filtration in the Wronskian application.

## 3. Supported norm and prime-power/zeta formulas

Assuming the stated supported odd-degree norm theorem, its kernel formula forces every nonzero mod-p supported kernel class to fail to lift to p^2. Its p-primary torsion is (Z/p)^(floor(log_p n)), not a cyclic group of order p^(floor(log_p n)). Thus the formula for T(n), its logarithm, and the distinction between multiplicity and lift depth are correct.

For odd n, the interval (n-2,n] contains only one possible odd prime power, namely n itself. Distinct primes cannot share a positive prime power. Hence T(n)/T(n-2) is p when n=p^r with p>=5, and is 1 otherwise, including n=3. The logarithm identity follows by separating the 2- and 3-power terms in psi(n). For Re(s)>1, absolute convergence of sum_p sum_(r>=1) log(p) p^(-rs) justifies rearrangement and gives exactly the displayed logarithmic zeta derivative after deleting p=2,3. There is no analytic estimate hidden in this identity.

The synthesis now phrases the change as a torsion-multiplicity increase and disclaims a comparison map between different degree cokernels. The numerical and abstract-group statements need no such map.

## 4. Remaining boundaries

The distinction between boundary-reconstruction kernels and norm-lifting obstructions is correct and essential. The supported theorem does not decide lifts with genuinely interior corrections. The explicit non-consequences and the caution about physical automaton admissibility are appropriate.

The displayed first-order scalar Frobenius correction identity follows directly from the binomial theorem, also at p=2. Its Fermat-quotient example now states **“for odd p”** for the identity E(1+X)(1)/p = 2 q_p(2), correctly respecting the requirement p not dividing 2. The all-precision nonlinear scalar theorem itself is outside this review's independent scope.

## Evidence and scope

- Audited synthesis: INTEGRAL_LIFTING_CONNECTIONS.md, SHA-256 9445101c59c25954173281fb566ac0b403073758e99294ce5adc0181204032b2.
- Integral tournament proof: frozen session6-automata-publish/PADDING_FROBENIUS_PROOF.md, SHA-256 7efc364d5bcf48ee1d3cb4cb10943119db4c0fccf5c99b36a0e850d46720bedd.
- No checkout or Git writes. No fresh computational census is claimed by this audit.
