# Parabolic scaling: integer clocks, arithmetic operations, and global invariants

19 October 2026. This session develops the proposed connection between near-cusp quadratic dynamics, base-independent real scaling, and statements about large integers. The results below are scoped theorems and exact transfers. They do not settle Collatz, twin primes, Goldbach, or the planar Jacobian conjecture; historical priority is not asserted.

The [preceding session map](liouville_session2_20261019_progress.md) contains the inherited bank-repair, prime-pair sieve, poset-correlation, and planar-flow results. This map records the new parabolic session and its precise remaining obligations.

## 1. The quadratic clock is more precise than a leading pi approximation

For x_0=0 and x_(k+1)=x_k²+1/4+epsilon, the strict escape count beyond a fixed R>1/2 satisfies N_R(epsilon)=pi/sqrt(epsilon)+O_R(1). This leading law and bounded-error refinement have established literature, credited in the [initial proof](parabolic_escape_precision_20261019/parabolic_escape_precision_20261019.md).

The exact corrected clock has a one-step defect bounded by a constant times the actual distance travelled. This gives convergent endpoint phases, a Cauchy occupation law through the bottleneck, and actual threshold widths asymptotic to 2pi²/n³. For epsilon=b^(-2k), every integer base b>=2 therefore gives N_R=pi b^k+O_R(1). Decimal notation has no special role.

The [next-order matching theorem](parabolic_escape_next_order_20261019/parabolic_escape_next_order_20261019.md) proves, at the exact parameter epsilon_n where the nth iterate equals R,

    n = pi/delta_n + C_R + (pi/3)delta_n
        + O_R(delta_n² log(1/delta_n)),   delta_n=sqrt(epsilon_n).

It resolves the first checkpoint's explicitly open next-order question. The coefficient suggested by the complex fixed-point multiplier is now supported by an actual endpoint/middle matching proof. The [independent certificate](parabolic_escape_next_order_20261019/review/README.md) also sharply encloses C_2, approximately −2.46383665339207711488559.

Across continuously varying epsilon, the additive escape-count error has cluster interval [C_R,C_R+1]. It cannot converge to one constant because the escape count is a staircase. No claim about the digit normality of pi, or the distribution of a fixed geometric parameter sequence, follows from this.

## 2. Smooth parameter maps can perform integer arithmetic

The [arithmetic-operations theorem](cusp_arithmetic_operations_20261019/cusp_arithmetic_operations_20261019.md) constructs closed real-analytic maps of epsilon implementing integer translations of the actual escape count, and an analytic involution implementing a prescribed complementary sum. They use the large inverse branch of the corrected phase; they do not query an escape-band table.

Their algebraic composition laws are exact. Their identities for integer counts hold outside rounding-boundary neighborhoods of relative measure O_R(log M/M²) at count scale M. The constants are proved to exist but are not numerically calibrated in this packet.

A common guard gives an exact integral-to-sum identity: each retained unit phase cell has the same length. Under the calibrated epsilon density

    pi/(2 epsilon^(3/2)) * (1−epsilon/3) d epsilon,

normalized guarded integrals count the corresponding integer pairs exactly. Translation maps are independent of C_R; the complementary reflection retains the endpoint phase. These are maps of parameter space, not an analytic embedding of the original quadratic recurrence.

This supplies exact guarded integral formulations for twin-prime counts and middle-third Goldbach counts. It does not estimate the prime predicates. Positivity of an unguarded integral can be a false positive; the guard is essential.

There is a positive unconditional application: for every sufficiently large target m, a set of ordinary epsilon-measure at least a constant times 1/(m² log⁴m) gives escape counts adding exactly to m, with each of the four neighboring forms 6a±1 and 6(m−a)±1 having at most 37 prime factors. The number-theoretic input is the inherited conservative four-form sieve, not a new prime theorem.

## 3. A finite jet determines every divergent clock term

The [higher-order proof](higher_order_bottleneck_clocks_20261019/HIGHER_ORDER_BOTTLENECKS.md) first treats even monomials and a quartic deformation. The [complete finite-jet classification](finite_jet_bottleneck_clocks_20261019/FINITE_JET_BOTTLENECK_CLOCKS.md) extends it to any real-analytic increment g(y)=a y^(2m)+... with a>0, for a fixed two-sided passage under y→y+epsilon+g(y).

There are exactly m possible divergent powers. With r=2m and epsilon=b^(-rk), they occur at scales

    b^((r−1)k), b^((r−3)k), ..., b^k.

Their coefficients are explicit finite Taylor coefficients obtained by Lagrange inversion. The jet of g through degree 2r−2 suffices to predict the count up to O(1), and that universal jet degree is sharp. Equality of the coefficient vector is exactly the criterion for bounded difference of passage counts. A fixed full passage has no logarithmically divergent term, even for an asymmetric germ; a one-sided passage can have one.

This is a positive use of finite local data to control arbitrarily large integer counts. The bounded remainder still contains the integer-selection phase. It also does not decide analytic flow embedding: an embedded flow and a non-embeddable polynomial germ can have identical divergent clocks. The classification concerns the specified constant-additive families in their chosen coordinates, not arbitrary coordinate changes of the unfolding.

## 4. Collatz has an exact parabolic limit with an arithmetic constraint

The [guarded compactification theorem](guarded_parabolic_compactification_20261019/GUARDED_PARABOLIC_COMPACTIFICATION.md) maps a bank branch sx+e to t/(s+et) under t=1/x. The near-neutral bank family approaches a genuine parabolic Möbius map. Its deformation and formal conjugacy invariant differ from those of the quadratic cusp, so the cusp's inverse-square-root clock cannot simply be imported.

Retaining the lawful parity guard gives a positive theorem: if a contracting parity word W repeats lawfully, its first copy ends below its source or returns exactly. For first-crossing bank words, the existing companion-cycle argument handles exact return. There is no episode cap under this explicit repetition condition. More generally, matching only enough initial parity bits to exceed the carry is sufficient.

The same congruence proves that lawful double repetitions of the near-neutral expansion banks lie outside their formal nonlinear crossover scale. The obstruction is arithmetic, not a defect in the real compactification.

A source-sensitive boundary criterion identifies the next target. For fixed offset d, coefficient stopping on the actual unbounded-source 2-adic limit set is equivalent to uniform finite-time original-height repair above some source height. The actual endpoint equation matters: an infinite symbolic bank path need not be a limit of sources realizing that offset. Arbitrary integer entry and coverage remain open.

## 5. Formal clocks and analytic first integrals are different data

The [cotangent integrability theorem](parabolic_cotangent_integrability_20261019/parabolic_cotangent_invariant_barrier.md) proves that a parabolic cotangent germ has a meromorphic first integral exactly when its normalized iterative logarithm converges. Its complete formal invariant ring is generated by the formal Hamiltonian. For every nonlinear polynomial scalar map, the global rational invariant field of its cotangent lift consists only of constants, extending the earlier polynomial-invariant obstruction.

The [unfolding calculation](parabolic_cotangent_integrability_20261019/formal_quadratic_unfolding_period.md) separates explicit convergent singular clock coefficients from a regular formal remainder that must diverge at the cusp. This explains how useful local passage periods coexist with failure of analytic embedding. The cited classical divergence theorem is distinguished from differential transcendence.

Cotangent lifts here are rational coverings with critical divisors. They are not polynomial Keller maps; finite-jet symplectic completions do not remove that global distinction.

## 6. Why 49 occurs at the neutral three-cycle

The [global period-three index proof](period_three_global_index_20261019/period_three_global_index.md) derives the iterative residue 51/49 at c=−7/4 from the other fixed points of the degree-eight return map. A global residue budget replaces a cubic local expansion.

The denominator has a proved structural explanation: a fixed-point/period-three collision resultant becomes a squared cyclic gap product, which equals the specified monic orbit polynomial's discriminant 49. The real cyclotomic-seven field and its Galois action are inherited facts, now connected to this index calculation. Polynomial discriminants depend on coordinate scaling; the iterative residue is conjugacy invariant.

The same 49 is the exact product of the local quadratic coefficient and parameter-unfolding coefficient. It gives local passage constant pi/7 per return, or 3pi/7 in original quadratic-step units. This is local intermittency, not escape of the real critical orbit: nearby real parameters remain in the real Mandelbrot interval and the critical orbit stays bounded.

## 7. Which arithmetic information scaling must retain

The [residue/correlation transfer note](scaling_residues_correlations_20261019/SCALING_RESIDUES_AND_CORRELATIONS.md) proves exact prefix-residue mixing for a uniformly sampled real input, its Fourier spectrum, and the affine local factors needed for twin-center patterns. This does not supply the same distribution for the fixed number pi. Full-wheel mixing and divisor-truncated sieving are kept separate.

It also gives a genuine poset connection: separate mantissa normalizations can reverse a deletion-count correlation inequality. All Boolean diamond ratios are preserved exactly when the rescaling factors factorize over coordinates; otherwise the exponent correction must be retained. The natural four-antichain example makes the lost information explicit.

## Next substantial proof obligations

- **Collatz:** characterize the realizable endpoint-offset boundary, or exhibit an actual finite guard cover; do not substitute all symbolic infinite paths for realizable sources.
- **Prime pairs:** estimate the guarded arithmetic integrals or the inherited four-form discrepancy. The smooth coordinate maps preserve the correlation problem rather than solve it by one-point density.
- **Clock precision:** turn existential phase-error constants into explicit finite-scale certificates if an algorithm for exact large escape integers is needed.
- **Planar Jacobian:** find a reduction that retains polynomial regularity and global inverse sheets; formal cotangent clocks and finite jets alone do not supply it.

The constructive theme is to keep scale, phase, and lawful residue data together. This session gives exact examples where that enables a transfer, and equally explicit examples where retaining only the leading real approximation loses the needed integer or global information.
