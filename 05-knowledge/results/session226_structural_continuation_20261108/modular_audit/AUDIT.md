# Independent adversarial audit: modular valuation correlation and sharp floor

Verdict: PASS, conditional on the explicitly inherited Session225 limiting law and kernels. No mathematical correction required. No claim about priority, height truncation, uniform-in-modulus convergence, or integer coverage follows or is needed.

## Scope and method

Read the proposed `session226-modular/PROOF.md` and `check.py`, and the inherited Session225 `PROOF.md` and `KERNELS.md` in `session226-work/05-knowledge/results/session225_joint_modular_law_recovery_20261108/`. This audit checks the continuation, not a fresh exhaustive audit of the inherited graph/Cesaro theorem. No producer files were changed and nothing was published.

Ran the producer program normally and with Python optimization. Both completed with the published numerical interval. The initially reviewed producer used Python `assert`, making its initial optimized run only a smoke test. The producer subsequently replaced every assertion with explicit `require` calls raising RuntimeError. I inspected that revision and reran both modes successfully. The final producer and the independent program now retain their checks under optimization; no outstanding qualification remains.

Independent program `independent.py` uses the inherited finite conditional sums over the U valuation, rather than copying the producer's closed forms. It verifies normalization, positivity, all adjacent likelihood ratios and the displayed increment quotients for seven primes, exponents 1 through 10 (385 adjacent steps). An independent integer triple enumeration checks every residue modulo 15 against branchwise CRT. It checks 324 valuation-vector pairs for strict MTP2 and density monotonicity and 20 strictly positive threshold covariances under every single-coordinate conditioning in modulus 3^2*5^2*7. A separate integer interval evaluation at scale 10^25 and N=10^6 certifies the published enclosure. Both normal and optimized runs pass.

## Mathematical attacks and resolutions

1. **Prime and exponent boundaries.** For p>=3, the unit A value is 1-2/p^2+1/p^3 > (1-p^-2)^2, and the unit B value is 1-1/p^2+1/p^3 > 8/9. The nonterminal increment ratio is t+1+p/(p^2-1)>2. The terminal ratio is k+1/(1-1/p)>2, including k=1. Since A/B<15/8, both give strict A/B increase. B/C is strict even beyond the C jump because B itself remains strict. No terminal plateau or p=3 exception was found.

2. **Zero residue and multiplicities.** The terminal state v=k correctly means residue zero with multiplicity one. Other multiplicities are (p-1)p^(k-v-1), positive everywhere. They cancel in likelihood ratios and cancel multiplicatively across meet/join. The analysis concerns valuation mass, not the residue density when asserting MTP2.

3. **Shared CRT branch and strictness.** The branch weights remain (9,6,1)/16 across all primes. Pairing each two branches in the expanded determinant gives exactly the positive E factor times (R_I^+-R_I^-)(R_J^+-R_J^-). Both differences are strict precisely when the vectors are incomparable. Comparable-vector equality is an identity. With one prime there are no incomparable pairs; with q=1 the empty product has no strict assertion. This is not a false product-of-mixtures argument.

4. **Association without overusing FKG.** The proof does not need an unquoted FKG theorem. Within each branch, product association follows by covariance decomposition. Strict monotone likelihood ratios and positive support imply strictly increasing expectations of every nonconstant increasing one-coordinate function, by the elementary two-copy covariance identity. Branch covariance therefore proves strict covariance for distinct-coordinate functions. Exact conditioning on any subset leaves all three posterior weights positive and leaves the remaining coordinates product conditional on the branch. Arbitrary event conditioning is expressly excluded. The general association claim is nonnegative, not an unjustified blanket strictness claim for arbitrary functions.

5. **Posterior reinforcement.** Increasing an observation multiplies every ordered branch odds ratio by a factor greater than one. Thus posterior weights move in strict likelihood-ratio order on the three branches. Combining that order with strictly increasing branch means yields stochastic increase of each unobserved coordinate. No reinforcement claim requires independence after mixing.

6. **Residue minima despite decreasing C.** C decreases only on the first valuation step. The proof correctly treats that step separately rather than claiming each branch is increasing. Other-prime A products exceed (8/pi^2)^2>1/2, including the empty product, while other-prime C products are below 2. The displayed k=1 and k>=2 A increments exceed 1/p^2. Consequently 9a Delta A-c/p^2>0 even before including the positive B increment. All subsequent steps are immediately strict. Unit residues are exactly the minimizers and zero is exactly the maximizer for q>1.

7. **Infimum, not an attained minimum.** Unit densities depend on prime support alone. Adding any missing odd prime strictly lowers the mixture, using the negative bracket displayed in the proof. All factors are positive and their deviations from one are absolutely summable. Every finite support has a strictly smaller extension, and exhausting all odd primes converges to the stated three-product constant. This proves the infimum over finite moduli and strict nonattainment, including q=1; no unjustified exchange of arbitrary infinities is needed.

8. **Rigorous tail and interval.** The A and B tail deficits are at most 2/N and 1/N. Their finite prefix products are at most one. The C tail contributes with the opposite sign in C(N)-Csharp and can be dropped for an upper bound. Strict positivity of the total difference follows separately from the prime-inclusion argument, not by ignoring that opposite sign. These facts yield 0<C(N)-Csharp<3/(2N). Directed integer rounding of each positive product encloses the finite products without floating-point assumptions. Independent arithmetic confirms

   0.771866451842 < Csharp < 0.771867951843.

## Receipts

- `check-normal.txt`, `check-optimized.txt`: producer stdout, successful exit codes.
- `independent-normal.txt`, `independent-optimized.txt`: independent stdout, successful exit codes.
- `independent.py`: independent, standard-library-only checker, safe under -O.
- `SHA256SUMS`: hashes of audited producer inputs and audit code/receipts.

Reproduction: run `python3 independent.py` and `python3 -O independent.py` from this directory. The finite tests supplement, rather than replace, the all-prime algebraic proof above.

Publication note: the producer proof status line was changed from audit requested
to audit PASS after review. Mathematical text is unchanged. The local input
hash record preserves reviewed originals; the package manifest covers published bytes.
