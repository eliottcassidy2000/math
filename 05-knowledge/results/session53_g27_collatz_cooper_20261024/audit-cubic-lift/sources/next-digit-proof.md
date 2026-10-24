# A fixed quadratic next-digit correction for the two-step root congruence

Status: FINITE-EXACT held-out evidence; general formula OPEN. This is a separate probe and does not change either frozen theorem or tower packet.

## Fixed model and calibration

For an odd prime p set

a_p = ((c(p²)−2)/p²)/2 modulo p.

The numerator's divisibility by p² follows conditionally from the base theorem and is independently checked in every calibration here. Division by2 is legitimate for odd p. No division by c(n) is used.

The proposed formula is

c(p²n) ≡ (1+p² a_p n²)c(n) modulo p³, n≥0.

Equivalently, with θ=t d/dt,

C_p²G ≡ G+p² a_p θ²G modulo p³.

The coefficient a_p is determined solely from n=1. Thus n=1 is a calibration identity, not evidence for the quadratic dependence. The shape n² and the calibration rule were supplied before the extended execution and are fixed in design.json. Every n≥2 in the declared range is a holdout. A linear n model with exactly the same coefficient is a negative control; no second fitting step is allowed.

## Exact holdout evidence

The verifier generates c(0)..c(40000) by the exact recurrence, checking every division. It tests all33 odd primes3..139 and every n≥2 with p²n≤40000. This ensures at least one true holdout for each prime, including101..139, which would have only calibration data at a20000 cutoff.

All7995 holdouts pass. Overlapping subgroups include1983 with p|n,7162 with p|c(n), and1864 with both. No such values are omitted or canceled; the displayed integer congruence is evaluated directly. For primes whose range contains no member of a subgroup, the count is zero, not a claim of subgroup coverage.

The first coefficients are (p,a_p)=(3,1),(5,2),(7,3),(11,10),(13,4),(17,2),(19,14),(23,18),(29,24),(31,5). Every remaining calibrated coefficient, subgroup count and first subgroup index appears in normal.json. Normal and optimized executions are identical.

The linear negative control fails for every tested prime on a holdout. At p=3,n=2 its residual is18 modulo27. The proposed p=3 formula is therefore

C_3²G ≡ G+9θ²G modulo27,

not G+9θG. In particular c(9)−c(1)≡18 modulo27 agrees with calibration a_3=1 because c(1)=2.

## Relation to prior derivative formulas

Straub, Gessel–Lucas congruences for sporadic sequences, Theorem1.3 and section4, proves a mod-p² digit formula involving a formal derivative of the recurrence sequence. That derivative is not the Euler operation producing n²c(n). The theorem concerns the15 original sporadic sequences; it is not directly a theorem for this square-root sequence. The current probe uses two Cartier steps, modulus p³, and an Euler second derivative, so no implication is asserted. Source: https://arxiv.org/html/2301.12248v1 . The related primary-author slides are https://arminstraub.com/downloads/slides/2025gessel-lucas.pdf .

A bounded search for “Gessel Lucas derivative sporadic Straub” and “A279619 derivative congruence” found the above related formulas but no exact match for this correction. This is not an external novelty claim. The base proof's inherited source already distinguishes excellent-Frobenius formulas from ordinary-lift Gessel–Lucas formulas.

## Structural target and unresolved step

Assuming the conditional base theorem, define the integral error series

E_p=(C_p²G−G)/p².

The proposed statement says E_p modulo p is in the one-dimensional span of θ²G, with scalar fixed by its t coefficient. This is a precise structural target, not an arbitrary coefficient fit.

The base theorem controls only the first-order excellent-lift error and kills its Cartier image modulo p. To identify E_p one order further, the exact p² Eisenstein U_p corrections, second-order lift terms and derivatives must be retained and shown to combine into a multiple of θ²G. The current proof supplies neither that cancellation nor a closed recurrence that uniquely determines E_p from its first coefficient. In particular reduction of a recurrence modulo p can become singular at indices divisible by p, so agreement of initial coefficients is not a substitute for this missing argument.

If the formula holds, the displayed leading correction vanishes whenever p divides n or c(n), matching the tested strata. This alone DOES NOT prove the higher-depth tower: an unspecified O(p³) remainder need not acquire additional p-adic divisibility when n is multiplied by p. The n² leading term suggests, but does not establish, the stronger n-dependent remainder bounds needed to explain two extra digits at every depth.

No statement for p=2 is made because the n=1 calibration divides by2 and must be reformulated there. No model coefficients were changed after inspecting holdouts.

## Reproduction

python replay.py, addressed by its path from any directory, checks exact recursive membership and hashes then reproduces normal and optimized fixtures. python hostile.py checks corrupted-proof, missing-fixture, and undeclared nested-file rejection. The verifier uses explicit exceptions, not removable assertions.
