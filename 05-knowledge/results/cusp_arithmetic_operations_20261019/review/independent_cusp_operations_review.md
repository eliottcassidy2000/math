# Independent audit: smooth arithmetic operations at the quadratic cusp

19 October 2026. **Publication recommendation: PASS.** The transfer proof,
including the universal translation formula and its parameter-flow
equation, requires no mathematical correction. This conclusion is from an
independent proof audit and separate controls; finite tests do not prove
the asymptotic estimates or supply a numerical onset.

## Scope and dependencies

The audited result is the elementary transfer from
tau(t_n)=n+O_R(log n/n^2) to analytic operations on the escape parameters,
guarded integral identities, and an ordinary-measure almost-prime
application. The next-order phase theorem and its strict escape convention
are inherited from the identified companion packets. The every-target
four-form rough-sieve theorem is an explicitly imported arithmetic input.
This review checks that both inputs imply exactly the conclusions used
here; it does not claim a new independent proof of the cited sieve.

## Proof audit

1. **Inverse branch and domains.** For T>sqrt(A), H(t)=t+A/t is strictly
   increasing on t>T. The stated large quadratic root is its inverse on
   the corresponding range. Xi is decreasing and analytic. The smaller
   root is below sqrt(A) and is inadmissible. F_q maps D_q bijectively onto
   D_-q. Composition is asserted only where each intermediate point
   exists; having F_(p+q) defined alone is insufficient. The reflection
   domain is preserved by u↦m−1−u. Both involution and the sign m−ell in
   R_m R_ell=F_(m−ell) are correct.

2. **Universal parameter translations.** In F_q, the C_R in U cancels
   the C_R subtracted by the inverse. Thus equation (10a) really is
   radius-independent on common cusp domains. Differentiating its phase
   coordinate gives d epsilon/ds=−2 epsilon^(3/2)/[pi(1−epsilon/3)]. This
   is a local flow of parameter transformations. It does not embed the
   original quadratic iteration into a global flow. Reflection retains
   the constant through H(t′)=m−1−2C_R−H(t).

3. **Uniform guard and strict endpoints.** From the imported threshold
   estimate, monotonicity of log x/x^2 supplies a single bound above
   L/4. The displayed choice K_R=128B_R+1 has ample slack. With strict
   escape, the staircase is n on [s_(n−1),s_n), and its value at s_n is
   n+1. The guard therefore gives the claimed ceiling identity without
   differentiating the threshold remainder. For translations, q≥−M/2
   keeps both arguments above the uniform lower scale; arbitrarily large
   positive q causes no new error. All constants and onsets may depend
   on fixed R, as stated.

4. **Reflection and full cells.** For nonintegral u=n−1+theta,
   ceil(m−1−u)=m−n. The offset m−1 is essential: using m instead gives
   m+1, while integer phases give m−1 under the ceiling model. The
   retained cells are exactly n=floor(m/3)+1 through floor(2m/3), which is
   the imported sieve's I_m. Endpoints do not introduce a missing cell.

5. **Integral identities and positivity.** Each guarded phase cell has
   length 1−2eta. This proves the exact weighted sums even for arbitrary
   bounded complex-valued statistics depending on the scale or shift.
   The full-interval error is at most 2||F|| times discarded measure,
   giving precisely 4||F|| M eta, or 4||F|| L_m eta for reflection.
   Positivity equivalence for nonnegative statistics is valid on the
   guarded integrals. It is correctly withheld for unguarded integrals.
   The hostile example is valid; an independent example below shows the
   same issue even with threshold perturbations O(n^-2).

6. **Parameter measure.** The Jacobian is
   −Xi′(u)=2pi^2/[t(u)^3(1−A/t(u)^2)], uniformly comparable to M^-3 on
   comparable-scale windows. Thus the discarded ordinary measure is
   O_R(log M/M^4), and its relative measure is O_R(log M/M^2). The exact
   inverse density is pi(1−epsilon/3)/(2epsilon^(3/2)); it is positive
   because epsilon<3 and independent of C_R. Window locations still
   depend on C_R. Ordinary parameter measure does not count cells
   uniformly.

7. **Every-target almost-prime consequence.** The imported u=38 result
   provides at least c m/(log m)^4 indices in exactly I_m for every
   sufficiently large integer m. Multiplication by the lower cell
   measure of order m^-3 proves the claimed c_R/[m^2(log m)^4]. Each
   individual form is greater than 1 and strictly less than 6m. With
   all prime factors at least (6m)^(1/38), 38 factors would force a
   value at least 6m, so Omega≤37 is correct, including multiplicity.

## Independent reproducible controls

The standalone `independent_cusp_operations_controls.py` imports no
author implementation or packet files. It uses explicit checks, retained
under Python optimization.

- Five exact symbolic phase/Jacobian/flow identities.
- 12,554 exact affine-domain and group-law controls, 180 rejected
  translation exits, a wrong-root example, and every middle-third
  endpoint class for 2≤m<500.
- 24 exact integral controls, computed by direct intersections of pairs
  of staircase cells. The original hostile integral is independently
  reproduced as 1/50. A second model s_n=n+(−1)^n/(100n^2) has a positive
  unguarded wrong-difference integral and zero guarded integral and
  arithmetic sum.
- 75 actual translation pairs and 54 actual reflection pairs at R=2,
  using an independent 75-digit outward-rounded Decimal implementation
  of x_(k+1)=x_k^2+1/4+epsilon. Pi is enclosed by the alternating-series
  Machin identity. The inherited coarse C_2 enclosure is used unchanged;
  every interval is propagated through the original recurrence and
  every earlier iterate is certified at most R before the strict crossing.
- Six exact-rational endpoint/interior recurrence tests at R=3/4.
- The author's verifier was replayed normally and under `-O`, with
  byte-identical output matching its published receipt. The independent
  verifier likewise passes in both modes with identical output.

## Boundaries of the conclusion

The exact value of C_R, the guard constant, and the asymptotic entry scale
are not numerically calibrated here. A fixed approximation to C_R needs
the extra guard stated in the proof. The twin-prime and balanced
Goldbach integrals are exact reformulations, not lower-bound estimates;
the balanced restriction is correctly distinguished from unrestricted
Goldbach. The rough four-form result does not assert actual primes.
No literature-priority claim is assessed. Reviewed-source hashes and all
replay receipts accompany this review.
