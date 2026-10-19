# Independent audit: next-order quadratic cusp phase

2026-10-19. **Publication recommendation: PASS.** No mathematical correction is required in the reviewed next-order theorem. This is an independent review of the supplied proof, with separate exact symbolic, directed-interval, exact-rational, and numerical controls. It is not a claim of worldwide novelty or a substitute for the analytic proof.

## Scope and verdict

The audited extension proves, for each fixed R>1/2 and the exact threshold x_n(epsilon_n)=R,

n = pi/delta_n + C_R + (pi/3) delta_n + O_R(delta_n^2 log(1/delta_n)),

where delta_n=sqrt(epsilon_n) and C_R is the initial note's endpoint phase. Its inverse threshold expansion and adjacent-width correction also pass. The initial note's exact-clock identity, phase convergence, adjacent-width argument, and strict threshold conventions were spot-checked as dependencies. Its coarse C_2 enclosure remains valid.

## Proof audit

1. **Exact defect and uniform regularity.** With w=1+y>=1/2, the integral formula for A_e is smooth on the claimed fixed rectangle. Direct differentiation gives E_y=(-2y^2-y^3-ey-e/3)/(w^2+e). The cancellations E_0(y)=O(y^3) and partial_e E_e(0)|_(e=0)=0 give both derivative estimates and the bound |E_e(y)|<=K_R(|y|+e)(y^2+e). No unproved analyticity at the apparent square root is needed.

2. **Endpoint comparisons.** The forward argument is a valid induction: the preceding point remains in [-1/2,0], where f_0 is increasing and 1-Lipschitz; the reciprocal lower bound then keeps the next point negative. On the reverse side, the explicit lower barrier L=v_j-je is positive and satisfies f_delta(L)<=y_(n-j+1). Because the actual orbit lies in [-1/2,b], where f_delta is increasing, its predecessor is at least L. This establishes positivity before relying on the positive inverse branch. The conditions hold uniformly for j,k<=floor(1/(10delta)), for sufficiently small delta depending on R. The already proved n delta->pi gives 2m<n.

3. **Matched summation.** The first and last blocks are paired with precisely u_0,...,u_(m-1) and v_1,...,v_m. Their perturbations sum to O_R(e log(m+1)+e^2m^2). The signed improved defect causes no problem: the middle is bounded in absolute value using exact telescoping of q=y^2+e, giving O_R(delta^2). The omitted cubic-defect tails are also O_R(delta^2). No differentiation of a limiting Fatou coordinate or unjustified exchange of divergent sums is used.

4. **Constant and universal correction.** At fixed endpoints, the improved clock contributes pi/delta+(pi/3)delta+h(b)-h(a)-(b-a)/2+O_R(delta^2). The two unperturbed improved-defect sums differ from the original sums by exactly -(b-a)/2. Thus the constant is the previously defined C_R, not a newly normalized phase.

5. **Inversion and spacing.** From M=t+pi^2/(3t)+O_R(log(t)/t^2), inversion gives epsilon_n=pi^2/M^2+2pi^4/(3M^4)+O_R(log(n)/n^5). Taking adjacent differences uses only a sum of the two remainder bounds. The n^-4 width coefficient is (6C_R-3)pi^2. An n^-5 spacing coefficient is not justified at this precision and is not claimed.

6. **Foundational phase.** The original nonnegative-defect estimate makes the fixed-m endpoint/middle argument valid; its two infinite sums converge by telescoping squared increments. Writing t_n=n-C_R+o(1) does imply t_(n+1)-t_n=1+o(1), so the original adjacent-width proof is sound. The inference would fail from a mere relative asymptotic; the note correctly avoids that shortcut.

7. **Equality and scope.** Equality x_n=R is not escape: N_R(epsilon_n)=n+1, and the n-th band is (epsilon_n,epsilon_(n-1)]. The next-order theorem concerns thresholds and does not assert a single additive limit for the continuous-parameter escape staircase. Constants depend on fixed R; no uniformity as R decreases to 1/2 is asserted. The logarithmic remainder is an upper bound, not an established nonzero logarithmic term.

## Independent controls

The standalone `independent_phase_controls.py` does not import the author's implementation or source files.

- SymPy verifies the exact derivative, cubic and parameter cancellations, formal inversion coefficient, and adjacent-width coefficient.
- A separately constructed degree-eight Abel normalization yields directed-interval phase certificates at four radii. Its exact residual identity and analytic tail estimate are fully supplied in `independent_phase_certificate.md`.
- The sharper certified interval is

  -2.463836653392077114885602 < C_2 < -2.463836653392077114885581.

- Certificates at cutoffs 256 and 512 overlap; the separately evaluated endpoint shift b=1.5 to b+b^2=3.75 is consistent with the exact phase increment 1.
- Exact rational calculations at R=3/4, epsilon_2=1/4, and epsilon_1=1/2 check both strict-band endpoints and an interior parameter.
- Independent safeguarded Newton shooting checks n=30,100,300,301,1000 for b=0.05,0.25,1.5,10. At every computed orbit, the summed improved clock identity holds to the prescribed high precision. Residual and adjacent-width tables are finite numerical controls, not asymptotic certificates. Large constants near b=0 are compatible with the theorem's R-dependent bounds.
- The author's next-order verifier was rerun; its output matches the packaged output byte-for-byte. Both supplied package manifests passed before the review-status-only amendments.

## Files and boundaries

`source_fingerprints.json` records the reviewed source and package fingerprints, including final review-status amendments when available. `independent_phase_controls.out` is the complete independent output; `reproduced_author_controls.out` is the reproduced author output. Python 3.12.14, SymPy 1.14.0, and mpmath 1.3.0 were used.

The analytic theorem passes independently of the finite numerical experiments. The review does not remove the logarithm from the error, establish a further phase coefficient, give an effective global asymptotic-entry threshold, or conduct a literature-priority audit. No Git state was changed by this review.
