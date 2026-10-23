# Independent runner-comb audit

Date: 2026-10-23. Status: PASS. No Git mutation.

Audited frozen producer deliverables `session36-runner/PROOF.md`, `NOISE.md`, and their accompanying checker. The exact kernel, sharp threshold, finite-sampling and boundary limitations, prime-only valuation extension, and quantitative constants all pass. The single wording correction from “entire function” to “whole phase-response function” was implemented. No mathematical correction was required.

## 5. Exact runner observation kernel

For D_n={t: ||nt||<1/(2b)}, its mean is 1/b and its nonzero Fourier multiplier at k=nm is sin(pi m/b)/(pi m); it is zero when n does not divide k. Probe n=b^j therefore sees exactly the nonzero band b^j Z minus b^(j+1) Z. The first J bands detect all frequencies outside qZ, q=b^J, and detect the mean. The exact kernel is therefore the zero-mean, 1/q-periodic L1 functions. Fourier uniqueness supplies both directions. Convolution here is continuous, so equality almost everywhere in phase is equality at every phase.

This divisibility-band argument works for composite b. The additional claim about *all* speeds with v_p(n)<J is correctly restricted to prime b=p. If q divides k and n divides k, then v_p(n)<J forces p to divide k/n, so those probes add no distinctions beyond the ladder. The extension is not justified by a composite-base valuation law and is not asserted for composite b.

## 6. Sharp bounded-component recovery

For c proper disjoint arcs, the signed product of the 2c endpoint sine factors is a real degree-c trigonometric polynomial, positive inside the carrier and negative outside. Equality of the observations when c<q kills its pairing with h=1_K-g. But 0<=g<=1 makes hP nonnegative pointwise and positive wherever h differs from zero away from finitely many endpoints. Thus h=0 almost everywhere. Regular-closedness upgrades equality almost everywhere to equality of the carriers. Empty and whole-circle cases follow from the mean.

The sharp hostile has q disjoint quarter-cell intervals, and a half-cell translate. Both have measure 1/4 and q components, and their difference is zero-mean q-periodic. All retained responses equal 1/(4b). The next q-speed probe gives 1/(2b) versus zero at phase zero, including b=2: the possible contact with the second carrier is an endpoint only and has zero measure. Thus the class of carriers with at most C components is injectively observed exactly when q>C.

`check_runner_independent.py` uses exact Fraction interval intersections, not bitmask convolution. It checks 1,368 actual breakpoint phases and intervening midpoints across bases 2,3,4,6,7, including q=49, together with independent rational ranks of finite-grid restrictions. The finite-grid tests are explicitly finite-dimensional evidence, not a replacement for the continuum Fourier proof.

The producer's all-grid-knot certification of its larger sharp hostiles is valid: retained comb walls and carrier endpoints lie on (1/(4q))Z; every overlap breakpoint is an endpoint difference on the same grid; and the response is affine between consecutive knots. Equality at all knots thus certifies the complete continuous response for those examples.

## 7. Finite phases and boundaries remain distinct

For any prescribed finite collection of scalar phase samples, the finitely many comb walls leave an open cell where every sampled indicator is constant. Two disjoint equal-length subintervals of that cell have identical observations. They each have just one component. Fixed point-membership samples can be added to the wall set without changing the argument. This does not contradict the component theorem: the latter has entire phase-response functions, not finitely many scalar samples.

Haar observations ignore isolated points and endpoint membership. A regular-closed finite-arc core and its isolated-point sidecar are separate data. The exact AP weak-safe set for speeds 1,...,13 is the six points {1,3,5,9,11,13}/14, and every Haar response is zero. The independent checker verifies every possible wall and every intervening open cell with rational arithmetic. This is not a demonstration that an empty weak-safe carrier is realized by another physical LRC row.

## 8. Quantitative constants

The main proof's bound eta >= e^(2c+1)/[2(4c)^(2c)] is correct. The sharper addendum is also correct. If delta=max_{|k|<=c}|hhat(k)|, then for radius r<=1/2 the endpoint neighborhoods give e <= 4cr+delta/(2r)^(2c). Choosing r=delta^(1/(2c+1))/2 is admissible because delta<=e<=1, and yields e<=min(1,(2c+1)delta^(1/(2c+1))). The delta=0 case follows by letting r decrease to zero.

For measured phase-response error rho, either in supremum norm or in L1, the ladder gives delta<=pi*c*rho/sin(pi/b). Indeed |sin(pi*m/b)|>=sin(pi/b) for integer m not divisible by b, and |m|<=c. The same constant covers the mean because pi*c/sin(pi/b)>=b. Also sin(pi/b)>=2/b and pi<22/7 imply the rational bound pi*c/sin(pi/b)<11bc/7. The addendum's N1--N4 constants and exponents are therefore correct. Empty and full carriers have the direct bound e<=b*rho and need no endpoint polynomial.

These are L1 estimates, not Hausdorff stability or endpoint localization, and make no optimal-exponent claim. A vanishingly small distant component explains the lack of a Hausdorff conclusion without a feature-size condition.

## Reproduction and exact file list

The audit uses `check_runner_independent.py`, which imports only Python's standard library and SymPy. It does not import producer code. Run it with `python check_runner_independent.py` or `python -O check_runner_independent.py` from this directory. The primary checks use explicit exceptions, not assertions.

Files to stage for this audit:

1. `RUNNER_AUDIT.md`
2. `check_runner_independent.py`
3. `runner_independent_receipt.json`
4. `runner_independent.out`
5. `runner_independent_optimized.out`

The receipt records 1,368 exact rational phase checks, six exact rational finite-grid rank/kernel controls, a distinct five-probe finite-sample hostile, and the complete six-point AP weak-safe skeleton. The proof, not finite experiments, establishes the continuum statements. The verdict does not certify literature priority, lawful owner-packet translations, an LRC row closure, or a topology/Hausdorff recovery claim.
