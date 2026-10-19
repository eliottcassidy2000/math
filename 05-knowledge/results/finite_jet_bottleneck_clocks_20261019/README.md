# Finite-jet determination of full-passage bottleneck clocks

19 October 2026. **PROVED, scoped:** for the constant-additive family

    y_(n+1) = y_n + epsilon + g(y_n),

where g is analytic with a positive even-order zero r, the full passage between fixed negative and positive endpoints has exactly r/2 possible divergent power terms and a bounded remainder. A coefficient formula computes them from the Taylor jet through degree 2r−2. There is no logarithmic divergence in this two-sided, fixed-endpoint setting.

The complete coefficient vector classifies bounded differences of these specified additive-family passage counts. It is not asserted to be a conjugacy invariant of the unperturbed germ. The jet bound is sharp in general; distinct jets can nevertheless give the same vector. An explicit analytic-flow/poly-nomial pair has identical divergent clocks but opposite analytic-embedding and cotangent-integrability behavior, using the cited classical divergence theorem.

Read [FINITE_JET_BOTTLENECK_CLOCKS.md](FINITE_JET_BOTTLENECK_CLOCKS.md) for the theorem, proof, asymmetric quartic and sextic formulas, and failure boundaries. The one-sided quadratic control explicitly shows a surviving logarithm when the full-passage hypothesis is removed.

Reproduce from this directory:

    python3 verify_finite_jet_clocks.py
    python3 -O verify_finite_jet_clocks.py

Dependencies: Python 3, SymPy, and mpmath. The verifier writes verification-results.json beside itself. It uses explicit checks rather than Python assertions, so its checks remain active under -O. The normal and optimized output receipts are included. Exact formal algebra and numerical diagnostics are separately labeled; neither finite computations nor uncertified floating-point orbit counts replace the analytic proof.

Companion results:

- [Higher-order polynomial bottleneck clocks](../higher_order_bottleneck_clocks_20261019/HIGHER_ORDER_BOTTLENECKS.md)
- [Parabolic cotangent integrability barrier](../parabolic_cotangent_integrability_20261019/parabolic_cotangent_invariant_barrier.md)

No novelty claim, prime claim, or global Collatz-entry claim is made.

## Independent review

[The independent review](review/REVIEW.md) passed with no unresolved correction. Its verifier independently recovers the coefficients through r=12 by differential recurrence and by scaled-denominator beta partitions; it also checks the finite-jet bound, quartic cancellation, and a nonlinear-coordinate hostile.

    python3 review/independent_controls.py
    python3 -O review/independent_controls.py

The author verifier passes 157 explicit checks, and the independent verifier passes 107 exact checks. Both normal/optimized stdout pairs and JSON replays agree. Review source pins and receipts are preserved in review/. SHA256SUMS at this directory's root covers the publication files; review/SHA256SUMS covers the unchanged review bundle.
