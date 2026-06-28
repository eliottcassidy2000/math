# k=8 Ordered-Tail Exchange

The useful new angle is not full distributional dominance.  Consecutive k=8
does not dominate the bounded bank in full stop-loss order: almost every
primitive row beats it somewhere in the low stop-loss functionals.  That makes
full convex order too strong and explains why raw distribution comparisons
keep turning into misleading scalar shadows.

The usable invariant is narrower.  Consecutive speeds are exact top in the
upper ordered tails from `tail_ge_4` onward and in stop-loss from `stop_ge_3`
onward.  The central layer `q3` can move the wrong way, but the exact scout
found a clean exchange rule over all primitive rows:

```text
(q3 - q3_consec)_+ <= (q0+q6)_consec - (q0+q6).
```

That is stronger than the current `L_y=q0+q6+q3/10` need.  It says the odd
Worpitzky/central mass may increase, but only by spending more bimodal mass
than it gains.  The proof idea becomes an insurance premium: price central
mass against the two ordered states, rather than trying to maximize every
moment or every coefficient.

Two dead ends are now explicit.  First, `tail_ge_3` is too low in the tail and
has `431` primitive beaters.  Second, raw `q3` has `2879` primitive beaters.
Both are useful diagnostics, but neither is a theorem target.  Naive
coordinate compression is also not a proof path: the one-step local maxima
counts are `432` for `L_y`, `420` for bimodality, and `363` for
`Sigma kappa_2`.

The tournament analysis should use proof angles as vertices:

```text
central_exchange_rate_lemma
-> upper_stoploss_barrier
-> q0_bimodality_atom
-> primitive_normal_form
-> full_convex_order_route
-> single_step_compression_gradient
-> raw_q3_maximization
-> raw_entropy_route
```

This preserves the LRC predicate through coefficient/tail sidecars and
destroys the runner/arc geometry on purpose.  The challenged assumption is
that tournament vertices must be runners, gaps, or arcs; here the live
vertices are proof obligations and functionals.

Next proof move: try to prove the exchange-rate lemma directly, then join it
to HYP-3200's covariance/ferromagnetic route for the `q0+q6` atom.  Keep
Lee-Yang/Joukowski radius, Worpitzky minority-edge, and root-curve data as
sidecars, not as replacements for the coefficient exchange.
