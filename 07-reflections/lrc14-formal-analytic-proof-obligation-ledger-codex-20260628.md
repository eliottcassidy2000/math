# LRC14 Formal/Analytic Obligation Ledger

The useful outcome of this pass is that the proof frontier is smaller than it
felt.

The Lean spine has already done a lot of work: concrete `Mreach` compactness,
denominator-sieve gates, cap RHS arithmetic, the small-k cover pigeonhole, the
concrete witness-floor reduction from the wide `p0` bound, and the gK8
per-shape Delsarte side are all closed or wired as theorem-level glue.  The
remaining hard work is not "formalize everything"; it is to prove a short list
of analytic and rigidity statements that the formal files already know how to
consume.

The most important analytic correction is S81/S82's pull-push pair: the clean
claim that a large speed removes exactly `1/7` of the seed lonely set is false
at the resonant apex, but the resonant danger sits on the 14-grid core while
the useful seed optimum lives off-grid.  So Part A should be formalized around
off-grid bulk survivor positivity and finite-ruler error budgets, not around a
generic `6/7` survival story.

The strongest current route is:

```text
hp0cap decorrelation / gK8 concentration
  -> concrete witness-floor margin
  -> witnessG2 floor + Part A off-grid survivor positivity
  -> concrete Mreach
  -> LRC14
```

The parallel route is tight-locus rigidity:

```text
finite equioscillation + contact graph + blind residue/height sidecars
  -> AP/GW/dilations only
  -> finite tight locus + uniform margin
```

What this ledger changes is workflow.  Future work should not add another
scalar without saying which obligation it reduces.  The best targets are now:

```text
O10 hp0cap decorrelation for binding k=8..12
O12 off-grid bulk survivor positivity / Part A
O15 full tight-locus rigidity
O13 gK8 concentration extremality
O14 doublet R-tail uniform bound
```

The Qsqrt(-7) idea is useful only if it attacks O10 or O12.  The test should be
concrete: pass floor packets through HYP-3267's `zeta_7` contact-holonomy lift,
then ask whether HYP-3254's signed floor terms become positive or finite-chamber
exact in a `Qsqrt(-7)` basis.  For O15, HYP-3257/HYP-3258/HYP-3259/HYP-3265 say
the unit core is a rank-3 index plus contact graph and real-manifold/census
sidecars, not the whole rigidity proof; the formal target must include the
blind residue/height ledger, the binding/covering split, and the contact graph.
