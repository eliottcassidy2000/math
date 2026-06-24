# LRC14 Counterexample Family And Sporadic Classifier

Date: 2026-06-24
Agent: codex-2026-06-24-S153
Related: HYP-2961, HYP-2960, HYP-2955, HYP-2953, HYP-2909, HYP-2905, OPEN-Q-108

The useful shift in this pass is to separate tight equality atoms from strict
counterexample candidates.  AP and Goddyn-Wong keep dominating the search, but
they are not possible strict counterexamples: `qdiv=14` already gives
`M>=1/14`.  They are equality boundary sporadics.  A strict bad row must enter
the covering branch `qdiv>14`.

The resulting classifier has five live families:

```text
L1 apex-multiple residual,
L2 genuine-wide zero-moment,
L3 bounded covering core,
L4 K33 zero-open state-lift,
L5 unnamed source-kernel.
```

Everything else is either already discharged by q-witness, scale peeling, Haar
open fronts, unit-petal/GW-strip migration, K33 state-lift labels, or known
wide concentration machinery.

The challenged assumption was that "sporadic" means "small."  In this proof
route, sporadic means "finite after the live family parameters have been
bounded."  AP/GW are equality sporadics now; future strict bad sporadics would
have to appear only after L1/L2/L3/L4/L5 are bounded into finite atlases.

Tournament Analysis again works best on proof carriers, not runners.  The
chosen vertices are classifier buckets.  The pairwise observable is whether a
bucket preserves strict-counterexample status, has a named proof exit, leaves
infinite parameters, and retains source-spectrum ownership.  Raw
residue/tournament shadows sit last because the recent gauntlets already found
loose impostors.

This is still not LRC14.  It is a cleaner target: prove the five live buckets
empty, or make one of them produce an actual counterexample.
