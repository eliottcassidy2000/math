# HYP-1841: LRC counterexamples require arithmetic endpoint cycles, not just endpoint cores

**Status:** OPEN; sharpened by THM-365 and S384 exact audits.

## Claim

After THM-357, THM-359, and THM-365, a Lonely Runner counterexample should be
viewed as a labelled arithmetic endpoint cycle:

```text
e_0 -> e_1 -> ... -> e_r = e_0
```

where each arrow records an owner speed `u`, a protector speed `p`, a center
`m`, a sign `eps`, and a strict integer protection inequality

```text
| p*(n*m + eps) - a*n*u | < u.
```

The conjectural obstruction is not that circular-arc protection cycles do not
exist.  They exist abundantly in the abstract.  The obstruction is that
primitive integer speed sets cannot realize such cycles while also satisfying
the full-measure balance needed for a counterexample.

## Evidence

`lonely_runner_endpoint_cycle_formal_s384.py` shows that bare circular-arc
topology is too weak:

```text
q=3: (0->2), (1->0), (2->1)
```

is already an all-protected full cover with a nonempty core and a directed
endpoint cycle.  Similar abstract mirages appear for `q=3,...,9`, including
short-arc restricted variants.

The same script audits sampled LRC systems.  All sampled systems peel to empty
terminal core:

```text
initial n=14:       coreE=0, unprotected first layer=6
n14 seven-ladder:  coreE=0, unprotected first layer=84
n14 single-gate:   coreE=0, unprotected first layer=24
initial n=15:       coreE=0, unprotected first layer=8
n15 3x5 ladder:    coreE=0, unprotected first layer=150
```

Thus even tiny-gap near-disproofs fail at the endpoint-cycle level after
peeling.

## Predictions

1. Speed-first near-disproofs will usually die by losing all directed cycles
   in the terminal endpoint core, even when their visible gap is very small.
2. Product-sum/torsion coordinates are short labelled-cycle attempts: they
   almost close a low quotient layer but expose higher endpoints.
3. A useful counterexample search should enumerate labelled endpoint cycles
   first, then solve the integer inequalities and full-measure constraints.
4. A useful proof should define a cycle slack potential and show that every
   primitive labelled cycle has positive leak: either a positive gap or a
   peelable endpoint.
5. The unit endpoint divisibility filter should be the first move in this
   cycle obstruction: cycles touching the unit layer must pass through speeds
   divisible by `n`.

## See

- THM-357
- THM-359
- THM-365
- `04-computation/lonely_runner_endpoint_cycle_formal_s384.py`
- `05-knowledge/results/lonely_runner_endpoint_cycle_formal_s384.out`
- `07-reflections/lonely-runner-endpoint-cycle-formal-s384.md`
- HYP-1811, HYP-1828, HYP-1836
