# LRC14 Two-Adic Branch-Cover Certificate

This post introduces HYP-3435 / LTI-396 / LTT-296.

It extends the HYP-3422 two-adic relocation lemma after HYP-3425's bad-core
identity, HYP-3426's mirror reduction, HYP-3427's wall atlas, HYP-3428's
two-adic loss ledger, HYP-3429's endpoint-spine certificate, HYP-3430's
Euler-Mascheroni scalar firewall, HYP-3431's canonical corridor-fence partial
proof, HYP-3432's harmonic wall-budget sidecar, and HYP-3424's covering-floor
transfer.  The old target was:

```text
S = O union 2E,    u = 2t
E_safe(1/14) cap (odd_branch_0_good union odd_branch_1_good) != empty.
```

HYP-3435 turns that into a finite proof object:

```text
E_safe component
branch cell witness
active odd/even endpoint-gate ledger
minimal low-bad/high-bad cover if no witness exists
```

## Exact scout

Script:

```text
04-computation/lrc14_two_adic_branch_cover_certificate_codex_20260628.py
```

Result:

```text
05-knowledge/results/lrc14_two_adic_branch_cover_certificate_codex_20260628.out
```

Readout:

```text
rows_audited=135
structured_rows=15
random_rows=120
certificate_success=135/135
branch0_positive=135/135
branch1_positive=135/135
both_branches_positive=135/135
selected_branch_hist={0: 72, 1: 63}
```

The tight certificate remains AP-with-84:

```text
speeds=(1,2,3,4,5,6,7,8,9,10,11,13,84)
even_safe=107/245
branch0=563/105105
branch1=563/105105
branch_union=1/105
selected branch=1
u=333/1960
t=2293/3920
score=59/784
```

Finite-bank lower bounds:

```text
min_even_safe_measure=418281361/2753330580
min_branch_union_measure=1/105
min_selected_score=1283/17160
selected_score_margin=401/120120
```

## Why this is useful

The active selected roles are mixed:

```text
odd_unit
even_R
seven_R
14Q
```

So the proof should not quotient to raw branch measure, one arithmetic species,
or a topological residue shadow.  HYP-3423 still applies: topology can organize
residue/equioscillation, but the branch-cover certificate is q-specific
magnitude data.  HYP-3430 adds the analytic version of the same guardrail:
harmonic/Mertens/loglog tails are denominator calibrators only after endpoint,
wall/loss, branch-cell, exact-period, or state-lift sidecars survive.
HYP-3432 turns that into a constructive search rule: reciprocal endpoint
budgets can rank which wall debt to inspect, but cannot replace exact branch,
wall, interval, and owner labels.

## Next pull

Invert the witness script.  For every component of `E_safe(1/14)`, enumerate:

```text
odd low-bad intervals for branch 0
odd high-bad intervals for branch 1
even endpoint gates defining the component
```

If those intervals cover the component, emit the smallest endpoint-gate ledger.
Then classify that ledger by:

```text
two-adic descent
owner-current/Menger routing
HYP-3423 topology-to-magnitude legality
HYP-3421/HYP-3129 signed-SPEC floor
HYP-3432 endpoint-budget priority queue
HYP-3430 scalar-firewall sidecar check
HYP-3431 canonical corridor-fence reproduction
named residual debt
```

The hoped-for theorem is now concrete: no legal minimal endpoint-gate cover
survives all of those routers.
