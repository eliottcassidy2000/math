---
id: HYP-3435
title: LRC14 two-adic branch-cover certificate
status: EVIDENCE / finite-ruler interval-certificate stress; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3422/HYP-3425/HYP-3426/HYP-3427, downstream of HYP-3428's loss ledger, HYP-3429's endpoint-spine certificate, HYP-3430's Euler-Mascheroni scalar firewall, HYP-3431's canonical corridor-fence partial proof, HYP-3432's harmonic wall-budget sidecar, and HYP-3424, using HYP-3423 as the quotient-legality guardrail
tangent: T1396
technique: LTI-396
tournament_technique: LTT-296
script: 04-computation/lrc14_two_adic_branch_cover_certificate_codex_20260628.py
result: 05-knowledge/results/lrc14_two_adic_branch_cover_certificate_codex_20260628.out
reflection: 07-reflections/lrc14-two-adic-branch-cover-certificate-codex-20260628.md
related:
  - HYP-3432
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3416
  - HYP-3415
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3435: LRC14 Two-Adic Branch-Cover Certificate

## Claim

HYP-3422 reduced the corrected covering-floor problem to the interval target,
incoming HYP-3425 rewrote the obstruction as the two-color bad core
`B0_odd cap B1_odd`, HYP-3426 mirrored the two branches into a one-branch
interval-piercing target, HYP-3427 put exact wall labels on survivor windows,
HYP-3428 named what the raw two-adic descent quotient is allowed to forget,
HYP-3429 compressed surviving windows to low-rank endpoint spines, and
HYP-3430 showed that Euler-Mascheroni/Mertens/loglog-style scalar tails cannot
replace endpoint sidecars.  HYP-3431 then proved the clean canonical
`{1..11,13,84m}` corridor-fence case for all `m`, and HYP-3432 sharpened the
harmonic side into an exact reciprocal endpoint-budget ranking sidecar.
HYP-3435 now refines that chain into a branch-cell and endpoint-gate
certificate target.

```text
S = O union 2E,    u = 2t,
E_safe(1/14) cap (odd_branch_0_good union odd_branch_1_good) != empty.
```

HYP-3435 sharpens that target into a proof object:

```text
finite E_safe component
+ branch-cell witness
+ active odd/even endpoint-gate ledger
+ minimal failure cover if no witness exists.
```

The promising theorem is not a scalar lower bound on branch measure.  It is a
finite-ruler obstruction statement layered over HYP-3425: an attempted cover of
`E_safe(1/14)` by the two odd branch-bad interval families should have to
expose a small active endpoint ledger.  Those ledgers can then be routed
through HYP-3417/HYP-3420 owner-current sidecars, HYP-3424's covering-floor
transfer, HYP-3425's component-gap Helly certificate, HYP-3426's mirror
reduction, HYP-3427's wall alphabet, HYP-3428's two-adic loss ledger,
HYP-3429's endpoint-spine certificate, HYP-3430's scalar-firewall rule, S302's
energy-sheet sidecar, HYP-3431's canonical corridor-fence base case,
HYP-3432's harmonic wall-budget ordering sidecar, or HYP-3423's
topology-to-magnitude guardrail.

Equivalently, a primitive covering counterexample to the branch lemma should
emit a finite certificate of the form

```text
component J of E_safe
odd_low blockers covering branch 0
odd_high blockers covering branch 1
even endpoint gates defining J
active equality contacts at the attempted witness
```

and the conjecture is that no such certificate can be legal after the known
LRC<=13, owner-current, and quotient-legality sidecars are attached.

## Executable Readout

Script:

```text
04-computation/lrc14_two_adic_branch_cover_certificate_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_two_adic_branch_cover_certificate_codex_20260628.out
```

Stress bank:

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

The smallest branch-union certificate in the bank is the familiar tight row

```text
covering_AP_with_84 = (1,2,3,4,5,6,7,8,9,10,11,13,84)
```

with exact data

```text
even_safe = 107/245 in 26 components
branch0 = 563/105105
branch1 = 563/105105
branch union = 1/105 in 4 components
selected branch = 1
u = 333/1960
t = 2293/3920
selected score = 59/784
active speed = 5
```

The lower bounds inside this finite stress bank are

```text
min_even_safe_measure = 418281361/2753330580
min_branch_union_measure = 1/105
min_selected_score = 1283/17160
threshold = 1/14
selected_score_margin = 401/120120
```

Active selected binders are not purely one arithmetic species:

```text
active_role_hist_top =
  odd_unit: 56
  even_R: 46
  seven_R: 17
  14Q: 16
```

That split is useful.  It says the certificate theorem should keep the active
endpoint-gate role, not collapse the witness to "odd", "even", "7-adic", or raw
branch measure.  It also gives the new HYP-3430 firewall a precise home:
harmonic, Mertens, or loglog tails may calibrate denominator entropy only after
the branch-cell, endpoint-owner, wall/loss, or exact-period sidecar has already
survived.  HYP-3431 supplies the first exact all-parameter success case of that
philosophy: the canonical corridor survives because its actual wall intervals
are retained.  HYP-3432 adds a useful priority queue for endpoint debt, but
still refuses to accept harmonic mass without exact branch/wall/interval
labels.

## The Branch-Cover Failure Form

For each `u in E_safe`, branch `0` fails if some odd speed is too close to
zero:

```text
||o*u/2|| < 1/14.
```

Branch `1` fails if some odd speed is too close to the antipodal side:

```text
||o*u/2|| > 3/7.
```

Thus a failure of HYP-3422 is a simultaneous cover:

```text
E_safe subset (union odd_low_bad) cap (union odd_high_bad).
```

HYP-3435 proposes that this cover cannot happen without a small endpoint-gate
certificate.  The proof should therefore stop asking for a global scalar
average first and instead classify the local finite covers on each component of
`E_safe`.

## Sensitivity Ledger

The script records two exact deletion derivatives for every audited row:

```text
odd_sensitivity(o)  = branch_union_measure(without o) - branch_union_measure
even_sensitivity(e) = branch_union_measure(without e) - branch_union_measure
```

In the tight `covering_AP_with_84` row, the largest listed odd blockers are

```text
7:207341/1261260
1:6/49
11:1543/12740
9:8675/84084
```

and the largest listed even gates are

```text
6:1201/21021
2:205/3822
8:13/294
10:373/9009
```

These are not proof by themselves.  They are a finite-ruler ledger for the next
lemma: a cover of an even-good component must identify which odd low blockers,
odd high blockers, and even endpoint gates are load-bearing.

## Tournament Analysis

Vertices are proof obligations and certificate interfaces, not runners, arcs,
residues, or raw rows.

Pairwise observable:

```text
certificate exactness
branch overlap payload
two-adic induction payload
active endpoint sensitivity
exception-routing compatibility
scalar-shadow risk
Euler-Mascheroni firewall compliance
```

Switch:

```text
higher weighted proof-facing score wins; ties use declared code order.
```

Fingerprint:

```text
vertices=8
score_hist={-5: 1, 51: 1, 69: 2, 77: 1, 90: 1, 93: 1, 97: 1}
directed_3cycles=0
hamiltonian_path =
  C00_finite_ruler_branch_cell_certificate
  -> C01_helly_interval_overlap_theorem
  -> C02_two_adic_descent_induction
  -> C03_active_constraint_sensitivity_ledger
  -> C04_owner_current_exception_router
  -> C05_signed_SPEC_constant_chase
  -> C06_topology_magnitude_guardrail
  -> C07_raw_branch_measure_scalar
```

Assumption challenge:

```text
considered vertices = runners, odd/even speeds, residues, E_safe components,
branch cells, endpoint gates, and proof obligations.
chosen vertices = proof obligations / certificate interfaces.
preserved predicate = exact branch-overlap certificate for the covering floor.
destroyed by this quotient = raw row identity, runner-level geometry, and some
individual endpoint coordinates.
repair = active odd/even gate ledger plus finite exception routing.
```

## Proof Route

1. Prove the finite-ruler normal form: every branch-cover failure has a minimal
   endpoint certificate whose intervals have endpoints in the explicit odd/even
   ruler generated by the row.
2. Classify minimal certificates by active role: `odd_unit`, `even_R`,
   `seven_R`, and `14Q`.
3. Show the `even_R` branch descends by HYP-3418/HYP-3422 two-adic induction or
   produces a smaller covering packet.
4. Route small mixed certificates through HYP-3417/HYP-3420 owner-current
   sidecars.
5. Use HYP-3423 to forbid a topology-only closure when the certificate still
   claims q-specific magnitude.
6. Use HYP-3430 to reject any analytic-tail shortcut that forgets endpoint
   owner, branch cell, wall/loss, exact-period, or state-lift sidecar data.
7. Treat HYP-3431's canonical corridor-fence proof as the base case the
   general minimal-cover extractor should reproduce when specialized to
   `{1..11,13,84m}`.
8. Use HYP-3432 reciprocal wall budgets only to order candidate endpoint debt;
   the accepted certificate still needs exact interval, branch, endpoint, and
   wall labels.
9. Feed any remaining wide nonlocal obstruction into HYP-3421/HYP-3129
   signed-SPEC/Rprime rather than into the failed `t=1/2` odd witness.

HYP-3436 completes the immediate next computation by extracting minimal
branch-cover certificates component-by-component and classifying the
complementary bad core.  The next finite theorem target is now the
minimal-cover no-gluing lemma: every emitted endpoint cover must either
reproduce the canonical corridor fence, descend two-adically, route through
owner-current/Menger, feed signed-SPEC/Rprime, or emit exact-period/state-lift
debt.
