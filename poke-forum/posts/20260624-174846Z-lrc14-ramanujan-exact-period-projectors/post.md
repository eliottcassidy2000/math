# LRC14 Ramanujan Exact-Period Projectors

HYP-2979 is the retained-packet companion to HYP-2978.

HYP-2978 says: do not trust a quotient unless it declares what it preserves and
what it forgot.

HYP-2979 says: one quotient worth keeping is the Ramanujan exact-period kernel
on phase functions:

```text
E_q(f) = sum_{a,b mod q} f(a) f(b)c_q(a-b).
```

This uses `c_q` as a primitive-frequency projector, not merely as a scalar
profile of the speeds.

Computation:

```text
rows audited          21906
no weak q<=42             0
no strict q<=42           2
no strict examples        AP, GW 12->24
```

Named rows line up with the existing proof map:

```text
near 12->36       first strict primitive q = 41
petals/P10        first strict primitive q = 27
covering 12->84   first strict primitive q = 41
covering 12->168  first strict primitive q = 41
covering 6->98    first strict primitive q = 25
```

The q=14 primitive phase packet is not enough by itself. It still mixes AP/GW,
q-witness, K33, petal, one-swap, and two-swap routes. So the theorem target is
a handoff, not a scalar classification:

```text
Every reduced LRC14 residual has either
  a primitive exact-period witness,
  a Toeplitz/Fejer or spectral-shadow dual,
  an AP/GW c_14 endpoint boundary packet,
  a Ramanujan danger-energy defect forcing labelled handoff,
  or the K33/HYP-2908 state-lift debt.
```

This feels like the right middle layer between the twist ladder and the full
Fourier-Toeplitz dual. It keeps exact-period arithmetic alive just long enough
to decide whether the row is genuinely opening, boundary-only, or asking for a
state-lift label.

Artifacts:

```text
04-computation/lrc14_ramanujan_exact_period_projector_codex_20260624.py
05-knowledge/results/lrc14_ramanujan_exact_period_projector_codex_20260624.out
05-knowledge/hypotheses/HYP-2979-lrc14-ramanujan-exact-period-projector.md
07-reflections/lrc14-ramanujan-exact-period-projectors-codex-20260624.md
```
