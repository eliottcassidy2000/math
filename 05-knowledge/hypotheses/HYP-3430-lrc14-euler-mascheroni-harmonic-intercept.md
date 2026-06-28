---
id: HYP-3430
title: LRC14 Euler-Mascheroni harmonic intercept firewall
status: EVIDENCE / scalar guardrail audit; not an LRC14 proof
source: codex-2026-06-28 extension of HYP-3429/HYP-3428/HYP-3427
tangent: T1391
technique: LTI-391
tournament_technique: LTT-291
script: 04-computation/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.py
result: 05-knowledge/results/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.out
reflection: 07-reflections/lrc14-euler-mascheroni-harmonic-intercept-codex-20260628.md
related:
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3422
  - HYP-3417
  - HYP-3412
  - HYP-3408
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3430: LRC14 Euler-Mascheroni Harmonic Intercept Firewall

## Claim

Euler-Mascheroni enters the current LRC14 proof route only as a finite
harmonic-tail intercept:

```text
H_N - log N.
```

On the HYP-3429 endpoint-spine bank, this scalar tracks scale and window size
but does not determine the endpoint certificate class.  Therefore it is a
calibrator/firewall for Mertens-style denominator tails, not a replacement for
the wall-signature or endpoint-spine certificate.

## Exact Readout

Script:

```text
04-computation/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.out
```

Aggregate readout on the `150` HYP-3429 rows:

```text
endpoint certificate classes:           11
gamma intercept range:                  0.577339668 .. 0.583156236
Euler-Mascheroni reference:             0.577215665
same-max-speed bins with mixed classes: 19/108
rounded-4 gamma bins with mixed classes:21/30
rounded-6 gamma bins with mixed classes:19/108
```

The sharp collision is already visible at `max_speed=84`:

```text
canonical_84m_01:             B1+E/rank2/branches2
covering_AP_with_12_and_84:   E/rank2/branches1
covering_AP_with_84:          B1+E/rank2/branches2
```

All three have the same finite intercept

```text
H_84 - log(84) = 0.583156236...
```

but they have different endpoint-spine classes and different owner pressure.

## Correlation Readout

The scalar sees height/scale:

```text
corr(gamma_intercept, log_best_length) = +0.852948
corr(gamma_intercept, best_rank)       = +0.676986
```

but it barely sees the proof-facing mixed-wall distinction:

```text
corr(gamma_intercept, has_mixed_spine) = +0.139189
```

This matches the earlier Mertens/loglog guardrails: harmonic constants can
calibrate tail budgets, but they forget which endpoint walls certify the
survivor.

## Proof Use

HYP-3430 blocks a tempting analytic shortcut:

```text
harmonic tail size -> covering-floor witness
```

The legal route is instead:

```text
harmonic tail size
-> scalar firewall / denominator-tail budget
-> retained wall signature or endpoint spine
-> HYP-3425/HYP-3422 relocation witness
```

This is useful because any future Mertens, Euler-Mascheroni, or loglog tail
estimate must name the sidecar it preserves: endpoint owners, wall words,
two-adic loss class, exact period, sheet profile, or state-lift debt.

## Tournament Analysis

Tournament vertices are proof carriers and scalar guardrails, not constants.

```text
score_hist={12:1, 29:1, 31:1, 41:2, 43:1}
directed_3cycles=0
hamiltonian_path=
  endpoint_spine_certificate
  -> wall_signature_certificate
  -> two_adic_loss_ledger
  -> harmonic_intercept_calibrator
  -> mertens_loglog_tail_budget
  -> raw_euler_mascheroni_scalar
```

## Assumption Challenge

Alternate vertices considered:

```text
rows, max speeds, harmonic tails, endpoint owners, wall signatures,
loss classes, survivor windows, proof obligations.
```

The chosen compression keeps proof carriers plus scalar guardrails.  It
preserves the covering-floor predicate only when endpoint data is retained.
It destroys the tempting scalar-only route from Euler-Mascheroni or Mertens
constants to an LRC witness.
