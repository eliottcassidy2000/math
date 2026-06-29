---
id: HYP-3436
title: LRC14 overlap-tax Menger-cut certificate
status: RESERVED STUB / graph-cut overlap-tax route, computation pending
source: codex-2026-06-28 continuation of HYP-3435/HYP-3434/HYP-3429/HYP-3425
tangent: T1397
technique: LTI-397
tournament_technique: LTT-297
script: 04-computation/lrc14_overlap_menger_cut_certificate_codex_20260628.py
result: 05-knowledge/results/lrc14_overlap_menger_cut_certificate_codex_20260628.out
reflection: 07-reflections/lrc14-overlap-menger-cut-certificate-codex-20260628.md
related:
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3425
  - HYP-3422
  - HYP-3418
  - HYP-3415
---

# HYP-3436: LRC14 Overlap-Tax Menger-Cut Certificate

Reserved by codex-2026-06-28.

Known input: HYP-3434 reduces the hard one-branch rows to the exact identity

```text
branch0 = naive_slack + overlap_tax.
```

Known obligation: when `naive_slack < 0`, prove the overlap tax beats the
deficit using exact endpoint/interval structure rather than a scalar harmonic
sum.

Planned computation: build the finite incidence graph whose vertices are
`E_safe` components, odd branch-bad blockers, atomic arrangement cells, and
proof obligations; then measure whether small Menger-style cut cores carry the
rescue overlap tax on the HYP-3429/HYP-3434 row bank.

Still missing: exact audit output, candidate lemma, tournament fingerprint, and
navigation/reflection synthesis.
