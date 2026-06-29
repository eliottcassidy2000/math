---
id: HYP-3438
title: LRC14 survivor-gate word audit
status: RESERVED STUB / survivor-gap gluing route, computation pending
source: codex-2026-06-29 continuation of HYP-3436 minimal bad-core cover extractor
tangent: T1399
technique: LTI-399
tournament_technique: LTT-299
script: 04-computation/lrc14_survivor_gate_word_audit_codex_20260629.py
result: 05-knowledge/results/lrc14_survivor_gate_word_audit_codex_20260629.out
reflection: 07-reflections/lrc14-survivor-gate-word-audit-codex-20260629.md
related:
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - THM-523
---

# HYP-3438: LRC14 Survivor-Gate Word Audit

Reserved by codex-2026-06-29 as the immediate local-to-global follow-up to
HYP-3436.

Known input: HYP-3436 classifies

```text
E_safe cap B0_odd cap B1_odd
```

component-by-component and shows that every audited primitive covering row
still has survivor gaps.  The tight row `{1,...,11,13,84}` has only four mixed
even-safe components, with survivor mass `1/105`, while most other even-safe
components are fully bad or clean.

Planned computation: decompose each mixed even-safe component into a gate word
whose letters are survivor gaps and adjacent bad-core blocks.  Each letter
should retain interval endpoints, endpoint wall labels, left/right bad-core
cover pairs, cover-owner deltas, branch mask, and the parent even wall.  The
audit should test whether the survivor gaps are forced by local endpoint
alternation, by owner-current debt, by corridor-fence geometry, or by a named
sidecar from the HYP-3427/HYP-3429/HYP-3431 chain.

Still missing: exact audit output, candidate lemma, tournament fingerprint, and
navigation/reflection synthesis.
