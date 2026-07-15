---
id: THM-846
title: Closed continued-fraction lift/face composites separate seam, theta sheet, and node-fibre data
status: RESERVED — exact coordinate and 58-cell operation-kernel audit in progress; no theorem claim yet
source: codex-2026-07-15-S13f
depends_on: [THM-828, THM-834, THM-838, THM-840, THM-842]
related: [THM-813, THM-829, THM-843, HYP-6880]
planned_verification:
  - 04-computation/n9_closed_cf_lift_face_kernels_codex_S13f.py
  - 05-knowledge/results/n9_closed_cf_lift_face_kernels_codex_S13f.out
  - 05-knowledge/results/n9_closed_cf_lift_face_kernels_codex_S13f.json
---

# THM-846 — reserved closed lift/face kernel audit

This namespace is reserved for the three closed words

```text
F_R = R_10 composed with Phi_(9->10),   R in {A,B,C},
```

where `Phi` is THM-838's centered continued-fraction coordinate copy and the
roles are THM-842's B3 faces.  The verifier will prove the literal coordinate
kernels before minimizing any quotient, then test seam, apex-relative theta
sheet, literal `Q`, bare merged-node marginals, and coupled `bar P` on all 58
THM-828 cells.

Preliminary agent algebra suggests that `F_B` is an idempotent rank-27 map
losing one seam bit, while ordered `(F_A,F_C)` has a coordinate left inverse.
These statements remain provisional until the exact stored replay is added.
