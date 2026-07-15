---
id: THM-831
title: Small half-frequency folded targets and two-ball symmetric radius budgets
status: RESERVED / PROOF AND EXACT REPLAY IN PROGRESS
source: codex-2026-07-15-S10 continuation
depends_on: [THM-774, THM-824]
related: [THM-821, HYP-6820]
---

# THM-831 — small half-frequency two-ball radius budgets

Namespace reservation after a live-main pull showed THM-830 assigned and
THM-831 free.

The target under audit is the reduced two-sheet folded set

```text
H_(alpha,beta)={t:||alpha t||+||beta t||>=11/13},
```

where `alpha>beta` are coprime and of opposite parity.  THM-774's exact
component formula implies that `4<=alpha<=9` is precisely its nonempty
one-positive-offset regime.  Preliminary exact arithmetic finds two
antipodal interval components for each of the sixteen admissible
`(alpha,beta)` pairs in this range, and the second-difference separation
needed by THM-824's no-switch argument holds in every row.

Still to prove and freeze: the exact centres and radii, a uniform or finite
exhaustive no-switch proof with all endpoint cases, the common-dilate
corollary, the sharp boundary at `alpha=10`, an independent exact replay, and
the preservation/tournament guardrail.  No claim in this stub may be used as
a theorem until the status is promoted.
