---
id: THM-2101
title: "An additive orbit-residue proof bypasses the small-root product in one-variable DvdK"
status: >
  RESERVED / PROVISIONAL PROOF UNDER AUDIT. Total constant-term vanishing
  makes the small-root residue sum equal to one; transitive Galois incidence
  conflicts with the zero full-root Lagrange sum. The abstract orbit and
  full-root lemmas are kernel-checked, but the analytic-germ-to-splitting-
  field subset bridge is still being written and audited. Do not cite as
  proved until promotion.
source: codex-2026-07-22-GMC2-additive-orbit-residue
depends_on:
  - THM-2067
related:
  - THM-1550
  - THM-1765
  - THM-2022
formalization:
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2LaurentShiftCheckA.lean
---

# THM-2101 -- reserved additive orbit-residue bypass

**RESERVED / UNPROVED PENDING FINAL AUDIT.** The intended proof replaces the
small-root product and its `t`-adic valuation by the residue weights
`alpha^(M-1)/Phi'(alpha)`. The exact analytic-to-algebraic bridge and formal
interface will be stated before promotion.
