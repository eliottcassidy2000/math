---
id: THM-2074
title: "Strict LRC(14) holds on a density-one set of thirteen-speed rows"
status: >
  RESERVED WITH COMPLETE PROOF UNDER AUDIT. THM-2051 confines every row
  without a positive-measure strict witness to one of finitely many
  support-three-to-five integer hyperplanes of coefficient height at most
  2^20. Their union has O(B^12) points in the thirteen-dimensional speed
  box, whereas primitive distinct rows have B^13/zeta(13)+O(B^12) points.
  A finite-field union bound also produces whole congruence packets of
  certified rows. This is an almost-everywhere theorem, not LRC(14).
source: codex-2026-07-21-LRC-density-one-hyperplane-sieve
depends_on:
  - THM-2051
related:
  - THM-934
  - THM-2065
script: 04-computation/lrc14_density_one_relation_hyperplanes_codex_20260721.py
result: 05-knowledge/results/lrc14_density_one_relation_hyperplanes_codex_20260721.out
---

# THM-2074 -- Density-one strict LRC(14)

This ID reserves the counting and finite-field consequences of THM-2051's
finite higher-relation trap. The exact constants and optimization-safe
referee will replace this stub after audit.
