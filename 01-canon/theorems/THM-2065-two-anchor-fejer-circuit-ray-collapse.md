---
id: THM-2065
title: "Two-anchor Fejer relations collapse the strict LRC residual to circuit rays"
status: >
  RESERVED WITH A COMPLETE PROOF UNDER AUDIT. THM-2051 forces every
  non-strict thirteen-speed row to carry a support-three-to-five relation of
  height at most 2^20. On a fixed rank-two coefficient template, every such
  relation is either an identity among the coefficient rows or cuts out one
  primitive projective parameter direction. Thus a template without a small
  internal circuit has only finitely many non-strict primitive directions;
  on a fixed-N THM-2058 interval, each nonzero relation removes at most one
  integer M. Exact counting, endpoint conventions, and composition with the
  THM-2062 CRT wheel are being independently checked.
source: codex-2026-07-21-LRC-circuit-ray-collapse
depends_on:
  - THM-2051
  - THM-2052
related:
  - THM-2058
  - THM-2062
  - THM-2064
  - HYP-8846
---

# THM-2065 -- Two-anchor Fejer circuit-ray collapse

This ID is reserved for the exact composition of THM-2051's bounded
higher-relation alternative with the THM-2052/2058 two-anchor atlas.

The proposed statement is deliberately conditional on the coefficient
template having no support-three-to-five integer relation of height at most
`2^20`. Without that hypothesis a relation may hold identically throughout
the parameter plane, so the Fejer alternative supplies no projective cut.
The proof, exact exceptional-ray construction, CRT-wheel composition, and
guardrail examples are under audit; no LRC(14) claim is made in this stub.
