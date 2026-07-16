# HYP-7085 — Farey divisor sheets and the Zarankiewicz guardrail

**Status:** CLAIMED / IN PROGRESS (codex-2026-07-16-S20).

This note reserves the exact denominator-layer refinement of the noncommon
owner-wall remainder in `HYP-7084`.  A slow wall at `x` has a unique reduced
coordinate `7x=h/d`; its simultaneous owner set should be exactly

```text
S_d(A)={a in A : d divides a}.
```

Thus the synchronized deck is the `d=1` layer and the remaining term is a
disjoint sum over Farey/divisor sheets `2<=d<=max(A)`.  The immediate task is
to prove and fraction-referee the resulting mechanical-word formula, then test
which quotient data control its signed extrema.

The cyclic Zarankiewicz book is a structural comparator, not yet a bound on
this signed functional.  Important status correction: the ordinary crossing
numbers are known exactly,

```text
cr(K_7,7)=81,  cr(K_7,8)=108,  cr(K_8,8)=144.
```

These equal the Zarankiewicz construction values; they are not open cases.
What remains open here is whether class parity supplies any information beyond
the unsigned quadratic incidence data that the owner-wall quotient destroys.

-> `HYP-7084`, `THM-887`, `THM-913`, cyclic Zarankiewicz `THM-922`, and the
parallel-class-circle reflection.
