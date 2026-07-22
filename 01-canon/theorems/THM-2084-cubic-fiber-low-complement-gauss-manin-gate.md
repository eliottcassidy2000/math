---
id: THM-2084
title: "Cubic-fiber low-complement Gauss--Manin gate"
status: >
  CLAIMED. For a planar complex Keller pair with one output-pencil member of
  source-fiber degree three, the exact depressed-cubic coefficient system is
  derived. The intended theorem is that reduced complementary fiber degrees
  four and five are impossible (degrees zero and three are removed directly,
  while degrees one and two are already tame by THM-2063/THM-2071), forcing
  any unresolved cubic-fiber pair to have reduced complementary degree at
  least seven. The denominator and degree-at-infinity audits are still being
  completed; do not cite the claimed exclusion until this file is promoted.
source: codex-2026-07-22-JC2-cubic-fiber
depends_on:
  - THM-2063
  - THM-2071
---

# THM-2084 -- cubic-fiber low-complement Gauss--Manin gate

Reserved for the exact paper proof and lightweight symbolic referee of the
source-fiber cubic recurrence.  The current working coordinates are

```text
P=z^3+p(x)z+q(x),
J(P,Q)=U(x)[(p'z+q') partial_z-(3z^2+p) partial_x]Q.
```

Writing `Q=A(x,P)+zB(x,P)+z^2C(x,P)` produces a rank-two connection.  The
degree-four and degree-five terminal equations yield, respectively, a conic
and a Weierstrass-cubic first integral together with an exact primitive.  The
remaining proof obligation is to show that polynomiality in the original
fiber coordinate removes every pole introduced by `U`, after which elementary
degree comparison eliminates every nonconstant polynomial branch.

