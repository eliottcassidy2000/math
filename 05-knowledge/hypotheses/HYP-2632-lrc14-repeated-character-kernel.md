---
id: HYP-2632
title: LRC(14) repeated-residue character kernel - the two-large tail reduces to a small chi_7 integer table
status: OPEN
source: codex-2026-06-19-S23
depends_on:
  - HYP-2630
  - HYP-2631
  - HYP-2626
  - HYP-2624
  - HYP-2617
related:
  - HYP-2614
  - HYP-2619
  - THM-538
  - OPEN-Q-108
---

# HYP-2632 - LRC(14) Repeated-Residue Character Kernel

## Claim Stub

This reserves the next sharp target after HYP-2630.  The two-large
repeated-residue tail should be rewritten as a finite character kernel before
any analytic reciprocal-tail estimate is attempted.

Scratch computation on the `d=9` support-six coimage coefficient shows:

```text
4+2 packet (1,1,1,1,a,a):
  a=2,4      chi_7(a)=+1   S_9 = -25 U
  a=3,5,6    chi_7(a)=-1   S_9 = -18 U
  a=0                       S_9 =  -4 U
  a=1                       S_9 =   0

4+1+1 packet (1,1,1,1,a,b), tail-side finite kernel:
  high surviving entries:  +8 U
  low surviving entries:   +1 U
  two candidate pairs:       0
```

Here `U` is the small `d=9` repeated-kernel unit
`0.009556483535904...`.  The observed HYP-2630 QR/NQR split is therefore not
just a real-valued census artifact.  It appears to be an integer packet table
over the quadratic character and the Jacobi-style signatures
`chi_7(a)`, `chi_7(b)`, `chi_7(ab)`, and `chi_7((a-1)(b-1))`.

## Proof Target

Build the exact finite transform:

```text
S_d(a_1,...,a_6)
 = sum_{r in (F_7^*)^6, a.r=0} C_d(r)
 = (1/7) sum_{t in F_7} C_hat(t a_1,...,t a_6),
```

then express the `4+2` and `4+1+1` repeated packets as a small
`chi_7`/Jacobi character table.  The analytic theorem should bound the
corresponding two-large reciprocal hyperplane sums using this signed table,
not the raw absolute support count.

## Status

Claim reserved from scratch computation only.  Next step is a stored script
and output that verify the additive-Fourier identity, list the integer packet
table after height-2 wall deletion, and record the Tournament Analysis
quotient.
