---
id: THM-2060
title: "CRT tail-coset saturation for scaled one-tail rows"
status: >
  RESERVED WITH A PROOF CANDIDATE UNDER ADVERSARIAL AUDIT. In the THM-2059
  fiber product, if the tail period has at least two lifts over the common
  reduction modulus, every reduction class contains a tail-safe residue.
  This would force any zero packet to have a|w and would close every primitive
  scaled-core one-tail row with a>1 once the core has one safe clock. The
  unit-coset and divisibility equivalences still require independent review.
source: codex-2026-07-21-LRC-CRT-coset-saturation
depends_on:
  - THM-2059
related:
  - THM-2057
  - THM-2058
  - HYP-8846
---

# THM-2060 -- CRT tail-coset saturation

In THM-2059 notation, put

```text
Q=Na,       g=gcd(w,Q),       h=Q/g,       d=gcd(N,h).
```

The reserved claim is:

```text
h/d>=2 and the core packet modulo N is nonempty
  => the full packet on the k/(Na) grid is nonempty.        (R1)
```

Candidate mechanism: after dividing `w` by `g`, tail safety is the central
packet modulo `h` multiplied by a unit. Every fiber over one residue modulo
`d` is an equally spaced orbit of size `h/d`; an orbit of size at least two
has a point at circle distance at least `1/4>1/14`. Thus the tail reduction
histogram should have full support, forcing its THM-2059 dot product with any
nonzero core histogram to be positive.

The proposed obstruction identity is

```text
h/d=1  iff  h|N  iff  a|w.                              (R2)
```

If (R1)--(R2) survive audit, settled lower-dimensional LRC supplies the core
clock and proves every row `aC union {w}` with `a>1` and `a` not dividing `w`;
the other branch has a common dilation and reduces to `C union {w/a}`. No
claim is made until the fiber, unit action, small moduli, and exact script have
been checked.
