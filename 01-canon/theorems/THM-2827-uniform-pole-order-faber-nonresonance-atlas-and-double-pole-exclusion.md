---
id: THM-2827
title: "Uniform pole-order Faber nonresonance atlas and double-pole exclusion"
status: RESERVED / UNPROVED EMPTY STUB
source: root/uniform-pole-order-faber-nonresonance-atlas-2026-07-28
depends_on: []
related:
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2826-uniform-triple-pole-faber-obstruction-and-simple-pole-chamber-exclusion
---

# THM-2827 -- the all-order Faber obstruction has one exact resonance ray

**RESERVED / UNPROVED EMPTY STUB.**

Proposed scope: for every normalized reduced Faber degree `m=4R-2`,
`R>=7`, and every integer `nu=ord(V)>=3` at a point with `ord(M)=1`,
prove the following nonresonance atlas in the nonsplit polynomial
exact-square-prefix chart:

1. an exact-prefix polar cone controlled by the constant term, leading term,
   or coprimality of THM-2760's `P_R,Q_R`;
2. a regular-`d`, pole-`q` cone controlled by the
   `Phi/H/Psi` mod-three pure face; and
3. a cone where `T,s,Td` are regular but `K` has a pole.

The intended exact conclusion is: every `nu` is excluded when
`R!=2 mod 3`; when `R=3k+2`, every `nu` is excluded except possibly

```text
nu=2+(4k+3)delta,                delta>=1,
a=1+(2k+1)delta.
```

This is only a resonance allowance, not existence.  In THM-2796,
`nu=p_j+2`, so every surviving pole part would have to be divisible by
`4k+3`.  In particular simple and double parts are always excluded.  The
first unresolved ray is `R=8,nu=13,a=6,delta=1`, corresponding to pole
part `p=11`.  MISTAKE-317 records why unique face support alone did not
prove more.

Until a proof, exact companion, and independent hostile audit are installed,
this reservation is not a proved result or proved dependency.
