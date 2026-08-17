---
id: THM-3542
title: "Fixed Keller depth-three newest-prime residue point/pair saturation"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.  Exact rational
  specialization presently gives point-factor degrees 1,2,6 and injective
  unordered-pair-sum factor degrees 1,2,6,6,9,12.  The specialization,
  no-extraneous-root, and generic-orbit transport gates are not yet promoted
  here.  This stub is not a proved dependency.
source: codex/depth-three-residue-saturation/2026-08-16
depends_on: []
related:
  - THM-3539-fixed-keller-newest-prime-decomposition-centralizer-and-lca-packet-floor
  - THM-3540-fixed-keller-depth-two-newest-prime-residue-saturation
---

# THM-3542 -- reserved depth-three residue-saturation namespace

**RESERVED / UNPROVED PROVISIONAL PROOF CANDIDATE UNDER AUDIT.**

The candidate specialization is the normalization point
`(tau,lambda)=(0,1)`, equivalently

```text
q0=(1/9,4/3,0) in V(L),   q1=F(q0),   q2=F^2(q0) in V(P2).
```

Exact exploratory elimination produces a squarefree degree-nine `x`
polynomial with irreducible-factor degrees `1,2,6`.  Its 36 unordered pair
sums are distinct, and their candidate resolvent has factor degrees
`1,2,6,6,9,12`.  These match the point/pair orbit sizes of the marked-leaf
stabilizer on the depth-two ternary tree.

Promotion requires a complete exact companion and prose audit proving that
the eliminant is the actual nine-point fibre, that this is a good
specialization of the `kappa(P2)` predecessor cover, and that specialization
orbits refine the generic decomposition action.  Until then this file proves
nothing and may not be cited as a dependency.
