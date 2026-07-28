---
id: THM-2830
title: "Disjoint positive adjacent-cone factorial moment-three detection"
status: RESERVED / UNPROVED EMPTY STUB
source: root/disjoint-adjacent-cone-factorial-orientation-2026-07-28
depends_on: []
related:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2828-lower-prefix-cone-factorial-moment-three-detection
---

# THM-2830 -- two disjoint factorial prefix cones

**RESERVED / UNPROVED EMPTY STUB.**

Proposed target: let

```text
U=sum_(i<b) lambda_i(f_(i+1)-f_i),
V=sum_(j>=b) mu_j(f_(j+1)-f_j),
```

where both finite coefficient families are nonnegative and nonzero.  Prove
the mixed orientation

```text
2L(V^3)L(UV)-3L(UV^2)L(V^2)>=0,
```

with equality only on the adjacent singleton boundary.  Together with the
strict positive cubic tensor of THM-2828, this would make factorial moments
one through three detect every complex plane spanned by two such disjoint
positive cones, and give a many-versus-many two-charge moment-six theorem
when the constant slot is absent.

The current evidence is only a hostile exact probe: `9,488` exhaustive
small cone pairs, `50,000` random integer cone pairs, and `4,116`
symmetrized quartic coefficient cells are nonnegative, with no negative
case.  No coefficient formula, universal proof, or independently audited
companion is yet installed.  This file supplies no proved dependency.
