---
id: THM-4428
title: "LRC14 two-direction network closure and sharp one-direction gap"
status: >
  PROVED by analytic tails and complete FINITE-EXACT heads, independently
  audited. Complete raw carrier supports with at most two primitive
  unoriented directions satisfy the degree-zero network target at every
  height. Three or more directions and LRC(14) remain OPEN here.
source: overnight-hexagon-sep05; rank-one overlap explicitly credited to concurrent THM-4425
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
  - THM-4422-lrc14-projection-deficit-and-beatty-row-reduction
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
proof: 05-knowledge/results/lrc14_two_ray_overnight_hexagon_sep05.md
companion_proof: 05-knowledge/results/lrc14_one_ray_overnight_hexagon_sep05.md
script: 04-computation/lrc14_two_ray_overnight_hexagon_sep05.py
output: 05-knowledge/results/lrc14_two_ray_overnight_hexagon_sep05.out
---

# THM-4428 -- LRC14 two-direction network closure and sharp one-direction gap

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** The complete
[one-ray proof and manifest](../../05-knowledge/results/lrc14_one_ray_overnight_hexagon_sep05.md)
and [two-ray proof and manifest](../../05-knowledge/results/lrc14_two_ray_overnight_hexagon_sep05.md)
are part of this theorem. The generic rank-one closure independently overlaps
concurrent [THM-4425, rank-one carrier closure](THM-4425-lrc14-all-height-rank-one-carrier-closure.md);
that namespace was reserved separately and is not used as a proved dependency.
The sharp non-norm-four gap and exactly-two-direction argument are additional.

Let w=(a,b,c) be primitive, positive, distinct, odd, ternary-unit, with a<b<c.
Retain every owner-live raw carrier

```text
Lambda={C in Z^3: w.C=0, C_i!=0 mod3,
                  |C_i|<3(sum(w)-w_i)/14 for all i}.
```

Let E_i be THM-4414's exact degree-zero network projections.
A direction is an unoriented primitive integer ray, not a real-span dimension;
all positive and negative eligible multiples remain separate carriers.

1. Empty support gives E_i=0. One direction gives min_i E_i<=6/77.
   If its l1 norm is not four, every E_i<=12/161<6/77, sharply and only
   at w=(1,19,23), i=2,3. Its physical mass is instead 12/437.
2. Exactly two directions give E_i<6/77 for every i, at every height.
3. In particular a full-support ternary-unit relation of l1 norm<=14 forces
   the complete support onto one direction and closes the network target.

For one ray v, its complete multiplier list is +/-k v, 1<=k<B_v,
3 not dividing k. Writing M_v=max|v_i| bounds each projection by
12/(49M_v)+4/(7c). In the non-norm-four case M_v>=4; c>=43 closes
analytically, and all 363 admissible triples with c<43 are independently
checked. Norm-four families use THM-4422 and require a selected projection:
(1,5,7) is a hostile to replacing min_i by every i in that branch.

For two rays u,v, short-relation rigidity gives M_u,M_v>=7. Owner residues
force u cross v=k w with 3 dividing nonzero k, hence M_u M_v>=3c/2.
The complete two-ray sum therefore satisfies

```text
E_i <12/343+16/(7c), every i.
```

This is below 6/77 for c>=55. The complete c<55 universe has 814 triples,
192 with exactly two rays, and every projection is strictly below the target.
Independent literal relation boxes and physical sheet networks audit the head.
The first two-ray hostile (17,23,25) forbids discarding an independent ray.
At least three live directions remain in a hypothetical network failure;
chart entry, synchronization, and LRC(14) do not follow.
