---
id: THM-4428
title: "LRC14 two-direction network closure and sharp one-direction gap"
status: >
  PROVED by analytic tails and complete FINITE-EXACT heads, independently
  audited. Complete raw carrier supports with at most two primitive
  unoriented directions satisfy the degree-zero network target at every
  height. Later three-direction and universal local-network proofs are
  routed below; chart entry, synchronization and LRC(14) remain OPEN.
source: overnight-hexagon-sep05; rank-one overlap explicitly credited to concurrent THM-4425
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
  - THM-4422-lrc14-projection-deficit-and-beatty-row-reduction
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
proof: 05-knowledge/results/lrc14_two_ray_overnight_hexagon_sep05.md
companion_proof: 05-knowledge/results/lrc14_one_ray_overnight_hexagon_sep05.md
script: 04-computation/lrc14_two_ray_overnight_hexagon_sep05.py
output: 05-knowledge/results/lrc14_two_ray_overnight_hexagon_sep05.out
script_sha256: 07cb6c797679335152eb850071456bf64d8fc7f1fc7da176b8f5f00790684d26
output_sha256: 4b387916df346e1b166aee081d78897bae5c089153a5f068d1b44e0ac817ea1c
one_ray_script: 04-computation/lrc14_one_ray_overnight_hexagon_sep05.py
one_ray_output: 05-knowledge/results/lrc14_one_ray_overnight_hexagon_sep05.out
one_ray_script_sha256: 6b41a879700632aa934650f27dafe9d076c051ddcee3262fabc07556a7aaf875
one_ray_output_sha256: c098ee8a43644f349d21f16257596c84991b732596059bb603db1769cbb73a2f
independent_script: 04-computation/lrc14_two_direction_network_closure_independent_referee_thm4428.py
independent_output: 05-knowledge/results/lrc14_two_direction_network_closure_independent_referee_thm4428.out
independent_script_sha256: c1e0964d82776c3dcf2c49b00955a95ca3a1b104e1e3b1bc43fc469a5d88ec52
independent_output_sha256: d41c3f2473a2dfff26e694b58a2fbb64c041f4e99ad7df2817a6c5eb6cb4950e
hash_basis: raw LF bytes
audit: >
  PASS. The import-free referee reconstructs all 814 finite-head rows with
  independent owner-interval, xy-box, and yz-box engines. On the 192 exact
  two-direction rows it also rebuilds literal sheet networks, verifies every
  determinant and reciprocal inequality, and matches all three projections.
  It audits the THM-4386 quantifier and the distinction between direction
  count and rank. There are 26,636 main and 195 literal-sheet gates; normal
  and optimized outputs byte-match the frozen LF artifact.
---

# THM-4428 -- LRC14 two-direction network closure and sharp one-direction gap

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** The complete
[one-ray proof and manifest](../../05-knowledge/results/lrc14_one_ray_overnight_hexagon_sep05.md)
and [two-ray proof and manifest](../../05-knowledge/results/lrc14_two_ray_overnight_hexagon_sep05.md)
are part of this theorem. The generic rank-one closure independently overlaps
concurrent [THM-4425, rank-one carrier closure](THM-4425-lrc14-all-height-rank-one-carrier-closure.md);
that independent theorem is now PROVED, but is not needed as a dependency here.
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
This theorem originally left at least three directions for a hypothetical
failure. [THM-4431, colored-basis/three-direction closure](THM-4431-colored-lattice-basis-and-three-direction-lrc-network-closure.md)
and the independently audited [universal local-network proof](../../05-knowledge/results/lrc14_global_slope_empty_core_certificate_sep06.md)
advance beyond that historical boundary. Chart entry, synchronization, and
LRC(14) do not follow. The universal proof's THM-4434 reservation is not
used as a proved dependency before its own audited promotion.

The independent referee additionally checks every one of the 814 head rows
through two unrelated lattice boxes, all 192 selected rows through literal
physical sheets, and the tail controls `(1,11,55)`, `(1,5,101)`, and
`(5,49,251)`. Its 26,636 arithmetic gates and 195 sheet gates remain active
under optimization. The semantic head and selected-row digests, exact hashes,
and notation identity `E_i(THM-4414)=S_i(THM-4422)` are frozen in the linked
output; this audit adds no entry or synchronization conclusion.

## Reproduction

Run from the repository root:

```powershell
python -B 04-computation/lrc14_one_ray_overnight_hexagon_sep05.py
python -B 04-computation/lrc14_two_ray_overnight_hexagon_sep05.py
python -B 04-computation/lrc14_two_direction_network_closure_independent_referee_thm4428.py
python -B -O 04-computation/lrc14_two_direction_network_closure_independent_referee_thm4428.py
```
