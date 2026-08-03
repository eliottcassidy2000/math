---
id: THM-3287
title: "Weighted-backbone dominance witness section and selector cut"
status: >
  RESERVED / PROVISIONAL CANDIDATE UNDER INDEPENDENT HOSTILE AUDIT.
  In the canonically weighted support-(1,3), bank-I2 response core, orient
  every maximizing clutch witness from its large-ratio trap state to its
  small-ratio trap state.  Every oriented simple path in THM-3277's
  seven-edge phase backbone, including all fourteen oriented phase
  minimizers, has a compatible witness lift, and the backbone has exactly
  four global sections.  The common state Q+{4} occurs on every backbone
  edge and realizes THM-3278's selector cut as the exact response-sign cut;
  an explicit positive integral blend is trapped there on every core edge
  and on all twelve rows.  This is a static dominance section, not a
  chronological response composition.  The statement remains outside the
  proved dependency graph until its independent audit lands.
source: root/creative-synthesis-recover/2026-08-03
audit: >
  The exact candidate companion rebuilds all eleven reset-link states and 22
  response rows directly from THM-3238's coefficient formulas; reconstructs
  every clutch maximum and oriented witness relation; reverses every edge as
  an orientation control; enumerates all simple paths in the backbone,
  selected tree and full core; exhausts every global section; independently
  replays THM-3277's fourteen phase minima; checks the common Q+{4} witness,
  selector-sign cut, primitive integral blends and canonical nonlift; and
  verifies that all physical Q-directed link transitions instead end at Q.
  Normal, optimized and stored outputs agree, and the source has zero
  assertion nodes and zero floating literals.  A separate hostile audit is
  still required before promotion.
depends_on:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization
  - THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas
  - THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary
related:
  - THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary
  - THM-3273-critical-group-c12-quotient-and-relative-c7-equivariance-boundary
script: 04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py
output: 05-knowledge/results/gmc_backbone_maximizing_witness_section_scout_20260803.out
script_sha256: aed89d67ec7acabfe5b4feae4a83f7c57b78053928be44a6cc319d81fa4a9cc6
output_sha256: 89200bc6cff7284dd33f636352f4f7f56294d90bcd2902fa96092d9f967f5fe3
hash_basis: LF-normalized bytes
---

# THM-3287 -- weighted-backbone dominance witness section and selector cut

**RESERVED / PROVISIONAL CANDIDATE UNDER INDEPENDENT HOSTILE AUDIT.**

The exact candidate and its primary replay are complete.  Nothing in this
file is a proved dependency until a separate hostile audit checks the witness
typing, path universe, common-state sign cut and chronological boundary.

## 1. Oriented maximizing-witness relation

Retain THM-3254's reset

```text
Q=(1,3,3,4,5,6,7,8)
```

and its eleven-state physical link.  For an ordered core edge `u -> v`, write

```text
F_lambda=lambda f_u+f_v,                 lambda>0.       (1)
```

Each link state has one trap interval for `(1)`.  A **small-ratio witness**
has interval `[0,U]`; a **large-ratio witness** has interval `[L,infinity)`.
THM-3269's exact clutch strength is

```text
kappa_{u,v}=max U/L.                                     (2)
```

Define the oriented dominance relation

```text
R_{u->v}={(tau,sigma):
  tau is a large-ratio witness,
  sigma is a small-ratio witness,
  U_sigma/L_tau=kappa_{u,v}}.                            (3)
```

Thus `(3)` follows dominance from row `u` at large `lambda` to row `v` at
small `lambda`.  Direct reconstruction proves

```text
R_{v->u}=R_{u->v}^{op}.                                  (4)
```

The complete relation bank has digest

```text
304204d695d4c7c23e63e1a809ccafc6e5bbbb4a4b97b8b95efeb0ce2c7870a3.
```

## 2. Path lifts and the exact backbone section

An oriented nonempty vertex-simple row path

```text
p=(v_0,...,v_m)                                         (5)
```

**lifts** when there are link states `s_0,...,s_m` with

```text
(s_i,s_{i+1}) in R_{v_i->v_{i+1}}                       (6)
```

for every edge of the path.  Exact dynamic enumeration gives

```text
THM-3277 backbone:  36/36 paths lift,
selected tree T_*:  68/132 paths lift,
full 22-edge core:  2442/21226 paths lift.               (7)
```

The 36 backbone paths have length histogram

```text
((1,14),(2,10),(3,6),(4,4),(5,2)).                      (8)
```

Independent replay of THM-3277's weighted phase ranking gives fourteen
oriented minimizing paths: one for each target except for the two reversal
ties at targets zero and six.  All `14/14` lift through `(3)`.

The stronger gluing statement is exact.  On row order

```text
(2,3,7,10,11,18,19,21,22),                              (9)
```

the backbone has exactly four global sections:

```text
(Q+4,Q-1,Q-7,Q-1,Q+4,X,Q-1,Q-1,Q+4),
X in {Q+1,Q+2,Q+4,Q+5}.                                (10)
```

The only free fibre is row 18.  Restricting any one of `(10)` to a path
explains the complete `36/36` census; it is not merely a collection of
unrelated pathwise choices.  The section-bank digest is

```text
5602fecdf065d428af010d12159f2c966d7127b0a99ffef93dec3ad0b3d5c99d.
```

## 3. A common physical witness and the selector-sign cut

The state

```text
B=Q+{4}=(1,3,3,4,4,5,6,7,8)                            (11)
```

occurs in a maximizing pair on every backbone edge.  With each numerical
edge oriented increasingly, its role is

```text
source: (2,19),(2,21),(18,19),
target: (3,11),(7,22),(10,11),(21,22).                  (12)
```

More strongly, the signs of all twelve exact response values at `(11)` are

```text
positive: {2,11,16,18,22},
negative: {3,7,10,13,17,19,21}.                         (13)
```

These are exactly THM-3278's small/full selector-origin classes.  Hence its
availability bipartition has a literal analytic representative: evaluation
at one physical reset-link state.

This match yields an explicit simultaneous obstruction.  In increasing row
order, take the primitive positive integral weights

```text
 2:      73474965566778464546130
 3:   26090324070418482091313896800
 7:      36737482783389232273065
10:   12696305040937586979075765600
11:    2821401120208352662016836800
13:      39300234402264634486662975
16: 1009443096386399316559081507200
17:     453439465747770963538420200
18:       1209237981608142599589630
19:        714549716404811536121145
21:     470233520034725443669472800
22:         11546066017636615857249.                     (14)
```

Let

```text
M=612782438251008238327605167443444079001600.            (15)
```

Then the weighted contribution `w_i f_i(B)` is `+2M` on every positive row
in `(13)` and `-M` on every negative row.  Every one of the 22 core edges
crosses the cut, so every two-row edge restriction has source value `M` at
`B`.  Its only Q-directed target has value zero, and the increment is
therefore `-M<0`: `B` is trapped on all 22 edge blends simultaneously.  The
full twelve-row blend has source value

```text
5(2M)-7M=3M
 =1838347314753024714982815502330332237004800.           (16)
```

Restricting to the nine backbone rows gives another primitive integral blend
with common edge margin

```text
355648542223452256719445831365899059200                   (17)
```

and full nine-row margin three times `(17)`.

## 4. Sharp failure beyond the backbone

The positive result stops at the exact first gluing obstruction.  In the
canonical rooted tree path

```text
2 -> 17 -> 16,                                           (18)
```

the first edge requires the middle witness `Q-1`, while the second requires
`Q-8`.  These fibres are disjoint, so `(18)` has no lift.  The complete
shortest nonlift census in `T_*` is

```text
2-17-16, 16-17-2, 11-3-16, 16-3-11.                    (19)
```

Thus the backbone property is sharp against a two-edge extension.  It does
not follow merely from every individual edge relation being nonempty.

## 5. Chronological-transition no-go

Every state appearing in `(3)` lies on the eleven-state physical reset link.
Its unique Q-directed physical transition ends at `Q`.  By contrast every
arrow in `(3)` joins two nonreset link states, and none is a physical
Q-directed transition.  Consequently

```text
dominance witness lift != chronological response-state path.              (20)
```

Even the global sections `(10)` say only that adjacent maximizing fibres
share a static witness at their common row.  They do not compose response
updates, transport a physical state, supply a positive current, or continue
after reaching `Q`.

## 6. Map contract and scope

The typed connection is

```text
source: oriented simple row paths in the weighted response core;
target: compatible sections of maximizing reset-link trap states;
map:    each row edge -> its relation R_{u->v};
preserved: exact clutch maximum and equality of the shared middle witness;
lost: ratio endpoint, response magnitude and chronological response state;
sidecars: reset Q and the ordered row-dominance gauge;
hostile: path 2-17-16 with disjoint Q-1/Q-8 middle fibres.                 (21)
```

The candidate is confined to the complete support-(1,3), bank-I2 reset-link
response core.  It is not a Gaussian-moment functional, response composition,
owner phase, same-ancestry carrier, physical walk, `GMC` consequence or
`LRC(14)` decrement.  MISTAKE-354's boundary applies: abstract compatible
realization is not original-response-compatible sequential dynamics.

## 7. Exact reproduction

Run

```text
python3 04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py
python3 -O 04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/gmc_backbone_maximizing_witness_section_scout_20260803.out.
```

The companion pins fifteen established artifacts, reconstructs the response
bank from coefficient formulas, has no assertion node or floating literal,
and prints every consequence object used above.

**End of provisional candidate; independent audit pending.**
