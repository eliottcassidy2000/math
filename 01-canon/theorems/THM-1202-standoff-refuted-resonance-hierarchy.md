---
id: THM-1202
title: Historical standoff refutation and measured resonance hierarchy; exact centre incidence and the continuum ceiling are supplied by THM-1211 and THM-1203
status: HISTORICAL MIXED RESULT.  The positive standoff conjecture is REFUTED by explicit witnesses and the scale-invariant resonance telemetry is retained.  The claimed relation-compatibility classification, 6.33-times next-best bound, and resulting maximiser proof were not established.  THM-1211 gives the exact mod-4 centre-hit criterion, THM-1203 independently proves the sharp 2/21 continuum ceiling and AP equality locus, and THM-1214 closes the clustered r=5 stratum without this resonance proposal
source: kind-pasteur-2026-07-18-S128 (cont.79; owner: prove the standoff bound for non-proportional directions)
depends_on:
  - THM-1181    # the corrected six-box/standoff proposal this note audits
  - THM-1149, THM-1148
related: [THM-1203, THM-1211, THM-1214, MISTAKE-180]
script: 04-computation/standoff_test_kps_S128c79.py (+ .out)
---

# THM-1202 — no standoff; a measured resonance hierarchy and its exact replacements

> **Supersession note (codex-S80).**  The standoff refutation below is valid,
> but its proposed resonance classification was never a proof.  THM-1211
> classifies centre hits exactly by primitive direction modulo four;
> THM-1203 proves the sharp continuum bad-measure ceiling by additive-triangle
> geometry; and THM-1214 closes the clustered `r=5` branch by owner windows.
> Consequently none of the unproved hierarchy claims below is an active LRC14
> dependency.

## (0) The conjecture is false

THM-1150 conjectured that non-proportional integer directions keep standoff > ρ ≈ 0.041 from
the six centres, so their geodesics miss B. They do not:

| direction | standoff |
|---|---|
| (41,73,131) | 0.01748 |
| (97,173,281) | 0.00494 |
| (211,367,593) | **0.00486** |

Far below ρ. This was predictable and I should have checked it before conjecturing: a
geodesic of direction d has length ~|d| in T³, so as |d| grows it equidistributes and comes
arbitrarily close to every point. The standoff route was never going to work.

## (I) They enter B — they just do not linger

The correct mechanism is equidistribution. Sojourn for large generic directions converges to
the volume |B| = 0.003367:

| direction | sojourn | ratio to \|B\| |
|---|---|---|
| (11,19,37) | 0.003470 | **1.03×** |
| (97,173,281) | 0.003560 | **1.06×** |
| (211,367,593) | 0.003000 | **0.89×** |

## (II) The proportional family is a scale-invariant resonance

| direction | sojourn | ratio |
|---|---|---|
| (1,2,3) | 0.09523 | **28.28×** |
| (5,10,15) | 0.09525 | **28.29×** |
| (20,40,60) | 0.09540 | **28.33×** |
| (97,194,291) | 0.09523 | **28.28×** |

Identical at every scale, because the geodesic is confined to the same 1-dimensional
subvariety {g₃ = 2g₂, g₄ = 3g₂}, which threads B. Scale-invariance is the signature of a
resonance rather than an accident.

## (III) Conjectural hierarchy telemetry — not a proved classifier

**(5,9,14) gives 6.33×**, an order between generic and full resonance. It satisfies
**d₂ + d₃ = d₄**. The historical note proposed the rule:

> a relation m·d = 0 concentrates the sojourn **iff** m·(1/4, 1/2, 3/4) ≡ 0 (mod 1).

For m = (1,1,−1) this is (1 + 2 − 3)/4 = 0 — compatible with the centre, hence partial
concentration. Neither my beat framing, my tube framing, nor my standoff framing predicted
this intermediate level; it only appears once the resonance structure is the object.

## (IV) Historical maximiser heuristic — superseded by THM-1203

It carries **two** independent compatible relations (d₃ − 2d₂ = 0 and d₄ − 3d₂ = 0), where a
partial resonance carries one. Two relations confine the geodesic to a curve that lies *in*
the bad configuration locus rather than merely crossing it — hence 28× against 6.33×.

The original text inferred a maximiser claim from these observations.  That
inference was not proved: compatible-relation rank alone does not bound exact
sojourn in the cyclic-gap polytope.  THM-1203 later proves the actual sharp
statement `mu(BAD)<=2/21`, with equality exactly for four-term AP frequencies,
by a non-AP additive-triangle reduction that does not use this hierarchy.

## Honest status

The standoff refutation and measured scale invariance are useful historical
telemetry.  The proposed iff relation classifier and 6.33-times next-best
bound are not theorems.  THM-1211 supplies the exact primitive mod-four centre
incidence law, THM-1203 supplies the continuum ceiling/equality theorem, and
THM-1214 independently closes the clustered `r=5` stratum.  No open LRC14
obligation now depends on completing this resonance heuristic.

## Archived next step

Classifying compatible relations may still illuminate the spectrum of
sub-extremal sojourns, but it is no longer a proof obligation for the ceiling
or clustered `r=5`.  Any revival must first prove that the proposed relation
data controls measure, rather than only centre incidence.
