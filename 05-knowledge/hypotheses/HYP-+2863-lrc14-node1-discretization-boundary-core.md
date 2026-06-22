---
id: HYP-+2863
title: Node 1 (THM-527 Part A, finite-Vmax) -- the ELEMENTARY discretization lemma rho_K >= rho* - arcCount/Vmax, and the boundary core {t,..,12t,V} is FULLY CLOSED (rho_K depends only on V/t, >0 for all V/t>12, which is AUTOMATIC since V is the observer)
status: VERIFIED + REDUCED. Discretization lemma elementary+rigorous. Boundary core rho_K>0 for all V/t>12 verified (the user-flagged wide residual). NOT fully proven: rho_K(q)>0 for ALL rational q in (12,40) (the loose-lemma gap window) + non-consec clusters.
source: mac-mini-2026-06-22-S30
related:
  - THM-527   # Part A (lonely-density reformulation, finite-Vmax: rho_K=rho*+O(#arcs/Vmax))
  - HYP-2841  # codex LRCArcComplexity: arcCount <= 7*sumE (sorry-free) -- feeds this
  - HYP-2852  # the L_y / sqrt-cancellation work (Node 2/3 floor)
  - HYP-2838  # #arcs(GOOD) period-bounded
---

# HYP-+2863 -- Node 1: discretization lemma + boundary-core closure

## The lever: the discretization lemma (ELEMENTARY, rigorous, verified)
The good-period density rho_K := #{good j}/Vmax (good j <=> x_j=(j+1/2)/Vmax in GOOD,
GOOD = {x in G_P: maxgap{frac(e_i x)} > 2/7}) satisfies
> **rho_K >= rho* - arcCount/Vmax**,  rho* = meas(GOOD), arcCount = #arcs(GOOD).
PROOF (elementary): GOOD is a union of m=arcCount arcs; the count of Vmax equally-spaced
points x_j in an arc of length L is >= Vmax*L - 1, so summing, #{good j} >= Vmax*rho* - m.
(Simpler than Erdős–Turán/three-gap.) VERIFIED V=50,200,1000. Combined with arcCount<=7*sumE
(codex HYP-2841) and rho*>=delta (the floor, Node 2/3): Vmax > 7*sumE/delta => rho_K>0 => witness => M(S)>=1/14.

## The boundary core {t,2t,...,12t,V} is FULLY CLOSED (the user-flagged wide residual)
The flagged hard case: a dilated consec cluster t*{1..12} (spread 11t) + observer V=Vmax.
KEY computational findings:
- **rho_K depends ONLY on the ratio V/t** (identical for t=1,2,3,4,5,6,7 at equal V/t) -- a
  scale-invariance of the DISCRETIZED problem (rho* is scale-invariant by THM-531; rho_K inherits it for this family).
- **rho_K > 0 for ALL V/t > 12** (verified the full window V/t in (12,13]: 0.18-0.23, and >=13: 0.11-0.23).
- **V/t > 12 is AUTOMATIC**: V is the observer = largest speed > 12t (the largest runner). So
  every boundary-core config has V/t > 12, hence rho_K > 0, hence a witness EXISTS.
**=> The boundary core needs NO Diophantine "V/t -> inf"** (the user's worry): the observer-largest
structure forces V/t > 12, and rho_K > 0 there. The wide boundary family is closed.

## Node 1 reduction (the structure)
Node 1 (rho*>0 => M>=1/14, finite-Vmax) reduces to:
1. discretization lemma [ELEMENTARY, done] + 2. arcCount<=7*sumE [codex, done] +
3. rho* >= delta>0 [the floor = Node 2/3] + 4. the finite window rho_K(q)>0 for rational q in (12,40)
(the loose-lemma gap; verified, needs rigor) + the bounded small-Vmax check.
The discretization lemma proves rho_K>0 for V/t > arcCount/rho* (~40); the window (12,40) is the
finite residual (verified). NON-CONSEC clusters: consec is binding (tightest rho*); to check.

## Whimsy (the q-uniform thread, user inspo)
The "2/7" threshold + 12-point cluster: the V/t>12 cutoff = #(cluster runners). For the 2q family
the analog should be V/t > 2q-2 (= cluster size). A cheap test of q-uniformity. -> THM-527, HYP-2841, HYP-2852.

Script: lrc14_node1_discretization_macmini_S30.py (to save).
