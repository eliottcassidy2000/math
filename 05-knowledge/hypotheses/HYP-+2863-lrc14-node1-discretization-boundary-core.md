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


## RIGOROUS closure of the boundary core (the widest-arc / s≈0 argument)
The good set in the orbit parameter s (good <=> maxgap{0, frac(1/2-is): i=1..12} > 2/7) has a
MAXIMAL arc of cyclic width ~0.130 around s=0: as s->0 the 12 cluster teeth frac(1/2-is) all
-> 1/2, so the teeth collapse to {0, 1/2} with maxgap = 1/2 > 2/7 (GOOD). The arc persists until
the teeth spread enough to drop maxgap to 2/7 (width ~0.065 each side, ~0.130 cyclic; computed
arcs: (0,.065),(.327,.339),(.476,.524),(.661,.673),(.935,1)).
PROOF that rho_K>0: the sample AP is s_j=(j+1/2)/q, q=V/t. The j=0 sample s_0=1/(2q). For q>7.7
(so 1/(2q)<0.065) s_0 lies in the s≈0 good arc => period j=0 is good => rho_K >= 1/V > 0.
Since V (observer) > 12t, q=V/t > 12 > 7.7 ALWAYS. Hence rho_K>0 for EVERY boundary-core config
{t,..,12t,V}. RIGOROUS (modulo the finite arc-width computation 0.065>1/(2·12)).
KEY INTUITION: at the first ruler period the cluster is maximally collapsed (all teeth near 1/2),
leaving a 1/2-gap -- a guaranteed lonely period -- as long as the observer outpaces the cluster
(V>12t, automatic). The user's "V/t -> inf" worry dissolves: V/t>12 is automatic AND sufficient.
NB: this boundary core has NO separate far-runner set P (all 12 runners are the cluster relative
to observer V); the criterion reduces to maxgap>2/7. GENERAL Node 1 (nontrivial G_P) needs the
G_P-aware version: the good arc must also be far-runner-safe (x_0≈0 fails G_P if a far p<V/7).


## Q-UNIFORMITY: the boundary-core closure ports to the whole 2q family (user's whimsy thread)
The s≈0 collapse argument is Q-UNIFORM. For LRC(2q) (apex prime q, q sectors, ruler maxgap
threshold 2/q, cluster of 2q-2 runners + observer V):
- s≈0 collapse: teeth -> {0,1/2}, maxgap = 1/2 > 2/q FOR ALL q>=5 (2/q <= 2/5 = 0.4 < 0.5). So the
  s≈0 good arc exists for every q.
- rho_K > 0 at V/t > 2q-2 (= cluster size), AUTOMATIC since V (observer) > (2q-2)t.
- VERIFIED q=5 (LRC10): rho_K(V/t=9)=0.333; q=7 (LRC14): 0.231; q=9 (LRC18): 0.176. s≈0 arc
  widths 0.086, 0.065, 0.052 (shrink with q but stay positive).
**=> the Node-1 boundary-core argument ports to the WHOLE composite-2q family** (proven THROUGH the
q-sector structure, not around it), supporting HYP-2846 (q-uniform witness route). The cutoff
"V/t > cluster-size" and "maxgap 1/2 > 2/q" are the q-uniform mechanism. 14=2·7 -> 18=2·9 etc.


## ⚠️ CORRECTION (mac-mini-2026-06-22-S30, same session): the boundary-core closure OMITTED G_P
The "rho_K>0 for all V/t>12" closure above is WRONG. THM-527 splits the runners into
**P = small part (speeds <=13), handled by G_P = {x: ||p x||>=1/14, p in P}** and
**L = large cluster (speeds >13), handled by maxgap{frac(e_L x)}>2/7**. The good-period
criterion is `x in G_P AND maxgap(L)>2/7` -- I computed maxgap over ALL co-offsets, omitting
that the small runners (<=13) live in G_P, NOT maxgap. The s≈0 period (where all teeth -> 1/2)
FAILS G_P: a small runner p has ||p*x_0|| = ||p/(2V)|| ~ p/(2V) < 1/14, so x_0 is unsafe.
VERIFIED: the TRUE rho_K (with G_P) for {1,..,12,V} is **0 for V=29,43,71**, only 0.04 at V=100
(sporadic, when the V-lattice happens to hit the small G_P set, meas(G_{1..12})~0.034). So the
boundary core does NOT close via the s≈0 maxgap collapse.

## The CORRECT structure (what's salvageable + the real residual)
- **Discretization lemma rho_K >= rho* - arcCount/Vmax: STILL CORRECT** (elementary). rho* = meas(GOOD), GOOD = G_P cap {maxgap(L)>2/7}.
- **rho* = meas(G_P cap {maxgap(L)>2/7})** -- the OVERLAP of the small-part safe set and the large-cluster gap event.
  - meas(G_P) > 0 via PROVEN LRC(|P|<=13) [the user's "route to LRC(<=13)"]. (meas(G_{1..12})~0.034.)
  - meas{maxgap(L)>2/7} > 0 via three-distance (mu-floor, THM-527 Part C).
  - rho* = their overlap > 0 = the DECORRELATION of the two scales (P small speeds vs L large speeds, frac decorrelate). This is the FLOOR = Node 3 (R' quasi-independence ~[0.66,1.27]) = my S29 sqrt-cancellation + kps Node-3 spectrum.
- So Node 1 reduces (correctly) to the OVERLAP floor rho*>0, NOT to the maxgap alone. The q-uniform s≈0 result is the maxgap(L) part only.

**Lesson:** the cluster/far split (L maxgap + P G_P) is essential; conflating them omits the
real constraint (G_P). The genuine residual is the overlap/decorrelation floor (Node 2/3).
