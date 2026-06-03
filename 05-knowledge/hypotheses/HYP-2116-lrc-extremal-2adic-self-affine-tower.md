---
id: HYP-2116
status: STRUCTURE -- extremal family = 2-adic self-affine tower; doubling self-similarity (verified exact); 2-adic stratification ratio 2; n=14 difficulty isolated to odd layer above prime 7
source: opus-2026-06-03-S579
related: [HYP-2115, HYP-2091, HYP-2063, HYP-2097, THM-369]
---

# HYP-2116: the extremal family is a 2-adic self-affine tower

**TWO RECURSIONS:** ladder A (n->n+2, append antipodal pair, flip in {0,1}; doubles worry-set 2^((n-2)/2); index=flip-set F) and doubling D (v->2v, t->t/2, n->2n; FRACTAL scaling; index=2-adic valuation v2). Product = integer affine action <x2,+2>.
**DOUBLING SELF-SIMILARITY (verified exact all n incl 14):** the even layer of AP_{2n-1} at t=1/(2n) equals AP_{n-1} at t=1/n (identical phase multiset). D: AP_{n-1} -> {2,4,..,2n-2} subset AP_{2n-1} is a phase-halving embedding.
**2-ADIC STRATIFICATION (ratio 2):** every speed v=2^r*odd; stratify by r=v2(v). At t=1/n the per-stratum MIN phase doubles with r: stratum r binds at ~2^r*delta. ODD stratum (r=0) is the UNIQUE binder ({1,n-1} at delta); higher strata 2x,4x.. safer. (n=14: r=0:0.071=delta, r=1:0.143, r=2:0.143, r=3:0.429.)
**SELF-AFFINE TOWER:** AP_{n-1} = disjoint union over r of 2^r * O_{ceil((n-1)/2^r)} (O_m=odd generators). Each sleeve measure exactly 2delta but arc geometry refines by depth (coarse wide arcs at r=0, fine thin at large r) = self-affine multifractal cover, contraction ratio 2, binding at the coarse top. D shifts the tower up one level (r->r+1, phases x1/2, fresh odd binder at bottom) = IFS fixed structure.
**PAYOFF n=14=2*7:** for n=2^a*q (q odd) even layer = D-scaled from level q with EXACTLY x2 slack; odd layer binds at delta (verified n=14,28,12,24: even slack x2.0 exactly). At n=14: even {2..12}=2*{1..6}=D(AP_6)=PRIME-7 AP (solved, phase>=1/7); odd {1..13} binds at 1/14 = the entire extremal constraint. WHOLE DIFFICULTY of n=14 is in the odd layer, ONE doubling above the solved prime 7.
**EXTENSIONS:** (1) inductive C'/loneliness up the doubling tower: D lifts level-q loneliness to even layer of 2q (x2 slack), residual=odd layer; reduce 14 to an odd-layer lemma off base prime 7. (2) apex 2q (HYP-2063) = the 2-adic seam (v=q largest odd binder, 2q=0). (3) stratify worry-set by (r,F) -- 2-adic depth x flip-set, self-affine Cantor boundary of the <x2,+2> tree. (4) saturation top-heavy: deficit Phi governed by r=0 stratum.
**See:** `07-reflections/lrc-extremal-family-2-adic-self-affine-tower-s579.md`, `04-computation/lrc_extremal_fractal_s579.py` (+.out); HYP-2115 (sleeves), HYP-2091 (+2 ladder), HYP-2063 (2q apex), HYP-2097 (flip-sets), THM-369.
