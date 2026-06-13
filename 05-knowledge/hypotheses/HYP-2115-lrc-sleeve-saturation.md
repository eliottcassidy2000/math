---
id: HYP-2115
status: FORMALIZATION (sleeve saturation = covering foundation) + Ext1 saturation recursion (rigorous given Lemma G) + Ext2 alignment principle (verified, = Lemma A, open)
source: opus-2026-06-03-S578
related: [HYP-2112, HYP-2114, THM-398, THM-369]
---

# HYP-2115: sleeve saturation -- the covering foundation

**DEFS:** k speeds, delta=1/(k+1). SLEEVE Sigma_v={t:||v t||<delta}, measure 2*delta exactly (v arcs width 2delta/v at j/v). Total sleeve U=union; safe=complement; DEFICIT Delta(S)=mu(safe)=1-mu(U). measure-saturated: Delta=0 (worry-set, safe measure-0 but nonempty). point-saturated: safe=empty (M<delta, COUNTEREXAMPLE). LRC <=> no point-saturation. Total sleeve mass 2k/(k+1)~2>1, so saturation is purely about OVERLAP.
**EXT1 SATURATION RECURSION (rigorous given Lemma G):** order S=(v_1..v_k), G_0=circle, G_j=G_{j-1} minus Sigma_{v_j}; mu(G_j)=mu(G_{j-1})-Phi_{v_j}(G_{j-1}) where Phi_{v_j} is the Lemma-G ramp functional with current safe set G_{j-1}. Delta(S)=mu(G_k) = telescoping sum of per-sleeve ramps. Single-multiple Lemma G = the last step. Verified by saturation curves.
**SATURATION CURVE / critical AP:** mu_j (safe measure after j sleeves, level delta) for AP {1..k}: mu_1=(k-1)/(k+1) (first sleeve removes exactly 2delta), descends, measure-saturates EXACTLY at sleeve #k -- every sleeve necessary, zero slack = the extremal tight signature.
**EXT2 ALIGNMENT PRINCIPLE (verified, =Lemma A, open):** deficit Delta decreases monotone with #3-term (k=8: 0.140->0.103 mean, 0.111->0.015 min; k=10 similar). Saturation needs COHERENT overlap = additive 3-term alignment (the fold makes 3 sleeves overlap coherently); circuit-free (no alignment) => Delta>=c~0.11 (Lemma A). measure-saturation REQUIRES folds; worry-set = maximally-aligned (AP-like).
**PROGRAM IN SATURATION LANGUAGE:** no-mult-of-(k+1) => delta-clock points in NO sleeve (permanent point-survivors, THM-369); mult-of-(k+1) => its thin periodic sleeve (Vitali cover) can't point-saturate the survivors, Delta=Phi>0 (C'); circuit-free => mis-aligned, Delta>=c (Lemma A); worry-set => measure-saturated, survives on clock points. LRC = point-survivors never vanish.
**OPEN:** Ext2 (alignment/Lemma A) proof = a discrepancy/overlap-coherence bound; the critical-saturation extremal characterization (tight <=> sleeves measure-saturate with zero slack at the final sleeve).
**See:** `07-reflections/lrc-sleeve-saturation-the-covering-foundation-s578.md`, `04-computation/lrc_sleeve_saturation_s578.py` (+.out); HYP-2112 (Phi=mu), HYP-2114 (folds), THM-398 (C'), THM-369, S551o (Vitali).
