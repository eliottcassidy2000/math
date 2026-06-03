---
id: HYP-2114
status: SYNTHESIS -- additive-circuit lemmas (A circuit-free, B 3-term fold) + multiple/Phi share ONE residual (C'); fold=delta-clock; measure depressed only by 3-term. Lemma A verified/open, Lemma B reduces to C'.
source: opus-2026-06-03-S577
related: [HYP-2112, HYP-2102, THM-398, THM-369, HYP-2063]
---

# HYP-2114: the additive-circuit lemmas meet C' -- the fold is the delta-clock

**SETUP:** k=#nonzero speeds, delta=1/(k+1). Lemma A (randomness): circuit-free (no 3-term v_c=v_a+v_b) => M>=delta. Lemma B (structure): a 3-term relation is a fold t v_c=t v_a+t v_b.
**R1 fold=shield:** v_c=v_a+v_b => at pinch t=m/v_c runner c at integer (shield, blocks (a,b) pinch). Circuit-free min M ~0.21-0.25 (>> delta, margin grows +0.05->+0.13 k=4..10); 4-term-rich circuit-free shares the min.
**R2 hardness<=>3-term:** M-delta drops monotone with #3-term (k=6: +0.088/+0.045/+0.034/+0.011 for #3t=0/1/2/>=3); near-tight configs are 0% circuit-free, mean #3term 4.7/7.7/12.5 (k=6/8/10).
**R3 CONVERGENCE (the fold IS the delta-clock):** near-tight splits exactly by multiple-of-(k+1): #(no multiple)==#(delta-clock j/(k+1) witness) EXACTLY (18=18 k=6, 4=4 k=8, 0=0 k=12). The folds strip every pair-pinch except the straddle pair (a,n-a) summing to n=k+1 = the delta-clock; witness iff no multiple of k+1. Had-multiple configs loose off-clock (C'/Phi>0).
**R4 measure depressed ONLY by 3-term:** circuit-free min safe-measure ~0.08-0.10 (4-term-rich-cf IDENTICAL: 0.0806=0.0806); 3-term-rich(>=3) ~0.004-0.018 ->0. Relevant energy = 3-term count, not 4-term.
**2x2 SYNTHESIS:** [circuit-free vs 3-term-rich] x [no-mult vs mult-of-(k+1)]. Worry-set (M=delta) lives ONLY in [3-term-rich AND no-multiple]; tightness REQUIRES folds. Open residual = [multiple of (k+1)] = C' = ker Phi cap {n|v}, SAME under both slicings.
**UNIFICATION:** Lemma B "fold delivers delta" = delta-clock (THM-369, PROVEN for no-mult) + C' (OPEN for mult) -- not separate. Lemma A verified (mu>=~0.08, depressed only by 3-term), proof open (discrepancy from 3-term count). All three frameworks (additive, multiple, ECCP/Phi) bottom out at C'.
**See:** `07-reflections/lrc-additive-circuit-lemmas-meet-Cprime-the-fold-is-the-delta-clock-s577.md`, `04-computation/lrc_additive_circuit_AB_round{1,2,2b,3,4}_s577.py` (+.out); HYP-2112 (Phi), HYP-2102/THM-398 (C'), THM-369, HYP-2063 (2q apex).
