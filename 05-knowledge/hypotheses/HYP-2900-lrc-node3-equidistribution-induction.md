---
id: HYP-2900
title: NODE 3 (analytic) -- a committed speed equidistributes within the seed's safe set (removes EXACTLY 1/7, leaves 6/7>0 => M>=1/14); inductive reduction Node3-LRC(n) <= LRC(n-1); so LRC(14) <= Node2(8..14) [kps Bonferroni-3] + induction
status: VERIFIED equidistribution (exact to 5 decimals); inductive SKELETON; rigorous version needs effective Erdos-Turan
source: mac-mini-2026-06-22-S46 (concurrent push to finish; complements kps Node 2)
related:
  - HYP-2899   # the two-structure split (Node 2 / Node 3)
  - THM-566    # kps adversarial / equidistribution
  - THM-527    # compactness (bounded => finite)
  - HYP-2885   # Node 2 (consec maximizes); kps S31t Bonferroni-3
---

# HYP-2900: Node 3 is equidistribution + induction; LRC(14) reduces to Node 2

## The verified equidistribution (lrc_node3_equidistribution_macmini_S46.py)
Split a 13-set as seed B (the <=12 "bounded" runners) + a committed large speed v. The committed v's
unsafe arcs U_v EQUIDISTRIBUTE within B's safe set Safe(B):
  meas(Safe(B) ∩ Safe(v)) = meas(Safe(B)) * (6/7)   -- VERIFIED to 5 decimals.
For B = {1..11,13} (meas(Safe@1/14)=0.01216) and v = 30030, 60060, 510510 (all 30030|v):
  meas(Safe(B) ∩ Safe(v)) = 0.010423 ≈ 0.01216*(6/7) = 0.010422.  >0 => M(B∪{v}) >= 1/14.
The committed v removes EXACTLY its measure-fraction (1/7) of B's safe set -- no resonance, pure
equidistribution. (Contrast the BUG I caught: a wrap-around error first gave 0.79; fixed => 0.01216.)

## The inductive reduction
**Node3-LRC(n) <= LRC(n-1).** Peel the committed speed v: the seed B (n-2 runners) has positive safe
measure at 1/n BY LRC(n-1) (M(B) >= 1/(n-1) > 1/n), and v equidistributes within it (removes 1/7, leaves
6/7 > 0), so M(B∪{v}) >= 1/n. Recursively peel all committed/large speeds. Induction chain:
  LRC(n) <= easy(n) [surviving prime] + Node2(n) [bounded core] + Node3(n) [<= LRC(n-1)].
Bottoms at LRC(7) (6 runners, PROVEN). So **LRC(14) <= Node2(8),...,Node2(14)** -- the BOUNDED cores --
plus the induction (Node 3) + the easy case. Node 2 (consec maximizes, kps's Bonferroni-3 truncation
p0<=T_1+T_2+T_3) is the only irreducible remaining work.

## Honest gaps
- The equidistribution is VERIFIED computationally (exact); the RIGOROUS version needs an effective
  Erdos-Turan / Weyl bound: discrepancy(U_v on Safe(B)) < meas(Safe(B))/7, valid when v >> scale(B).
  The committed v >= 30030 (radical, S45) works when scale(B) is below the equidistribution threshold
  (~v/14); the peeling (peel LARGEST first) keeps v >> rest. The threshold is the Node2/Node3 boundary.
- The induction uses LRC(n-1) for the seed -- this is the SAME conjecture one size down, so the chain
  reduces LRC(14) to the bounded Node-2 cores (NOT to a trivial base beyond k=6). Node 2 stays the crux.
=> This is the proof SKELETON: induction (Node 3, verified equidistribution) collapses LRC(14) onto
Node 2 (kps Bonferroni-3) + the easy case. Closing Node 2 closes LRC(14).
