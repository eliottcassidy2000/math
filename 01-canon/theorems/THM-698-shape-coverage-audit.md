---
id: THM-698
title: The shape-coverage audit — THE SUPPLY'S DOMAIN IS UNIFORMLY TWO-SCALE: shapeOf v = (speeds <= 13, co-offsets Vmax - a of speeds > 13) by definition, so every family in THM527ACertificateSupply's domain is literally a (P, E, V)-two-scale family (middle speeds are just large co-offsets — no exotic shapes exist, the taxonomy-gap question DISSOLVES); and witnessG2(P,E) = slowmu(goodSet E ∩ safeSet P) = meas(G_P ∖ D(E)) — EXACTLY the dead-zone program's object, whose positivity (the supply's HYPOTHESIS) is what THM-690/691/692 proved unconditionally for every class; hence the supply's remaining content is PRECISELY the realization step [witnessG2 > 0 ⟹ finite certificates] on spread clusters at moderate V — per-family decidable, packed clusters already constructive (THM-693/696/697)
status: AUDIT COMPLETE (definitional reading of LRCFourteenSkeleton.shapeOf/witnessG2 + LRCReachDecitation.THM527ACertificateSupply; the identifications are definitional, not conjectural). CONSEQUENCES: (i) the "non-taxonomy shapes" gap named at S247/S248 is EMPTY — the domain is two-scale by construction; (ii) the dead-zone theorem (THM-690/691/692) = "witnessG2 > 0 for every class" — the supply's hypothesis is a THEOREM of this program; (iii) the supply's remaining content = spread-cluster realization at moderate V; packed clusters have kernel-pure constructive witnesses via the S248/S249 dichotomy (test-point route for q*-admitting P; first-window route for top-blocks, which force pmin >= 9).
source: klein-2026-07-11-S250 (HYP-5975; the critical open mapping named at S247/S248, executed)
depends_on:
  - LRCFourteenSkeleton (shapeOf, witnessG2 — the definitions audited)
  - THM-690/691/692 (the positivity of the supply's hypothesis)
  - THM-693/696/697 (the packed constructive witnesses)
related:
  - mac-mini LRCReachDecitation (THM527ACertificateSupply — the citation this audit scopes)
---

# THM-698 — the shape-coverage audit

## The three identifications (definitional)

1. **shapeOf v = (P, E)**: P = the |speeds| <= 13, E = Vmax − a for |speeds|
   a > 13. Every family in the supply's domain is a two-scale family
   P ∪ {Vmax − e : e ∈ E} — the "cluster" is everything above 13, middle
   speeds appearing as large co-offsets. THERE ARE NO EXOTIC SHAPES: the
   multi-scale/ray taxonomy (THM-688/694/695) refines the WITNESS
   CONSTRUCTIONS by co-offset structure, but the DOMAIN is uniformly
   two-scale.

2. **witnessG2(P, E) = slowmu(goodSet E ∩ safeSet P) = meas(G_P ∖ D(E))**:
   safeSet P is G_P; goodSet E is "the cluster phases leave a gap > 1/7",
   i.e. the 1/14-arcs do not cover, i.e. m_E > 0 — the complement of the
   dead zone. The supply's hypothesis witnessG2 > 0 is exactly
   meas(G_P ∖ D(E)) > 0 ⟺ μ∞(P,E) > 0.

3. **The dead-zone theorem = the hypothesis's unconditional truth**:
   THM-690/691/692 proved μ∞(P,E) > 0 for EVERY class — the supply's
   hypothesis is never vacuous, and the supply is demanded on every
   two-scale shape.

## The verdict

The supply's remaining content, precisely: **the realization step
[witnessG2 > 0 ⟹ an interval with thirteen strict certificates] for spread
clusters at moderate V.** Packed clusters (co-offsets ≲ 10·pmin) have
kernel-pure constructive witnesses through the S248/S249 dichotomy:

> for |P| ≤ 5, EITHER some q* ∈ [8,13] avoids P with |E| < q* (the
> test-point route, THM-693/696), OR P = [k+1, 13] exactly, forcing
> pmin ≥ 9 and an enormous first window (the 6/7-phase route, THM-697).

Spread clusters at moderate V remain per-family decidable (the banks'
territory), with the S239/S240 census showing all extremal shapes positive
with 2×+ margins. Bounded V belongs to THM-686. The other two citations
(≤7-arcs pigeonhole; THM-661's interior) are independent of the shape
question.

## Files

The audit is definitional — the readings are LRCFourteenSkeleton.lean
(shapeOf ~line 148, witnessG2 ~line 223) and
LRCReachDecitation.THM527ACertificateSupply. The packed dichotomy's Lean:
LRCTestDataSupply.lean (S248) + LRCFirstWindowWitness.lean (S249).
