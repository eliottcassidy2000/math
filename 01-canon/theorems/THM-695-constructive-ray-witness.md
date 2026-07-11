---
id: THM-695
title: The constructive ray witness — ray-coupled families made explicit: for S = P ∪ {V_mid − f : f ∈ E_mid} ∪ {V − e : e ∈ E_top} with V_mid = (a_r·V − c)/b, the witness is THM-693's TIME at the b-scaled modulus (q, w) = (q*·b·V, b·w₆₉₃): small and top-cluster speeds INHERIT their bands by pure scaling ((bx) % (bQ) = b(x % Q)), and the ray-mid speeds get residue ρ'_f·V − δ(c + bf) with ρ'_f = (a_r(aV + δ) − a(c + bf)) % (q*b), strictly in-band under a per-f DECIDABLE digit hypothesis + threshold 14q*(c + bf) < V. FORMALIZED SORRY-FREE (LRCRayWitness.lean, kernel-pure); with THM-693/694 the constructive certificate supply covers two-scale, multi-scale, AND ray-coupled residual shapes — feeding mac-mini's cont.27 de-citation (THM-527-A narrowed to certificate supply)
status: PROVED AND FORMALIZED (scaled_band, ray_speed_band, rayFamily_strictlyLive, rayFamily_strictWitness, demoRay_strictWitness — all [propext, Classical.choice, Quot.sound]; one staged-nlinarith fix from first compile; demo P = {1,2,3}, ray (1,2,0) at V_mid = 5000, E_mid = {0,1}, E_top = {0,1,2,3,5,8,13,21} at V = 10000 witnessed THROUGH the theorem at (q, w) = (260000, 20002) = (2·130000, 2·10001); Python cross-check all 13 strictly in-band, mid digits ρ' = 17, 15).
source: klein-2026-07-11-S247 (HYP-5960; the named ray follow-up of THM-694/S246)
depends_on:
  - THM-693 / LRCTwoScaleWitness (the inherited small/top bands and the witness time)
  - kps LRCStrictRuler
related:
  - THM-688/S246 (the ray taxonomy + evaluator this makes constructive)
  - mac-mini cont.26/27 (lrc14_from_citations_only; the certificate-supply de-citation this supply feeds)
---

# THM-695 — the constructive ray witness

## Statement

Let S contain small speeds P (test residues nonzero at (q*, a)), a top
cluster {V − e : e ∈ E_top} (coarse missed class c₁ mod q*), and a RAY
cluster {V_mid − f : f ∈ E_mid} with b·(V_mid − f) = a_r·V − c − b·f. Put
δ = (c₁ − aV) mod q*, w₆₉₃ = aV + δ. Then at

> **(q, w) = (q*·b·V, b·w₆₉₃)** — THM-693's witness time, b-scaled modulus —

every speed is strictly in-band provided V > 2184, 168e < V per top
co-offset, and per ray index f: the DECIDABLE digit hypothesis
q*b < 14·ρ'_f < 13·q*b with ρ'_f = (a_r(aV + δ) − a(c + bf)) mod (q*b),
plus 14q*(c + bf) < V. The residues: small/top speeds inherit THM-693's
values scaled by b (since (bx) mod (bQ) = b·(x mod Q)); the ray-mid speed u
gets (u·w) mod q = ρ'_f·V − δ(c + bf). ∎ (The Lean file is the proof; the
key algebra: u·b·(aV + δ) = (a_rV − c − bf)(aV + δ) = ρ'V − δ(c+bf) + qm·V
after the ediv/emod split — one linear_combination.)

The digit hypothesis is a condition on V mod q*b — computable per family;
its wholesale satisfiability over V-residues is the (a,b)-level dead-zone
game whose continuum version S246's census answered (all positive).

## What this closes

The constructive certificate supply now covers **two-scale (THM-693),
multi-scale (THM-694), and ray-coupled (THIS)** residual shapes — every
witness an explicit mixed-radix number, every proof kernel-pure. This is
exactly the supply that mac-mini's cont.27 de-citation consumes (THM-527-A
narrowed to certificate existence): the citations' arithmetic content is
now constructively realized on all the taxonomy's class shapes.

## Formalization & files

`04-computation/lean/TournamentH7/TournamentH7/LRCRayWitness.lean`
(kernel-pure, root-wired): scaled_band, ray_speed_band,
rayFamily_strictlyLive/strictWitness, demoRay_strictWitness. Demo
cross-check: `05-knowledge/results/lrc14_ray_witness_demo_klein_S247.out`.
