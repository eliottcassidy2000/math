---
id: THM-883
title: THE BODY-WEIGHTED RAMANUJAN EXPANSION (THM-873 for ARBITRARY bodies) + THE WITNESS-TORSION LAW OF PEELS — (A) for any finite body E and any λ, the bad set is ONE effective arc per primitive fraction a/l, l ∈ D(E) (divisor closure), of half-width λ/m_E(l) with m_E(l) = min{v ∈ E : l | v}; #arcs = Σ_{l∈D(E)} φ(l). (B) sweep-exact transform: −ĝ_E(h) = Σ_{l∈D(E)} c_l(h)·sin(2πhλ/m_E(l))/(πh) − Σ_r T(J_r) with EXPLICIT prefix-max shadow intervals J_r (cascades/swallowed arcs cancel exactly; interval cores recover THM-873's closed-gap corrections). (C) PEEL FORM: ĝ_{P∪{v}}(h) = ĝ_P(h) − Σ_{new arcs A} T(A ∩ G_P); grouping whole survivors by level gives RESTRICTED RAMANUJAN SUMS c_l^S(h) = Σ_{a∈S_l} e(−ha/l) over the SURVIVING torsion S_l ⊆ (Z/l)^*; ≤ 2r_P straddle lenses each |T| ≤ 2λ/v. (D) WITNESS-TORSION LAW (census, exact): for the covering extremals the surviving new torsion is EXACTLY the core's witness classes, FULL — deep well {1..12}+182: S₁₃ = all 12 primitive 13ths, levels 14/26/91/182 DEAD (level 14: the four abutment points 3,5,9,11/14 swallowed = THM-882's isolated flat points; ±1/14 straddle); {1..11,13}+84: S₁₂ full, else dead. So ĝ_deepwell = ĝ_core − c₁₃(h)·sin(2πh/2548)/(πh) − (two half-lenses at ±1/14): a THREE-TERM exact spectrum for the covering-min extremal. (E) disc₁₃({1..12}; 1/14) = 104999803919/6363107150400 ≈ 1.650134e-2 EXACT (autocorrelation ℚ-arithmetic + monotone h-sum referee)
status: PROVED ((A) nesting argument; (B) sweep identity + level grouping, elementary; (C) measure additivity + classification; all machine-refereed: 7 bodies × set-identity EXACT + transform identity ≤ 3e-15 over h ≤ 40, peel identities ≤ 5e-17, censuses exact, disc exact-ℚ with h-sums converging from below)
source: boxeph-2026-07-16-S23; executes kind-pasteur cont.28's named critical task (a) "body-weighted Ramanujan spectra — the bridge from the closed interval-core case to open stratum [B]"
depends_on:
  - THM-873   # the interval-core case (recovered as the m_E(l) = l specialization)
related: [THM-731/732 (disc certificates — (E) supplies exact-ℚ disc), THM-826 (effective arcs), THM-879 (Möbius-sinc closure at k=13 — (B) is its arbitrary-body front end), THM-882 (the abutment points reappear as isolated flat points), klein-S287 (isolation inverts difficulty — quantified by the 1/m_E(l) weights), klein-S296 (witnesses = torsion points — (D) is its peel-spectral face)
script: 04-computation/lrc14_body_weighted_ramanujan_boxeph_S23.py -> 05-knowledge/results/lrc14_body_weighted_ramanujan_boxeph_S23.out
---

# THM-883 — body-weighted Ramanujan expansion; the witness-torsion law

**(A) Concentric collapse.** For a finite body E ⊂ Z⁺ and λ ∈ (0, 1/2), the bad set
B(E; λ) = {t : ∃v ∈ E, ‖vt‖ < λ} is the union, over levels l in the divisor closure
D(E) and primitive residues a ∈ (Z/l)^*, of the arcs A_l(a) = (a/l − W, a/l + W),
W = W_E(l) = λ/m_E(l), m_E(l) = min{v ∈ E : l | v}. *Proof:* the speed-v arcs at a fixed
reduced center a/l are concentric with half-widths λ/v, v ∈ E ∩ lZ; their union is the
widest. ∎ (#arcs = Σ_{l∈D(E)} φ(l); for E = {v} this is Σ_{d|v} φ(d) = v.)

**(B) Sweep-exact expansion.** Rotate so a good point is the origin, sort arcs by left
edge, let M_r = max_{s<r}(right edges) and J_r = (left_r, min(right_r, M_r)] (empty when
M_r ≤ left_r). Then 1_B = Σ_r 1_{I_r} − Σ_r 1_{J_r}, so for h ≠ 0

> −ĝ_E(h) = b̂_E(h) = Σ_{l ∈ D(E)} c_l(h)·sin(2πhλ/m_E(l))/(πh) − Σ_r T(J_r),

with c_l(h) the Ramanujan sums and T(interval) its transform. The first sum is the
**body-weighted Ramanujan term** — THM-873's expansion with the level weight deformed
from λ/l to λ/m_E(l); the shadows are the exact correction in EVERY regime: pairwise
overlaps give THM-873's closed-gap lenses; a fully swallowed arc has J_r = I_r and
cancels exactly (cascades are handled, which is where the interval-core proof's
"arcs interact only within their own gap" fails for covering bodies — e.g. the deep
well swallows 156 of its 220 arcs).

**(C) Peel form.** For E = P ∪ {v} (core widths unchanged, i.e. m_E = m_P on D(P)):
B_E = B_P ⊔ (new arcs ∩ G_P), hence exactly

> ĝ_{P∪{v}}(h) = ĝ_P(h) − Σ_{new arcs A} T(A ∩ G_P).

New arcs (levels l ∈ D(v)∖D(P), all of width λ/v) classify as: **swallowed** (⊂ B_P,
contribute 0), **whole survivors** (⊂ G_P, contribute their full thin-arc transform),
**straddlers** (≤ 2r_P of them, one per crossed good-interval endpoint, each
|T(lens)| ≤ 2λ/v). Grouping whole survivors by level gives the **restricted Ramanujan
sums** c_l^S(h) = Σ_{a ∈ S_l} e(−2πiha/l), S_l = {a : A_l(a) ⊂ G_P}:

> ĝ_{P∪{v}}(h) = ĝ_P(h) − Σ_{l new} c_l^{S_l}(h)·sin(2πhλ/v)/(πh) − Σ_{straddle lenses}.

**(D) The witness-torsion law (exact census at λ = 1/14).**

| peel | new levels | surviving torsion |
|---|---|---|
| {1..12} + 182 | 13, 14, 26, 91, 182 | S₁₃ = ALL 12 primitive 13ths; others DEAD |
| {1..11,13} + 84 | 12, 14, 21, 28, 42, 84 | S₁₂ = ALL 4 primitive 12ths; others DEAD |
| {1..10,13,22} + 84 | 12, 14, 21, 28, 42, 84 | S₁₂ full; 2 straddlers at ±1/14; else DEAD |

The surviving new torsion is exactly the core's **witness set** (the primitive classes
of the core's missing level — klein-S296's clean witnesses), and for these extremal
peels it survives as a FULL class. Level-14 anatomy at the deep well: 3/14, 5/14, 9/14,
11/14 are swallowed as ABUTMENT points of two core arcs — precisely THM-882's isolated
flat points — while 1/14, 13/14 straddle ∂G_P (the good set's extreme endpoints). So the
covering-min extremal's spectrum is THREE terms:

> ĝ_{deep well}(h) = ĝ_{{1..12}}(h) − c₁₃(h)·sin(2πh/2548)/(πh) − (two half-lenses at ±1/14).

**(E) Exact disc.** disc₁₃({1..12}; 1/14) = Σ_{h≠0}|ĝ(13h)|² =
**104999803919/6363107150400** = 1.650134e-2, computed in exact ℚ via the v-grid
autocorrelation (1/13)Σ_j λ(G ∩ (G − j/13)) − λ(G)², with the spectral h-sums
converging to it from below (1.649774/1.650045/1.650112e-2 at |h| ≤ 2000/8000/32000).
The fleet-quoted 1.648e-2 (opus-S271 via kps THM-873) agrees only to ~1e-4 — the exact
rational above supersedes it for certificate use.

## Consequences for stratum [B]

1. **Quantitative isolation.** Every far level is 1/m_E(l)-suppressed in the weight —
   klein-S287's "isolation inverts difficulty" with constants: the far element protects
   its torsion by arcs of width λ/v only (deep well: 1/2548).
2. **Where the cancellation is NOT.** For the extremal peels the restricted sums are
   FULL primitive classes (c^S = c_l, max |c_l^S| = φ(l) at l | h): no equidistribution
   gain from the restriction at the surviving level. The disc gain must come from the
   sinc aliasing in h (the Möbius-sinc mechanics of THM-879), not from partial torsion.
   Partial-torsion cancellation becomes available only for NON-extremal cores whose
   witness classes are incomplete — a concrete dichotomy for the remaining [B] work.
3. **Certificate-grade exact disc.** The autocorrelation route computes disc_v in exact
   ℚ for any body (Part D pattern) — THM-731/732 certificates can drop float verifieds.

## Evidence log
- [x] (A) set identity exact on 7 bodies incl. cascade/wild regimes; #arcs = Σφ(l)
- [x] (B) transform identity ≤ 3.02e-15 (h ≤ 40, all bodies); shadows censused
- [x] (C) peel identity ≤ 4.8e-17 on three peels; censuses exact
- [x] (E) exact-ℚ disc with two-sided referee
- [ ] extend the witness-torsion census across the multikiller bank (mechanical)
- [ ] the non-extremal partial-torsion cancellation (consequence 2) — the named opening
