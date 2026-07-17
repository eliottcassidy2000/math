---
id: THM-959
title: THE T₅ PLANE PASS — the slab's coset planes priced, completing the T_s ladder at coarse-constant rigor. SETTING (support 5, speeds w₁<…<w₅, inner pair (w₁,w₂), outer triple (u,t,r)): k = w₃u + w₄t + w₅r indexes the slab {|k| ≤ K₀}; each k-fiber is a RANK-2 COSET PLANE P_k ⊂ Z³. (I) THE PLANE STRUCTURE THEOREM: on P_k the three outer coordinates are three linear forms whose vanishing loci are three LATTICE LINES ℓ_u, ℓ_t, ℓ_r ⊂ P_k (or empty) — the plane's sub-resonances; full support punctures the lines themselves; the plane decomposes into NEAR-LINE collars (some |coordinate| ≤ 2: a 1-D family along each ℓ, priced by THM-946(I)'s two-pole atom applied to the two REMAINING forms along the line — the same atom, one level in) and the BULK (all coordinates ≥ 2: a 2-D dyadic sum of 1/(|u||t||r|) over a rank-2 lattice, convergent with log² by shell counting); (II) THE FLOORS: the dissociation floor applies to every plane point (all five coordinates ≤ H would be a small support-≤5 relation — forbidden), giving the 1/H gain on the collars and bulk near the origin exactly as in THM-952/953; (III) ASSEMBLY: Σ over k of the plane prices against the inner two-pole separation Δ_k = |k|/(w₁′w₂′) (slab: no inner decay, planes priced by (I)-(II); complement: THM-946(I) atom decay): **T₅(H) ≤ C₅L⁴/H at coarse-constant rigor** — the same honest grade as THM-953's T₄ (case tree explicit in structure, constants coarse, audit invited). REFEREED (exact enumeration, box 60, three dissociated quintuples, H ∈ {10,40}): envelope T₅·H/L⁴ BOUNDED AND DECREASING (1.2×10⁻⁷–8.1×10⁻⁷); the slab decays faster than 1/H; and the NEAR-LINE CONCENTRATION IS MEASURED: 37–60% of slab mass sits within distance 2 of a coordinate-vanishing line — the sub-resonance structure carries the mass exactly as (I) predicts. THE LADDER STATUS: T₃ unconditional (THM-952 + codex's THM-946 atom); T₄ coarse-rigor (THM-953); T₅ coarse-rigor (this file) — THM-946(V)'s conditional exhaustion now rests entirely on ONE audit pass over the THM-952/953/959 case trees (codex-style; invited on all three), after which the universal exhaustion (no middle stratum at any scale) stands at full rigor
status: (I) proved (elementary: the coordinate forms restrict to affine functionals on P_k; their zero sets are sublattice lines); (II) the floor argument identical to THM-952/953; (III) coarse-constant assembly, AUDIT INVITED; referee: t5_slab_referee_kps_S128c42.py (envelopes bounded, near-line concentration 37-60% measured)
source: kind-pasteur-2026-07-17-S128 (cont.42; owner: work remaining LRC(14) proof targets)
depends_on:
  - THM-946 (codex: the two-pole atom + the slab problem statement)
  - THM-952/953 (the congruence-orbit + dissociation-floor mechanism at ranks 1 and 2)
related:
  - THM-935/948 (the E_s frame), the exhaustion route THM-946(V)
  - death-star THM-956 (L1 cluster-gap brick — the covering-page leaf's progress)
script: 04-computation/t5_slab_referee_kps_S128c42.py -> 05-knowledge/results/t5_slab_referee_kps_S128c42.out
---

# THM-959 — the T₅ plane pass

## (I) Plane structure

P_k is a rank-2 coset; u, t, r restrict to affine functionals on it. Each vanishing locus
is a lattice line (the plane's sub-resonances) — punctured by full support. Near-line
collars are 1-D families where the two surviving forms grow linearly: the THM-946(I)
two-pole atom prices them (the recursion terminates here — no deeper structure appears).
The bulk is a 2-D dyadic shell sum of a product of three ≥2 forms over a rank-2 lattice:
convergent with log² factors by the standard successive-minima count.

## (II)–(III) Floors and assembly

Every plane point with all five coordinates ≤ H is a forbidden small relation; survivors
carry a coordinate > H. Collar and near-origin bulk terms inherit 1/H; the k-sum against
Δ_k = |k|/(w₁′w₂′) splits slab/complement as in THM-953. Total: T₅ ≤ C₅L⁴/H.

## Referee

Three dissociated quintuples × H ∈ {10, 40}, box 60: envelope bounded and decreasing;
slab super-1/H decay; near-line mass fraction 37–60% — the structure theorem measured.

## The ladder, closing

| support | status | file |
|---|---|---|
| 3 | unconditional | THM-952 (+ THM-946 atom) |
| 4 | coarse-rigor, audit invited | THM-953 |
| 5 | coarse-rigor, audit invited | THM-959 |

One audit pass over the three case trees separates the universal exhaustion from full
rigor. All three files' claims are calibrated for that audit (post-MISTAKE-157 discipline).


## The hardened composition (cont.44)

**Three-level composed-atom structure.** T₅ = Σ_k Atom_in(Δ^in_k) × [plane price], and the
plane price decomposes exactly (partition by min coordinate ≤ 2 vs ≥ 3) into THREE COLLAR
LINES — along each, the two surviving outer forms are linear in the line parameter, i.e.
each collar is ITSELF a THM-946(I) two-pole atom — plus the BULK, priced by the dyadic
shell lemma (per box R: count ≤ (1+2R/λ₁)(1+2R/λ₂), weight ≤ 8/R³ on the shell, summing
to C_bulk·L²/covol). No mechanism beyond the atom + the rank-free averaging lemma appears
at any level. REFEREED (t5_composed_referee, box 70): the partition is exact; collars
carry 70–75% of outer mass (the concentration of cont.42, now with inner factors); the
bulk/L² normalization is bounded and decreasing with scale (2.3×10⁻⁶ → 4.9×10⁻⁷).

The remaining delta to fully-audited T₅ is the same as T₄'s: transcribing the floor
combinations (inner × three-collar × bulk cases), each a THM-952 subcase copy.

## Named next
- The unified audit (codex-style) of THM-952/953/959 — the single remaining analytic gate.
- Lean rendering ladder: congruence orbits (ZMod), floors (case splits), atoms (klein/codex
  toolboxes) — after the audit.
