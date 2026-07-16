---
id: THM-903
title: THE REFLECTION FRAME FOR THE RESONANT CROSS-SECTION — Error_t(a) is EXACTLY invariant under x ↦ −x (verified 2e-16 at residues 1 and 6; five-line substitution proof: missed-sets and g-kernels twist coherently by s ↦ 6−s), so NO residue transfer 6 ↔ 1 exists (the S116 proposal corrected); the stationary offset PINS sector 0 and breaks sector-measure reflection except on pin-clean pairs (A₁₃ = A₃₅ exact; A₂₄ reflection-FIXED); the fixed pair-sectors are exactly {s, 6−s} — the locus where THM-891's negative-residue-6 mass concentrates: the hard mass sits on the involution's fixed locus (the Rédei/locker/SC template), and certificates may be built on the half-domain with symmetrized kernels
status: invariance PROVED (substitution; pin bookkeeping via the joint R/g twist) + verified numerically (2e-16 at a = 1, 6, 13, 20, t = 501, core {0,1,2,3,4,6}); sector-measure laws exact on the fine grid; transfer verdict = corrected refutation of my own S116 proposal
source: mac-mini-2026-07-16-S117 (owner: investigate equivariance concepts; work toward finishing LRC(14))
depends_on: [codex THM-891 (definitions; the open negative side), THM-727 (centered sections)]
related: [THM-466(iii) (OCF reversal equivariance), S112 negation-reversal law, klein ±pairs/2I descent, kps coprime-v invariance, THM-819 inverse pairs — the equivariance ledger]
script: 04-computation/fixed_locus_residue6_macmini_S117.py -> 05-knowledge/results/fixed_locus_residue6_macmini_S117.out
---

# THM-903 — the reflection frame

**(i) Invariance.** Under x ↦ −x every nonstationary point's sector reflects s ↦ 6−s, and
g_s(−y) = g_{6−s}(y) a.e.; the missed-set and kernel indices twist coherently, so
Error_t(a) = Σ_s ∫_{R_s} g_s(atx) dx is EXACTLY invariant (machine: |δ| ≤ 2·10⁻¹⁶ at
a = 1, 6, 13, 20). Consequence: the reflection acts WITHIN each residue class of a — the
naive neg-6 ↔ pos-1 transfer (my S116 letter) does NOT exist. Corrected in place.

**(ii) The pin.** The stationary offset occupies sector 0 at BOTH x and −x (0 does not
reflect to 6), so individual sector measures do NOT all reflect: A₁₂ ↛ A₄₅ (0.0238 vs 0);
the pin-clean identities hold exactly (A₁₃ = A₃₅ = 0.0119; A₂₄ = A₂₄ fixed). The
reflection-fixed pair-sectors are exactly {s, c} with s + c = 6: **{1,5} and {2,4}** —
THM-891's concentration set for the negative residue-6 mass. The hard mass sits on the
involution's fixed locus, as everywhere in this project (SC classes under T ↦ T^op,
locker squares under d ↦ n/d, tight corners of the corridor, a² ≡ 1 witnesses).

**(iii) Actionable for the residue-6 closure** (@codex/death-star): (a) build certificates
on the half-domain x ∈ [0, ½] with symmetrized kernels — the odd modes integrate to zero
by (i), halving the search space; (b) on the fixed sectors the pair data is intrinsically
symmetric, so "pair marginals insufficient" is structural: the next certificate must use
reflection-EVEN cubic (triple) data — matching death-star's relation-type refinement and
opus's triple-beat T1/T2 exactly; the three routes are one route.
