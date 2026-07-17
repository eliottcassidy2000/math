---
id: THM-892
title: THE INVARIANT MEAN ⟨Q_s⟩ — CLOSED FORM (the frame-average of the second moment over all dilation frames; the reflection's named statement, proved) — (K) k̂_P(n) = −csc²(πn/P)/(2P²) EXACTLY (upgrading THM-886(III)/888 from bound to IDENTITY), so Q_s(w) = (π²/P²)Σ_{n≠0}csc²(πn/P)|S(nw)|² exactly; (C\*) Σ_{u∈(Z/q)\*}csc²(πu/q) = J₂(q)/3 (Jordan totient; Möbius + MI0); (P) Σ_{m∈dZ_P}|S(m)|² = (P/d)·N(P/d) with N(h) = Σ_{k,k'}ε_kε_k'[p_k ≡ p_k' mod h] the SIGNED ENDPOINT COINCIDENCE SPECTRUM; hence ⟨Q_s⟩ = (π²/P²)·Σ_{g|P,g<P}(J₂(P/g)/3)·φ(P/g)⁻¹·Σ_{g|d|P}μ(d/g)(P/d)N(P/d) — π² times an explicit rational functional of the coincidence spectrum on the divisor lattice. REFEREED EXACT (rel. err ≤ 5e-14 vs direct φ(P)-frame averages, four cluster types). EMPIRICAL LAW: ⟨Q_s⟩/M = 3.18/3.00/2.83/3.42 across family/two-owner/balanced/family-120 — THE INVARIANT MEAN IS UNIFORMLY ~3M while single-frame values fluctuate [1.2, 3300]: the objective part of Q_s IS the O(M) law. COROLLARY (generic-frame theorem, Markov): for every cluster, at least (1−1/λ) of all coprime dilation frames have Q_s(w) ≤ λ·⟨Q_s⟩ ≈ 3λM — the measure-theoretic complement to the Diophantine classifier
status: PROVED + LEAN-PARTIAL (klein-S318: (K)'s rational heart tent_second_difference [general P] and (C*)'s Möbius/Jordan layer J2_divisor_sum are kernel-checked in TournamentH7.Thm892Shadows/.Thm892Jordan; (P) open) ((K) discrete second difference, 3 lines; (C\*) Möbius over MI0; (P) subgroup orthogonality; the assembly is bookkeeping) + REFEREED EXACT on four clusters (every lemma separately + the final identity at 1e-14–1e-16)
source: boxeph-2026-07-16-S29 (owner directive: prove the unit-average closed form; the-intersubjective-object reflection's named statement)
depends_on: [THM-880/881 (frame), THM-886/887/888 (the spectrum machinery; (K) upgrades their kernel bound to an identity), the reflection 07-reflections/the-intersubjective-object-boxeph-S28.md]
related: [THM-884(E)/THM-879(i) (the N(h) are signed versions of the v-grid coincidence moments — one object family), klein cont.8 (relation-as-object: composition defect = 3c₃ — the tournament face of the same owner prompt)]
script: 04-computation/lrc14_invariant_mean_boxeph_S29.py -> 05-knowledge/results/lrc14_invariant_mean_boxeph_S29.out
---

# THM-892 — the invariant mean of the second moment

**Namespace note.** Renumbered from `THM-891` after rebase exposed the earlier-pushed
`THM-891` cross-section-cancellation claim (`cb88cc992`).  The mathematics is unchanged.

**(K) The kernel is exactly the cosecant.** For n ≠ 0: k̂_P(n) = −csc²(πn/P)/(2P²).
*Proof.* K(j/P) = (j/P)(1−j/P) has discrete second difference Δ²K = −2/P² + (2/P)δ_0
on Z_P; taking DFT, −4sin²(πn/P)·k̂(n) = 2/P². ∎ (THM-886(III)'s "bound" and THM-888's
collapse masses are therefore IDENTITIES; the measured ratio 2.000 was equality.)
Hence **Q_s(w) = (π²/P²)·Σ_{n≠0} csc²(πn/P)·|S(nw)|² exactly** (refereed vs the
THM-880 bilinear form).

**(C\*) Unit cosecant sums are Jordan totients.** Σ_{u∈(Z/q)\*} csc²(πu/q) = J₂(q)/3.
*Proof.* Möbius over divisors reduces to full sums; Σ_{v mod q', v≠0} csc²(πv/q') =
(q'²−1)/3 (MI0); Σ_{d|q}μ(d)((q/d)²−1) = J₂(q) − [q=1]. ∎

**(P) Subgroup Parseval = coincidence counts.** Σ_{m ∈ dZ_P} |S(m)|² = (P/d)·N(P/d),
N(h) := Σ_{k,k'} ε_kε_k'[p_k ≡ p_k' mod h] — with N(1) = (Σε)² = 0 and N(P) = M.

**Theorem (the invariant mean).** Averaging over the frame group (Z/P)\*:
> ⟨Q_s⟩ = (π²/P²) · Σ_{g|P, g<P} (J₂(P/g)/3) · φ(P/g)⁻¹ · Σ_{g|d|P} μ(d/g)·(P/d)·N(P/d).
*Proof.* The orbit of n under units is its gcd-class; each class average of |S|² is the
Möbius resolution of the subgroup sums (P); the class k̂-masses are (C\*) via (K). ∎
Refereed exact against the direct average over all φ(P) frames on four cluster types
(rel. err ≤ 5·10⁻¹⁴).

## Readings

1. **The objective part of Q_s is the O(M) law.** ⟨Q_s⟩/M = 3.18, 3.00, 2.83, 3.42 on
   the battery — uniformly ~3M — while individual frames fluctuate across three orders
   of magnitude ([1.2, 3300] at {1..6,120}). The week's arc lands exactly where the
   relativity frame said it would: the frame-average is the invariant "acceleration",
   and it obeys the clean law the uniform-frame conjectures kept failing to be.
2. **Generic-frame corollary (Markov, free):** for every cluster and every λ > 1, at
   least (1 − 1/λ)·φ(P) of the coprime dilation frames satisfy Q_s(w) ≤ λ⟨Q_s⟩ ≈ 3λM.
   The Diophantine classifier says WHICH frames are bad; the mean says HOW FEW.
3. **The arithmetic content is the coincidence spectrum.** ⟨Q_s⟩ is a divisor-lattice
   functional of N(h) — the signed congruence structure of the endpoint multiset. On
   the battery N(7) dominates (218/180/88/800): the seven-section skeleton is the
   leading torsion of the mean. The N(h) are signed cousins of the v-grid coincidence
   moments (THM-879(i), THM-884(E)) — the disc machinery and the frame-average are one
   object family: torsion sieves of endpoint sets.
4. **The isotypic program (named next):** the mean is the trivial character of the
   frame action. The full character expansion Σ_w χ̄(w)Q_s(w) should have a twisted
   closed form (J₂ → Jacobi-type sums, N(h) → χ-twisted coincidences) — "the character
   theory of resonance"; the comb diagonal (THM-888(A), already factoring through
   w mod 7) is its first nontrivial isotype.

## Evidence log
- [x] (K)/(Q)/(C\*)/(P) each refereed; final identity exact at 1e-14–1e-16, four clusters
- [x] fluctuation ranges + coincidence spectra recorded
- [ ] the character expansion (isotypic decomposition of Q_s under (Z/P)\*)
- [x] ⟨Q_s⟩ ≈ 3M: THE CONSTANT IS π²/3 = 3.2899… — and the deviation identity is now
      PROVED (S30 addendum below).

## ADDENDUM (boxeph-S30): THE DEVIATION IDENTITY (N2, proved)

Using the Jordan identity Σ_{q|P} J₂(q) = P²:
> **⟨Q_s⟩ = (π²/3)·M·(1 + 1/P) + (π²/(3P²))·Σ_{q|P, q>1} J₂(q)·δ_{P/q}**,
> δ_g := A_g − PM/(P−1) (the gcd-class imbalance of the spectrum; Σ_g φ(P/g)δ_g = 0).
*Proof.* Split A_g = Ā + δ_g in THM-892's closed form; the Ā-part telescopes by the
Jordan identity to (π²/3)M(P+1)/P exactly. ∎
Referee: EXACT on all four battery clusters (identity error ≤ 4·10⁻¹⁶); measured
deviations −5.06/−13.56/−16.56/+8.64 against universal terms 138.5/151.5/118.5/224.0 —
the universal term carries 88–116% of the mean, with the arithmetic correction signed
either way. **The invariant mean is a TWO-TERM LAW: a universal (π²/3)M(1+1/P), exact,
cluster-independent — plus a J₂-weighted coincidence-spectrum imbalance.** The
uniform-frame program's correct target, in final form.
