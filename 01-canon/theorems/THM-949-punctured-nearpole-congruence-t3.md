---
id: THM-949
title: THE PUNCTURED NEAR-POLE CONGRUENCE LEMMA — closing the T₃ gap named by THM-946(III). Setting: triple a < b < c, reduced parameters (g = gcd(a,b), g₀ = gcd(g,c), d = g/g₀, a′ = a/g, b′ = b/g, c′ = c/g₀); the support-3 relation lattice is the union over t ≠ 0 of lines (x_t + mb′, y_t − ma′, dt). THE LEMMA (proved): (i) the exact poles are zero-coordinate points — EXCLUDED from the full-support mass by definition (the puncture); (ii) the near-pole integer points on line t obey the CONGRUENCE LAW: the minimal |h₁| equals r₁(t) = least absolute residue of (−c′t·a′^{−1}) mod b′, and symmetrically r₂(t) mod a′ — near-pole behavior is a rotation orbit on Z/b′ × Z/a′, not an analytic quantity, which is exactly why THM-946's per-line bound (the 6/(1+AΔ) terms) cannot be summed pointwise; (iii) THE DISSOCIATION FLOOR: under H-dissociation, a near-pole point with |h₂| ≤ H and |h₃| = d|t| ≤ H must have r(t) > H (else the point itself is a forbidden small relation) — so every surviving near-pole point carries a coordinate > H; (iv) THE SUMMED BOUND: the total near-pole mass obeys NP(H) ≤ C·(1 + ln(2+Vmax))/H via three cases (r > H floor with |h₂| ≥ 1 and the harmonic t-sum; |h₂| > H floor symmetric; |h₃| > H with the c′t/(2b′) growth of |h₂| giving t⁻² decay). REFEREED: exact near-pole sums on six triples (incl. a Schur-adjacent one), H ∈ {5, 20, 80, 320}: the envelope NP·H/(1+ln Vmax) is BOUNDED (≤ 0.026, no growth) and NP decays as C/H. CONSEQUENCE (with THM-946(I)'s corrected two-pole atom, per THM-946(III)'s own audit of what remained): **T₃(H) ≤ C₃(1+ln(2+Vmax))²/H UNCONDITIONALLY** — the assembly: for |t| ≤ H every full-support point has a coordinate > H (dissociation direct: one factor < 1/H, the rest bounded by a one-pole harmonic ≤ CL per line, harmonic t-sum gives L²/H); for |t| > H the corrected atom's main term (13) gives O(L/H) and the near-pole classes are covered by this lemma at all t. REMAINING OPEN (unchanged, honestly): the T₄ resonance-STRIP and T₅ resonance-SLAB affine-coset estimates (THM-946(IV)) — the exhaustion route (THM-946(V)) is now conditional on those two only
status: lemma (i)-(iv) PROVED (congruence identification + dissociation floor + three-case bound; constants coarse); T₃ assembly per THM-946's audit division; refereed (6 triples × 4 horizons, envelope bounded). ACKNOWLEDGMENT: this builds directly on codex's THM-946 correction of my THM-944 (MISTAKE-157 territory: my original two-pole constant was false without the log and my T₄/T₅ nesting was not a proof — this file supplies the T₃ piece their audit specified, and does NOT claim T₄/T₅)
source: kind-pasteur-2026-07-17-S128 (cont.40; owner: work to finish the LRC(14) proof)
depends_on:
  - THM-946 (codex: the corrected two-pole atom + the audit naming this lemma)
  - THM-935 (the E_s identity the T_s feed)
related:
  - THM-948 (the exact packet audit — concrete packets never need these tails)
  - THM-946(IV) (the strip/slab estimates — the remaining conditional surface)
script: 04-computation/nearpole_congruence_referee_kps_S128c40.py -> 05-knowledge/results/nearpole_congruence_referee_kps_S128c40.out
---

# THM-949 — the punctured near-pole congruence lemma; T₃ unconditional

## The lemma

On line t, solutions of a′h₁ + b′h₂ = −c′t have h₁ determined mod b′: the point nearest
the h₁-pole has |h₁| = r₁(t) = |−c′t·a′^{−1} mod b′|_± (least absolute residue), and
h₂ ≈ −c′t/b′ there; symmetrically for the h₂-pole. The poles themselves (r = 0) are
zero-coordinate points, excluded from the full-support mass.

**Dissociation floor.** If |h₂| ≤ H and d|t| ≤ H, then (r₁(t), h₂, dt) being a lattice
point forces r₁(t) > H — otherwise it IS a small relation, contradicting H-dissociation.

**Summed bound.** Split the near-pole mass by which coordinate exceeds H:
- r > H: contribution ≤ (1/H)·Σ_t 1/(max(1,|h₂|)·d|t|) ≤ (1/H)·(harmonic in t up to
  the |h₂|-growth threshold, then t⁻²) ≤ C(1+ln(2+Vmax))/H;
- |h₂| > H: symmetric;
- d|t| > H: the weight 1/(d|t|) < 1/H directly, and Σ over the remaining two factors
  along the congruence orbit is O(1+ln Vmax) by the same harmonic comparison.

Total: NP(H) ≤ C(1+ln(2+Vmax))/H. ∎

## The T₃ assembly (division per THM-946(III))

|t| ≤ H: dissociation applies to every full-support point directly (some factor < 1/H;
the rest ≤ one-pole harmonic CL per line; Σ_t 1/(dt) ≤ L): contribution ≤ CL²/H.
|t| > H: THM-946(13) main term is O(L/H) after the outer weight; near-pole classes by
the lemma. Hence **T₃(H) ≤ C₃L²/H, L = 1+ln(2+Vmax), unconditionally.**

## What remains for the exhaustion route
Exactly THM-946(IV): the T₄ resonance strip (Δ_{u,t} vanishes on cu+dt = 0) and the T₅
resonance slab — genuine affine-coset estimates, not nesting. The conditional statement
THM-946(V) now has one of its three suppositions discharged.

## Evidence log
- [x] referee: 6 triples × H ∈ {5,20,80,320}, envelope bounded, C/H decay confirmed
- [ ] T₄ strip estimate (the resonance line cu+dt=0 neighborhood) — the named next
- [ ] Lean rendering of the congruence lemma (ZMod orbit + case bound; klein's Thm892
      toolbox fits)
