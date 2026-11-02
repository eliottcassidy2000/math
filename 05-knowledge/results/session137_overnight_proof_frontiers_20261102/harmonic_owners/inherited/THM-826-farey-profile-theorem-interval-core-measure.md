---
id: THM-826
title: THE FAREY PROFILE THEOREM — the COMPLETE measure profile of the interval core. For 0 ≤ λ ≤ 1/(k+1), m({1,…,k}; λ) = Σ over consecutive pairs a/i < b/j of the Farey sequence F_k (on the circle) of max(0, (1 − λ(i+j))/(ij)). The profile is piecewise linear and convex, with breakpoints exactly at λ = 1/s for s = k+1, …, 2k−1, except that s = 2k−2 is omitted when k is even; on [1/(s+1), 1/s] it is m = A_s − λ·B_s with A_s = Σ_{i+j≤s} 1/(ij), B_s = Σ_{i+j≤s} (i+j)/(ij) (sums over consecutive-F_k pairs). ENDPOINT IDENTITIES: A_{2k−1} = Σ 1/(ij) = 1 (the classical Farey telescope: m(0)=1); initial slope B_{2k−1} = 2·Σ_{l≤k} φ(l)/l (the HYP-2856 totient constant, ~12k/π²); THM-819 (the primitive harmonic law) is the FIRST SEGMENT (s = k+1); the area under the profile is the Franel-type sum Σ 1/(2ij(i+j))
status: PROVED (two-lemma proof below: gap-locality + no-intrusion; both elementary) + REFEREED exact (k = 2..12, dense rational λ grid including all breakpoints; every identity checked in ℚ)
source: kind-pasteur-2026-07-15-S128 (cont.10; owner: prove the full measure profile with Farey breakpoints)
depends_on:
  - THM-819   # the s = k+1 segment (the primitive harmonic law); this theorem is its completion
related:
  - HYP-2856 / OPEN-Q-108 / HYP-2899 (the product-Möbius packet ledger): the profile's initial slope IS the totient floor constant 2Σφ(l)/l; the profile shows the Div(D)-axis of the Div×B_r ledger scalarizes CLEANLY on interval cores (disjoint Farey-gap atoms, no Boolean inclusion-exclusion needed) — locating all B_r (far-element/Venn A..G) complexity in the multi-scale bodies
  - THM-805 (the Tutte/Dirichlet-kernel bridge), HYP-6865 (the electrical staircase), mac-mini THM-736
  - codex THM-801 (Möbius–Čech descent lifting the A+B+C−D−E−F+G corner recursion — the staircase-side Venn whose measure-side triviality this theorem establishes)
---

# THM-826 — the Farey profile theorem

**Setup.** G(λ) = {t ∈ [0,1) : ‖jt‖ ≥ λ, j = 1..k}, m(λ) = |G(λ)|. F_k = the Farey fractions of
order k on the circle (all primitive a/i, i ≤ k, including 0/1 once). Between consecutive fractions
a/i < b/j: the GAP has length b/j − a/i = 1/(ij) (unimodularity bi − aj = 1).

**Lemma 1 (nesting).** The bad arc of speed j at an imprimitive center a/j (gcd d > 1) is contained
in the arc of speed j/d at the same point (radius λ/j < λ·d/j). Hence the bad set is the union over
PRIMITIVE c/l ∈ F_k of arcs of radius λ/l — one effective arc per Farey fraction.

**Lemma 2 (no intrusion).** For λ < 1/(k+1), the arc at c/l never crosses a neighbouring Farey
fraction's arc INTO the gap beyond it: crossing a/i (i ≠ l) from outside would need
λ(1/l − 1/i) > |a/i − c/l| ≥ 1/(il), i.e. λ > 1/(i−l) ≥ 1/(k−1) > 1/(k+1) — impossible. So each gap
is consumed exactly by its two endpoint arcs, at rates 1/i (left) and 1/j (right).

**Theorem.** The good set is the disjoint union of gap remnants, so
m(λ) = Σ_gaps max(0, 1/(ij) − λ(1/i + 1/j)) = Σ_gaps max(0, (1 − λ(i+j))/(ij)). ∎

**Structure.** A gap dies at λ = 1/(i+j); consecutive pairs have k+1 ≤ i+j ≤ 2k−1 (mediant property:
i+j > k; max at j=k−1,i=k).  The realized values are k+1..2k−1, except that 2k−2 is absent when
k is even, giving the corresponding breakpoint set {1/s}.  Indeed, sum 2k−2 would force the
denominator pair (k−2,k), (k−1,k−1), or (k,k−2); for even k these are non-coprime or equal. On
[1/(s+1), 1/s] the live gaps are those with i+j ≤ s: m = A_s − λB_s. Endpoints: A_{2k−1} = Σ1/(ij) =
Σ gap lengths = 1; B_{2k−1} = Σ(1/i + 1/j) = 2Σ_{l≤k}φ(l)/l (each fraction borders two gaps).
THM-819 is the s = k+1 segment: the only pairs with i+j = k+1 are the witnesses-adjacent gaps, and
A_{k+1} − λB_{k+1} at λ = 1/(k+2) reproduces 2δ·H^prim. ∫₀^{1/(k+1)} m = Σ 1/(2ij(i+j)).

**Dictionary notes (owner's mandate: metagraph ⟷ tournament ⟷ measure).**
- The gap atoms are the measure-side analogue of the staircase corner cells: the A..G corner
  recursions (A+B−C = B2 even half; A+B+C−D−E−F+G = B3 full; the χ₇ odd mode) do Boolean
  inclusion–exclusion because staircase corners OVERLAP; the Farey gaps DON'T — the interval core's
  Div-axis is inclusion-exclusion-free. That asymmetry is exactly why interval cores are the
  tractable stratum of route [B] and multi-scale bodies (overlapping far arcs) are not.
- Electrical: the profile is the staircase Smith network's response curve under threshold load; the
  s-th segment's slope B_s counts "1/denominator current" through the live gaps.

## Evidence log

- [x] referee exact ℚ: profile formula == direct arc-union measure at all λ ∈ {c/720 : c} ∩ [0,1/(k+1)]
      plus every breakpoint 1/s, k = 2..12; segment form A_s − λB_s checked; A=1, B=2Σφ/l, THM-819
      segment, and the area identity checked exactly (thm826_farey_profile_referee_kps_S128c10.py)
- [x] even-k omission `2k−2` and the corrected realized-sum range checked through k=64
      (thm841_no_dyadic_breakpoint_codex_S14.py)
