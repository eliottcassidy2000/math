# HYP-3953: closing hpartA via the c-averaged ruler — exact identity, Fubini gap integral, and the average inhomogeneous lonely measure

**Status:** CLAIMED (stub) — kind-pasteur-2026-07-01-S30, in progress. kps block (3950+).

## Target
**hpartA** (Lean skeleton axiom): `0 < witnessG2(shape) → Mreach(v) ≥ 1/14`, where
witnessG2 = meas{x ∈ G_P : maxgap{frac(u_i x)} > 1/7} (semantics currently opaque in Lean).
Pure-math session (owner directive: no Lean).

## The three claimed moves (to be developed + stress-tested this session)
1. **EXACT c-RULER IDENTITY (no drift):** for cluster elements ℓ_i = V − u_i, at ruler times
   τ_j = (j+c)/V one has ‖ℓ_i τ_j‖ = ‖u_i τ_j − c‖ EXACTLY. So M(S) ≥ 1/14 iff some ruler point
   lies in G_c := {x ∈ G_P : ∀i ‖u_i x − c‖ ≥ 1/14} — replacing the drift/window argument of
   THM-527-A with a lattice count: hpartA ⟸ ∃c: V·meas(G_c) > arcCount(G_c) (THM-565 part 2).
2. **FUBINI GAP IDENTITY:** ∫₀¹ meas(G_c) dc = ∫_{G_P} F(x) dx, F(x) = Σ_gaps(gap − 1/7)⁺.
   Positivity of witnessG2 ⟺ positivity of the integral ⟹ sup_c meas(G_c) > 0. Closes the
   TIGHT-cluster case (Σu_i small vs V·∫F) outright.
3. **THE c-AVERAGE KILLS Σm≠0 RELATIONS:** E_c[L^c(U)] (average inhomogeneous lonely measure)
   has Fourier support only on the sum-zero sublattice {m : Σm_i u_i = 0, Σm_i = 0} — the
   covering adversary cannot cover ALL targets c. Conjecture: min_U E_c[L^c(U)] has a robust
   positive floor (unlike the homogeneous inf, which vanishes on covering sets), giving the
   wide-case level-cut induction a uniform Λ* and feeding the k=8..13 witness floor
   (rhoGlobFloorRat ledger). k_T = 7 critical backstop = THM-594-E Parseval/DMNR floor.

## Artifacts (to fill)
- 04-computation/lrc14_hpartA_cruler_fubini_kps.py (+ .out)

## Depends on / relates to
THM-527, THM-565 (three-gap sampling), LRCGapReach, HYP-2827 (pigeonhole), THM-594-C/E
(DMNR + Parseval floor), HYP-3950 (arc-count budget), opus HYP-3901 (renormalization),
rhoGlobFloorRat ledger, OPEN-Q-108.
