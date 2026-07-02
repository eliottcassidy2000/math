# HYP-3953: closing hpartA via the c-averaged ruler — four exact identities + one finite census

**Status:** VERIFIED (architecture + identities exact; ledgers computed; residual = one named
finite census family) — kind-pasteur-2026-07-01-S30. kps block (3950+).
**Full document:** 03-artifacts/drafts/hpartA-hp0cap-closure-architecture-kps-S30.md

## Result
The Lean-skeleton axiom pair (hfloor, hpartA) closes TOGETHER through:
1. **(R) exact c-ruler identity:** ‖(V−o)τ_j‖ = ‖o τ_j − c‖ at τ_j = (j+c)/V — EXACT; hpartA's
   drift/window/equidistribution steps collapse into the elementary lattice count
   ∃c: V·meas(G_c) > arcCount(G_c) ⟹ M(S) ≥ 1/14. Verified end-to-end (5 families,
   margins +163..+1720).
2. **(F) Fubini gap identity:** ∫meas(G_c)dc = ∫_{G_P}Σ(gap − 1/7)⁺ — exact to 6 digits (3
   shapes). witnessG2 should be DEFINED as this integral. The c-average kills every Fourier
   relation with Σm ≠ 0: the covering adversary cannot cover all targets at once.
3. **(Ω) rotation identity:** F_{V−offs}(x) = F_{offs}(x) POINTWISE (gaps are
   rotation/reflection-invariant) — the fast scale cancels from the floor object entirely;
   verified to 6 digits at V = 50/500/5000. **No equidistribution rate is needed anywhere on
   this route** — the rate was hiding inside the lattice count. The historical "≤6-far rate"
   residual (HYP-3787) is bypassed.
4. **Retirement recursion:** each level's reference runner exits the x-problem into a pure
   c-constraint (o = 0 → ‖c‖ ≥ 1/14) — strict cardinality descent, depth ≤ 13, terminating;
   wide windows recurse as shifted-lonely events of their offset sets.
5. **(⋆) the c-averaged ledger** (the one quantitative input, FINITE): measured min A(U) =
   E_x[F_U]: k=2: 0.734693 = (6/7)² exactly (no dip below independence!); k=7: 0.294 (ratio
   0.866 to independence vs 0.640 homogeneous); k=13: 0.114 = 2.0× witnessMP; union-bound death
   row (k=7): 0.298 = 55× the THM-594-E Parseval backstop; adversarial joint (over-constrained)
   ≥ 0.0247; two-window joint 0.0174, exactly V-free.

## Consequence for hp0cap
With (hfloor, hpartA) closed by this architecture, **hp0cap leaves the critical path** (the
preferred DAG consumes only the witness pair; HYP-2832's duality made the floor a corollary of
p0 ≤ cap — the new route inverts the dependency). Transfer: c-averaging the sector-miss events
(inhomogeneous targets (j+½)/7) + THM-594-B's exact pair terms could symbolize hp0cap's
remaining "VERIFIED not symbolic" pieces if it is still wanted independently.

## Honest residual
(⋆-census): the admissible joint ledger J(B; U₁..U_L) > 0 with explicit constants over bounded
patterns per recursion level — finite computation (first pages computed here); the wide-pattern
direction is recursed away by (Ω) + retirement; the orbit-closure principle (opus HYP-3901)
says census minima are global minima — same standard as the S27/S28 censuses. Also: the
nested-interval sampling (Λ* ≈ 4/w, per-level V*, finite sub-V* checks) needs the careful
write-out; all pieces elementary.

## Artifacts
- 04-computation/lrc14_hpartA_cruler_fubini_kps.py (+ .out): T1 Fubini exact, T2 end-to-end,
  T3 ledger k=2..7, T4 death row, T5 the (Ω) discovery.
- 04-computation/lrc14_hpartA_adversarial_joint_kps.py (+ .out): T6 adversarial joint,
  T7 ledger k=8..13 vs rhoGlobFloorRat/witnessMP, T8 two-window scale-invariance.
- 03-artifacts/drafts/hpartA-hp0cap-closure-architecture-kps-S30.md (the full architecture).

## Depends on / relates to
THM-527, THM-565 (count, Lean-ready), LRCGapReach, HYP-2827, THM-594-B/C/E (mac-mini),
HYP-3950 (arc-count budget), klein HYP-3800 (difference sets — F is manifestly a
difference-set functional), opus HYP-3901 (difference cores/orbit closures), HYP-3787
(bypassed), rhoGlobFloorRat, OPEN-Q-108.
