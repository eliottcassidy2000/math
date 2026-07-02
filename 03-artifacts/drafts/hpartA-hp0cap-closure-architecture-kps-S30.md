# Closing hpartA (and the witness floor with it): the c-averaged ruler architecture

**kind-pasteur-2026-07-01-S30 (HYP-3953).** Pure-math session on the two Lean-skeleton axioms.
Numerics: `lrc14_hpartA_cruler_fubini_kps.py`, `lrc14_hpartA_adversarial_joint_kps.py` (+ .out).

## 0. The targets

- **hpartA**: `0 < witnessG2(shape) → Mreach(v) ≥ 1/14`, where witnessG2 = meas{x ∈ G_P :
  maxgap{frac(u_i x)} > 1/7} for the cluster co-offsets u_i.
- **hfloor** (its sibling; the pair is what the preferred DAG consumes): witnessG2 ≥ m_P =
  14249/252252, split as k ≤ 7 (pigeonhole, done) + k = 8..13 (the rhoGlobFloorRat ledger).
- **hp0cap**: p0(E) ≤ cap_k (the sector route; §6 below).

The honest conclusion of this session: **the right proof does not factor through the literal
hpartA interface.** The pair (hfloor, hpartA) closes TOGETHER through four exact identities and
one finite-census floor. The architecture below replaces both axioms; each piece is either
proved here, elementary, or a named finite computation.

## 1. The exact c-ruler identity (kills the drift argument)

Let V be a reference speed and write a cluster element as ℓ = V − o (offset o = V − ℓ). At the
ruler times τ_j = (j + c)/V, j ∈ ℤ:

    ‖ℓ τ_j‖ = ‖(j + c) − o τ_j‖ = ‖o τ_j − c‖        (EXACT — no approximation, no drift).

So with G_c := {x ∈ G_P : ∀i ‖o_i x − c‖ ≥ 1/14} (note: the reference runner V itself has o = 0,
contributing the pure c-constraint ‖c‖ ≥ 1/14 — see §4 "retirement"):

    (∃ c, ∃ j : (j+c)/V ∈ G_c)  ⟹  M(S) ≥ 1/14.

THM-527-A's drift/window argument, the "Vmax-ruler embedding", and the "equidistribution
ρ_K → ρ*" step all collapse into the elementary lattice count (THM-565 part 2, already
sorry-free in Lean): #{j : (j+c)/V ∈ G_c} ≥ V·meas(G_c) − arcCount(G_c). Hence

    **(R)  ∃c : V·meas(G_c) > arcCount(G_c)  ⟹  M(S) ≥ 1/14.**

Verified end-to-end (T2): five 13-set families, margins +163..+1720, all match exact M ≥ 1/14.

## 2. The Fubini gap identity (what G2 really is)

    **(F)  ∫₀¹ meas(G_c) dc = ∫_{G_P} F(x) dx,   F(x) := Σ_gaps (gap{frac(o_i x)} − 1/7)⁺.**

Proof: Fubini on the joint (x, c)-event; at fixed x the free-c set is exactly the union of the
gap surpluses. Verified exact to 6 decimals on three shapes (T1). Consequences:
- witnessG2 > 0 ⟺ ∫_{G_P} F > 0 ⟹ sup_c meas(G_c) > 0. (The skeleton's opaque witnessG2 should
  be DEFINED as ∫_{G_P}F — the integral form is the one the proof consumes, and positivity is
  equivalent.)
- The c-average diagonalizes the arithmetic: in Fourier, E_c kills every relation with Σm ≠ 0;
  only the sum-zero sublattice {m : Σ m_i o_i = 0, Σ m_i = 0} survives. **The covering adversary
  can cover the target c = 0 but cannot cover all targets at once** — this is why the c-averaged
  ledger (§5) is robust exactly where the homogeneous lonely ledger collapses.

## 3. The rotation identity (the fast scale cancels; T5/T8's discovery)

    **(Ω)  F_{V−offs}(x) = F_{offs}(x) for every x and every V.**

Proof: the phases {frac((V−o_i)x)} are a rotation (by Vx) of the reflected offset phases
{frac(−o_i x)}; circular gaps are invariant under rotation and reflection. ∎

This is the structural heart of the session. Consequences:
- The free-gap functional of a narrow (bounded-offset) window is a SLOW, BOUNDED-COMPLEXITY
  function of x — the scale V appears nowhere in it. There is NO two-scale equidistribution
  lemma to prove at the top level; the fast scale enters only through the pointer frac(Vx),
  i.e. through the elementary count (R). ("The rate was hiding inside the lattice count.")
- The tower functional is exactly scale-free: E[1_{G_B}·F_{U₁}·F_{U₂(V₂)}] is numerically
  IDENTICAL for V₂ = 50, 500, 5000 (T8) — because it is identical, by (Ω), as a function.
- The historical "≤ 6-far equidistribution rate" residual (HYP-3787) is NOT needed on this
  route. Where prior routes estimated |meas − limit| ≤ C/V, here the limit identity is exact
  and only the count (R) carries V.

## 4. Retirement, and the recursion that terminates

In (R), the reference runner V contributes no x-constraint — only ‖c‖ ≥ 1/14. **Each level of
the scale recursion retires its reference runner into a c-condition**, so the x-problem
strictly loses one runner per level: depth ≤ 13, termination structural.

The wide case (offsets not bounded relative to V — the multi-scale sets) recurses: the window's
avoidance event {x : ∀i ‖o_i x − c‖ ≥ 1/14} is itself a (shifted-target) lonely event of the
offset set — the same species one level down, with its own reference and its own c. The
multi-target Fubini (F) applies verbatim per level (independent c_ℓ averages). The sampling is
the standard nested-interval walk: each level needs its pointer to wrap ≥ 2 times inside the
current good interval, which pins the uniform window-ratio constant Λ* ≈ 4/w (w = the working
gap-surplus width) and per-level V* thresholds with FINITE sub-V* checks — all elementary,
same species as THM-565's V* atlas (worst V* ≈ 234 there).

## 5. The floor: the c-averaged ledger (finite, and robust where it must be)

The one genuinely quantitative input left is the per-level joint floor

    **(⋆)  J(B; U₁,…,U_L) := ∫ 1_{G_B}(x) · Π_ℓ F_{U_ℓ}(x) dx ≥ (explicit positive ledger value)**

over the FINITE data (bounded bottom part B, bounded offset-difference patterns U_ℓ — bounded
after (Ω) + recursion). Measured ledgers (T3, T6, T7, T8):

| object | k | min found | at | comparison |
|---|---|---|---|---|
| A(U) = E_x[F_U] | 2 | 0.734693 | (5,37) | = (6/7)² to 6 digits — **no dip below independence** |
| A(U) | 7 | 0.294221 | (3,5,6,7,8,9,11) | homogeneous L⁰ min = 0.217; ratio to (6/7)^k: 0.866 vs 0.640 |
| A(U) | 13 | 0.114168 | (odd-heavy) | clears witnessMP = 0.0565 ×2.0 |
| union-bound death row (k=7) | 7 | 0.2980 | consecutive | 55× above THM-594-E's Parseval floor |
| adversarial joint J/meas(G_B) | — | 0.0247 | B={1..6}, U={0..9} (inadmissible, over-constrained) | positive with room |
| two-window joint | — | 0.0174 | doubly-consecutive worst | exactly V-free |

Reading: the c-average is the smooth face of the lonely measure. The homogeneous inf collapses
on covering sets (the whole difficulty of LRC); the c-averaged inf stays within 13% of
independence through k = 7 and clears the m_P target by 2–4.4× through k = 13. The T4 row is
the union-bound death (j = 7 = 1/(2r), where mac-mini's THM-594-E Parseval floor 0.00545 is the
unconditional backstop): observed 0.298 — the floor's slack is enormous; the Parseval bound is
the in-principle guarantee, the census the working constant.

Honest normalization note: A(U) is the P-free layer. The admissible objects (|P| + k = 13)
restrict to G_P, which SHRINKS J (T6's 0.0247 is an over-constrained stress, |P|+k = 16); the
admissible census is the finite computation the architecture names as its base case, exactly
parallel to the S27/S28 census standard.

## 6. What this gives hp0cap

Less, honestly — hp0cap's assembled route (bounded check + THM-563 periodicity + dichotomy
HYP-2788 + R-tail/gK8) stands as the plan of record. Two transfers from this session:
- (F)-type c-averaging applied to the SECTOR side smooths the same arithmetic: the sector-miss
  events are inhomogeneous lonely events (target (j+½)/7), and their c-smoothed versions have
  the same sum-zero-lattice support — a route to replace the "VERIFIED not symbolic" pieces of
  the dichotomy by exact second-moment arithmetic (THM-594-B makes pair terms exact rationals).
- Via HYP-2832 (cap = min meas(G_P); floor = corollary of p0 ≤ cap), any strengthening of the
  witness floor reduces the LOAD on hp0cap: with (⋆) closing the witness route, hp0cap becomes
  REDUNDANT for LRC(14) (the preferred DAG consumes (hfloor, hpartA) only). Closing the pair
  via §1–§5 retires hp0cap from the critical path entirely.

## 7. The irreducible residual, stated once

Everything above is exact/elementary except one family of numbers:

    **(⋆-census) the admissible joint ledger J(B; U₁,…,U_L) > 0 with explicit constants,
    over bounded B and bounded difference-patterns U_ℓ, per level of the recursion.**

For bounded patterns this is a finite computation (the T3/T6/T7 ledgers are its first pages).
The only infinite direction is pattern WIDTH inside a single window — and (Ω) + retirement
recurse it down with strict cardinality descent, so the census is over patterns of ≤ 12
elements with entries bounded by the level's Λ*-window — finite at every level. The
orbit-closure/difference-core heuristics (opus HYP-3901: rational directions = closed orbits =
census minima; irrational = larger measure) say the census minima ARE the global minima — the
same principle already verified for the S27/S28 censuses.

## 7.5 [S31 addendum] The torus-band theorem: the (⋆)-census is SYMBOLIC (HYP-3954)

The c-averaged ledger lives on the (x,c)-torus as the complement of a union of integer-slope bands
B_i = {‖c − u_i x‖ < h}. Elementary theorem (verified V1-V3, 6 digits everywhere):
(i) vol(B_i) = w = 2h; (ii) **pairs are EXACTLY independent** — (x,c) ↦ (c−u_i x, c−u_j x) has
determinant u_j − u_i ≠ 0, so vol(B_i ∩ B_j) = w² with NO arithmetic correction (THM-594-B's
fluctuations live on the x-circle at fixed c; the c-average removes them); (iii) d-fold overlaps are
2-dim SUBTORUS BOX VOLUMES with Fourier series over the saturated sum-zero relation lattice
Λ_S = {m : Σm_i = 0, Σm_i u_i = 0} — exact rationals at h = 1/14; the AP relation (1,−2,1) gives
2h², the maximal triple weight, which is WHY AP patterns minimize A. Checks: A(pair) = 36/49,
A(AP-triple) = 61/98, A(k=4 argmin) = 11/21 — inclusion-exclusion equals the sampled ledger to six
digits on every case. Bonferroni-3 certifies the floor ≥ witnessMP through k = 8; beyond, run the
census at full depth or by exact c-breakpoint integration (L^c is piecewise-linear in c with rational
breakpoints — exact rational integration per pattern). The bands are the 2-dimensional LR-ZONOTOPE
picture (generators (u_i, 1)) — the c-ruler recursion is an LR-zonotope recursion.

## 9. [S31 addendum] The strategic layer: MSS finite checking and the sLRC-false filter

- **MSS (arXiv:2411.06903, Forum Math. Sigma 2025):** speeds up to C(n+1,2)^{n−1} suffice to verify
  LRC for n+1 runners — for LRC(14): **91^12 ≈ 3.2·10^23, unconditional**. Every "unbounded/wide"
  residual in this project is FINITE-IN-PRINCIPLE by citation; the project's censuses, cutoffs
  (W* ≤ 513, V* ≈ 234, Λ = 10^4) are effective-bound improvements of ~20 orders of magnitude on the
  studied branches. The honest framing of the whole program: **making MSS effective at n = 14.**
- **BCS (arXiv:2603.24784):** the SHIFTED LRC is FALSE from n = 5, and the Lonely Vector Property
  fails from n = 12. DESIGN FILTER for every open lemma: any statement quantified over ALL shifts
  has sLRC strength and is at risk of being false — state floors as AVERAGES over shifts (Fubini) or
  EXISTENCE of good shifts, never ∀-shift. (This architecture complies: (F) is an average, (R) is an
  existence. Audit flagged for: the odd-covering bridge (covering systems carry free offsets = the
  shifted side; Erdős–Selfridge lives there) and any windowed-tiling lemma quantified over windows.)
- **Guo–Sun (math/0412217):** odd covering systems with distinct squarefree moduli need lcm with
  ≥ 22 prime divisors — the arithmetic-budget species of our THM-522/2-adic-tower arguments; note
  the LRC-pinned analogue of "odd covering" is trivial (all-odd sets miss q = 2 ⇒ q-witness at 1/2),
  which locates the depth of Erdős–Selfridge strictly on the shifted side of the dictionary.

## 8. Summary for the DAG

Replace the axiom pair (hfloor, hpartA) by:
  (D1) definition witnessG2 := ∫_{G_P} F (§2) — makes the opaque symbol concrete;
  (D2) the count reduction (R) (§1) — elementary, mostly Lean-ready (THM-565 part 2 exists);
  (D3) the identities (F), (Ω) (§2, §3) — one-paragraph proofs;
  (D4) the nested sampling with Λ*, V* (§4) — elementary, finite sub-V* checks;
  (D5) the (⋆)-census (§5, §7) — finite, first pages computed, 2–4.4× margins.
No equidistribution rate appears. hp0cap leaves the critical path (§6).
