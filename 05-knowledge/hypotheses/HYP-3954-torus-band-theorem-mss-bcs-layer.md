# HYP-3954: the torus-band theorem — the c-averaged ledger is symbolic; + the MSS/BCS strategic layer

**Status:** VERIFIED — kind-pasteur-2026-07-01-S31. kps block (3950+).
**Companion to:** HYP-3953 (the c-ruler architecture); architecture doc §§7.5, 9.

## The theorem (elementary; the session's core)
Write the c-averaged ledger as a torus volume: A(U) = 1 − vol(∪_i B_i) on the (x,c)-torus, with
B_i = {(x,c) : ‖c − u_i x‖ < h} the width-w band of integer slope u_i (h = 1/14, w = 2h).
1. vol(B_i) = w.
2. **Pair independence, exact:** (x,c) ↦ (c − u_i x, c − u_j x) has determinant u_j − u_i ≠ 0, hence
   is measure-preserving T² → T²; vol(B_i ∩ B_j) = w² with NO arithmetic correction. THM-594-B's
   two-branch fluctuations live on the x-circle at fixed target; the c-average removes them.
3. **d-fold = subtorus box volume:** the joint law of (c − u_i x)_{i∈S} is Haar on the 2-dim subtorus
   cut by the saturated sum-zero relation lattice Λ_S = {m : Σm_i = 0, Σm_i u_i = 0};
   vol(∩_S B_i) = Σ_{m∈Λ_S} Π_i ŵ_h(m_i) (Poisson), exact rationals at h = 1/14.
   AP relation (1,−2,1): closed form 2h² — the maximal triple weight — WHY AP minimizes A.
4. A is translation- and dilation-invariant (a difference-pattern functional).

## Verification (lrc14_torus_band_theorem_kps.out)
- V1: pair overlap = 0.020408 = 1/49 at all slope pairs incl. resonant (7,21) and adjacent (100,101).
- V2: Fourier = grid to 6 digits on all triples; AP = 2h² exactly; (2,−3,1) → 0.006803,
  (6,−7,1) → 0.002915 (decay with coefficient size).
- V3: full inclusion-exclusion = the S30 sampled ledger to SIX DIGITS on every case:
  A(pair) = 36/49; A(AP-triple) = 61/98 (both (1,2,3) and (10,24,38) — pattern invariance);
  A(1,2,4) = 0.625850; A(k=4) = 11/21 at both (1,2,3,4) and (5,6,7,8) (translation invariance).
- V4: Bonferroni-3 (exact terms) certifies A ≥ witnessMP for k ≤ 8 (0.517/.412/.312/.201/.099 vs
  0.0565); dies at k = 9 (alternating blowup). Beyond: full-depth IE (2^k terms, all exact) or the
  **c-breakpoint engine** (L^c(U) is piecewise-linear in c with rational breakpoints — exact rational
  integration per pattern). **HYP-3953's (⋆)-census is SYMBOLIC.**

## The strategic layer (owner's references, ingested)
- **MSS (arXiv:2411.06903):** speeds ≤ C(n+1,2)^{n−1} suffice — LRC(14) ⟸ finite check ≤ 91^12,
  UNCONDITIONAL. Every wide/unbounded residual in the repo is finite-in-principle; the program's
  honest framing: making MSS effective at n = 14 (current effective cutoffs W* ≤ 513, V* ≈ 234).
- **BCS (arXiv:2603.24784):** shifted LRC FALSE from n = 5; Lonely Vector Property false from n = 12.
  **Design filter:** ∀-shift lemmas have sLRC strength = false-risk; state floors as shift-averages
  (Fubini) or shift-existence only. This architecture complies; audit flagged for the odd-covering
  bridge (free offsets = the shifted side) and windowed-tiling lemmas.
- **Guo–Sun (math/0412217):** odd covering ⟹ ≥ 22 prime divisors — arithmetic-budget species of
  THM-522/2-adic-tower; the pinned analogue of odd-covering is TRIVIAL (all-odd sets miss q = 2 ⇒
  q-witness at t = 1/2) — Erdős–Selfridge's depth lives strictly on the shifted side.
- The bands are the 2-dim **LR-zonotope** picture (generators (u_i, 1)) — the c-ruler recursion is an
  LR-zonotope recursion; the covering-radius school (MSS/BCS/Schymura) is its theorem backend.

## Artifacts
- 04-computation/lrc14_torus_band_theorem_kps.py (+ .out)
- 03-artifacts/drafts/hpartA-hp0cap-closure-architecture-kps-S30.md (§7.5, §9 addenda)
- TANGENTS: T-crossing-numbers (Zarankiewicz/Harary-Hill ↔ equidistribution extremality; annulus
  crossing profile of the winding tournament — logged, unexplored).

## Depends on / relates to
HYP-3953, THM-594-B, mac-mini HYP-3852/3853 (the d-fold Bernoulli program — the torus version is the
clean home; their MT-slice c_AP identity is a band-triple layer), klein HYP-3847/3848/3849,
MSS arXiv:2411.06903, BCS arXiv:2603.24784, Guo–Sun math/0412217, HYP-2764, OPEN-Q-108.
