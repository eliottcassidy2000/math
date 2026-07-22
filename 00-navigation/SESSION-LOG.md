## kind-pasteur-2026-07-21-S128c146 -- GMC(2) orthogonality CORE extracted as a standalone Mathlib-PR-ready lemma (ThreeTermRecCoprime.lean)

**Owner directive:** get the GMC(2) formalization complete + Mathlib-PR-ready; be creative bypassing
hard pieces; poke for ideas.

- **STATE (honest):** the full GMC2NC2Capstone is sorry-free but CONDITIONAL on two premises --
  `DvdK1` (external DvdK analytic input) + `HeightWitnessSupplier`. `DvdK1` is being actively
  completed by death-star/boxeph/codex via THM-2067/THM-1550/Henselian (two routes, no sorries in
  the spine). That whole artifact is NOT a Mathlib PR: research-grade, repo-specific, gated on
  unformalized analysis. So I did not fight the analytic endpoint (well-covered).
- **CREATIVE BYPASS -- extract the piece that IS Mathlib-worthy.** New file `ThreeTermRecCoprime.lean`
  (builds green, 8475 jobs; axioms `[propext, Classical.choice, Quot.sound]`): for a monic three-term
  recurrence `p(n+2)=(X-C a_{n+1})·p(n+1)-C b_{n+1}·p n` over ANY `CommRing`, if every `b_{n+1}` is a
  UNIT then **`IsCoprime (p n) (p (n+1))`** (`isCoprime_succ`). This is the estimate-free "consecutive
  orthogonal polynomials are coprime" (Favard; Hermite/Legendre/Laguerre/Chebyshev uniform) -- the
  replacement for the REFUTED domination step (THM-1585 -> THM-1660). Plus field version, `noCommonRoot`
  `[Nontrivial R]`, `monic` + `natDegree = n`, probabilists' Hermite instance.
- **STRENGTHENS prior art:** the repo's `GMC2HermiteNoCommonRoot.lean` (functions `ℝ→ℝ`) and
  `ThreeTermRecurrence.lean` (`no_common_root_poly [IsDomain R]`) only give NO-COMMON-ROOT; this gives
  `IsCoprime` (algebra-level, stronger) over general `CommRing` (drops `IsDomain`) + monic/degree.
  Mathlib GAP confirmed: has `Polynomial.hermite`/`Chebyshev` but NO general three-term-recurrence
  coprimality. Proof = one Bézout update closed by `linear_combination` (no analysis).
- **IDEAS worth reusing (reflection):** `IsCoprime` beats "no common root" for Mathlib-readiness;
  bypass the coprime-API scavenger hunt with an EXPLICIT Bézout witness; `Ring.inverse` +
  `Ring.inverse_mul_cancel` give the unit inverse over a general ring; two bugs -- `noCommonRoot`
  needs `[Nontrivial R]`, and the recurrence consumes `b(n+1)` so Hermite is `b k = k` (giving
  `He₂=X²-1`), not `b k=k+1`.
- **Deliverable:** `ThreeTermRecCoprime.lean` + `ThreeTermRecCoprime-MATHLIB-PR-NOTES.md` (statement,
  gap justification, provenance, pre-submission checklist: narrow imports via `shake`, optional
  Legendre/Chebyshev instances). Full GMC(2) is one analytic theorem (`DvdK1`) from unconditional --
  correctly owned by the others, not blocked by anything here. Reflection
  gmc2-orthogonality-core-is-mathlib-pr-ready-kps-S128c146. HYP-8970. Complementary to HYP-8965/8951/8960.

## death-star-2026-07-22-S112 -- GMC2 formalization: the Vieta valuation-0 coeff-ratio input for the THM-2067 wrapper (kernel-pure) + critical-path coordination

**Owner directive:** keep improving the GMC2 formalization, working substantial pieces.

- **SURVEYED** the crowded, fast-moving GMC2/DvdK1 formalization: codex's additive-orbit-residue route (THM-2101: Check A, additive orbit-sum, full-root Lagrange, small-root existence via reciprocal-monic Hensel -- which IMPORTS AND USES my S111 GMC2Henselian.HenselianLocalRing, confirming it on the critical path) + boxeph's multiplicative orbit-product wrapper (S235 GMC2Thm2067Wrapper, premising THM-1550 + Vieta). No sorries in the spine; DvdK1 is being completed via two routes.
- **DELIVERED (kernel-pure, lake-built):** GMC2PhiVieta.lean **coeff_ratio_Phi_eq_const** -- for Phi = X^M - t*R over F(t), Phi.coeff 0 / Phi.leadingCoeff = algebraMap(R.coeff 0 / R.leadingCoeff), i.e. the ratio is t-FREE (a constant in the image of F): coeff 0 = -t*r0, leadingCoeff = -t*lc(R), the t cancels. This is the number-theoretic content of the wrapper's hOmega input (Vieta: prod_roots = (-1)^d*(coeff0/leadingCoeff) has valuation 0).
- **CRITICAL COORDINATION:** boxeph S235 asked me to derive HenselianLocalRing (PowerSeries F) + the degree-dropping factorization -- unaware BOTH are done: HenselianLocalRing is my S111 (GMC2Henselian, HYP-8960), and the degree-drop is SOLVED by codex's reciprocal/reversed-monic trick (GMC2ReciprocalSmallRoots, importing my instance). Messaged boxeph to prevent duplicated work + offered the Vieta coeff-ratio for hOmega.
- **STATE:** THM-2067 is a kernel-pure SKELETON down to (a) boxeph's concrete Gal instantiation (equivariance via rootsEquivRoots) + Check A + Vieta (core supplied here) and (b) THM-1550 piece (iii), the Wiener-Hopf D_m=0 forall m <=> prod_small=c*t bridge (mine, the one deep analytic gap; obstacles (i) Henselian + (ii) degree-drop are now handled).
- **HONEST:** a modest kernel-pure lemma + high-value coordination (prevented boxeph re-deriving my Henselian). The DvdK1 core is being competently completed by codex/boxeph; my clear contribution (Henselian foundation) is on the critical path. HYP-8965.
## boxeph-2026-07-22-S235 -- THM-2067 wrapper assembled kernel-pure (full orbit-product argument as a Lean reduction, HYP-8951)

**Owner:** work the Galois wrapper and the deep analytic gap.

**DELIVERED kernel-pure (GMC2Thm2067Wrapper):**
- pow_monomial_eq_const_absurd: (c*t)^N = d^K impossible over F(t) (via S233 monomial_pow_ne_const).
- thm2067_contradiction: the COMPLETE THM-2067 argument as a reduction -- transitive Galois-type action on Omega + equivariant f:Omega->E over F(t) + THM-1550 (prod_S f = c*t, Gal-fixed) + Vieta (prod_Omega f = const d) => False. Via prod_pow_card_group_eq (S232) => pullback along injective algebraMap F(t)->E => monomial closing.
So THM-2067 is now a kernel-pure Lean SKELETON: GMC2OrbitProduct -> GMC2RatFuncClosing -> GMC2PhiIrreducible (A) -> GMC2Thm2067Wrapper (13 thms).

**REMAINING (2 hard pieces, scoped):**
1. Concrete Gal instantiation: instances + transitivity (galAction_isPretransitive <- A) free; EQUIVARIANCE routes through Gal.smul_def = rootsEquivRoots (needs characterizing, E=SplittingField) + Vieta + Check A. Mine, bounded.
2. THM-1550 deep gap (prod_small=c*t Gal-fixed) = unramified Hensel. Obstacles: IsAdicComplete(PowerSeries) not synthesizable, degree-dropping Weierstrass factorization, Wiener-Hopf bridge. Multi-session, death-star owns.

**Honest:** wrapper/argument assembled kernel-pure (genuine milestone -- orbit product provably yields the contradiction, isolating THM-1550 + Vieta + Gal-instantiation). NOT full GMC(2). Artifacts: reflection the-thm2067-wrapper-assembled-...-boxeph-S235.md, HYP-8951, GMC2Thm2067Wrapper.lean.

