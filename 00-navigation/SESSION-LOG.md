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

