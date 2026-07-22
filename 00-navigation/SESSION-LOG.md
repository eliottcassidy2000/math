## boxeph-2026-07-22-S236 -- concrete Galois instantiation of THM-2067, kernel-pure (HYP-8956)

**Owner:** work the two remaining pieces (Gal wrapper + Check A; deep Hensel gap); push/pull; mine threads.

**DELIVERED kernel-pure:**
- GMC2GalRootAction: direct action sigma.x=sigma(x) on Phi.rootSet L. mem_rootSet_smul, rootAction (MulAction), coe_smul (TAUTOLOGICAL equivariance, sidesteps galAction/rootsEquivRoots), isPretransitive_rootAction (transitivity for irreducible Phi over Normal L via IsConjRoot.exists_algEquiv).
- GMC2Thm2067Concrete.thm2067_contradiction_concrete: instantiates the wrapper at Phi.Gal on Phi.rootSet SplittingField. Irreducible Phi over F(t) + small-root product=c*t Gal-fixed (THM-1550) + full product=const d (Vieta) => False. Instances resolve via Phi.Gal type.
=> ENTIRE THM-2067 = concrete kernel-pure reduction (6-module chain, 16 thms), open inputs = Vieta + THM-1550.

**COORDINATION (split w/ death-star):** death-star drove THM-1550 -- obstacle (i) HenselianLocalRing(PowerSeries F) DONE kernel-pure (GMC2Henselian.lean) + monic M-th-root Hensel + a_j*Y_j reparametrization (no degree-dropping factorization theorem needed). Split: death-star = fixed-point convergence + Vieta-for-Pi + Wiener-Hopf; me = Gal-instantiation (DONE) + Vieta(hOmega) + Check A. Wrapper takes THM-1550 as hS/hfix.

**REMAINING (mine, bounded):** Vieta (prod roots=(-1)^d r0/rd; non-monic Vieta + separability) + Check A (CT(Lambda^m)=[u^Mm]R^m; no coeff_pow in Mathlib => custom lemma).

**Honest:** the Galois wrapper/instantiation piece is COMPLETE kernel-pure (transitive action + tautological equivariance + concrete reduction). NOT full GMC(2): Vieta + Check A (mine) + THM-1550 (death-star, obstacle (i) done) remain. Artifacts: reflection the-concrete-galois-instantiation-...-boxeph-S236.md, HYP-8956, GMC2GalRootAction.lean, GMC2Thm2067Concrete.lean.

