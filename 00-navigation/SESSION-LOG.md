## kind-pasteur-2026-07-21-S128c149 -- additive route: hfull discharged (parallel to boxeph) + assembled to the b=1 wrapper + the honest Abel-duality analysis

**Owner:** finish the additive route; many push/pulls; look at prior additive/multiplicative duality
for inspiration.

**DELIVERED (GMC2FullRootConcrete.lean, kernel-pure [propext,Classical.choice,Quot.sound], pushed):**
- `fullRootSum_eq_zero` (hfull): `∑_{α:Φ.rootSet L} α^{M-1}/Φ'(α) = 0`, via Φ separable ⇒
  `Φ.map = C(lc)·nodal(roots)` ⇒ `Φ' = lc·nodal'`, the rootSet↔Finset bridge, and codex's nodal
  Lagrange lemma. (Mathlib `Splits` is now single-arg `Splits (f)`; the rootSet-subtype bridge +
  `conv_lhs` to avoid rewriting `Ψ.leadingCoeff` were the sharp edges.)
- `additive_dvdk_reduces_to_smallSum`: `fullRootSum_eq_zero` + boxeph's `root_packet_eq_zero` ⇒ the
  additive one-variable DvdK contradiction is COMPLETE MODULO the single b=1 wrapper `hb`, with NO
  THM-1550 / product / Hensel / Wiener-Hopf.

**COLLISION (honest):** boxeph-S239 (HYP-8985) discharged hfull CONCURRENTLY via the SAME Weierstrass
method (GMC2FullRootPhi). My `fullRootSum_eq_zero` PARALLELS it (duplicate; defer to first-pusher).
My net-new: the `additive_dvdk_reduces_to_smallSum` assembly + the honest DUALITY ANALYSIS.

**DUALITY ANALYSIS (reflection).** The additive/multiplicative duality is THM-2101 §6's ABEL OPERATOR
`A(G-1)=∫(G-1)ds/s`: integrating the residue SUM produces `log(Π)`, the `1/m`, Hensel factorization,
Wiener-Hopf. So the additive b=1 wrapper (`∑_S w = 1`) and the multiplicative THM-1550 (`Π = c·t`) are
THE SAME valuation content seen through the Abel operator. This TEMPERS my S128c148 over-claim
("additive fully escapes analysis"): the CLOSING ALGEBRA escapes THM-1550 (proved kernel-pure:
orbit-sum + hfull + root-packet), but the SMALL-PACKET SELECTION (b=1) remains -- as a SUM
(partial-fraction residue, a cleaner Mathlib target than the multiplicative product/Hensel) not a
product.

**HONEST.** Did NOT close the b=1 wrapper (the Newton-polygon small-root-packet-residue-sum identity,
plausibly multi-session). The additive route now: irreducibility (mac-mini) + Galois (boxeph) + hfull
(boxeph/me) + root-packet (boxeph), all kernel-pure ⇒ DvdK1 ⟸ b=1 wrapper alone. NEXT: the formal
partial-fraction `∑_S α^{M-1}/Φ'(α) = [u⁰]u^{M-1}/Φ` (mac-mini-S163 root-free F(t) = RHS); Newton-polygon
small-root selection; `F(t)=1` under vanishing. Files GMC2FullRootConcrete.lean + reflection
additive-route-hfull-discharged-reduced-to-b1-wrapper-kps-S128c149. HYP-8990. cf boxeph HYP-8985.

## death-star-2026-07-22-S113 -- GMC2: hOmega (Vieta rootSet form) DISCHARGED + concrete multiplicative THM-2067 REDUCED to THM-1550 alone; hB independently derived then deferred to boxeph (convergence)

**Owner:** finish BOTH (a) Gal-instantiation and (b) THM-1550 myself, pulling in concurrent pieces; push/pull often.

- **DELIVERED kernel-pure [propext,Classical.choice,Quot.sound] (both pushed, in TournamentH7 root):**
  1. **GMC2PhiVieta.prod_rootSet_Phi** -- the rootSet-form hOmega EXACTLY as boxeph's concrete wrapper needs: prod over the rootSet SUBTYPE of Phi R M = algebraMap(RatFunc.C d), d=(-1)^(deg R)*(R.coeff 0/R.leadingCoeff). Bridge: separable=>nodup, Finset.prod_attach + prod_eq_multiset_prod + toFinset_val + Nodup.dedup collapses subtype-product to multiset-root-product, then prod_roots_Phi + fold sign/nested-algebraMap into RatFunc.C. Discharges the hOmega that was WIP at compaction.
  2. **GMC2Thm2067Reduced.thm2067_reduced_to_thm1550** (new file) -- instantiates boxeph's thm2067_contradiction_concrete at Phi R M with hPhi=GMC2DvdKAssembly.irreducible_Phi (boxeph) + hOmega=prod_rootSet_Phi (mine): takes ONLY the THM-1550 data (S,x0,c,hc,hfix,hS)+separability, returns False. Concrete MULTIPLICATIVE THM-2067 for Phi R M now formally reduced to THM-1550 ALONE (needed maxHeartbeats 1000000 + open Polynomial).
- **CONVERGENCE (honest):** independently derived the additive-route hB (full-root Lagrange sum=0, Weierstrass: Phi.map=C(lc)*nodal => Phi'(a)=lc*nodal'(a), factor 1/lc, close with codex's monic Lagrange) kernel-pure -- but boxeph S239 landed the SAME Weierstrass proof ~1 min earlier (GMC2FullRootPhi + self-contained wrapper). DEFERRED to boxeph, removed my duplicate GMC2FullRootVanishing. Two independent kernel-pure derivations = strong confirmation of boxeph's hfull.
- **PULLED IN / connected:** boxeph S237 (Phi_eq_map_swap + irreducible_Phi = hPhi for my exact Phi; offered me the rootSet Vieta corollary -> I took+delivered it); boxeph S236 (concrete Gal instantiation); mac-mini S162 (degree-(M+N) irreducibility); kps S128c147->c148 (bypass verdict then its correction: additive route stays before integration); boxeph S238 (root_packet_eq_zero uses MY GMC2GalRootAction.coe_smul).
- **NET STATE:** both routes' algebraic cores are now kernel-pure and reduced to the SAME shared valuation-selection input -- multiplicative=THM-1550 (small-root product=c*t), additive=boxeph's b=1 packet-sum wrapper. That shared valuation/Newton-polygon core is the one genuine remaining analytic gap (kps's verdict), multi-session. HONEST: did NOT close GMC(2); delivered the multiplicative hOmega + reduced theorem (non-dup) and confirmed boxeph's hfull. HYP-8976.

## boxeph-2026-07-22-S239 -- hfull DISCHARGED via Weierstrass; additive root-packet lemma self-contained (kernel-pure, HYP-8985)

**Owner:** keep working hfull (full-root Lagrange sum), think Weierstrass.

**DELIVERED kernel-pure (GMC2FullRootPhi):** discharged hfull = sum_{roots} a^k/Phi'(a) = 0 for the NON-MONIC Phi via the Weierstrass product form (codex's Lagrange was for the monic nodal).
- phi_map_eq: Phi.map = C(lc)*Lagrange.nodal (Splits.eq_prod_roots + dedup).
- aeval_deriv_eq: aeval a Phi' = lc*nodal'(a) (derivative_map + derivative_C_mul).
- full_root_sum_eq_zero: full-root sum for non-monic Phi = 0 (reindex rootSet<->aroots.toFinset via Finset.sum_subtype + mem bridge; factor 1/lc; codex Lagrange + card_rootSet_eq_natDegree).
- root_packet_eq_zero_selfcontained: root_packet_eq_zero with hfull discharged internally.

**STATE:** additive-route ALGEBRAIC CORE now FULLY self-contained kernel-pure -- b_k(S) in F => 0, no hfull hyp, no THM-1550/product/Hensel/Vieta. Whole additive skeleton kernel-pure (GMC2GalRootAction + codex GMC2LaurentShiftCheckA + GMC2RootPacketConcrete + GMC2FullRootPhi). ONLY remaining: the b=1 wrapper (sum_{S_+} a^(M-1)/Phi'(a) = F(t) = sum D_m t^m => b=1; shared valuation input, additive form; cf mac-mini-S163).

**Honest:** hfull DISCHARGED (last hypothesis-level gap in the additive core removed). Not full GMC(2): b=1 residue-sum wrapper remains. Artifacts: reflection hfull-discharged-via-weierstrass-...-boxeph-S239.md, HYP-8985, GMC2FullRootPhi.lean.

