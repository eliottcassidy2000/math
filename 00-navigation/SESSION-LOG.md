## boxeph-2026-07-22-S239 -- hfull DISCHARGED via Weierstrass; additive root-packet lemma self-contained (kernel-pure, HYP-8985)

**Owner:** keep working hfull (full-root Lagrange sum), think Weierstrass.

**DELIVERED kernel-pure (GMC2FullRootPhi):** discharged hfull = sum_{roots} a^k/Phi'(a) = 0 for the NON-MONIC Phi via the Weierstrass product form (codex's Lagrange was for the monic nodal).
- phi_map_eq: Phi.map = C(lc)*Lagrange.nodal (Splits.eq_prod_roots + dedup).
- aeval_deriv_eq: aeval a Phi' = lc*nodal'(a) (derivative_map + derivative_C_mul).
- full_root_sum_eq_zero: full-root sum for non-monic Phi = 0 (reindex rootSet<->aroots.toFinset via Finset.sum_subtype + mem bridge; factor 1/lc; codex Lagrange + card_rootSet_eq_natDegree).
- root_packet_eq_zero_selfcontained: root_packet_eq_zero with hfull discharged internally.

**STATE:** additive-route ALGEBRAIC CORE now FULLY self-contained kernel-pure -- b_k(S) in F => 0, no hfull hyp, no THM-1550/product/Hensel/Vieta. Whole additive skeleton kernel-pure (GMC2GalRootAction + codex GMC2LaurentShiftCheckA + GMC2RootPacketConcrete + GMC2FullRootPhi). ONLY remaining: the b=1 wrapper (sum_{S_+} a^(M-1)/Phi'(a) = F(t) = sum D_m t^m => b=1; shared valuation input, additive form; cf mac-mini-S163).

**Honest:** hfull DISCHARGED (last hypothesis-level gap in the additive core removed). Not full GMC(2): b=1 residue-sum wrapper remains. Artifacts: reflection hfull-discharged-via-weierstrass-...-boxeph-S239.md, HYP-8985, GMC2FullRootPhi.lean.

