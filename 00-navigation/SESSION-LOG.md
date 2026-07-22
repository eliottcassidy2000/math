## kind-pasteur-2026-07-21-S128c147 -- DvdK1 elementary reductions (zero-charge + both-signs + gcd) + the honest bypass verdict

**Owner:** prove GMC(2) unconditionally, own remaining tasks, explore repo for a DvdK bypass /
alternative route. Long session, frequent push/pull.

**VERIFIED STATE.** GMC(2)/NC2 is a kernel-pure reduction to a SINGLE external input DvdK1:
`heightWitnessSupplier_holds` is PROVED (GMC2NC2, green), so `nc2_of_dvdK1 : DvdK1 -> NC2`. DvdK1
route = THM-2067 concrete 16-thm skeleton; the two open inputs are Vieta `hΩ` (core done by
death-star, assembly bounded/boxeph) and THM-1550 `hS/hfix` (the deep Wiener-Hopf analytic gap,
death-star; Henselian obstacle done). Check A / Galois action / irreducibility all done+green.

**MY CONTRIBUTION (pushed, kernel-pure [propext,Classical.choice,Quot.sound]) -- GMC2DvdKZeroCharge.lean:**
- `dvdK1_of_zero_mem`: if some `q i0 = 0`, DvdK1 holds at m=1 (the single mass `delta_{i0}` is the
  unique balanced channel, CT = c i0 != 0).
- `dvdK1_of_bothSigns : DvdK1BothSigns -> DvdK1`: since `ChargesStraddleZero = (zero charge) OR
  (both signs)` and the zero disjunct is now discharged, the ENTIRE DvdK1 reduces to the both-signs
  (no-zero) case -- formally isolating exactly what the Galois route must prove.
- `constantTermRelation_scale`: rescaling all charges by a nonzero common factor leaves the
  constant-term relation unchanged => the hard case may assume gcd of charges = 1.
Combined with prior elementary cases (two-charge/unique-channel/positive), the NO-CANCELLATION
territory of DvdK1 is now elementary+formal; the hard case is pinned to {both signs, no zero, gcd 1,
>=2 coincident balanced channels at every m}.

**BYPASS VERDICT (honest; the explored task).** NO elementary bypass of the analytic core exists.
THM-2067 (multiplicative orbit-product) and ALL THREE THM-2101 additive proofs (monodromy /
transcendental specialization + zero Lagrange sum / t-adic Newton-polygon + partial fractions) bottom
out at the SAME small-root-product / Wiener-Hopf content (THM-1550). The difficulty is intrinsic:
complex-coefficient phase cancellation across >=2 channels for ALL m simultaneously (positive-real
coeffs => DvdK1 trivial). So GMC(2) is exactly ONE theorem -- THM-1550 -- from unconditional, a genuine
multi-session analytic formalization owned by death-star, not sidesteppable.

**HONEST.** I did NOT make GMC(2) unconditional; THM-1550 is untouched and is the blocker. I sharpened
and formalized the perimeter (zero/both-signs/gcd reductions), isolated the hard case, and delivered
the bypass verdict. Files GMC2DvdKZeroCharge.lean (3 lemmas, green) + reflection
dvdk1-elementary-reductions-and-the-bypass-verdict-kps-S128c147. HYP-8975. Complementary to
mac-mini-S162 / boxeph-S236 / death-star (no collision: zero/both-signs/gcd reductions were unowned).

## death-star-2026-07-22-S113 -- GMC2: full Vieta (prod_roots_Phi) kernel-pure = the hOmega input for the concrete THM-2067; fleet advanced (a) concurrently

**Owner directive:** finish BOTH (a) Gal-instantiation and (b) THM-1550 myself, pulling in concurrent pieces; push/pull often.

- **DELIVERED (kernel-pure, GMC2PhiVieta.lean, lake-built):** **prod_roots_Phi** -- over any E where Phi = X^M - tR splits, (Phi.map E).roots.prod = (-1)^(deg R) * algebraMap(algebraMap(R.coeff 0 / R.leadingCoeff)), a CONSTANT (t cancels) = valuation 0. Built on my coeff_ratio_Phi_eq_const (S112) + Splits.coeff_zero_eq_leadingCoeff_mul_prod_roots. This is the number-theoretic CONTENT of boxeph's hOmega hypothesis.
- **PULLED IN (fleet advanced (a) while I worked):** boxeph S236 (HYP-8956) -- the CONCRETE Gal instantiation, GMC2Thm2067Concrete.thm2067_contradiction_concrete + GMC2GalRootAction (incl. the equivariance/rootsEquivRoots piece); mac-mini S162 -- irreducibility of the degree-(M+N) Phi over F(t) (GMC2DvdKParameterIrreducible.irreducible_map_ratfunc, discharges hPhi). So thm2067_contradiction_concrete now needs only: hPhi [mac-mini, DONE], hOmega [my Vieta content, DONE], hS [THM-1550, the DEEP gap].
- **ATTEMPTED (a) equivariance myself (~14 iterations):** localized it precisely to a Mathlib instance-diamond -- the goal's `algebraMap SF SF` (from IsScalarTower.toAlgHom F SF SF) vs the canonical Algebra.id SF; rfl/simp/algebraMap_self_apply/Subsingleton all fail on the mismatch. Handed the precise blocker to boxeph -- who had already solved the full instantiation (S236), superseding my attempt.
- **hOmega final wiring WIP:** connecting my prod_roots_Phi (roots multiset) to boxeph's hOmega (product over the rootSet subtype) hit dependent-type/Fintype plumbing (rw rootSet_def breaks the rootSetFintype motive); ~90% there, not landed this session.
- **HONEST:** Vieta landed kernel-pure (the hOmega content); (a) is near-complete via the fleet's pieces + mine (only the hOmega subtype-product wiring remains); (b) THM-1550 (Wiener-Hopf product / additive germ identity) is the deep multi-session crux, being worked by codex. Did NOT fully finish -- the deep bridge (b) is genuinely multi-session. Pushed/pulled/coordinated throughout as directed. HYP-8976 (kps took HYP-8975 concurrently). NOTE: kps-S128c147 (above) independently reached the SAME honest verdict from the additive side -- the bypass verdict: THM-2067 AND all 3 THM-2101 additive routes bottom out at THM-1550; the analytic core is intrinsic, confirming (b) is the sole remaining crux.
## mac-mini-2026-07-22-S162 -- GMC(2)/DvdK Galois route: FORMALIZED death-star-S106 target (2) -- the degree-(M+N) parameter polynomial Phi = Y^M - t*R(Y) is IRREDUCIBLE over F(t) (kernel-pure, self-contained, general field F). Discharges the open `Irreducible Phi` hypothesis of boxeph's GMC2Thm2067Concrete.
## boxeph-2026-07-22-S237 -- assembly bridge: Irreducible(death-star's Phi R M) from my (A); irreducibility + Vieta on one Phi (kernel-pure, HYP-8975)

**Owner:** concurrent GMC(2) push; work my remaining pieces; pull; coordinate.

**DELIVERED kernel-pure (GMC2DvdKAssembly):**
- algebraMap_comp_C: composite F->F[t]->F(t) = algebraMap F (RatFunc F).
- Phi_eq_map_swap: death-star's GMC2PhiVieta.Phi R M = my map(algebraMap)(Bivariate.swap(C(X^M)-C R X)). Via map_map + algebraMap_comp_C + map_C/algebraMap_X/map_X + ring.
- irreducible_Phi: Irreducible(Phi R M) from phi_irreducible_ratfunc.
SETTLES mac-mini's degree flag (my A = degree-(M+N), not degree-1; mac-mini conflated with phi_t_irreducible). Irreducibility + Vieta now on the SAME Phi.

**REMAINING assembly:** instantiate thm2067_contradiction_concrete at Phi R M: hPhi=irreducible_Phi [DONE], hOmega=Vieta (needs rootSet<->multiset bridge: separable=>nodup=>Finset-prod=multiset-prod, + sign fold), hS/hfix=THM-1550 (death-star). SPLIT: death-star adds rootSet-form Vieta OR I take it; I start Check A (no coeff_pow => MvPolynomial.coeff_add_pow/custom).

**Honest:** clean kernel-pure integration bridge; not the final assembly (hOmega bridge + Check A remain, split with death-star) or full GMC(2) (THM-1550, in progress). Artifacts: reflection the-assembly-bridge-...-boxeph-S237.md, HYP-8975, GMC2DvdKAssembly.lean.

