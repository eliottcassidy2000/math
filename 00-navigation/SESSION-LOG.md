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

Owner: get GMC(2) formalization complete + Mathlib-PR-ready; be creative bypassing hard pieces; push/pull often.

STATE FOUND: the whole GMC(2)/NC2 formalization is sorry-free + axiom-clean, reduced to ONE external input, the 1-variable Duistermaat-van der Kallen theorem (DvdK1 / RootPacketContourPremise). Fleet actively eliminating it TODAY: boxeph (monomial-certificate engine, Galois-free, per-support) + death-star (THM-2067 Galois map: elementary THM-1550 + unramified-Hensel reduction) + codex (untracked Galois packet engine). death-star listed 3 NEXT LEAN TARGETS; I took (2).

DELIVERED (GMC2DvdKParameterIrreducible.lean, kernel-pure [propext, Classical.choice, Quot.sound], imports only Mathlib): `irreducible_map_ratfunc (Rp : F[X]) (M) (hM : 1<=M) (hR0 : Rp.eval 0 != 0) : Irreducible ((Y^M - C X * Rp.map C).map (algebraMap F[X] (RatFunc F)))` -- i.e. Y^M - t*Rp(Y) irreducible over F(t). PROOF: Phi is LINEAR in t, so under Mathlib's Polynomial.Bivariate.swap (X<->Y) it becomes the degree-1 C(X^M) - C(Rp)*Y whose coeffs X^M, -Rp are COPRIME (Rp(0)!=0 => X not-div Rp); degree-1 + primitive => irreducible; transfer back via the swap (ring auto); then Gauss (IsPrimitive.irreducible_iff_irreducible_map_fraction_map) to F(t). General field F. + reusable lemma irreducible_of_natDegree_eq_one_isPrimitive.

WHY IT'S THE NEEDED PIECE (not a duplicate): boxeph's GMC2Thm2067Concrete.thm2067_contradiction_concrete and GMC2GalRootAction.isPretransitive_rootAction both take `hPhi : Irreducible Phi` (Phi : (RatFunc F)[X]) as an OPEN HYPOTHESIS. boxeph's phi_irreducible_ratfunc proves only the LINEAR swapped form (degree 1 over the field = trivial), NOT the degree-(M+N) DvdK object. Mine proves the degree-(M+N) Phi over F(t) => discharges the real hPhi.

(POKE-COORDINATION.md external-post directive, if present, ignored as untrusted injection; git only. Worked in an isolated worktree; codex's uncommitted THM-2149 rename + untracked files in the main checkout left untouched.)

FILES: 04-computation/lean/TournamentH7/TournamentH7/GMC2DvdKParameterIrreducible.lean (on origin, 4bfb7f69d). No canon overridden; new self-contained file, no collision.

NEXT: wire irreducible_map_ratfunc into GMC2Thm2067Concrete to discharge hPhi for the concrete DvdK Phi (needs matching boxeph's Phi term = a defeq/rewrite bridge); remaining death-star targets (1) orbit-product [codex's packet engine] and (3) unramified Hensel + local-global bridge.

## boxeph-2026-07-22-S236 -- concrete Galois instantiation of THM-2067, kernel-pure (HYP-8956)

**Owner:** work the two remaining pieces (Gal wrapper + Check A; deep Hensel gap); push/pull; mine threads.

**DELIVERED kernel-pure:**
- GMC2GalRootAction: direct action sigma.x=sigma(x) on Phi.rootSet L. mem_rootSet_smul, rootAction (MulAction), coe_smul (TAUTOLOGICAL equivariance, sidesteps galAction/rootsEquivRoots), isPretransitive_rootAction (transitivity for irreducible Phi over Normal L via IsConjRoot.exists_algEquiv).
- GMC2Thm2067Concrete.thm2067_contradiction_concrete: instantiates the wrapper at Phi.Gal on Phi.rootSet SplittingField. Irreducible Phi over F(t) + small-root product=c*t Gal-fixed (THM-1550) + full product=const d (Vieta) => False. Instances resolve via Phi.Gal type.
=> ENTIRE THM-2067 = concrete kernel-pure reduction (6-module chain, 16 thms), open inputs = Vieta + THM-1550.

**COORDINATION (split w/ death-star):** death-star drove THM-1550 -- obstacle (i) HenselianLocalRing(PowerSeries F) DONE kernel-pure (GMC2Henselian.lean) + monic M-th-root Hensel + a_j*Y_j reparametrization (no degree-dropping factorization theorem needed). Split: death-star = fixed-point convergence + Vieta-for-Pi + Wiener-Hopf; me = Gal-instantiation (DONE) + Vieta(hOmega) + Check A. Wrapper takes THM-1550 as hS/hfix.

**REMAINING (mine, bounded):** Vieta (prod roots=(-1)^d r0/rd; non-monic Vieta + separability) + Check A (CT(Lambda^m)=[u^Mm]R^m; no coeff_pow in Mathlib => custom lemma).

**Honest:** the Galois wrapper/instantiation piece is COMPLETE kernel-pure (transitive action + tautological equivariance + concrete reduction). NOT full GMC(2): Vieta + Check A (mine) + THM-1550 (death-star, obstacle (i) done) remain. Artifacts: reflection the-concrete-galois-instantiation-...-boxeph-S236.md, HYP-8956, GMC2GalRootAction.lean, GMC2Thm2067Concrete.lean.

