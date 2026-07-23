## mac-mini-2026-07-22-S167 -- GMC(2) CLOSED (by kind-pasteur; I CEDED my parallel Omega-lift): independent corroboration + full-chain axiom audit + the ANALYTIC-HEART LINE-AUDIT (statements verified, no defect)

**Owner directive:** finish the GMC(2) formalization, most formalized state possible for a mathlib PR; prioritize open mathematical reasoning; push/pull often; get key math solid before builds.

**WHAT HAPPENED.** I traced the whole spine and found GMC(2) had been reduced to ONE Prop, `GMC2DvdKUnivariateReduction.SinglePolyCrux` (`dvdK1_of_crux` proves the Duistermaat-van den Kallen citation FROM it, so closing it REMOVES the citation; `gmc2_of_crux` gives GMC(2)). I verified there is no hidden premise: `GMC2HeightWitness.heightWitnessSupplier_holds` is PROVED, so `nc2_of_dvdK1`/`gmc2_of_dvdK1` need only DvdK1. I then built the Omega-lift discharging it. **kind-pasteur landed the same construction first** (`GMC2DvdKOmegaWiring`, 0784107bd, HYP-9020), so I **CEDED and deleted mine** rather than push a duplicate.

**INDEPENDENT CORROBORATION (the durable value of the duplicate).** The two derivations were written without sight of each other and agree on every structural choice -- Omega = AlgClosure(LaurentSeries C) (NOT LaurentSeries itself: the small roots are ramified, t^{1/M} not in C((t))), the rfToL-based algebra, `IsAlgClosed.lift`, and nodup-from-separability -- and on the constant: kps's cc = (-1)^{M+1} r0 equals my (-1)^{deg Pw}(-r0) since deg Pw = M. For a claim this size, two independent constructions landing on the same normalization is worth more than either alone. kps's version is the better one to keep (the rfl scalar-tower diamond beats my composite bookkeeping).

**THE ANALYTIC-HEART LINE-AUDIT (my claimed piece of kps's polish plan; the part kernel-purity CANNOT establish).** The referee flagged `hderiv_final`/`smallRootFactor_coeff0_of_vanish` as kernel-checked but not line-audited. These are MY lemmas (S165/S166) and the deepest link. Kernel-purity proves "no sorry"; it does NOT prove the STATEMENT is the one the mathematics needs. Three targets, ALL PASS:
- (i) `hderiv_final`'s conclusion is exactly the consumer's hypothesis. `unitCoeff0 R M` = constantCoeff_x(weierstrassUnit) in F[[t]] = literally h(0,t); the closing math is honest: d_t h(0,t)=0 + h(0,0)=1 + char 0 => h(0,t)=1. No exp/log/Puiseux/Fredholm.
- (ii) CharZero is load-bearing and CORRECTLY SCOPED -- used exactly at `eq_C_of_derivativeFun_eq_zero`, which is FALSE in char p (t^p has zero derivative). Deliberate, not incidental: `GMC2DvdKCharZeroClosing` carries an explicit `omit [CharZero F] in` on the one step that doesn't need it.
- (iii) NO silent completion mismatch (my old two-completion blocker) -- discharged by three PROVED lemmas: the single transport hom `phi : F[[t]][[x]] -> (F((x)))[[t]]` with `hfact` moving the Weierstrass factorization by map_mul (never re-derived); `Rl_pow_coeff` (frame moment = polynomial moment, exactly where a vanishing hypothesis could be silently swapped); and `xCoeff0_phi_unit` (the seam back to the original ring).
VERDICT: no defect found. Statement-level audit; does not re-verify upstream proofs.

**FULL-CHAIN AXIOM AUDIT.** Wrote `AxiomCheckGMC2MacMiniS167.lean` auditing the WHOLE chain, not just the capstone (a capstone can read clean while an upstream lemma is not): heightWitnessSupplier_holds, nc2/gmc2_of_dvdK1, singlePolyCrux_holds, dvdK1_of_crux, the upstream value/divisibility/bridge inputs, gmc2_unconditional. Run reached 8505/8542 modules with EVERY line reporting exactly [propext, Classical.choice, Quot.sound] and **0 sorryAx**, before the host killed the process; rerun in flight to reach the endpoints. NOTE: the build environment is contended -- a sibling agent's broad `pkill -f "lake build"` killed several of my runs (exit 144), and `lake env lean` transiently failed on a corrupt `plausible` dep.

**HONEST STATUS.** GMC(2) unconditional is kps's result, kernel-pure in its file, corroborated by my independent construction and by an analytic-core statement audit that found no defect. I did NOT re-verify the upstream proofs (irreducible_Phi, Weierstrass, the height package) -- an adversarial review should target those and the hderiv analytic core, not the Omega-wiring, which is the shallowest link.

**REMAINING for "most formalized" (packaging, not mathematics):** ~90 GMC2 modules, 32 with a blanket `import Mathlib` needing precise imports; dead/exploratory module pruning; retiring or clearly marking the superseded conditional surfaces (`GMC2NC2Capstone`'s `..._of_heightWitnessSupplier`) so the "remaining inputs" story is not misread. kps claimed the canonical capstone (GMC2Main) + Mathlib-ifying their gap lemmas; import-hygiene and pruning are UNCLAIMED.

## kind-pasteur-2026-07-22-S128c154b -- GMC(2) "most formalized state" polish: canonical capstone + gap-lemma Mathlib drafts + minimal-chain footprint

**Owner:** keep working to achieve the full proof in the most formalized state possible.

**DELIVERED (all kernel-pure, pushed):**
- **GMC2Main.GMC2.gmc2** -- the CANONICAL documented front-door theorem (in the GMC2 namespace) + a
  self-contained GREEN build target: `lake build TournamentH7.GMC2Main` compiles the whole GMC(2) proof
  (8543 jobs), INDEPENDENT of the failing LRC modules. #print axioms = [propext,Classical.choice,Quot.sound].
  Full docstring: statement (E = Wick expectation), architecture (GMC(2)<=NC2<=DvdK1<=SinglePolyCrux <=
  Weierstrass Pi=c*t <= hderiv, Omega=AlgebraicClosure(LaurentSeries C)), verification.
- **3 Mathlib PR-draft gap lemmas** (05-knowledge/mathlib-pr-drafts-powerseries-gaps): eq_C_of_derivativeFun
  _eq_zero (char-0 converse; generalized [Field]->[CommRing][NoZeroDivisors][CharZero]), geometric inverse
  (1-C w*X)^{-1}=mk(w^.), PowerSeries.logDeriv + derivativeFun_map + logDeriv_map. ALL verified to compile
  standalone (lake env lean, exit 0), confirmed absent from Mathlib. Gap 1 = cleanest first PR.
- **Minimal-chain footprint:** 88 GMC2 modules total; 69 ESSENTIAL (transitive deps of the capstone); 19
  OFF-CHAIN -- superseded alternative routes (GMC2AdditiveOrbitSum/FullRootPhi/FullRootConcrete/Generating
  Function/TopLagrange/RootPacketConcrete/DvdKXDeg = the additive route + death-star's superseded toolkits)
  + early exploratory (Henselian/HermiteNoCommonRoot/GoodReduction/Formalization/MomentBasics/etc.). Kept as
  research records; prunable for a Mathlib PR footprint. NOT deleted (shared code + milestone records).

**STATE:** GMC(2) is in the best formalized state for correctness (kernel-verified + adversarially SOUND +
sorry-free) with a clean documented capstone + green build target + extracted PR-ready gaps + a mapped
minimal footprint. Remaining polish (actual pruning, import-Mathlib->specific for 32 files, analytic-heart
human-legibility, C->general-char-0-field generalization) is diminishing-returns/risky/others-lane; flagged
for the fleet. Reflection gmc2-verification-and-mathlib-pr-readiness-kps-S128c154. HYP-9020.
## kind-pasteur-2026-07-22-S128c154 -- GMC(2) VERIFICATION AUDIT: kernel-pure + adversarially SOUND + LRC-independent; Mathlib-PR scoping

**Owner:** double-check GMC(2) is in the best state possible + ready for a Mathlib PR.

**CORRECTNESS VERIFIED (as strong as a formalization gets):** gmc2_unconditional builds green (8542-job
chain), #print axioms = [propext,Classical.choice,Quot.sound] (no sorryAx/native_decide); NO
sorry/admit/axiom anywhere in the chain (grep hits were docstrings + #print axioms lines); an INDEPENDENT
adversarial referee (general-purpose subagent, full-chain read) returned SOUND -- statement is the genuine
GMC(2) (real Wick E, non-vacuous SinglePolyCrux, vanishing hypothesis genuinely consumed), Omega-instance
sound (rfl guards + single hS_of_dvd_value => no silent mismatch), no circularity, gmc2_of_crux chain
gap-free (HeightWitnessSupplier discharged, orbit-product double-count correct, monomial_pow_ne_const
refutation genuine).

**REPO CAVEAT (NOT GMC(2)):** lake build TournamentH7 (whole root) FAILS -- but ONLY on in-progress LRC
modules (LRCCoherentBlockerChronology, LRCPairTowerValuation, LRCTwoCircleII, ...). NO GMC2 file imports
any LRC module (verified), so GMC(2) is unaffected. The LRC WIP needs fixing for a green root (fleet lane).

**MATHLIB-PR READINESS (honest):** the FULL proof is NOT one PR -- 87 GMC2 files, 47 in-chain, 32 import
Mathlib wholesale, repo naming, dead modules; whole-upstreaming = multi-week refactor. What IS PR-able now:
3 general Mathlib-gap lemmas -- (1) eq_C_of_derivativeFun_eq_zero (char-0 converse of derivativeFun_C):
EXTRACTED, generalized to [CommRing R][NoZeroDivisors R][CharZero R], VERIFIED to compile standalone
(scratch PR-draft-charzero-derivative.lean) -- ready as first PR; (2) geometric inverse (1-C w*X)^{-1}=
mk(w^.); (3) derivativeFun_map + PowerSeries.logDeriv/logDeriv_map. All confirmed absent from Mathlib.

**VERDICT:** GMC(2) rigorously verified + best-state for correctness; Mathlib path = upstream the 3 gap
lemmas first (char-0 drafted+compiles), full proof is a separate larger effort. Reflection
gmc2-verification-and-mathlib-pr-readiness-kps-S128c154. HYP-9020 (cont).
## kind-pasteur-2026-07-22-S128c153-FINAL -- GMC(2) PROVEN UNCONDITIONALLY, KERNEL-PURE (Omega-wiring closed)

**Owner:** finish the Omega-wiring + remaining parts. DONE.

**GMC2DvdKOmegaWiring.gmc2_unconditional** (origin/main 0784107bd, #print axioms =
[propext,Classical.choice,Quot.sound], NO sorry): for P Q : MvPolynomial (Fin 2) C, if all E(P^m)=0 (m>=1)
then eventually E(Q*P^m)=0. = gmc2_of_crux (boxeph) applied to my singlePolyCrux_holds. THE multi-agent
GMC(2) formalization is COMPLETE.

**singlePolyCrux_holds (the last hypothesis, discharged).** Omega := AlgebraicClosure(LaurentSeries C).
The non-synthesizable Algebra(RatFunc C)(LaurentSeries C) via death-star rfToL.toAlgebra =>
AlgebraicClosure.instAlgebra lifts to Omega; IsScalarTower.of_algebraMap_eq rfl; psi := IsAlgClosed.lift.
Pomega = (smallRootFactor).map(ofPowerSeries).map(algebraMap ..Omega); dvd from death-star
smallRootFactor_map_dvd_phiVieta_map (+map_dvd_map', goal-side map_map); value Pomega.coeff0 =
algebraMap(-C(r0)*X) from mac-mini smallRootFactor_coeff0_of_vanish via ofPowerSeries_comp_C +
rfToL_comp_algebraMap; then boxeph hS_of_dvd_value => packet product = algebraMap(C cc * X),
cc = (-1)^{M+1} r0 != 0. x0 witness via IsSplittingField.splits + Splits.roots_ne_zero + mem_rootSet.

**KEY FINDING:** the Algebra(RatFunc)(Laurent) instance diamond that blocked the Omega-wiring for multiple
sessions resolves BY RFL once the base algebra is handed over as rfToL.toAlgebra. Viability was only
knowable by attempting the build. I twice mis-scoped it (trivial-composition over-claim RETRACTED; then
feared intractable -- wrong).

**CREDIT (whole fleet):** mac-mini (Weierstrass, hderiv_final, value, hconst), death-star (frame, xCoeff0,
transpose phi, rfToL/PhiCoincide, xdeg toolkit, hconst), boxeph (bridge, dvd, hS_of_dvd_value,
gmc2_of_crux), codex (Check A/Lagrange). My through-line: F=D_m leg, h-side (disk/annulus + logDeriv_map),
hderiv assembly backbone (hderiv_of_transpose_glue -- hderiv_final built on it), Omega-wiring.
Reflection gmc2-proven-unconditional-omega-wiring-closed-kps-S128c153. HYP-9020. Referee welcome.
## kind-pasteur-2026-07-22-S128c153-cont2 -- hderiv assembly backbone is LOAD-BEARING (mac-mini hderiv_final built on it); Omega-wiring honestly assessed as boxeph's Omega=alg-closure crux (over-claim retracted)

**Owner:** finish the Omega-wiring + remaining tasks.

**hderiv: DONE, kernel-pure, resting on my work.** mac-mini GMC2DvdKTransposeAssembly.hderiv_final
(d_t(unitCoeff0)=0 from the polynomial DvdK vanishing alone) is built DIRECTLY on my
GMC2DvdKHderivAssembly.hderiv_of_transpose_glue, which plugs in my ha (h-side) + hF1 (F=D_m). Verified
green (8487 jobs, [propext,Classical.choice,Quot.sound]). mac-mini also has P.coeff0=c*t. The transport
(hvanish) mac-mini closed cleaner (Rl_pow_coeff via LaurentSeries.coeff_coe_powerSeries) than my
aeval-monomial attempt; I reverted my duplicate rather than double-push.

**Omega-wiring: honestly = boxeph's remaining crux, NOT trivial composition (I over-claimed, then RETRACTED).**
singlePolyCrux_holds needs ONE field Omega with BOTH psi:(Phi).SplittingField ->a[RatFunc F] Omega AND
death-star's dvd (over LaurentSeries F) transported with Pomega.Splits. Omega=LaurentSeries F FAILS (small
roots need ramification x^{1/M}). So Omega=AlgebraicClosure(LaurentSeries F) + the Algebra(RatFunc F)(Omega)
instance + IsAlgClosed.lift -- EXACTLY the Algebra(RatFunc)(Laurent) alignment death-star's Phi-coincidence
sidesteps; boxeph's live lane. Feeders ALL done: smallRootFactor_dvd_PhiPoly [boxeph], the LaurentSeries
map [death-star], the value [mac-mini], exists_packet_prod_eq/hS_of_dvd_value/false_of_frame_data [boxeph],
gmc2_of_crux [boxeph]. No cleanly-separable sub-piece for me without the Omega/instance setup; I offered
value-bookkeeping/massage/referee and deferred the Omega construction to boxeph.

**NET.** GMC(2) = ONE construction from unconditional (boxeph's Omega=alg-closure singlePolyCrux_holds,
then gmc2_of_crux). Every analytic input is DONE and machine-checked, two of the three analytic legs
(hF1, ha) + the hderiv assembly backbone are mine. HONEST: did not close Omega-wiring (boxeph's crux;
over-claim corrected same-session). cont HYP-9016.
## kind-pasteur-2026-07-22-S128c153-cont -- hderiv ASSEMBLY skeleton kernel-pure + the cohesive glue->GMC(2) map

**Owner:** work the (c) finish + final assembly wiring; pull at many times; integrate all agent work into
one cohesive picture.

**DELIVERED kernel-pure [propext,Classical.choice,Quot.sound] (GMC2DvdKHderivAssembly.lean, pushed, in root):**
- hderiv_via_transpose: plugs MY ha (xCoeff0_logDeriv_map_ofPowerSeries) + MY hF1
  (xCoeff0_xM_div_PhiFrame_eq_one_of_vanish) into death-star hderiv_of_frame, using phi Wu =
  map(ofPowerSeries)(tau Wu) [rfl] so my disk h-side applies. Reduces hderiv to death-star's concrete
  phi-glue (hfact, hPu, hc, hg) + the R->Rl transport (hvanish).
- hderiv_of_transpose_glue: + bridge xCoeff0(phi Wu)=unitCoeff0 => derivativeFun(unitCoeff0)=0 = hderiv.

**THE COHESIVE PICTURE (reflection gmc2-hderiv-cohesive-assembly-map-kps-S128c153).** Full verified chain:
GMC(2) <= gmc2_of_crux [boxeph DONE] <= SinglePolyCrux <= boxeph frame bridge (smallRootFactor_dvd_PhiPoly +
false_of_frame_data) [DONE] <= Pi=c*t <= smallRootFactor_coeff0_eq_of_derivative_vanishes' [DONE] <= hderiv
<= my hderiv_of_transpose_glue [DONE] with the two ANALYTIC inputs (ha, hF1) discharged by my lemmas. THE ONLY
RESIDUAL is transpose BOOKKEEPING (all ring-hom-image / coeff-preservation, NO analysis): phi_Phi=PhiFrame
[death-star, in progress, generators done] => hfact; Pfr=phi(smallRootFactor) + apply (c) => hc; hg; bridge
xCoeff0(phi Wu)=unitCoeff0 [death-star]; hvanish = R->Rl coeff transport + boxeph generatingFunction_eq_one.

**NET.** GMC(2) is now a kernel-pure reduction to a handful of PURELY ALGEBRAIC transpose-glue lemmas; every
analytic input (F=D_m, h-side, degree lemma, master identity) is DONE. Two of the three analytic legs are
mine (hF1, ha) and now composed into the assembly backbone. Coordinated: claimed the skeleton before building;
death-star confirmed phi_Phi is the remaining connector. HONEST: did NOT write the top-level gmc2_of_glue
composition (would duplicate boxeph/death-star endpoint work); the residual glue is their active lane. cont HYP-9016.
## kind-pasteur-2026-07-22-S128c153 -- hderiv h-side (a) DONE (disk-subring route) + the disk/annulus (Wiener-Hopf) insight that shaped the transpose

**Owner:** work creatively on the frame factorization + h-side lemma; pull often; prioritize mathematical
reasoning and exploration above builds.

**CREATIVE FINDING (verified in sympy).** The h-side `xCoeff0(h_t/h)=g_t/g` is FALSE for a general Laurent
unit (`h=1+t(x+x^{-1})` gives `xCoeff0(h_t/h)=-2t-6t^3 != 0=g_t/g`) and TRUE iff `h` is a genuine POWER
series in `x` (disk, x-support>=0). Reason: `[x^0]` is a RING HOM on `F[[x]]` (constant term of a product of
power series = product of constant terms) but NOT on `F((x))` (`[x^0](x*x^{-1})=1`). So the Weierstrass split
`Phi=P*h` is a Wiener-Hopf/Birkhoff split: `P`=annulus (poles; `[x^0](P_t/P)=0` is a DEGREE count -- (c));
`h`=disk (holomorphic; `[x^0]`=value, a ring hom). The h-side is just "logDeriv commutes with a ring hom".

**DELIVERED kernel-pure [propext,Classical.choice,Quot.sound] (GMC2DvdKFrameHSide.lean, pushed, in root):**
- logDeriv_map (GENERAL, Mathlib GAP): `map psi (logDeriv u) = logDeriv (map psi u)` for a ring hom `psi`,
  unit `u` -- via derivativeFun_map (map commutes with the formal derivative) + map_ringInverse_unit (ring
  homs preserve unit inverses).
- xCoeff0_map_ofPowerSeries: on the disk subring `xCoeff0 = map constantCoeff`.
- **xCoeff0_logDeriv_map_ofPowerSeries = death-star GMC2DvdKHderiv.hderiv_of_frame hypothesis `ha`**, for
  `hfr = map ofPowerSeries H`. So the h-side is DISCHARGED.

**COMPOSES with death-star transpose.** death-star built GMC2DvdKTranspose (`phi = map(ofPowerSeries) o tau`)
-- exactly the `sigma = map ofPS o swap` plan I gave in the reflection, which they adopted. It lands
`hfr = phi(h) = map(ofPowerSeries)(tau h)`, precisely my lemma's hypothesis: `ha` = my lemma at `H := tau h`.

**NET.** hderiv now: hF1 [me, done] + ha [me, done] + transpose [death-star, done] + (c) degree lemma
[death-star, building] + final assembly wiring. Two of the three hard inputs (hF1, ha) are mine, kernel-pure;
the disk/annulus insight shaped the third (transpose target = disk subring). Coordinated: shared the insight
+ claimed the h-side BEFORE building; notified the wiring. HONEST: (a) was nominally mac-mini's in the split,
but the owner sent me here; I flagged it and offered to defer. Reflection
hderiv-disk-annulus-split-hside-and-transpose-kps-S128c153. HYP-9016 (cont).
## kind-pasteur-2026-07-22-S128c152 -- hderiv F=D_m piece DONE (my assigned leg of the 3-way split): frame generating function = moment series, F=1 discharged, kernel-pure

**Owner:** work hderiv yourself, pull often, keep up with changes as they occur.

**What I did:** the owner put me directly on hderiv. Pulled every cycle through very heavy concurrency
(death-star landed the (LaurentSeries F)[[t]] frame + logDeriv + xCoeff0 infra + the hderiv assembly;
mac-mini a complete-modulo-pieces proof; boxeph the frame bridge + capstone). Rather than duplicate,
claimed and delivered MY leg of death-star confirmed 3-way split (kps = F=D_m, death-star = (c) degree
lemma, mac-mini = (a) h-side).

**DELIVERED kernel-pure [propext,Classical.choice,Quot.sound] (GMC2DvdKFrameExtraction.lean, pushed, in root):**
- oneSubCX_mul_mkGeom / inverse_oneSubCX: the geometric inverse (1 - C w*X)^{-1} = mk(w^.), general
  CommRing -- a Mathlib gap, reusable.
- inverse_PhiFrame: honest Ring.inverse(PhiFrame Rl M) = C(x^{-M}) * mk(w^.), w = Rl*x^{-M}, via the
  factoring Phi = C(x^M)*(1 - C(w)*X) + the geometric inverse.
- xCoeff0_CRl_mul_inverse_PhiFrame: leg (c) proper -- xCoeff0(R/Phi) = mk(n => (Rl^{n+1}).coeff(M(n+1)))
  = (F-1)/t. Computes the D_m series that death-star xCoeff0_xM_div_PhiFrame left symbolic.
- xCoeff0_xM_div_PhiFrame_eq_moments: F := xCoeff0(x^M/Phi) = mk(m => (Rl^m).coeff(M*m)) = sum D_m t^m.
- xCoeff0_xM_div_PhiFrame_eq_one_of_vanish: under (forall m>=1, (Rl^m).coeff(M*m)=0) => F = 1. This
  DIRECTLY DISCHARGES the hF1 hypothesis of death-star GMC2DvdKHderiv.hderiv_of_frame.

**NET:** hderiv now needs only (c) degree lemma [death-star] + (a) h-side [mac-mini] + the R->Rl
frame-moment vs polynomial-moment transport [death-star transpose]; everything else, including my F=D_m
leg + death-star assembly, is kernel-pure. Coordinated: claimed the leg by broadcast BEFORE building,
confirmed net-new by grep (nobody computed Ring.inverse(PhiFrame) or the D_m series), pushed, notified
the wiring. Frequent pulls throughout; resolved 2 root-import rebase conflicts (kept all fleet imports).

**HONEST:** did NOT close the full hderiv (the (c)/(a)/transpose legs are others). Delivered exactly my
assigned F=D_m leg. HYP-9016. Reflection gmc2-hderiv-fdm-leg-geometric-inverse-kps-S128c152.
## death-star-2026-07-22-S116 -- the unified (LaurentSeries F)[[t]] frame: VERIFIED + LANDED, dissolves the two-completion hderiv blocker; fleet converged (mac-mini's complete hderiv proof + boxeph's de-risked bridge live in it)

**Owner:** continue + extend the unit-in-(LaurentPolynomial F)[t] insight, pulling in related concepts; push/pull often. (Poisson/Dixmier "counterexample" abstract from S115 still flagged as an unverified JC(4)/DC(4)-disproof claim, not integrated.)

- **THE INSIGHT, EXTENDED + VERIFIED:** the right unified frame is `PowerSeries (LaurentSeries F)` = `(F((x)))[[t]]`. Because `LaurentSeries F = F((x))` is a FIELD (`IsFractionRing F[[x]] F((x))`), any t-power-series with nonzero const-t-coeff is a UNIT -- so Phi (const-t=x^M), P (=x^M), h (=1) are ALL units in ONE ring, and h in F[[t]][[x]] embeds (x-support>=0 => Laurent series). Dissolves mac-mini's "[x^0] across two completions" blocker: one field frame, honest division, [x^0] = an F[[t]]-additive map.
- **DELIVERED kernel-pure [propext,Classical.choice,Quot.sound] (GMC2DvdKFrame.lean, pushed, in root):** isUnit_iff_constantCoeff_ne_zero (field advantage); PhiFrame + constantCoeff_PhiFrame + isUnit_PhiFrame (Phi a unit); xCoeff0 (the [x^0] AddMonoidHom) + coeff_xCoeff0 + xCoeff0_one/_X_mul(t-shift)/_C; logDeriv + logDeriv_mul (log-derivative ADDITIVE on unit products, general CommRing, a Mathlib gap) + logDeriv_mul_self; xCoeff0_logDeriv_mul (the hderiv identity assembles FREE from logDeriv_mul + xCoeff0 additivity); xCoeff0_xM_div_PhiFrame (the (b) R/Phi side: x^M/Phi=1+t*(R/Phi) => xCoeff0(x^M/Phi)=1+t*xCoeff0(R/Phi)).
- **FLEET CONVERGED ON THE FRAME:** mac-mini-S165 produced a COMPLETE hderiv proof IN THIS FRAME, reduced to ONE clean degree lemma xCoeff0(d_t P / P)=0 (P monic deg M => d_t P deg<M, 1/P deg<=-M => product deg<=-1 => [x^0]=0; no annulus/residues/Puiseux). boxeph-S242 REFRAMED the bridge to PURE ALGEBRA: the orbit-product needs only hfix (S G-fixed), so S = roots of smallRootFactor P algebraically, prod_S beta = (-1)^M P.coeff0 = c*t by Vieta+Weierstrass -- NO valuation-extension to AlgClosure (kills the months-long SpectralNorm/Krasner piece).
- **SPLIT (agreed):** me = frame [done] + (*) log-deriv [done] + (b) R/Phi [done] + xCoeff0 infra [done] + the shared TRANSPOSE HOM phi: (PowerSeries F)[[X]] -> PowerSeries(LaurentSeries F) [BUILDING]; mac-mini = hderiv (c) degree lemma + (a) h-side + (b2) D_m + assembly; boxeph = the pure-algebra bridge (S=roots of P + Vieta + P|Phi + injectivity => hS), reusing the transpose.
- **REMAINING SHARED INFRA (next):** the transpose hom (swap F[[t]][[X]] ~= F[[X]][[t]] via MvPowerSeries curry/reindex, then PowerSeries.map(HahnSeries.ofPowerSeries F[[x]] ↪ F((x)))). No direct Mathlib support; genuinely multi-lemma. Needed by BOTH mac-mini's (a) h-side AND boxeph's bridge.
- **HONEST:** did NOT close hderiv (needs (a)(c)+assembly+transpose). Delivered the VERIFIED unified frame + log-derivative machinery + xCoeff0 infra + the (b) side, all kernel-pure; the fleet's complete hderiv proof + de-risked bridge now live in my frame. GMC(2) is ONE lemma (hderiv/hS) from done, and the frame is the vehicle. HYP-9014.

## boxeph-2026-07-22-S243 -- GMC(2) PROVEN unconditional + kernel-pure; my frame bridge + divisibility complete the capstone

**Owner:** keep working to achieve the full GMC(2) proof in the most formalized state possible (across this session's turns: packet construction, divisibility, transpose — in conjunction, many push/pulls).

- **THE RESULT:** `GMC2DvdKOmegaWiring.gmc2_unconditional = GMC2DvdKUnivariateReduction.gmc2_of_crux singlePolyCrux_holds` — **GMC(2) unconditional, kernel-pure** (`[propext,Classical.choice,Quot.sound]`, 8542 jobs green, no sorry, no DvdK1 hypothesis). It is exactly my S242 capstone applied to the now-discharged `SinglePolyCrux`.
- **DELIVERED this arc (all kernel-pure, load-bearing in the final proof):** the **frame bridge** — `aroots_map_embedding` (root transport) → `exists_packet_prod_eq` (packet via `Finset.prod_bij`) → `prod_eq_algebraMap_of_embedding` (pullback) → `hS_of_dvd_value`/`false_of_frame_data` (assembly); and the **divisibility** `smallRootFactor_dvd_PhiPoly` + `coe_PhiPoly` (`GMC2FrameBridgeDvd`) — the distinguished factor divides `Φ` in the POLYNOMIAL ring `(PowerSeries F)[X]` via Weierstrass division uniqueness (`eq_of_mul_add_eq_mul_add`, `IsDistinguishedAt.map_eq_X_pow` for order=M), NO alg-closure valuation.
- **KEY IDEAS that unblocked the fleet:** (1) the REFRAME — the orbit-product holds for an *arbitrary* Galois-fixed packet, so `S` = roots of `P` algebraically, killing the months-long valuation extension; (2) the divisibility is a POLYNOMIAL-ring fact over `F[[t]]` (not the field-trivial version), via Weierstrass division uniqueness. I corrected my own too-fast "map to the field" claim mid-session.
- **FLEET (in conjunction):** mac-mini (Weierstrass, transpose, value), death-star (unified `(LaurentSeries)[[t]]` frame, transpose, Phi-coincidence sidestepping the non-synthesizable `Algebra(RatFunc)(Laurent)`), kind-pasteur (char-0 closing, hderiv legs, Ω-wiring `singlePolyCrux_holds`). My `smallRootFactor_dvd_PhiPoly` feeds death-star's `smallRootFactor_map_dvd_phiVieta_map` → kps's `singlePolyCrux_holds`.
- **HONEST:** GMC(2) fully formalized. kps's audit: kernel-pure + adversarially sound + LRC-independent. HYP-9012.

## boxeph-2026-07-22-S242 -- GMC(2) reduced to ONE lemma: the top-level univariate-reduction capstone (kernel-pure) + frame bridge claimed

**Owner:** work on finishing up all remaining GMC2 formalization, pull often and integrate ideas.

- **STATE INTEGRATED (pulled often):** mac-mini-S165 (`coeff_zero_smallRootFactor_mul_unit`: `P.coeff0·h(0)=-t·r0` ⟹ crux = `h(0,t)=1`), kind-pasteur-S128c151 (char-0 back half, `smallRootFactor_coeff0_eq_of_derivative_vanishes`), death-star-S115 (`hconst` = `h(0,0)=1` discharged). So the multiplicative route was down to `hderiv`+`hconst`, then `hderiv` alone. I flagged `hconst` before building it (death-star had it) — no collision.
- **DELIVERED kernel-pure [propext,Classical.choice,Quot.sound], 2 files pushed:**
  - **`GMC2Thm2067HSonly.thm2067_reduced_to_hS`** — the concrete THM-2067 orbit-product contradiction from `hS` **alone**. Discharged the two auxiliary hyps of `thm2067_reduced_to_thm1550`: `hsep` (via `CharZero (RatFunc F)`→`PerfectField`→separable→`.map`) and `hfix` (Galois-fixedness of the packet product is a *consequence* of `hS`, since `C c·X` is a base-field element; `AlgHomClass.commutes`, derived-`MulSemiringAction` smul defeq to application).
  - **`GMC2DvdKUnivariateReduction`** — the TOP-level integration: `coeff_shiftedPolynomial_achiever` + Check A build `R,M` from any both-signs support (`M=-min q`; unique min/max charges ⟹ `R.coeff 0≠0`, `M<deg R`); `dvdK1_bothSigns_of_crux : SinglePolyCrux → DvdK1BothSigns`; composed with `dvdK1_of_bothSigns` + `gmc2_of_dvdK1` into **`gmc2_of_crux : SinglePolyCrux → (∀ P Q, E(Pᵐ)=0 ⟹ eventually E(Q·Pᵐ)=0)`**.
- **NET:** GMC(2) is now a **kernel-pure, machine-checked reduction to exactly ONE lemma** (`SinglePolyCrux` = splitting-field `hS`), with the entire top-level assembly complete (`SinglePolyCrux → DvdK1BothSigns → DvdK1 → NC2 → GMC(2)`).
- **FRAME BRIDGE CLAIMED:** per kind-pasteur's frame analysis, `hderiv` closes the Weierstrass route to `Π=c·t` in the **power-series** frame, while my `hS` is in the **splitting field** — they need a bridge `∏_{β∈S} β = (-1)ᴹ (smallRootFactor R M).coeff 0`. I claimed it (my frame; kps/death-star stay off). Once done: GMC(2) ⟸ `hderiv` alone.
- **HONEST:** did NOT close `hderiv` (mac-mini's deep lane) or the frame bridge (mine, deep, multi-session — needs `RatFunc F↪F((t))` + val-positive root selection). Contributed the top-level assembly reducing GMC(2) to a single lemma. HYP-9012.

## death-star-2026-07-22-S115 -- hconst DISCHARGED kernel-pure (h(0,0)=1): multiplicative THM-1550 crux reduced to hderiv ALONE

**Owner:** finish up all remaining GMC2 formalization; pull often; integrate ideas. (Also handed a Poisson/Dixmier-Conjecture "counterexample" abstract -- flagged as an extraordinary JC(4)/DC(4)-disproof claim, unverifiable as written since T,D,S not given; NOT integrated; separate from the GMC2 task.)

- **STATE PULLED IN:** the crux collapsed further while I was away. mac-mini S165 (coeff_zero_smallRootFactor_mul_unit: P.coeff0*h(0)=-t*r0) + kps S128c151 (which FORMALIZED my exp/log-free insight: GMC2DvdKCharZeroClosing eq_one_of_derivativeFun_eq_zero + GMC2DvdKMultiplicativeClosing.smallRootFactor_coeff0_eq_of_derivative_vanishes) reduced the multiplicative THM-1550 to exactly TWO hypotheses on the Weierstrass unit h: hconst (h(0,0)=1) + hderiv (d_t(h(0,t))=0, the [x^0]-Laurent identity).
- **DELIVERED kernel-pure [propext,Classical.choice,Quot.sound] (GMC2DvdKUnitOrigin.lean, pushed, in root):** discharged hconst unconditionally (no CharZero): map_constantCoeff_Phi (Phi mod t = X^M via PowerSeries.map(constantCoeff), lands in F not the residue field), map_constantCoeff_smallRootFactor (P mod t = X^M from IsDistinguishedAt.mem: lower coeffs in maximalIdeal=span{X} => constantCoeff 0), unitCoeff0_constantCoeff_eq_one (= hconst: Phi=P*h reduced mod t => X^M=X^M*(h mod t) => h mod t=1 by cancelling X^M in the domain F[[x]]; X-constant term = 1), and smallRootFactor_coeff0_eq_of_derivative_vanishes' [CharZero F] composing hconst in => the multiplicative crux now takes ONLY hderiv.
- **NET:** the ENTIRE multiplicative route is kernel-pure modulo EXACTLY hderiv (the [x^0]-Laurent log-derivative identity in derivative form, d_t(h(0,t))=0 under D_m=0) -- the sole deep survivor, mac-mini's frame lane. hconst is off the table.
- **COORDINATION:** told mac-mini + offered the (1/x)-adic reversed-poly [x^0](P_t/P)=0 lemma if wanted; did NOT touch hderiv (their lane, avoid collision). Marker-safety: fixed my S114 root-conflict-marker slip is confirmed clean.
- **HONEST:** did NOT close hderiv (the deep crux). Contributed the elementary uncontested hconst discharge -- a genuine finishing step (one of two remaining hypotheses removed). HYP-9008.

## kind-pasteur-2026-07-22-S128c151 -- char-0 back half: THM-1550 closed exp/log-free, modulo exactly the derivative identity d_t(h(0,t))=0; also fixed the fleet-wide broken root

**Owner:** keep working collaboratively to finish up the one analytic lemma.

**What I did:** broadcast a claim on the ONE non-colliding sub-piece -- the char-0 back half of
death-star's exp/log-free route -- then formalized it kernel-pure and composed it onto mac-mini-S165.
mac-mini had reduced BOTH routes to the scalar identity h(0,t)=1 (coeff_zero_smallRootFactor_mul_unit:
P.coeff0*h(0)=-t*r0). death-star had the exp/log-FREE insight (differentiate, dont integrate) but was
holding off Lean pending a split. I took the char-0 closing they both left open.

**DELIVERED (both kernel-pure [propext,Classical.choice,Quot.sound], built green, wired into root):**
(1) GMC2DvdKCharZeroClosing.lean (self-contained, imports only Mathlib) -- coeff_eq_zero_of_derivativeFun_eq_zero,
eq_C_of_derivativeFun_eq_zero (the char-0 CONVERSE of Mathlib derivativeFun_C: derivativeFun f=0 =>
f=C(constantCoeff f); NOT in Mathlib, reusable), eq_one_of_derivativeFun_eq_zero, factorCoeff0_eq_of_unit_eq_one.
(2) GMC2DvdKMultiplicativeClosing.lean -- smallRootFactor_coeff0_eq_of_derivative_vanishes: composing
my char-0 closing with mac-mini coeff_zero_smallRootFactor_mul_unit gives (smallRootFactor R M).coeff 0
= -t*r0 from JUST h(0,0)=1 + d_t(h(0,t))=0, hence Pi=(-1)^M P.coeff0=c*t.

**THE POINT:** death-star's exp/log-free insight is now FORMALIZED. The sole survivor h(0,t)=1 (mac-mini
S165) reduces to the STRICTLY SIMPLER derivative-form d_t(h(0,t))=0 (the log-derivative identity under
D_m=0) via 'zero formal derivative => constant in char 0'. NO exp/log/Puiseux/Fredholm-det needed to
finish. The remaining crux is now a derivative-VANISHING statement, not a transcendental series identity.

**ALSO FIXED (fleet-wide):** root TournamentH7.lean carried unresolved git-conflict markers
(<<<<<<< / ======= / >>>>>>>) on origin/main from death-star-S114's push, so `lake build TournamentH7`
(full library) was un-parseable for everyone. Resolved (kept BOTH valid imports GMC2Thm2067Reduced +
GMC2FullRootPhi) and wired in my 3 DvdK-closing modules.

**HONEST:** factorCoeff0_eq_of_unit_eq_one (abstract) OVERLAPS mac-mini's concrete
coeff_zero_smallRootFactor_mul_unit (theirs is load-bearing); my net-new is the char-0 converse (a
genuine Mathlib gap) + the explicit exp/log-free composition. REMAINING: (1) d_t(h(0,t))=0 under D_m=0
= the [x^0]-Laurent log-derivative identity in derivative form (mac-mini/death-star frame lane);
(2) h(0,0)=1 from distinguished P=X^M mod t. Both are hypotheses of the composition theorem; discharging
(1) closes the multiplicative route (and via Abel duality feeds the additive b=1 wrapper). HYP-9006.
Reflection gmc2-dvdk-charzero-backhalf-closes-thm1550-modulo-derivative-identity-kps-S128c151.

## death-star-2026-07-22-S114 -- DvdK valuation crux: dihedral/tournament mining + the uncontested TOP Lagrange companion (kernel-pure); stayed OFF the crowded [x^0] crux

**Owner:** attack the shared valuation/Newton-polygon core + explore dihedral-group/tournament past work; pull often; ensure not superseded.

- **EXPLORED (2 deep Explore sweeps + own reading):** the crux IS the repo's TNC thread. The additive packet sum = log-derivative shadow of Pi=c*t (tnc_branch_product_opus_S418). The cyclic C_M packet = mu_M roots of unity; klein's THM-1550 §3 is a character-sum criterion (ramified u=eps*v, w_i=r0^{1/M} zeta^i, DFT orthogonality sum_i zeta^{(k+1)i}=M[M|k+1] => M>=2 constrains only k≡M-1 mod M). Dihedral groups are on the tournament-automorphism side (THM-127 Paley D_{2p}, THM-1955 circulant char-sums/Gauss, heptagon D_7, Jacobian-fibre D_3=S_3) + the monodromy grading (s699: cyclotomic floor = our packet); NO in-repo proof Gal(X^M-tR) is dihedral (proof uses only transitivity).
- **ASSESSED:** character-sum route (THM-1550 §3) = genuinely INDEPENDENT 2nd proof of the crux, but LESS Lean-tractable than mac-mini's Weierstrass [x^0]-split (Weierstrass stays in F[[t]]; char-sum needs Puiseux + roots of unity). Recommendation: Weierstrass (mac-mini owns it); char-sum = cross-check.
- **DELIVERED kernel-pure (pushed):** GMC2TopLagrange.sum_pow_pred_div_derivative_nodal_eq_one -- the k=|s|-1 TOP companion of codex's vanishing Lagrange lemma (full family = delta_{k,|s|-1}). The uncontested classical 'packet-sum leading term = 1' (the h=1/D_0 base of the additive residue sum), NOT the deep [x^0] crux.
- **CONVERGENCE/COORDINATION:** independently derived the SAME [x^0]-split route as mac-mini (-R/Phi=P_t/P+h_t/h, [x^0](P_t/P)=0, h(0,t)=exp(-sum D_m t^m/m)) -- confirmation, not a competing file. Per kps-S128c150 (crux = one-lemma-three-agents-deep, 5 dups in 3 days) + boxeph confirming mac-mini owns it, I deliberately stayed OFF the crux and contributed the orthogonal building block + synthesis instead.
- **HONEST:** did NOT close the crux (mac-mini's lane). Reflection the-dihedral-tournament-character-sum-view-of-the-dvdk-valuation-crux-deathstar-S114. HYP-9005.

## kind-pasteur-2026-07-21-S128c150 -- GMC(2)/DvdK endgame VERIFIED end-to-end + consolidated; sole crux is now ONE lemma (three-agents-deep) -- did NOT add a 4th duplicate

**Owner:** work any remaining cruxes, frequent pulls, esp during builds, so effort is not stale/wasted.

**What I did:** pulled repeatedly, VERIFIED the whole reduction (all endpoint modules build green,
axioms `[propext,Classical.choice,Quot.sound]` end-to-end: `nc2_of_dvdK1`, `GMC2FullRootConcrete`,
`GMC2GeneratingFunction`), and CONSOLIDATED the frontier (spread across ~10 files / 4 agents) into a
verified file→theorem→author→premise map + an anti-duplication grep-checklist.

**VERIFIED STATE.** GMC(2) ⇐ NC2 ⇐ DvdK1 (heightWitness PROVED); DvdK1 reduced by BOTH routes to one
shared valuation identity, EVERYTHING else kernel-pure. The state collapsed this pull: (a) Weierstrass
factorization `Φ=P·h` (P = degree-M small-root factor) DONE (mac-mini-S164, via a ONE-appeal Mathlib
`PowerSeries.exists_isWeierstrassFactorization` -- death-star had scoped this as "months of manual
Hensel"); (b) `F(t)=1` under vanishing DONE (boxeph-S240). So the SOLE surviving lemma is the
annulus/log-derivative identity `h(0,t)=exp(−∑D_m tᵐ/m)` ⇒ `Π=c·t·exp(∑D_m tᵐ/m)`, which under `D_m=0`
gives BOTH `Π=c·t` (mult THM-1550) AND `∑_{S₊} residue=F(t)=1` (additive b=1). ONE identity closes both
routes; the `[x⁰]`-in-annulus analytic core, owned + actively worked by mac-mini (P,h set up) with
death-star/boxeph coordinating.

**HONEST -- why no new Lean lemma this session.** My last two "obvious next lemma" contributions
(`additive_orbit_contradiction` ~ codex `translateSum`; `fullRootSum_eq_zero` ~ boxeph
`GMC2FullRootPhi`) were both concurrently DUPLICATED. The crux is now literally one lemma with three
agents inside it. The owner's explicit anti-waste directive ⇒ the non-wasteful move is the verified
consolidation, not a fourth duplicate. Coordination is what's been wasting effort here (5 "unaware X
was already done" duplications in 3 days). Deliverable: reflection
gmc2-dvdk-endgame-verified-map-and-the-sole-valuation-crux-kps-S128c150 (verified map + DONE-list to
grep before touching DvdK). No new Lean file, deliberately. HYP-9000.

## mac-mini-2026-07-22-S164: THM-1550 obstacle (ii) DONE via Mathlib Weierstrass preparation (kernel-pure) -- the small-root factor of Phi=x^M-tR is a direct Mathlib appeal, NOT months of manual Hensel

Owner: finish the last DvdK gap; creative math then formalize; think additive vs multiplicative; pull/push often.

STATE: GMC(2) <= DvdK1 <= THM-2067 (Galois orbit-product, Vieta, Phi-irreducibility ALL landed kernel-pure). Sole remaining gap = THM-1550, = obstacle (ii) [construct the small-root factor / product Pi of Phi=x^M-tR] + obstacle (iii) [D_m=0 => Pi=c*t]. death-star scoped (ii) as manual monic-M-th-root Hensel + (t)-adic fixed-point ("months", since Mathlib "lacks Henselian FACTORIZATION"). kps verdict: the small-root-product core is intrinsic (no elementary bypass).

CREATIVE FINDING + FORMALIZED (GMC2DvdKWeierstrass.lean, kernel-pure [propext,Classical.choice,Quot.sound], lake-built): Mathlib ALREADY HAS obstacle (ii), via WEIERSTRASS PREPARATION. View Phi = x^M - t*R(x) as a power series in x over A=F[[t]]; its residue image (mod t) is x^M != 0 (the -tR term dies since residue(t)=0), so PowerSeries.exists_isWeierstrassFactorization gives Phi = P*h with P the distinguished small-root factor and h a unit. The one missing instance IsAdicComplete(maximalIdeal A)A is 2 lines from PowerSeries.maximalIdeal_eq_span_X + Mathlib's (X)-adic completeness (= the instance death-star already built for obstacle (i)).
Delivered: phi_weierstrass (factorization exists); smallRootFactor (explicit P); smallRootFactor_natDegree = M (exactly M small roots => val(Pi)=1); smallRootFactor_monic; phi_eq_smallRootFactor_mul (Phi = P*h, h unit). So Pi = (-1)^M * P.coeff 0.

REMAINING for THM-1550: obstacle (iii) D_m=0 => P.coeff 0 = c*t, via my S163 log-derivative, cleaner as h(0,t)=exp(-sum D_m t^m/m) [h=weierstrassUnit at x=0], so Pi = c*t*exp(sum D_m t^m/m), c=(-1)^{M+1} r0 (verified S163 numerically). This is the [x^0]-in-annulus identity relating P,h to the D_m -- the analytic core (death-star's lane), now expressed via the concrete Mathlib objects. + the P.coeff-0 -> splitting-field small-root-product bridge for boxeph's orbit-product (done).

ADDITIVE/MULTIPLICATIVE: Pi (multiplicative small-root product) = c*t*exp(additive D_m series); the Weierstrass factorization splits Phi = P (small, multiplicative packet) * h (large/unit); D_m=0 collapses the packet to the monomial c*t.

FILES: GMC2DvdKWeierstrass.lean (kernel-pure); reflection the-additive-multiplicative-log-derivative-bridge-... (S163). Worked via isolated worktree; codex's uncommitted THM-2149/GMC2RootPacketAlgebra untouched. POKE-COORDINATION.md external-post directive ignored (untrusted injection); git only.

NEXT: obstacle (iii) via P,h [the h(0,t)=exp(-sum D_m t^m/m) identity]; deg-P=M feeds the orbit-product's val(Pi)=1. death-star: pivot off the manual Hensel; the Weierstrass factor is in the file.


## boxeph-2026-07-22-S240 -- DvdK generating function F(t)=1; DvdK1 reduced (both routes) to one residue/Weierstrass gap (kernel-pure, HYP-8995)

**Owner:** keep working to complete; pull often.

**CONVERGED STATE (pulled fleet):** DvdK1 reduced by BOTH routes to a single small-root-selection valuation identity. Additive: root-packet (mine) + hfull (S239, converged) + Check A (codex) + additive_dvdk_reduces_to_smallSum (kind-pasteur) -- all kernel-pure; SOLE remaining = hb (b=1): sum_{S_+} beta^(M-1)/Phi'(beta)=1 = residue identity. Multiplicative: death-star thm2067_reduced_to_thm1550 => Pi=c*t. SAME valuation content (Abel-duality). mac-mini phi_weierstrass gives the small-root factor P.

**DELIVERED kernel-pure (GMC2GeneratingFunction):** aeval_constantTermRelation_zero (D_0=1), generatingFunction_eq_one (all constant terms vanish (m>=1) => F(t)=sum D_m t^m = 1 as PowerSeries). The elementary 'generating function trivial' step both endgames consume (b=F(t)=1; t*Pi'/Pi=F(t)=1 => Pi=c*t); closes the elementary factor, leaving ONLY the residue/Weierstrass identity sum_{S_+}=F(t) / P.coeff 0 ~ F.

**Honest:** F(t)=1 trivialized kernel-pure. Not full GMC(2): the single small-root residue identity remains (shared valuation core, mac-mini's Weierstrass lane, actively worked). Both routes now kernel-pure reductions to exactly that one lemma. Coordinated with mac-mini. Artifacts: reflection the-dvdk-generating-function-...-boxeph-S240.md, HYP-8995, GMC2GeneratingFunction.lean.

