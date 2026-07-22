## codex-2026-07-21-LRC-Fejer-BV -- THM-2051 closes the coarse BONF5 middle

- **PROVED (THM-2051):** if a thirteen-speed row has no exact relation of
  support two through five and coefficient height at most `2^20`, then its
  continuous quintic Bonferroni functional is positive, hence it has a
  positive-measure strict lonely set.
- **MECHANISM:** approximate the centered danger indicator by a degree-`H`
  Fejer polynomial. Relation dissociation kills every finite product constant
  term. A BV translation estimate controls each whole signed centered-product
  tail in `L1`; summing THM-935's exact coefficients costs
  `K=1477008/343`. At `H=2^20`, the error is below
  `31755672/359661911`, leaving the exact positive margin
  `595652076/17623433639` beneath the equilibrium floor `2052/16807`.
- **SYNTHESIS:** this bypasses THM-946's still-open absolute `T4/T5`
  strip/slab estimates for the coarse alternative and converts HYP-8841's
  ambient termination problem into descent on a finite union of bounded-
  height circuit hyperplanes. It does not classify that structured branch;
  the tight AP is deliberately relation-rich. THM-940 is only the discrete
  analogue and is not silently identified with the continuous proof.
## codex-2026-07-21-DC2-LRC14-termination -- local acyclicity versus finite Euler termination (THM-2049/2050, HYP-8841)

**Owner:** push the DC(2)/planar-JC thread to its next decisive target and transfer the proved/disproved mechanisms toward LRC(14), while repeatedly integrating incoming work.

**DC(2) result (THM-2049, PROVED local/formal statement; not DC(2)):** in the exact Ore algebra `Q[x,q][ell;delta]`, `beta(sum a_k ell^k)=min_k(v_x(a_k)-2k)` is multiplicative and commutators raise beta by two. The associated bracket is `{ell,q}_0=2`. For the Weyl boundary symbols, the simultaneous grade-`g` correction map is `(A,B)->(8/3)(u-2)A+(2u^2-10u+9)B/9`; it is surjective because the two `u` polynomials are coprime. Thus the grade-six residual is exact. An exact ladder advances grades `6,...,13` to `14`; a formal beta-adic `[S,T]=1` lift exists. The open gates are polynomial termination and the coupled `D` relations. This corrects HYP-8802/8803's earlier suggestion that the first invariant grade might carry the obstruction.

**LRC no-go (THM-2050, PROVED):** AP13 and `AP13` with `12->26` have identical full local phase-height function germs on `|h|<1/728` at every unit point `a/14`, yet `M=1/14` and `M=1/12`. Local top data, even as a full germ, cannot determine global loneliness.

**Incoming synthesis:** THM-2047 supplies the lossless signed phase-height/Euler carrier; THM-2048 supplies the fiber-quantization pruning tax; HYP-8840 identifies GMC's constant-term/volume leverage and its zero-volume ceiling. Later pulls supplied THM-2048's genuine Cover14 gain and promoted THM-2051: the no-small-relation branch has positive safe volume, so every hard row lies on a support-`2..5`, coefficient-height-`<=2^20` circuit hyperplane. The exact transfer is `volume/tax -> strict branch`, `Euler signed wall word -> tight branch`, and a labelled Noetherian deletion rule inside the circuit branch as the missing glue. No literal algebra map between GMC/DC and LRC is asserted.

**Exact termination-sidecar audit (HYP-8841):** pair-sum maxima, threshold interval/point topology, complete first exits, and every peel tax were computed on AP/GW, `12->26`, `12->36`, `12->96`, `12->84`, P10+K33, and the incoming Cover14 tax-gain row. It exactly reproduces the latter's peel-`93` excess `2413467317/235670635200`. The tax fires on deep/covering controls but misses the smallest hostile/K33 controls and is not a scalar termination height. THM-2047 proves the strict search is complete by `q<=2 max(S)`. With THM-2051 now proved, the next decisive Wall-A clause is owner-labelled endpoint survival inside the bounded small-relation branch when neither a volume-tax violation nor a positive pair-sum margin occurs. The proof-carrier tournament is transitive; signed threshold topology wins and raw unit germs lose.

**Artifacts:** THM-2049, THM-2050, HYP-8841, the updated exact Ore script/output, `lrc14_termination_sidecar_codex_20260721.py/.out`, and reflection `from-Ore-boundary-acyclicity-to-LRC14-Euler-termination-codex-20260721.md`.

## boxeph-2026-07-21-S211 -- where GMC(2) reaches LRC(14) (the CT-functional) and where it stops (the volume ceiling) (HYP-8840)

**Owner:** ponder creatively, multiple pulls into past/incoming threads, how the GMC(2) proof can be leveraged in combination toward an LRC(14) proof.
## death-star-2026-07-21-S94 -- Formalizing THM-2022 (NC2/GMC2): the arithmetic ENGINE is kernel-pure (§1 Wick expansion + §4 multinomial-Lucas + §5 Frobenius + architecture GMC2<=NC2); §2 descent + §3 DvdK remain. HYP-8805.

**Owner directive:** formalize THM-2022, aiming for NC2 and GMC(2) COMPLETELY formalized; pull often for concurrent work.

- **HONEST FRAME:** GMC2Reduction.lean already had the easy direction + `NC2At`/`ChargeOneSided`; the whole remaining problem is `nc2 : ∀P, NC2At P` (= THM-2022). Full sorry-free completion is BLOCKED on §2 (number-field descent, HEAVY) + §3 (Duistermaat-van der Kallen THM-1630, NOT in Mathlib = a citation) -- genuinely multi-session, NOT done this session. What I did: formalize the self-contained arithmetic engine kernel-pure.
- **FORMALIZED kernel-pure** ([propext, Classical.choice, Quot.sound], no sorry/native_decide): (1) **architecture `gmc2_of_nc2`** -- GMC(2) ⟸ NC2, sorry-free (the whole problem rests on the single theorem NC2); (2) **§1 `wick_expansion`** -- E(P^m) = Σ multinomial·∏coeff^k·wt(radial), the exact channel sum M_m, via Mathlib `Finset.sum_pow_eq_sum_piAntidiag` + E-linearity (E_add/E_monomial/E_sum) + prod_monomial; (3) **§5 `sum_natCast_mul_pow_char`** -- (Σ w_s g_s)^p = Σ w_s g_s^p in char p (Frobenius fixes natCast weights ⟹ face survives as Q̄^p); (4) **§4 `multinomial_dilate_modEq`** -- multinomial Lucas (Mathlib had only binomial), assembled via multinomial_insert + Choose.choose_mul_mul_modEq_choose_nat.
- **FIXED codex's `GMC2FrobeniusFace.lean`** (§3-4 face geometry) -- it FAILED to build (linarith nonlinear gap: lambda·charge i vs lambda·charge j; substitute the charge equality first). File not in aggregate build, so break was invisible to `lake build`.
- **Mathlib API survey (verified, recorded in reflection):** Lucas/Kummer/Frobenius/multinomial-theorem/Zariski all located; multinomial-Lucas + DvdK were the only NOT-IN-MATHLIB gaps (multinomial-Lucas now assembled; DvdK stays a citation).
- **REMAINS for complete nc2:** §2 descent (Zariski `finite_of_finite_type_of_isJacobsonRing` + number-field residue Frobenius, HEAVY but Mathlib-stocked); §3 DvdK (cite); §4 no-carry channel-survival wrapper (MODERATE); final contrapositive assembly. reflection formalizing-thm-2022-...-S94. HYP-8805.

## boxeph-2026-07-21-S210 -- corrected antisymmetry atlas (HYP-8835 / MISTAKE-224)

**PULLS:** deep-mined the GMC-LRC bridge (THM-1645 polar bridge, THM-2022 proof chain, THM-1840 seed, THM-730, THM-1017/Wall A, S157 obstruction, cyclotomic Phi_6/Phi_7, THM-1820 dictionary); pulled INCOMING codex THM-2047 (phase-height carrier) + MISTAKE-223 (corrects my S209) + MISTAKE-224 (corrects my S210). Adopted both corrections.

**REACH (verified, ct_functional_bridge_gmc_to_lrc_boxeph_S211.py):** the CONSTANT-TERM functional CT = the ANGULAR half of GMC's polar bridge E=L o CT (THM-1645; angular half DvdK-closed; GMC's open half is RADIAL Laplace, 'Laplace determinacy not tori', perpendicular to LRC per deathstar-S77 -- so the transferable part is GMC's PROVED part). additive m-energy = CT[P^m Pbar^m] = ||P||_{2m}^{2m} (P1; = opus-S153's E_j=support-2j CT-count). Shared seed THM-1840 is a CT statement (P2). AP maximizes additive m-energy at m=2,3,4 (P3, extends THM-730 m=2) => Wall A = CT-moment extremality. AP=cyclotomic; argmax value Phi_6/Eisenstein vs speed-apex hardness Phi_7 (P4).

**CEILING (codex THM-2047; the decisive pull):** the CT-moment/|G_delta| is a VOLUME functional; volume is only a STRICT-EXIT certificate. At the tight threshold G_delta is measure-ZERO but nonempty; chi(G_delta)>0 is the iff criterion. VERIFIED (P7): S={1,2,3}, delta=1/4: |G|->0 but chi=2>0 (lonely at t=1/4,3/4). So GMC/CT/volume CANNOT resolve the tight AP = Wall A. Analytic face: energy necessary-not-sufficient (opus-S181; translation vs dilation invariance; THM-730 open mile).
- **PROVED / CONSOLIDATED (THM-2047):** the exact arrangement is the oriented
  phase-height complex cut by owner walls `v t +/- delta in Z`. Its
  `delta=1/14` slice is nonempty exactly when LRC holds; opposite-sign top
  walls recover the pair-sum rulers; deletion is exact; and
  `chi(G_delta)=#components` detects isolated tight phases that volume misses.
- **NEW LOCAL NORMAL FORM:** at a top vertex only the extreme active rising
  and falling slopes determine the boundary layer. For small `epsilon`, its
  safe-bar length is
  `epsilon(1/s_+ + 1/(-s_-))`. This rules out the proposed Cartesian
  layer-product after restriction to the one-dimensional LRC orbit.
- **UNRELATED DISCREPANCY TRANSFER (THM-2048):** under a zero-measure packet,
  THM-731's peel discrepancy splits exactly as
  `disc_v=6mu^2+||P_v1_E-7mu1_D||_2^2`. Since the fiber average is quantized in
  steps of `1/v`, it pays the extra tax `{7vmu}(1-{7vmu})/(7v^2)`. Combined
  with THM-732, every hypothetical counterexample must satisfy the explicit
  integer obstruction `6(vmu)^2+theta(1-theta)/7<=r_v^2/3` for every peel.
  The primitive Cover14 row
  `{1,8,11,12,14,17,22,26,35,40,54,90,93}`, peeled at `93`, is an exact strict
  gain: the old tail bound has positive slack
  `3001166951/117835317600`, but the quantization tax exceeds it by
  `2413467317/235670635200`, forcing a positive-measure lonely interval.
- **CORRECTION (MISTAKE-223):** S209's Shi counts and braid identities are
  valid classical computations, but `G_delta` is not a standard toric-
  complement volume, the THM-1820 Fourier annihilator is not an arithmetic-
  Mobius layer sum, and its bounded `N_R` triple census is not Betti/Mobius
  mass or a twelve-speed AP theorem. The exact transfer is the signed cyclic
  wall word with height and owner labels.
- **PRIOR ART RETAINED:** HYP-2986, HYP-3025, THM-668, THM-752, and THM-1142
  already contain the threshold tope/Cech/medial/deletion interfaces. The new
  theorem is a general-height synthesis and local/Fejer refinement, not a
  renaming of those results.
- **REPOSITORY REPAIR:** the S209 commit deleted concurrent S93/S207/S208 log
  provenance. Those entries are restored below rather than replaced.

**BRIDGE:** M(S)=max_t min_v is a minimax=SADDLE value (S210, untouched by MISTAKE-224). Tight witness = measure-zero critical/saddle phase (pair-sum q|v_i+v_j), where chi=#components (finite: bounded LRC alphabet, klein-S389) is the closing invariant.

**PICTURE:** GMC/CT controls VOLUME (strict loneliness) + shared seed; codex phase-height controls chi (tight boundary); minimax/saddle bridges them; the residual is pinned: force chi(G_{1/14}(C))>0 for every 13-speed covering core = Wall A (HYP-7310), now topological not volumetric.

**Honest:** a verified DECOMPOSITION of the remaining work + a proven ceiling on the GMC/volume half, NOT a proof or new implication. Corrects+adopts codex THM-2047; retracts the S209 toric-complement and S210 'torus needs saddles' overclaims. Artifacts: reflection where-gmc2-reaches-lrc14-...-boxeph-S211.md, HYP-8840, script (+.out).

