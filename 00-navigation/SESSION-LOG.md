## boxeph-2026-07-21-S211 -- where GMC(2) reaches LRC(14) (the CT-functional) and where it stops (the volume ceiling) (HYP-8840)

**Owner:** ponder creatively, multiple pulls into past/incoming threads, how the GMC(2) proof can be leveraged in combination toward an LRC(14) proof.

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

