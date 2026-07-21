## boxeph-2026-07-21-S209 -- Orlik-Solomon is a repo-wide pattern; toric arrangements are the native LRC lens (HYP-8830)
> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).
## boxeph-2026-07-21-S206 -- what an LRC(14) disproof MUST be, and why Fibonacci is the foil (HYP-8815)
## codex-2026-07-21-LRC-Hasse -- THM-2043 local completeness, global no-go, and resolved q=41 exit

- **PROVED (THM-2043):** `F_7[C_14]` splits into the two length-seven local
  rings at `+1,-1`; the fourteen Hasse coordinates are complete for a reduced
  period-14 function. This is local completeness only.
- **SHARP NO-GO:** `T_n={1,...,11,13,96+3444n}` has exactly the AP's full
  mod-14 owner-residue/phase/Hasse packet, blockedness for every `q<=13`, and
  `q_threshold=14`, but has the uniform strict certificate
  `(q,a,integer margin)=(41,17,1)`, i.e. minimum `3/41>1/14`. For every fixed
  `k`, a CRT family additionally matches the AP lift height modulo `7^k` while
  retaining the same exit. Raw jets, q-threshold, and any fixed finite height
  precision are therefore insufficient.
- **POSITIVE CARRIER:** exact owner height or an adaptive resolved
  `(q,a,margin)` certificate is faithful. Owner-labelled residues mod 13 and
  14 recover every bounded height `<=181`. All 156 one-lift aliases have exact
  primitive strict certificates by `q<=42`; the eleven named HYP-2979 rows
  keep a seven-row route-mixed fiber at every Hasse depth.
- **UNRELATED-WORK SYNTHESIS:** THM-2000's support-measure/pushforward
  discipline explains the lost lift kernel; THM-1605 orbit products suggest a
  possible global packet transport; THM-2033 confluence suggests a labelled
  boundary-to-jet step; THM-1775 recurrence, THM-1960 regular-seed
  decomposition, THM-1745 leaf filtrations, THM-856 positive initial forms,
  and THM-346 cube transport were audited with their exact blockers. General
  bad-prime template: for `N=Mp^a`, `F_p[C_N]=F_p[X]/((X^M-1)^(p^a))`.
- **HONEST FRONTIER:** LRC(14) remains open. The best routes remain THM-671
  resolved-modulus/B5 supply and familywise Fejer/Toeplitz rigidity. The new
  Hasse result redirects the characteristic-seven lane to an adaptive
  owner-labelled phase sheaf, especially the `{14,27,41}` atlas.
- **MERGE RECONCILIATION:** historical S88--S92 entries below that call
  NC2/GMC(2) open predate or did not contain THM-2022's later whole-face
  Frobenius proof. Canonical status: NC2/GMC(2) is proved on paper by THM-2022;
  full Lean formalization remains incomplete.

## death-star-2026-07-21-S92 -- GMC(2)/NC2 formalization → Mathlib-submission readiness: verified kernel-pure + extracted the GENERAL three-term no-common-root (new to Mathlib, PR-ready). HYP-8805.

**Owner:** explore creatively other repo areas where Orlik-Solomon could be leveraged; think abstractly, find similar structural patterns.

**ABSTRACT SIGNATURE:** Mobius lattice of flats/layers + inclusion-exclusion characteristic polynomial (finite-field point-count) + complement cohomology (OS algebra) + localization-at-a-flat product factorization. Applies wherever an object is "count/measure configurations avoiding a lattice of coincidences."

**FOUR arrangement types the repo meets (all anchors verified, orlik_solomon_across_the_repo_boxeph_S209.py):**
1. BRAID A_{n-1} (tournaments/NC2, S208): OS Poincare pi(t)=prod_{k=1}^{n-1}(1+kt), Betti = unsigned Stirling-1st c(n,n-i), top (n-1)! -- a GRADED cohomology lens on tournament/ordering space, finer than char_A. Per-tournament refinement via graphic sub-arrangements (chromatic poly / acyclic orientations = THM-805 Tutte).
2. SHI / deformed braid {x_i-x_j in {0,1}} (LRC resonance x_i-x_j=integer): chi(q)=q(q-n)^{n-1}, regions (n+1)^{n-1}=parking functions, by finite-field count. LRC resonances ARE a deformed braid arrangement; finite-field point-count = Diophantine = LRC's native flavor.
3. TORIC / De Concini-Procesi (LRC relation lattice {k.v=0}, THM-1820): |G_delta| = int prod 1[||v_j t||>=delta] = sum_{k.v=0} prod ghat(k_j) = toric-COMPLEMENT VOLUME = arithmetic-Mobius sum over LAYERS (verified match <3e-3). The repo's LRCMod-mod-q ladders ARE the finite-field method on the torus. DCP theory is the correct, under-named home.
4. GENERIC (figurate cutting cake/bagel, S207): regions = sum C(n,k), deficit = mu of degenerate flats.

**TIGHT = RELATION-RICHEST:** AP (1,2,3) uniquely maximizes toric-arrangement Betti/Mobius mass N_R (verified; matches THM-1820 B2 up to +/-k convention). The LRC tight extremal = the config whose toric arrangement is richest in low-height layers = the deep degenerate flat where Mobius mass concentrates (= reify-ladder cold vertex).

**LEVERAGE:** (a) tournament cohomology (graded, beyond char_A); (b) finite-field/Ehrhart engine as the systematic tool for LRCMod ladders (Wall A = AP is unique max-layer residue pattern, a point-count); (c) TORIC LAYER-LOCALIZATION: near a resonance layer |G_delta| should factor into lower-rank toric complement x transverse braid factor -- the toric analog of the S208 braid flat-factorization that gave HYP-8775 -- a candidate Wall-A tool.

**Honest:** anchors 1-4 are verified identities/point-counts; the toric layer-localization (the actual new LRC lever) is proposed BY ANALOGY with S208, NOT verified (only |G|=complement-volume=arithmetic-Mobius verified). Value = the unification + naming DCP toric arrangements as the LRC home. Artifacts: reflection orlik-solomon-is-a-repo-wide-pattern-...-boxeph-S209.md, HYP-8830, script (+.out).

