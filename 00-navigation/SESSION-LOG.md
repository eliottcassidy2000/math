## boxeph-2026-07-21-S209 -- Orlik-Solomon is a repo-wide pattern; toric arrangements are the native LRC lens (HYP-8830)

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

