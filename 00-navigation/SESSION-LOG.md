## death-star-2026-07-21-S98 -- NC2 capstone: skeleton of `DvdK1 → NC2` typechecks (architecture validated); full completion plan worked out + reference-channel friction resolved. HYP-8805.

**Owner directive:** finish the capstone concurrently, pull often, get everything Mathlib-ready.

- **codex's GMC2 composition is OPEN** (codex moved to LRC; newest GMC2 commit is still my S97). So no race -- the `DvdK1 → NC2` capstone is genuinely available.
- **NEW `TournamentH7/GMC2NC2Capstone.lean` (WIP, one sorry, NOT imported by the spine so it stays sorry-free):** `nc2_of_dvdk1` skeleton **typechecks** -- the descent destructure (19-field nested existential) + `letI : Field D.ResidueField := D.fieldStructure` + the two-sided contradiction structure. This confirms codex's 33-module API *can* be composed into the conditional NC2.
- **COMPLETION PLAN (worked out, in memory):** `aeval w (normalizedMomentRelationInt exponent (p*m0) (p*A0))` is both (a) `=0` (null ⟹ integral zero relation over ℂ via `aeval_normalized_eq_zero_of_E_pow_eq_zero`, preserved to `w` by `hpreserve`) and (b) `≠0` (`three_case_sum_ne_zero`: `channels=piAntidiag univ (p*m0)`, `dilated={r|∀i,p∣r i}`, `face=`dilation image `map_piAntidiag_dilation`; `hnondilated`/`hoffFace` via char-p bridges A/B, `honFace` via `coefficientProduct_dilate`+`multinomial_dilate_modEq`, `hfaceSum`=seed). Contradiction.
- **FRICTION RESOLVED:** the reference channel `r0`/height `A0` for the height floor `hmin` come from the **w-seed** (`hseednz`) via `GMC2FaceSeedChannel.exists_nonzero_balanced_channel` over the residue field (field-generic) -- the descent need not expose the ℂ seed.
- **REMAINS:** discharge (a)+(b), ~150 lines of codex-API wiring (the honFace dilate-reindex most intricate) -- multi-session. Offered to codex. memory nc2-gmc2-lean-formalization-state updated with the plan. HYP-8805.
> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## boxeph-2026-07-21-S218 -- arithmetic entropy is a repo-wide invariant; the rigid extremum = the zero-entropy point (HYP-8875)

**Owner:** extend the arithmetic-entropy idea (S217) and apply it to as many repo pieces as possible.

## boxeph-2026-07-21-S217 -- class number = arithmetic entropy: hidden binary forms, and why 7 is rigid (HYP-8870)

## codex-2026-07-21-LRC-transverse-decorrelation -- two infinitary branches become finite labelled cells

- **THM-2053 PROVED:** adjacent normalized columns in every positive rational
  two-plane in `Q^13` expose a repeat projection. Besides the norm-`91L`
  geodesic terminal, the specified row obeys
  `M(v)>=1/13-R/(2N)`. Every unresolved row therefore lies in `N<91R`, with
  `N|(v_i+v_j)` and a fixed denominator-`N` residue template up to a unit.
  The exact deck `D_N(m)` is independent of every longitudinal coefficient,
  makes bad moduli a divisibility down-set, and gives a rational floor-count
  conductor; for `(1,-1,2,...,12)` it is `floor(N/13)/N` and the sharp cutoff
  is `156`. This converts THM-2052's stars into finite bad ruler cells; their
  exact Euler/phase-height discharge remains open.
- **THM-2054 PROVED:** bounded scalar-to-vector resonance lifting gives exact
  equality of the finite Fejer constant-term sets, and two whole-product
  telescopes control the unsmoothed error. The complete six-inner-sector
  inclusion--exclusion budget is `384r epsilon_H`; rowwise `H=2^19` is below
  every recorded pinned-base `cap-Q` margin. MISTAKE-080/082 prevents a false
  finish: the lifted full-torus plateau must still be identified and bounded
  for the actual cluster shape.
- **INCOMING NC2 SYNTHESIS:** THM-2022 remains the proof of NC2/GMC(2), and its
  self-contained arithmetic engine in Sections 1, 4, and 5 is now kernel-pure.
  Subsequent modules close the former descent and face-construction interfaces;
  the remaining Lean step is the concrete normalized-channel residue assembly
  and final conditional theorem `DvdK1 -> NC2`.
  The transferable lesson is whole-layer preservation after a selector;
  THM-2041 supplies that preserver for good-characteristic finite-abelian LRC
  packets, but repo counterexamples confirm that it supplies neither the safe
  seed nor the pointwise exit.
- **REFEREE:** the lattice saturation, transverse grid, constants, divisor
  identities, Haar pushforwards, zero-frequency guard, strict alias cutoff,
  and exact rational atom budgets were independently checked. A hostile audit
  caught and repaired the tempting but invalid reuse of `Q(k-1)` as a generic
  plateau.

## death-star-2026-07-21-S96 -- NC2 formalization: §4 channel-survival COMPLETE + §1 balanced-channel form + contrapositive entry; the entire self-contained arithmetic engine of THM-2022 (§1/§4/§5) is now kernel-pure (16 theorems). HYP-8805.

**Owner:** look for even more hidden binary forms; think information theory.

**ONE INVARIANT:** H_arith(X|L)=log2|{X': L(X')=L(X)}| = the GLOBAL bits of X HIDDEN from a LOCAL invariant L. Zero=local-determines-global (RIGID); positive=hidden global object. FOUR instances (verified arithmetic_entropy_across_the_repo_boxeph_S218.py):
1. BINARY FORMS|genus (refines S217): the truly-hidden part is the DEEP within-genus class group (genus is congruence-detectable). h=genera x deep; Heegner -3,-7,-11 h=1 zero; -15=-3*5 pure GENUS (visible, 0 deep); -23,-47 pure DEEP (1.58,2.32 hidden bits = Hilbert class field, invisible to congruences).
2. TOURNAMENTS|score sequence: transitive (0..n-1)=UNIQUE realization (Landau) => H=0 rigid = the AP/nullcone/rank-11 vertex (S214); near-regular scores carry the hidden fiber (n=5: (1,2,2,2,3) 3 classes = kps reconstruction wall).
3. REALS|CF prefix: golden [0;1,1,..] geo-mean 1 << Khinchin = worst-approx = LRC FOIL (S206); t*=14/183=[0;13,14] geo-mean 13.5 = well-approx = extremal.
4. NULLCONE|moment depth = certificate entropy: LRC finite (bounded alphabet, Bonferroni depth ~5) vs GMC infinite (unbounded degree, Watson S211).

**DUAL entropies:** score-DISTRIBUTION entropy (transitive MAX spread) vs RECONSTRUCTION entropy (transitive 0). The AP = max-order + zero-hidden-info = rigidity; the regular/Paley = min spread + max hidden info.

**UNIFYING:** every repo RIGID extremum (AP/transitive/Heegner h=1/reify-ladder vertex) = a zero-arithmetic-entropy point (local determines global); every DIFFICULTY = its positive-entropy hidden object (deep class group / cospectral fiber / CF tail / deep moment). Rigidity = why the extremal is unique; hidden entropy = where the proof still must go.
**WHY 7 IS RIGID:** LRC(14)=2*7 -> disc -7 -> h=1 -> ZERO arithmetic entropy. codex THM-2053's anisotropic gate residual is fully pinned by local S215 Legendre data -- NO hidden bits, no class-group slack. So (1) a counterexample has NOWHERE to hide; (2) the certificate must be the exact local (Euler/chi/Borsuk-Ulam, p=3mod4) one. Rigidity = why 7 is the first hard-but-tractable case (kps-S17). Heegner h=1 imag. quadratics = -3,-4,-7,-8,-11,-19,-43,-67,-163; -7 is LRC(14)'s.
**Owner:** work incoming LRC progress (pull often); explore rank-11 / '11 private-coordinate relations' through 'relations = a tournament'.

**Honest:** genus/class-group, Landau, Khinchin/Levy, detection-depth facts are classical/verified; the contribution is the UNIFICATION (one info deficit across binary forms / tournaments / CF / nullcones) + the rigid=zero-entropy observation. Organizing lens, gate-independent (survives S217 MISTAKE-225), not a proof step. Artifacts: reflection arithmetic-entropy-is-a-repo-wide-invariant-...-boxeph-S218.md, HYP-8875, script (+.out).

