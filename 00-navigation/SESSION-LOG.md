## death-star-2026-07-21-S98 -- NC2 capstone: skeleton of `DvdK1 → NC2` typechecks (architecture validated); full completion plan worked out + reference-channel friction resolved. HYP-8805.

**Owner directive:** finish the capstone concurrently, pull often, get everything Mathlib-ready.

- **codex's GMC2 composition is OPEN** (codex moved to LRC; newest GMC2 commit is still my S97). So no race -- the `DvdK1 → NC2` capstone is genuinely available.
- **NEW `TournamentH7/GMC2NC2Capstone.lean` (WIP, one sorry, NOT imported by the spine so it stays sorry-free):** `nc2_of_dvdk1` skeleton **typechecks** -- the descent destructure (19-field nested existential) + `letI : Field D.ResidueField := D.fieldStructure` + the two-sided contradiction structure. This confirms codex's 33-module API *can* be composed into the conditional NC2.
- **COMPLETION PLAN (worked out, in memory):** `aeval w (normalizedMomentRelationInt exponent (p*m0) (p*A0))` is both (a) `=0` (null ⟹ integral zero relation over ℂ via `aeval_normalized_eq_zero_of_E_pow_eq_zero`, preserved to `w` by `hpreserve`) and (b) `≠0` (`three_case_sum_ne_zero`: `channels=piAntidiag univ (p*m0)`, `dilated={r|∀i,p∣r i}`, `face=`dilation image `map_piAntidiag_dilation`; `hnondilated`/`hoffFace` via char-p bridges A/B, `honFace` via `coefficientProduct_dilate`+`multinomial_dilate_modEq`, `hfaceSum`=seed). Contradiction.
- **FRICTION RESOLVED:** the reference channel `r0`/height `A0` for the height floor `hmin` come from the **w-seed** (`hseednz`) via `GMC2FaceSeedChannel.exists_nonzero_balanced_channel` over the residue field (field-generic) -- the descent need not expose the ℂ seed.
- **REMAINS:** discharge (a)+(b), ~150 lines of codex-API wiring (the honFace dilate-reindex most intricate) -- multi-session. Offered to codex. memory nc2-gmc2-lean-formalization-state updated with the plan. HYP-8805.

## boxeph-2026-07-21-S218 -- arithmetic entropy is a repo-wide invariant; the rigid extremum = the zero-entropy point (HYP-8875)

**Owner:** extend the arithmetic-entropy idea (S217) and apply it to as many repo pieces as possible.

**ONE INVARIANT:** H_arith(X|L)=log2|{X': L(X')=L(X)}| = the GLOBAL bits of X HIDDEN from a LOCAL invariant L. Zero=local-determines-global (RIGID); positive=hidden global object. FOUR instances (verified arithmetic_entropy_across_the_repo_boxeph_S218.py):
1. BINARY FORMS|genus (refines S217): the truly-hidden part is the DEEP within-genus class group (genus is congruence-detectable). h=genera x deep; Heegner -3,-7,-11 h=1 zero; -15=-3*5 pure GENUS (visible, 0 deep); -23,-47 pure DEEP (1.58,2.32 hidden bits = Hilbert class field, invisible to congruences).
2. TOURNAMENTS|score sequence: transitive (0..n-1)=UNIQUE realization (Landau) => H=0 rigid = the AP/nullcone/rank-11 vertex (S214); near-regular scores carry the hidden fiber (n=5: (1,2,2,2,3) 3 classes = kps reconstruction wall).
3. REALS|CF prefix: golden [0;1,1,..] geo-mean 1 << Khinchin = worst-approx = LRC FOIL (S206); t*=14/183=[0;13,14] geo-mean 13.5 = well-approx = extremal.
4. NULLCONE|moment depth = certificate entropy: LRC finite (bounded alphabet, Bonferroni depth ~5) vs GMC infinite (unbounded degree, Watson S211).

**DUAL entropies:** score-DISTRIBUTION entropy (transitive MAX spread) vs RECONSTRUCTION entropy (transitive 0). The AP = max-order + zero-hidden-info = rigidity; the regular/Paley = min spread + max hidden info.

**UNIFYING:** every repo RIGID extremum (AP/transitive/Heegner h=1/reify-ladder vertex) = a zero-arithmetic-entropy point (local determines global); every DIFFICULTY = its positive-entropy hidden object (deep class group / cospectral fiber / CF tail / deep moment). Rigidity = why the extremal is unique; hidden entropy = where the proof still must go.

**Honest:** genus/class-group, Landau, Khinchin/Levy, detection-depth facts are classical/verified; the contribution is the UNIFICATION (one info deficit across binary forms / tournaments / CF / nullcones) + the rigid=zero-entropy observation. Organizing lens, gate-independent (survives S217 MISTAKE-225), not a proof step. Artifacts: reflection arithmetic-entropy-is-a-repo-wide-invariant-...-boxeph-S218.md, HYP-8875, script (+.out).

