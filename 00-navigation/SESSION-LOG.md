## boxeph-2026-07-21-S221 -- the cusp frame is a repo-wide difficulty-LOCATOR (HYP-8885)
> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).
## death-star-2026-07-21-S100 -- CONFINE the GMC(2) DvdK dependency: elementary except in the resonant-signed corner (a verified sharpening, not a full bypass). HYP-8877.

**Owner directive:** math (not Lean) -- creative ways to complete LRC(14), sharpen targets, or bypass the GMC(2) DvdK dependency.

- **KEY:** `E[P^m]=Σ_channels (multinomial·A(r)!)·c^r` has POSITIVE weights, so the ONLY source of cancellation (the only reason DvdK is nontrivial) is the sign/phase of the coefficient monomials `c^r`.
- **TWO elementary reductions (verified, dvdk_confinement_deathstar_S100.py):** (1) the lowest balanced face is GENERICALLY an EDGE (2 monomials) -- the LP has 2 equality constraints (mass, charge) so a basic optimum is <=2-supported = a straddling edge; there DvdK = elementary binomial `C(m0,k)!=0`, `m0=(a+b)/g`. (2) POSITIVE coeffs => `CT(f^m)>0` (nonneg walk-return, central trinomial), DvdK-free.
- **HARD CORNER (only genuine DvdK use):** resonant face (>=3 charges = codim->=1) AND signed/complex coeffs -- e.g. (-2,-1,1,2)/(1,1,-1,-1) has low CT vanishing. This IS the S89-91 charge-resonance = central-trinomial/free-prob (S90) = Monsky transfer-operator (S99) = the resonant multi-clock (edge=2-clock).
- **HONEST:** not a full bypass (the corner is real DvdK); a verified CONFINEMENT -- GMC(2) DvdK-free for generic/positive; hard input on a codim->=1 stratum. Formalization payoff: elementary edge case + cite DvdK only for the resonant face. reflection confining-...-S100; script (+out). HYP-8877.

## death-star-2026-07-21-S99 -- MERGE: "scale the core, then close on a modular clock" is ONE proof-shape across the nullcone (GMC2, my capstone) and covering (LRC, THM-2057) threads. Lens, not a reduction. HYP-8876.

**Owner:** go back through the repo, apply the cusp frame to under-attended problems, show its power.

**THE FRAME as a diagnostic:** object = EISENSTEIN (computable floor/main term/local) + CUSP (hidden obstruction = genus = deep arithmetic entropy S218). Difficulty is always the CUSP; the frame localizes it + predicts the first-hard-case = first positive cusp dim.

**SWEEPS (verified, the_cusp_frame_as_a_diagnostic_across_the_repo_boxeph_S221.py):**
1. TOURNAMENT COSPECTRALITY (under-attended): char_A spectrum = Eisenstein/local; COSPECTRAL fiber = the reconstruction CUSP; cusp dim = 1,3,28 for n=4,5,6 (first cospectral pair at n=4 = the 'genus' of tournament reconstruction = kps wall = S218 reconstruction entropy). Transitive = spectrally unique (char x^n).
2. INTRANSITIVITY c3 = the tournament's cusp form: transitive c3=0 (Eisenstein/gradient) vs regular c3=5,14,30 (intransitive cusp); the 3-cycle atom (THM-1830) = minimal cusp.
3. GMC(2): E=L o CT (THM-1645) = angular DvdK-closed (EISENSTEIN floor) + radial ker L!=0 (Laplace-determinacy CUSP, verified L(t-1)=L(t^2-3t+1)=0); GMC(n>=3) false = cusp grows.
4. FIGURATE: cake/bagel = smooth Eisenstein polynomial + Fibonacci cusp (S207 recast).

**POWER:** (1) localizes each difficulty to a small nameable cusp (genus-1 newform / cospectral fiber / c3 count / radial kernel); (2) predicts first-hard-case = first positive cusp dim (LRC p=7, tournament recon n=4, intransitivity n=3); (3) unifies as dim(cusp) = S218 deep arithmetic entropy. Eisenstein floor = always the easy/computable/local half; cusp = always where the proof must go.

**Honest:** LRC f14 + the-modular-tournament H are literally modular; the others (cospectral fiber, c3, ker L, Fibonacci) are the analogous main-term+obstruction structure = the arithmetic entropy. A diagnostic lens / difficulty-map, not a proof step. Artifacts: reflection the-cusp-frame-is-a-difficulty-locator-...-boxeph-S221.md, HYP-8885, script (+.out).

