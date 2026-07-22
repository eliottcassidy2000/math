## boxeph-2026-07-21-S222 -- bypass GMC(2)'s DvdK dependency via the saddle-point/Watson method (HYP-8890)

**Owner:** creative ways to bypass the GMC(2) dependency on DvdK.

**TARGET:** DvdK (THM-1630) is the SOLE imported premise of THM-2022 -- used to get a nonzero face constant term Q; residues+Liouville, NON-effective.

**BYPASS (verified, bypass_dvdk_via_saddle_point_watson_boxeph_S222.py):** the needed direction 'f two-sided => CT(f^m)!=0 for some/all large m' is a SADDLE-POINT/WATSON (Laplace) integral. CT(f^m)=[z^0]f^m = (1/2pi) int f(r* e^{i th})^m d th on the saddle circle |z|=r* (r* = mean-exponent-zero radius, exists IFF 0 in int Newton polytope = two-sided). Dominant-saddle asymptotic CT(f^m) ~ rho^m c/sqrt(m), rho=dominant modulus>0 => NONZERO for large m, EFFECTIVE, no residues/Liouville/DvdK.
## codex-2026-07-21 -- modular-form bridge audited; second scaled AP-tail plane closed

- **MISTAKE-226:** HYP-8880/S220 conflated divisor-indexed LRC clocks with
  modular cusps, dilation with `Gamma_0` level, and Hecke coefficients with
  cusp values. The level-14 weight-two newform is a rational non-CM eta
  product, but no transform connects it or its symmetric square to the signed
  phase-height predicate. The modular proposal is retained only as a sidecar
  search prompt.
- **THM-2057 extension:** the missing-clock sieve also closes every
  `{a,2a,...,12a,w}`. The `13a` and `14a` clocks cover all cases except
  `182a|w`; on that ray `t=14m/[a(182m+1)]` gives exact strict margin
  `14m/(182m+1)>1/14`. The exact audit passed all `800000` rows with
  `a<=80,w<=10000`.
- **Next decisive target:** join THM-2058's primitive phase-packet interval
  carrier to THM-2057's missing-clock lcm tax inside each THM-2053 transverse
  deck. Use modular coefficients only if they can be pulled back to signed
  owner-channel sums; the eta-product factorization suggests an
  inclusion--exclusion sidecar, not an obstruction theorem.
## death-star-2026-07-21-S99 -- MERGE: "scale the core, then close on a modular clock" is ONE proof-shape across the nullcone (GMC2, my capstone) and covering (LRC, THM-2057) threads. Lens, not a reduction. HYP-8876.

**VERIFIED cases:** (A) two-sided<=>saddle<=>CT eventually nonzero; one-sided=>CT==0 (the DvdK conclusion, trivial). (B) positive-coeff f=2z+3/z+1: CT(f^m)~f(r*)^m/sqrt(2pi m sig^2), ratio->1. (C) MIXED-sign f=z^2+1/z-1 (the real DvdK case): CT!=0 all m, growth rate |CT|^(1/m)->rho~2.3>0. (D) periodicity (equal-modulus saddles cancel) = the coprime m0 = THM-1840 (elementary); DEGENERATE f(r*)=0 (f=z+1/z-2 => CT=(-1)^m C(2m,m)) = the coalescing/confluent saddle = my S208/HYP-8775 hyper-Bessel cusp (in hand).

**REDUCTION:** DvdK -> standard analytic combinatorics (dominant-saddle nonvanishing) + THM-1840 periodicity + S208 confluent cusp -- all effective, none DvdK's machinery. Makes GMC(2)'s angular/Eisenstein floor (S221) DvdK-free + effective (yields the open effective-DvdK bound m0).

**Honest:** a bypass ROUTE verified in parts, not a complete replacement theorem; full write-up needs steepest-descent through the general (complex/off-axis) dominant saddle + the aperiodicity=>unique-dominant-saddle lemma (both standard Hayman/Pemantle-Wilson, neither needing DvdK). Creative core: DvdK's angular non-vanishing IS a Watson/Laplace saddle count, whose only hard residue is the confluent cusp the repo already resolved. Artifacts: reflection bypassing-dvdk-the-saddle-point-watson-route-...-boxeph-S222.md, HYP-8890, script (+.out).
**SWEEPS (verified, the_cusp_frame_as_a_diagnostic_across_the_repo_boxeph_S221.py):**
1. TOURNAMENT COSPECTRALITY (under-attended): char_A spectrum = Eisenstein/local; COSPECTRAL fiber = the reconstruction CUSP; cusp dim = 1,3,28 for n=4,5,6 (first cospectral pair at n=4 = the 'genus' of tournament reconstruction = kps wall = S218 reconstruction entropy). Transitive = spectrally unique (char x^n).
2. INTRANSITIVITY c3 = the tournament's cusp form: transitive c3=0 (Eisenstein/gradient) vs regular c3=5,14,30 (intransitive cusp); the 3-cycle atom (THM-1830) = minimal cusp.
3. GMC(2): E=L o CT (THM-1645) = angular DvdK-closed (EISENSTEIN floor) + radial ker L!=0 (Laplace-determinacy CUSP, verified L(t-1)=L(t^2-3t+1)=0); GMC(n>=3) false = cusp grows.
4. FIGURATE: cake/bagel = smooth Eisenstein polynomial + Fibonacci cusp (S207 recast).

**POWER:** (1) localizes each difficulty to a small nameable cusp (genus-1 newform / cospectral fiber / c3 count / radial kernel); (2) predicts first-hard-case = first positive cusp dim (LRC p=7, tournament recon n=4, intransitivity n=3); (3) unifies as dim(cusp) = S218 deep arithmetic entropy. Eisenstein floor = always the easy/computable/local half; cusp = always where the proof must go.

**Honest:** LRC f14 + the-modular-tournament H are literally modular; the others (cospectral fiber, c3, ker L, Fibonacci) are the analogous main-term+obstruction structure = the arithmetic entropy. A diagnostic lens / difficulty-map, not a proof step. Artifacts: reflection the-cusp-frame-is-a-difficulty-locator-...-boxeph-S221.md, HYP-8885, script (+.out).
**CLOCKS ARE CUSPS (the merge):** cusps of X0(N) = divisors of N = the sub-clocks (verified #cusps=Sum phi(gcd(d,N/d))). 14-clock=X0(14) cusps {1,2,7,14} (primes 2,7, apex Paley-7); 12-clock=X0(12) cusps {1..12 divs} (primes 2,3, Eisenstein/argmax Phi_6). gcd 2 (chirality), lcm 84 (double-kill 84a|w). SCALING = the Gamma0 level: M(cS)=M(S) verified, scaled zeta-core (codex THM-2057) = same object on a refined clock.

**THE OBSTRUCTION = FIRST CUSP FORM:** genus X0(2p) = 0,0,1,2,2 (p=3,5,7,11,13) => FIRST cusp form at p=7 (X0(14) genus 1, f14=14a) = the first hard case = apex 7. 12-clock genus 0 (CUSPLESS, argmax, rigid) vs 14-clock genus 1 (f14 obstruction). f14 spells 2*7: a_2=-1 (2-cusp), a_7=+1 (7-cusp), w_2=+1,w_7=-1, rank 0 (L(14a,1)>0), period field Q(sqrt-7) = S215 apex disc -7 (as PERIOD field, not weight-1 theta); sym^2 f14 = GL(3) 2nd-moment obstruction.

**ARITHMETIC ENTROPY HOME (S218):** the genus (dim cusp forms) IS the deep/hidden entropy -- genus 0 (argmax) zero deep entropy rigid, genus 1 (LRC) the cuspidal obstruction. General backing (S219 script): theta of binary form = Eisenstein(+)cusp; disc -7 PURE EISENSTEIN (h=1) = GL(2) shadow.
**PROVED TWO FULL ONE-TAIL PLANES:** every row
`{a,2a,...,11a,13a,w}` is LRC(14)-safe. A central-unit orbit gives a witness on
the `12a` clock unless killed; then the `14a` clock works unless killed; double
killing forces `84a|w` and scales HYP-2896's affine binding phase with strict
margin `7m/(84m+5)`. More generally, a missing clock `N<=14` in a core forces
`Na|w`, yielding an lcm divisibility tax over all missing clocks. Exact audit
passed on all `a<=120,w<=12000`. The same sieve closes every
`{a,2a,...,12a,w}`: the `13a` and `14a` clocks leave only `182a|w`, where the
explicit phase has margin `14m/(182m+1)`. This second plane passed all
`a<=80,w<=10000`, alongside a general missing-clock box and named exact
pair-sum controls.

**HONEST:** split is weight-2 X0(14) NOT L(S); VALUE modular / EXISTENCE+constant COMBINATORIAL (floor constant NOT L(14a,1)/L(sym^2)/period -- in the descent); MISTAKE-087 (Phi_6 non-extremal); no weight-1 dihedral construction. Artifacts: reflection the-lrc14-obstruction-is-the-first-cusp-form-scaled-cores-clocks-are-cusps-boxeph-S220.md; HYP-8880; scripts S220 + S219 (+.out).

