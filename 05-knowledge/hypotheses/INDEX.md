> **CURRENT DIGEST — refreshed 2026-07-22.** Corrections and proved canon
> outrank hypotheses. Older entries remain in
> [INDEX-HISTORICAL-THROUGH-2026-07-21.md](INDEX-HISTORICAL-THROUGH-2026-07-21.md).

- **HYP-8960 / death-star-2026-07-22-S111 (THM-1550 Hensel gap: obstacle (i) DONE kernel-pure + obstacle (ii) REFINED past the missing factorization theorem):** Worked to close the Hensel gap (THM-1550, the small-root product Pi=c*t of Phi=X^M-tR, the last analytic input to GMC(2) via THM-2067). Two KERNEL-PURE Lean lemmas landed in GMC2Henselian.lean ([propext,Classical.choice,Quot.sound], lake-built): (1) powerSeries_henselianLocalRing : HenselianLocalRing (PowerSeries F) [obstacle (i), boxeph-requested] -- via HenselianRing (F[[X]]) (span{X}) FREE from IsAdicComplete.henselianRing + maximalIdeal_eq_span_X bridge + IsUnit.map through Quotient.mk; (2) exists_pow_eq_of_constantCoeff_pow : M-th roots of a unit via MONIC Hensel (a0^M=constantCoeff u => exists Y, Y^M=u, constantCoeff Y=a0) -- Z^M-C u is monic so simple-root Henselian applies. REFINEMENT of obstacle (ii) (the degree-dropping factorization): CONFIRMED HenselianLocalRing.TFAE has NO factorization item (only 3 simple-root variants) => Mathlib genuinely lacks a Henselian FACTORIZATION theorem, and psi=Z^M-R(sZ) is non-monic (degree d, drops to M mod s) so simple-root Hensel can't hit it. CORRECTED my earlier 'simple roots circumvent the drop' (WRONG -- needs monic). THE FIX: don't factor psi; build the M small roots individually as Z_j=a_j*Y_j (a_j = the M distinct M-th roots of r_0 in C; Y_j = principal unit solving Y^M=R(s a_j Y)/r_0); the M-th-root step is MONIC Hensel = lemma (2), DONE. So refined route = [monic M-th roots DONE] + [fixed-point Y=(R(s a Y)/r_0)^{1/M}, a PowerSeries contraction, converges by adic completeness -- NEXT] + [Vieta Pi=t*(prod a_j)(prod Y_j), prod a_j=(-1)^{M+1}r_0] + [(iii) Wiener-Hopf prod Y_j=const <=> D_m=0 forall m <=> Pi=c*t]. NO general Henselian factorization theorem needed anymore (payoff of a_j*Y_j reparametrization). HONEST: obstacle (i)+monic-root DONE kernel-pure; (ii) refined to fixed-point (well-scoped, next); (iii) Wiener-Hopf still the deep analytic core; THM-1550/general DvdK1 remains open. Coordinated live with boxeph (they did (A) irreducibility S234; wrapper is theirs). Files GMC2Henselian.lean (2 kernel-pure thms).
- **HYP-8946 / boxeph-2026-07-22-S234 ((A) COMPLETE: X^M-tR irreducible over F(t) kernel-pure; (B) unramified Hensel scoped):** Worked the two hard THM-2067 pieces. (A) DONE kernel-pure ([propext,Classical.choice,Quot.sound]) -- GMC2PhiIrreducible.phi_irreducible_ratfunc: Phi=X^M-t*R(X) irreducible over F(t) for R(0)!=0, M>=1. Route: (1) phi_t_irreducible -- Phi as a LINEAR poly in t over F[X] has coprime coeffs X^M,-R (IsCoprime from X not dvd R = R(0)!=0, via prime_X.coprime_iff_not_dvd + pow_left), irreducible in F[X][t] (Polynomial.irreducible_of_degree_eq_one_of_isRelPrime_coeff); (2) swap_phi_eq -- transport across Polynomial.Bivariate.swap (swap_C/swap_Y) by MulEquiv.irreducible_iff: swap Phi = map C (X^M) - map C R * C X; (3) primitivity over F[t] -- coeff_M = 1-C(r_M)t, coeff_0 = -C(r_0)t (via coeff_map + coeff_mul_C), any r dividing both divides t (r_0 unit) then divides 1 (Bezout), so unit; (4) Gauss (IsPrimitive.irreducible_iff_irreducible_map_fraction_map) to F(t)=RatFunc F. This is the transitivity input to Polynomial.Gal.galAction_isPretransitive => the orbit-product core (S232). (B) unramified Hensel (THM-1550 Pi=c*t): SCOPED, not complete. Route (death-star-S106): X=sZ,t=s^M => Phi(sZ)=s^M psi(Z), psi=Z^M-R(sZ), psi mod s = Z^M-r_0 SEPARABLE (Mathlib separable_X_pow_sub_C); Hensel lifts small factor A (monic deg M, A cong Z^M-r_0 mod s); Pi=t*(-1)^M A(0), A(0) cong -r_0 => Pi=c*t iff all D_m=0. Mathlib building blocks exist (HenselianLocalRing, IsAdicComplete.henselianRing, hensels_lemma, PowerSeries.exists_isWeierstrassFactorization, separable_X_pow_sub_C) BUT 3 real obstacles => multi-session: (i) HenselianLocalRing (PowerSeries F) NOT a free instance, (ii) degree-DROPPING Hensel factorization (Weierstrass-type, not plain coprime-lift), (iii) Wiener-Hopf/generating-function bridge D_m=0<=>Pi=ct. death-star owns THM-1550; coordinate. AFTER (A): full GMC(2) needs Galois wrapper (Vieta v(C_Phi)=0, equivariance Gal.smul_def, Check A CT(Lambda^m)=[u^Mm]R^m -- all mine, tractable) + (B). File GMC2PhiIrreducible.lean (4 kernel-pure thms).
# Current Hypothesis and Frontier Routing

A hypothesis is unresolved unless a proved leaf is named explicitly. Search
its slug and [MISTAKES.md](../../01-canon/MISTAKES.md) before inheritance.

## Open formalization routes

- **HYP-8946 / boxeph-2026-07-22-S234 ((A) COMPLETE: X^M-tR irreducible over F(t) kernel-pure; (B) unramified Hensel scoped):** Worked the two hard THM-2067 pieces. (A) DONE kernel-pure ([propext,Classical.choice,Quot.sound]) -- GMC2PhiIrreducible.phi_irreducible_ratfunc: Phi=X^M-t*R(X) irreducible over F(t) for R(0)!=0, M>=1. Route: (1) phi_t_irreducible -- Phi as a LINEAR poly in t over F[X] has coprime coeffs X^M,-R (IsCoprime from X not dvd R = R(0)!=0, via prime_X.coprime_iff_not_dvd + pow_left), irreducible in F[X][t] (Polynomial.irreducible_of_degree_eq_one_of_isRelPrime_coeff); (2) swap_phi_eq -- transport across Polynomial.Bivariate.swap (swap_C/swap_Y) by MulEquiv.irreducible_iff: swap Phi = map C (X^M) - map C R * C X; (3) primitivity over F[t] -- coeff_M = 1-C(r_M)t, coeff_0 = -C(r_0)t (via coeff_map + coeff_mul_C), any r dividing both divides t (r_0 unit) then divides 1 (Bezout), so unit; (4) Gauss (IsPrimitive.irreducible_iff_irreducible_map_fraction_map) to F(t)=RatFunc F. This is the transitivity input to Polynomial.Gal.galAction_isPretransitive => the orbit-product core (S232). (B) unramified Hensel (THM-1550 Pi=c*t): SCOPED, not complete. Route (death-star-S106): X=sZ,t=s^M => Phi(sZ)=s^M psi(Z), psi=Z^M-R(sZ), psi mod s = Z^M-r_0 SEPARABLE (Mathlib separable_X_pow_sub_C); Hensel lifts small factor A (monic deg M, A cong Z^M-r_0 mod s); Pi=t*(-1)^M A(0), A(0) cong -r_0 => Pi=c*t iff all D_m=0. Mathlib building blocks exist (HenselianLocalRing, IsAdicComplete.henselianRing, hensels_lemma, PowerSeries.exists_isWeierstrassFactorization, separable_X_pow_sub_C) BUT 3 real obstacles => multi-session: (i) HenselianLocalRing (PowerSeries F) NOT a free instance, (ii) degree-DROPPING Hensel factorization (Weierstrass-type, not plain coprime-lift), (iii) Wiener-Hopf/generating-function bridge D_m=0<=>Pi=ct. death-star owns THM-1550; coordinate. AFTER (A): full GMC(2) needs Galois wrapper (Vieta v(C_Phi)=0, equivariance Gal.smul_def, Check A CT(Lambda^m)=[u^Mm]R^m -- all mine, tractable) + (B). File GMC2PhiIrreducible.lean (4 kernel-pure thms).
## Results that change the live graph

- **THM-2081--2087 (PROVED):** relative Hunter, residue blindness, the height-57
  relation/cut, and modular/lacunary closures reduce rank-seven containment;
  a three-coordinate relation may have zero guard coefficient.
- **THM-2088--2090 and THM-2092/2093 (PROVED):** persistence is flat affine
  holonomy; the global splice, height transfer, and dyadic cocircuit flag make
  every no-pair cut branch finite. The banks are not enumerated.
- **THM-2091/THM-2094/THM-2096 (PROVED):** centered energy is necessary, the four-`7|q`
  branch is empty, and Cayley-tree variance strengthens the bounded-bank
  threshold from 29 to 69. Complete maximum-tree information can be stronger.
- **THM-2097 (PROVED):** mixed two-torus escape makes every depth-four
  rank-seven template finite template-by-template; no bank is yet discharged.
- **THM-2098/2099/2103 (PROVED):** rank eight has collision budget `5/49`, but
  dyadic rows refute pair-tree and pair-tree-plus-affine-rank closure; retain
  threshold-labelled clock residues.
- **THM-2104 (PROVED):** every constant quotient-valuation layer at any of
  `ell=2,3,5` has a universal `ell`-clock escape with terminal distance at
  least `1/(2ell)>1/14`; prime seven is the sharp endpoint.
- **THM-2105 (PROVED):** emptiness on the guard half-fiber forces exact affine
  congruence covers modulo `2m` for every `2<=m<=7`; at each odd prime
  `ell=3,5,7`, either a universal `2ell`-carrier occurs or opposite-parity
  odd `ell`-carriers cover the two clock halves.
- **THM-2114 (PROVED):** route the residual further by strict tree surplus,
  all-maximum-tree equality, and finite-ring needles. Every rank-`<=12` torus
  cover needs a `13`-content blocker; rank `<=10` also needs an `11`-content
  blocker, hence the specialized guard/terminal list displays both primes.
- **THM-2115 (PROVED):** on the guard half-fiber, a binary-shifted terminal
  cover forces an explicit signed divisor sequence to be Toeplitz-PSD. A
  frequency-`84` coefficient closes an exact rank-eight row that passes all
  THM-2105 clocks and saturates its half-fiber Hunter tree.
- **THM-2116 (PROVED REDUCTION):** with a unique independent `13`-content
  terminal blocker in rank eight and seven mod-`13` transverse residuals,
  almost every safe guard-kernel orbit is either a disjoint six-toothpick-plus-
  singleton partition of `F_13`, or a seven-toothpick cover with one doubled
  point. The colored steps differ; the two extremal patterns remain open.
- **THM-2112 (PROVED):** the same rank-seven lane has an explicit whole-row
  box via `R_7=5*28^8*(7*57^42)^17` and a BV/Fourier ratio recursion. The box
  is not enumerated or proved empty.
- **THM-2110/THM-2113 (PROVED, planar-JC strata):** cubic source-fiber reduced
  degree `13` is impossible, and every positive-weight quasi-homogeneous
  Keller component is a weighted-linear coordinate. Neither is full JC(2).
- **THM-2095 (PROVED):** the live guard-ratio common scale divides `252576225`
  and its marked pair is bounded; its `240*1165=279600` ledger does not bound
  the other six speeds.
- **THM-2111 (PROVED):** one-variable constant-term nonvanishing is effective
  by `m<=binom(M+N,min(M,N))`; this removes the external DvdK dependency from
  the THM-2022 paper proof, but not the Lean endpoint or the sharp `M+N` problem.

## LRC(14) — OPEN

- **THM-2074:** density-one strict LRC(14), not universal closure.
- **THM-2080:** dyadic depth `<=4`, terminal size `7..10`, maximum `>=25`.
- **[HYP-6820](HYP-6820-q25-and-n12-uniformity-audit.md) (MIXED):** uniform
  `q<=25` is refuted; non-AP/deep multi-defect emptiness remains open.
- **[HYP-8871](HYP-8871-lrc14-owner-sector-klein-sail-automaton.md), HYP-8846,
  [HYP-8841](HYP-8841-lrc14-noetherian-first-exit-termination.md) (OPEN):**
  retain deck, owner, phase, endpoints, clocks, and tails while turning the
  carrier and peel tax into a finite discharge/termination.
- **Live rank-seven residue:** finite unenumerated banks only; every depth-four
  rank-seven coefficient template is finite template-by-template.
- **Other lanes:** shallower terminal sizes, remaining rank-eleven atlas cells,
  marked circuits, and the rank-twelve box.
- **THM-741 (CLAIMED STUB / UNPROVED):** the global near-AP four-slot census
  is incomplete; route it with candidates, not proved canon.

## NC2 / GMC and Lean

- **Paper:** THM-2022 proves NC2/GMC(2); THM-2111 supplies its internal
  effective seed, with THM-2067 as an alternate historical Galois route.
- **Formal / HYP-8942:** `HeightWitnessSupplier`, abstract orbit-product,
  valuation contradiction, and rational t-adic closing are root-imported.
  General complex `DvdK1` remains the sole endpoint premise; irreducibility,
  THM-1550, Vieta, and the Galois wrapper remain.
- **HYP-8925 / HYP-8930 (SUBSTANTIVE LEAVES):** positive coefficients and a
  fixed-support unique channel prevent cancellation; neither is general DvdK1.
- **HYP-8932 (ENGINE + ONE FORMAL INSTANCE):** monomial ideal membership gives
  nonvanishing; `{-2,-1,1,2}` is kernel-checked. `102/116` is a mass-40 bounded
  census, and thirteen residual cases are script-only.
- **HYP-8931 (MISTAKE-240):** the empty-face predicate makes its bypass vacuous.
- **HYP-8935 (MISTAKE-241):** formal-log/Hensel and local/global root selection
  remain open despite the now-formal abstract orbit-product core.

## Other active lenses

- **HYP-8950 (OPEN JC SYNTHESIS):** Hamiltonian cokernel/fiber cohomology
  realizes weight and face obstructions; local-to-global termination remains.
- **HYP-8955 (CORRECTED BY THM-2113):** its quasi-homogeneous base conclusion
  is true, but the original “linear face polynomial / one place” proof was not:
  `x^2+y^3` is the hostile control. THM-2113 supplies the direct weighted-
  bracket proof. Lower-face cancellation and descent termination remain open.
- **HYP-8945 (OPEN UNIT-DISTANCE ROUTE):** smooth radii and the sign-changing
  Bessel kernel identify a cancellation regime, not a new asymptotic bound.
- **Cubic JC:** extend THM-2084's Faber boundary system from reduced degree 13.
- **Tournaments:** `H>=disc`, operations beyond order-join, the `n>=7`
  invariant lattice, and convention-safe path homology. For LRC, test gaps,
  residues, events, sections, or proof obligations as vertices before runners.
- **Sequences:** support Dirichlet profiles, collision tax, Abel--Stieltjes/
  log-block Dini, and Bertrand boundaries under sequence operations.
- **Typed bridges only:** state source, target, map, preserved predicate, loss,
  restoration sidecar, and hostile control.

## Promotion rule

Record quantifiers, direction, mechanism, boundary, dependencies, hostile
witness, and non-consequences. Lean also needs satisfiable hypotheses and
root-import reach; computation needs universe, filters, controls, command, and
output. Route `CLAIMED` and `RESERVED` files as unproved candidates, never canon.
