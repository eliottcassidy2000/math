- **HYP-8946 / boxeph-2026-07-22-S234 ((A) COMPLETE: X^M-tR irreducible over F(t) kernel-pure; (B) unramified Hensel scoped):** Worked the two hard THM-2067 pieces. (A) DONE kernel-pure ([propext,Classical.choice,Quot.sound]) -- GMC2PhiIrreducible.phi_irreducible_ratfunc: Phi=X^M-t*R(X) irreducible over F(t) for R(0)!=0, M>=1. Route: (1) phi_t_irreducible -- Phi as a LINEAR poly in t over F[X] has coprime coeffs X^M,-R (IsCoprime from X not dvd R = R(0)!=0, via prime_X.coprime_iff_not_dvd + pow_left), irreducible in F[X][t] (Polynomial.irreducible_of_degree_eq_one_of_isRelPrime_coeff); (2) swap_phi_eq -- transport across Polynomial.Bivariate.swap (swap_C/swap_Y) by MulEquiv.irreducible_iff: swap Phi = map C (X^M) - map C R * C X; (3) primitivity over F[t] -- coeff_M = 1-C(r_M)t, coeff_0 = -C(r_0)t (via coeff_map + coeff_mul_C), any r dividing both divides t (r_0 unit) then divides 1 (Bezout), so unit; (4) Gauss (IsPrimitive.irreducible_iff_irreducible_map_fraction_map) to F(t)=RatFunc F. This is the transitivity input to Polynomial.Gal.galAction_isPretransitive => the orbit-product core (S232). (B) unramified Hensel (THM-1550 Pi=c*t): SCOPED, not complete. Route (death-star-S106): X=sZ,t=s^M => Phi(sZ)=s^M psi(Z), psi=Z^M-R(sZ), psi mod s = Z^M-r_0 SEPARABLE (Mathlib separable_X_pow_sub_C); Hensel lifts small factor A (monic deg M, A cong Z^M-r_0 mod s); Pi=t*(-1)^M A(0), A(0) cong -r_0 => Pi=c*t iff all D_m=0. Mathlib building blocks exist (HenselianLocalRing, IsAdicComplete.henselianRing, hensels_lemma, PowerSeries.exists_isWeierstrassFactorization, separable_X_pow_sub_C) BUT 3 real obstacles => multi-session: (i) HenselianLocalRing (PowerSeries F) NOT a free instance, (ii) degree-DROPPING Hensel factorization (Weierstrass-type, not plain coprime-lift), (iii) Wiener-Hopf/generating-function bridge D_m=0<=>Pi=ct. death-star owns THM-1550; coordinate. AFTER (A): full GMC(2) needs Galois wrapper (Vieta v(C_Phi)=0, equivariance Gal.smul_def, Check A CT(Lambda^m)=[u^Mm]R^m -- all mine, tractable) + (B). File GMC2PhiIrreducible.lean (4 kernel-pure thms).
> **CURRENT DIGEST — refreshed 2026-07-22.** Corrections and proved canon
> outrank hypotheses. Older entries remain in
> [INDEX-HISTORICAL-THROUGH-2026-07-21.md](INDEX-HISTORICAL-THROUGH-2026-07-21.md).

# Current Hypothesis and Frontier Routing

A hypothesis is unresolved unless a proved leaf is named explicitly. Search
its slug and [MISTAKES.md](../../01-canon/MISTAKES.md) before inheritance.

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
- **THM-2104 (PROVED):** every constant `3`-adic quotient-valuation layer has
  one universal sixth-clock escape, with terminal distances `1/6` or `1/3`.
- **THM-2084 (PROVED, cubic source-fiber scope):** a non-tame cubic stratum
  begins, if at all, at unresolved reduced degree `13`; this is not JC(2).
- **THM-2095 (PROVED):** the live guard-ratio common scale divides `252576225`
  and its marked pair is bounded; this does not bound the other six speeds.

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

- **Paper:** THM-2022 proves NC2/GMC(2); THM-2067 supplies its internal seed.
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
