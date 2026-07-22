- **HYP-8975 / kind-pasteur-2026-07-21-S128c147 (DvdK1 elementary reductions: zero-charge discharged + DvdK1 reduced to the both-signs case + gcd-invariance; the bypass VERDICT):** Owner: prove GMC(2) unconditionally, own remaining tasks, explore repo for a DvdK bypass / alternative route. VERIFIED STATE: GMC(2)/NC2 is a kernel-pure reduction to a SINGLE external input DvdK1 -- heightWitnessSupplier_holds is PROVED (GMC2NC2, builds green), so nc2_of_dvdK1 : DvdK1 -> NC2. DvdK1 route = THM-2067 concrete 16-thm skeleton; open inputs = (1) Vieta hΩ (core done by death-star coeff_ratio_Phi_eq_const; assembling prod-over-rootSet bounded, boxeph) + (2) THM-1550 hS/hfix (the deep Wiener-Hopf/unramified-Hensel analytic gap, death-star; obstacle (i) Henselian DONE). Check A (GMC2LaurentShiftCheckA), Galois action, irreducibility all DONE/green. MY CONTRIBUTION (GMC2DvdKZeroCharge.lean, kernel-pure [propext,Classical.choice,Quot.sound], pushed): dvdK1_of_zero_mem (0 in support => CT nonzero at m=1, delta_{i0} unique balanced channel); dvdK1_of_bothSigns : DvdK1BothSigns -> DvdK1 (FORMALLY reduces all of DvdK1 to the both-signs no-zero case, since ChargesStraddleZero = zero-charge OR both-signs); constantTermRelation_scale (gcd-rescaling invariance => hard case may assume gcd 1). Combined with prior elementary cases (two-charge/unique-channel/positive-coeff), the NO-CANCELLATION territory of DvdK1 is now elementary+formal; hard case pinned to both-signs, no-zero, gcd 1, >=2 coincident balanced channels at every m. BYPASS VERDICT (honest, task done): NO elementary bypass of the analytic core exists -- THM-2067 (multiplicative orbit-product) and all THREE THM-2101 additive proofs (monodromy / transcendental specialization / t-adic Newton-polygon+partial-fractions) bottom out at the SAME small-root-product/Wiener-Hopf content (THM-1550). The difficulty is intrinsic: complex-coefficient phase cancellation across >=2 channels for ALL m (positive coeffs => trivial). GMC(2) is ONE theorem (THM-1550) from unconditional, a genuine multi-session analytic formalization owned by death-star, NOT sidesteppable. HONEST: I did NOT make GMC(2) unconditional; I sharpened the perimeter + isolated the hard case + gave the bypass verdict. Files GMC2DvdKZeroCharge.lean + reflection dvdk1-elementary-reductions-and-the-bypass-verdict-kps-S128c147. Complementary to HYP-8956/8965/8951/8960/8970.
- **HYP-8976 / death-star-2026-07-22-S113 (GMC2: full Vieta prod_roots_Phi kernel-pure -- the hOmega input for the concrete THM-2067; hOmega discharge WIP):** Owner: finish (a) Gal-instantiation + (b) THM-1550, pulling in concurrent pieces. DELIVERED (kernel-pure, GMC2PhiVieta.lean, lake-built): **prod_roots_Phi** -- over any E where Phi=X^M-tR splits, (Phi.map E).roots.prod = (-1)^(deg R) * algebraMap(algebraMap(R.coeff 0/R.leadingCoeff)), a CONSTANT (the t cancels) = valuation-0. Built on coeff_ratio_Phi_eq_const (S112) + Splits.coeff_zero_eq_leadingCoeff_mul_prod_roots. This is the number-theoretic CONTENT of boxeph's hOmega hypothesis (thm2067_contradiction_concrete: prod over rootSet = algebraMap(RatFunc.C d)). PULLED IN (the fleet advanced (a) concurrently): boxeph S236 (HYP-8956) did the CONCRETE Gal instantiation (GMC2Thm2067Concrete + GMC2GalRootAction -- incl. the equivariance/rootsEquivRoots that I independently localized to a tower Algebra-SF-SF instance diamond, then superseded); mac-mini S162 proved irreducibility of the degree-(M+N) Phi over F(t) (GMC2DvdKParameterIrreducible.irreducible_map_ratfunc, discharging hPhi). So thm2067_contradiction_concrete now needs only hPhi [mac-mini, done], hOmega [my Vieta content, done -- final subtype-product/rootSet-Fintype wiring WIP], hS [THM-1550, the DEEP gap]. ATTEMPTED (a) equivariance myself (~14 iters), localized precisely to `algebraMap SF SF = id` for the IsScalarTower.toAlgHom F SF SF instance vs canonical Algebra.id (a Mathlib instance-diamond) -- handed to boxeph, who'd already solved the full instantiation. HONEST: Vieta landed kernel-pure (the hOmega content); (a) is near-complete via fleet pieces + mine; the hOmega final wiring (rootSet subtype product = aroots.prod, dependent-type-finicky) + (b) THM-1550 (Wiener-Hopf/additive germ identity, multi-session crux) remain. Pushed/pulled/coordinated throughout. Files GMC2PhiVieta.lean. NOTE: kps-S128c147 (HYP-8975) independently reached the same honest verdict from the additive side (the bypass verdict): confirms (b) THM-1550 is the sole remaining crux.
- **HYP-8956 / boxeph-2026-07-22-S236 (CONCRETE Galois instantiation of THM-2067 -- kernel-pure; the whole THM-2067 argument is now a concrete reduction to Vieta + THM-1550):** Worked the Galois wrapper piece to completion. NEW kernel-pure ([propext,Classical.choice,Quot.sound]): GMC2GalRootAction -- the DIRECT action sigma.x = sigma(x) of L~=aK L on Phi.rootSet L: mem_rootSet_smul (autos preserve roots), rootAction (MulAction instance), coe_smul (equivariance coe(sigma.x)=sigma(coe x), TAUTOLOGICAL/rfl -- sidesteps Mathlib galAction/rootsEquivRoots), isPretransitive_rootAction (transitivity for irreducible Phi over Normal L via IsConjRoot.exists_algEquiv + minpoly K y associate of irreducible Phi). GMC2Thm2067Concrete.thm2067_contradiction_concrete -- instantiates the abstract wrapper at Phi.Gal on Phi.rootSet Phi.SplittingField: irreducible Phi over F(t) + small-root product = c*t and Gal-fixed (THM-1550) + full product = const d (Vieta) => False. All Galois instances resolve via the Phi.Gal type (Fintype, MulDistribMulAction, IsSplittingField.finiteDimensional). So the ENTIRE THM-2067 argument is now a CONCRETE kernel-pure reduction: GMC2OrbitProduct -> GMC2RatFuncClosing -> GMC2PhiIrreducible (A) -> GMC2Thm2067Wrapper -> GMC2GalRootAction -> GMC2Thm2067Concrete (16 thms), open inputs = Vieta + THM-1550. COORDINATION w/ death-star: they DROVE THM-1550 -- landed obstacle (i) HenselianLocalRing(PowerSeries F) kernel-pure (GMC2Henselian.lean) + monic M-th-root Hensel + a_j*Y_j reparametrization (killing the need for a degree-dropping factorization theorem); split agreed: death-star = fixed-point convergence + Vieta-for-Pi + Wiener-Hopf (D_m=0<=>Pi=ct), me = Gal-instantiation (DONE) + Vieta(hOmega full product) + Check A. My wrapper takes THM-1550 as exactly hS/hfix so their output slots in. REMAINING (mine, bounded): Vieta (prod roots=(-1)^d r0/rd, non-monic Vieta + separability rootSet<->multiset) + Check A (CT(Lambda^m)=[u^Mm]R^m; no Polynomial.coeff_pow in Mathlib => custom coeff-of-power lemma). Files GMC2GalRootAction.lean, GMC2Thm2067Concrete.lean.
- **HYP-8970 / kind-pasteur-2026-07-21-S128c146 (GMC(2) orthogonality CORE extracted as a standalone Mathlib-PR-ready lemma; full GMC(2) stays conditional on DvdK1):** Owner: get the GMC(2) formalization complete + Mathlib-PR-ready; be creative bypassing hard pieces. The full GMC2NC2Capstone is sorry-free but CONDITIONAL on DvdK1 (external DvdK analytic input, which death-star/boxeph/codex are actively completing via THM-2067/THM-1550/Henselian) + HeightWitnessSupplier -- NOT a Mathlib PR (research-grade, gated on unformalized analysis). CREATIVE BYPASS: extracted the clean general algebraic CORE that IS Mathlib-worthy. NEW FILE ThreeTermRecCoprime.lean (builds green, axioms [propext,Classical.choice,Quot.sound]): for a monic three-term recurrence p(n+2)=(X-C a_{n+1})p(n+1)-C b_{n+1} p n over ANY CommRing, if every b_{n+1} is a UNIT then IsCoprime (p n) (p (n+1)) (isCoprime_succ) -- the estimate-free 'consecutive orthogonal polynomials are coprime' (Favard; Hermite/Legendre/Laguerre/Chebyshev uniform), the replacement for the REFUTED domination step (THM-1585->THM-1660). Plus: field version (b ne 0), noCommonRoot [Nontrivial R], monic + natDegree=n, probabilists' Hermite instance. STRENGTHENS the repo's prior GMC2HermiteNoCommonRoot.lean (functions R->R) + ThreeTermRecurrence.lean (no_common_root_poly [IsDomain R]): IsCoprime (algebra-level, stronger than no-common-root) over general CommRing (not IsDomain) + monic/degree. Mathlib GAP confirmed: has Polynomial.hermite/Chebyshev but no general three-term-recurrence coprimality. Proof = one Bezout update closed by linear_combination (no analysis). Deliverable: ThreeTermRecCoprime.lean + ThreeTermRecCoprime-MATHLIB-PR-NOTES.md (statement, gap justification, provenance, pre-submission checklist: narrow imports via shake, optional Legendre/Chebyshev instances). Reflection gmc2-orthogonality-core-is-mathlib-pr-ready-kps-S128c146. Full GMC(2) is one analytic theorem (DvdK1) from unconditional -- correctly owned by death-star/boxeph/codex, not blocked by anything here. Complementary to HYP-8965/8951/8960.
- **HYP-8965 / death-star-2026-07-22-S112 (GMC2 formalization: the Vieta valuation-0 coeff-ratio input for the THM-2067 wrapper, kernel-pure + critical-path coordination):** Continued improving the GMC2/DvdK1 Lean formalization. Surveyed the (crowded, fast-moving) state: codex's additive-orbit-residue route (THM-2101: Check A, additive orbit-sum, full-root Lagrange identity, small-root existence via reciprocal-monic Hensel -- which IMPORTS AND USES my S111 GMC2Henselian.HenselianLocalRing) + boxeph's multiplicative orbit-product wrapper (S235 GMC2Thm2067Wrapper, premising THM-1550 + Vieta). DELIVERED (kernel-pure [propext,Classical.choice,Quot.sound], lake-built): GMC2PhiVieta.lean coeff_ratio_Phi_eq_const -- for Phi=X^M - t*R over F(t), Phi.coeff 0 / Phi.leadingCoeff = algebraMap(R.coeff 0 / R.leadingCoeff), i.e. the ratio is t-FREE (a constant in the image of F), because coeff 0 = -t*r0 and leadingCoeff = -t*lc(R) and the t cancels. This is the number-theoretic content of the wrapper's hOmega input (Vieta: prod_roots = (-1)^d*(coeff0/leadingCoeff) is a valuation-0 constant). CRITICAL COORDINATION: boxeph S235 asked me to derive HenselianLocalRing (PowerSeries F) + the degree-dropping factorization, unaware both are DONE -- HenselianLocalRing is my S111 (GMC2Henselian, HYP-8960), and the degree-drop is SOLVED by codex's reciprocal/reversed-monic trick (GMC2ReciprocalSmallRoots, which imports my instance); messaged boxeph to prevent duplication. So THM-2067 is now a kernel-pure SKELETON down to: (a) boxeph's concrete Gal instantiation (equivariance via rootsEquivRoots) + Check A + Vieta (its core supplied here), and (b) THM-1550 piece (iii) the Wiener-Hopf D_m=0 forall m <=> prod_small=c*t bridge (mine, the one deep analytic gap; obstacles (i)+(ii) handled). HONEST: modest lemma + high-value coordination; my Henselian is confirmed on the critical path (used by codex). Files GMC2PhiVieta.lean (kernel-pure).
- **HYP-8951 / boxeph-2026-07-22-S235 (THM-2067 WRAPPER assembled kernel-pure -- the full orbit-product argument as a Lean reduction; isolates THM-1550 + Vieta + Gal-instantiation):** Worked the Galois wrapper + the deep gap. DELIVERED kernel-pure ([propext,Classical.choice,Quot.sound]) GMC2Thm2067Wrapper: pow_monomial_eq_const_absurd ((c*t)^N=d^K impossible over F(t), via monomial_pow_ne_const) and thm2067_contradiction (the COMPLETE THM-2067 argument as a reduction: transitive Galois-type action on Omega + equivariant embedding f:Omega->E over F(t) + THM-1550 [small-root product prod_S f = c*t and Gal-fixed] + Vieta [prod_Omega f = const d != 0] => False; via prod_pow_card_group_eq (S232) => pullback along injective F(t)->E => monomial closing). So the WHOLE THM-2067 chain is now a kernel-pure Lean SKELETON: GMC2OrbitProduct (core) -> GMC2RatFuncClosing -> GMC2PhiIrreducible (A, irreducibility) -> GMC2Thm2067Wrapper (13 kernel-pure thms). REMAINING (2 hard pieces, scoped): (1) concrete Gal instantiation -- plug G=Phi.Gal, Omega=Phi.rootSet E; instances + transitivity (galAction_isPretransitive <- my A) free, but EQUIVARIANCE coe(phi.x)=phi(coe x) routes through Gal.smul_def = rootsEquivRoots(phi . rootsEquivRoots.symm x), needs characterizing Gal.rootsEquivRoots (finicky, needs E=SplittingField) + Vieta (prod_roots, t cancels => (-1)^d r0/rd) + Check A (CT(Lambda^m)=[u^Mm]R^m); (2) THM-1550 deep gap (prod_small=c*t Gal-fixed) = unramified Hensel; obstacles confirmed: IsAdicComplete(maximalIdeal(PowerSeries F)) NOT synthesizable => HenselianLocalRing must be derived, degree-dropping Weierstrass-type factorization, Wiener-Hopf D_m=0<=>Pi=ct bridge; multi-session, death-star owns. Full GMC(2) <= DvdK1 <= THM-2067 = kernel-pure skeleton, open inputs = Gal-plumbing (mine, bounded) + THM-1550 (deep gap). File GMC2Thm2067Wrapper.lean.
> **CURRENT DIGEST — refreshed 2026-07-22.** Corrections and proved canon
> outrank hypotheses. Older entries remain in
> [INDEX-HISTORICAL-THROUGH-2026-07-21.md](INDEX-HISTORICAL-THROUGH-2026-07-21.md).

- **HYP-8960 / death-star-2026-07-22-S111 (RECIPROCAL-MONIC REPAIR, kernel-pure):** `GMC2Henselian` proves `F⟦s⟧` Henselian and lifts roots of units. Root-imported `GMC2ReciprocalSmallRoots` monicizes the genuine-degree reverse and now proves the concrete specialization: if `M<deg R`, `R(0)≠0`, `a^M=R(0)`, `a≠0`, and the characteristic does not divide `M`, then `Z^M-R(sZ)` has a power-series root reducing to `a`. All public theorems use only `[propext, Classical.choice, Quot.sound]`. This bypasses both the self-referential fixed-point proposal and unavailable degree-dropping factorization. Enumerating all branches and proving their product/Wiener--Hopf identity remain open on historical THM-1550; paper-proved THM-2101 bypasses that product, but Lean still lacks analytic block continuation and the final concrete `DvdK1` wrapper.
- **HYP-8946 / boxeph-2026-07-22-S234 (ALGEBRAIC CORE COMPLETE; LOCAL ROUTE SUPERSEDED):** `GMC2PhiIrreducible.phi_irreducible_ratfunc` kernel-checks `X^M-tR(X)` over `F(t)` and supplies transitivity. Its original Hensel plan listed a missing local-ring instance and degree-dropping factorization; HYP-8960 now bypasses both by reciprocal monicization and directly constructs every chosen ramified root. Only branch enumeration/product descent and the Wiener--Hopf identity remain on historical THM-1550; neither is needed by THM-2101's additive paper proofs.
# Current Hypothesis and Frontier Routing

A hypothesis is unresolved unless a proved leaf is named explicitly. Search
its slug and [MISTAKES.md](../../01-canon/MISTAKES.md) before inheritance.

## Open formalization routes

- **HYP-8946:** `X^M-tR` irreducibility over `F(t)` is kernel-checked and used
  by THM-2101; degree-dropping Hensel factorization is now optional.
- **HYP-8960:** reciprocal monicization now constructs each root of
  `Z^M-R(sZ)` from a simple nonzero residue root. Simultaneous branch/product
  control and the Wiener--Hopf identity remain open.
- **THM-2101:** strict two-sided DvdK has three product-free paper proofs:
  monodromy, transcendental specialization, and a purely t-adic packet proof.
  Their global constructions and final `DvdK1` wrapper remain Lean assembly.

## Results that change the live graph

- **THM-2081--2087 (PROVED):** relative Hunter, residue blindness, the
  height-57 cut, and modular/lacunary closures reduce rank-seven containment;
  a three-coordinate relation may have zero guard coefficient.
- **THM-2088--2090 and THM-2092/2093 (PROVED):** flat persistence, the global
  splice, height transfer, and dyadic cocircuit flags make every no-pair cut
  branch finite. The banks are not enumerated.
- **THM-2091/THM-2094/THM-2096 (PROVED):** centered energy is necessary, the
  four-`7|q` branch is empty, and Cayley-tree variance raises bounded-bank
  thresholds; complete maximum-tree information can be stronger.
- **THM-2097/2100/2112 (PROVED):** every depth-four rank-seven template has an
  explicit whole-row finite box, but no box is enumerated or proved empty.
- **THM-2098/2099/2103/2119 (PROVED):** guarded sizes `8..10` split into pure,
  low-mixed, and high-vertical lanes; only pure rows inherit the collision
  budget. Pair trees miss dyadic rows; pure-transverse rank eight is
  projectively three-sparse under THM-2119's exact hypotheses.
- **THM-2104/2105 (PROVED):** quotient valuations cross `2,3,5` walls, and
  emptiness forces exact affine carriers through denominator fourteen.
- **THM-2114/2115 (PROVED):** finite-ring needles force `13`- and `11`-content
  blockers; a joint Toeplitz certificate closes a row missed by scalar clocks.
- **THM-2116/2120/2122/2123/2124/2125/2128/2131/2133/2135 (PROVED):** a rank-eight cover needs a guard
  `13`-blocker or at least five nonblocker terminals projectively parallel to
  the guard modulo thirteen. THM-2124 proves the complementary finite-plane seven-pencil theorem; THM-2128 kills `(7,1)` and THM-2131's digit lift kills
  `(8)` when all terminals are nonblockers. THM-2133 leaves scalar `6+1`/`5+2` tails; THM-2135 forces `13^4|v` in `6+1` and excludes `5+2` depths `(1,1),(1,2)`.
- **THM-2117/2121 (PROVED):** clocks, the maximum Hunter tree, and all scalar
  minors can miss an open safe cell; every strict safe cell has a full
  Toeplitz/Fejer certificate of order at most `14nV^2+1`. Boundary-only points
  are outside that statement.
- **THM-2126 (PROVED):** scalar pair filters have a neutral wall and no uniform
  positive margin. Relation labels, common fibers, and higher multiplicity are
  essential sidecars; its examples are excluded by other gates.
- **THM-2130 (PROVED):** root capacity at `11/13` extends selected higher ranks
  and forces a rank-eight sparse mod-11 alternative or `143|det(g,c_i)`.
- **THM-2095 (PROVED):** the live guard-ratio scale divides `252576225`; the
  `240*1165=279600` marked-pair ledger does not bound the other six speeds.
- **THM-2101/2111 (PROVED PAPER):** three additive DvdK proofs avoid root products, including a purely t-adic Newton-packet proof; the
  effective first return is a compound-determinant order at most
  `binom(M+N,min(M,N))`. Sharpening and Lean assembly are separate.
- **THM-2102/2110/2113/2118/2127/2129/2132/2134/2136 (PROVED planar-JC strata):** cubic source
  fibers and exact two-face trains close. A nonradial factor edge is a coarsened
  toric scalar power or terminally short; global Hermite compatibility remains.

## LRC(14) — OPEN

- **THM-2074:** density-one strict LRC(14), not universal closure.
- **THM-2080:** dyadic depth at most four, terminal size `7..10`, maximum at
  least `25`; depth-zero eleven and rank twelve remain separate.
- **HYP-6820 (MIXED):** uniform `q<=25` is refuted; non-AP/deep multi-defect
  emptiness remains open.
- **HYP-8871/HYP-8846/HYP-8841 (OPEN):** retain deck, owner, phase, endpoints,
  clocks, and tails through any discharge.
- **Live residues:** rank-seven finite banks; the guard-`13` pencil and
  nontransverse guarded lanes; depth-zero eleven, marked circuits, and the
  rank-twelve box.
- **THM-741 (CLAIMED STUB / UNPROVED):** the global near-AP four-slot census is
  incomplete; route it as a candidate, not proved canon.

## NC2 / GMC and Lean

- **Paper:** THM-2022 proves NC2/GMC(2); THM-2111 supplies an effective seed,
  THM-2101 two independent additive proofs, and THM-2067 the historical route.
- **Formal / HYP-8942:** `HeightWitnessSupplier`, orbit-product,
  irreducibility, additive incidence, and full-root Lagrange identities are
  root-imported; reciprocal monicization also constructs each ramified small
  root. General complex `DvdK1` is the sole endpoint premise; analytic wrappers remain.
- **HYP-8925/HYP-8930:** positive coefficients and a fixed-support unique
  channel prevent cancellation; neither is general `DvdK1`.
- **HYP-8932:** monomial membership gives nonvanishing; `{-2,-1,1,2}` is
  kernel-checked. `102/116` is bounded evidence with thirteen script-only rows.
- **HYP-8931 (MISTAKE-240):** its empty-face predicate makes the bypass vacuous.
- **HYP-8935 (MISTAKE-241):** reciprocal-monic Hensel now lifts individual
  roots concretely; branch products, analytic block transport, and the
  `DvdK1` wrapper remain.

## Other active lenses

- **HYP-8950 (OPEN JC SYNTHESIS):** Hamiltonian cokernel/fiber cohomology
  realizes weight and face obstructions; local-to-global termination remains.
- **HYP-8955 (CORRECTED BY THM-2102/2113/2127):** power-free faces, exact
  coprime two-face trains, and affine-root families close; general repeated
  proper-power descent remains open.
- **HYP-8945 (OPEN UNIT-DISTANCE ROUTE):** the sign-changing Bessel kernel
  identifies cancellation, not a new asymptotic bound.
- **Cross-domain incidence lesson:** HYP-2448's full barycentric root packet is
  invertible, whereas HYP-2430's design-incidence map has enormous nullity.
  Transfer must carry the rank defect and the lost-coordinate sidecar.
- **Tournaments/sequences:** test gaps, residues, events, proof obligations, or
  support-collision profiles before choosing runners as vertices.
- **Typed bridges only:** state source, target, map, preserved predicate, loss,
  restoration sidecar, and hostile control.

## Promotion rule

Record quantifiers, direction, mechanism, boundary, dependencies, hostile
witness, and non-consequences. Lean needs satisfiable hypotheses and root reach;
computation needs universe, controls, command, and output. Route `CLAIMED` and
`RESERVED` files as unproved candidates, never canon.
