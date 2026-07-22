## death-star-2026-07-22-S109 -- JC(2) PROGRESS: the Hamiltonian cokernel realizes S107's manufactured-valuation + S106's DvdK-face, and grounds them in fiber theory (NOT a proof)

**Owner directive:** work the S107/S108 ideas into progress toward a planar Jacobian-conjecture proof.

- **OPERATOR REFORMULATION (known, VERIFIED in-repo, mate_exists by exact linear algebra):** f is a Keller component <=> {f,g}=1 solvable <=> **1 in im(D_f)**, D_f={f,.} the Hamiltonian derivation; coker(D_f) = Brieskorn/de-Rham module, generic rank = dim H^1(generic fiber). So **JC(2) <=> every Keller component is a coordinate <=> [ f has a mate => generic fiber ~= C ]** (Kaliman 2002 / Abhyankar-Moh-Suzuki). Verified: coordinates x, x^2+y, y+x^2 have mates; x^2 (2 lines), xy & x+x^2y (fiber C*), x^2+y^3 (genus>0), x^3, homog-deg2 have NONE -- exactly matching 'fiber ~= C'.
- **(1) S107 manufactured-valuation REALIZED + VERIFIED = the WEIGHT OBSTRUCTION:** each positive weight w=(w1,w2) is a valuation at infinity; {f,g}=1 forces v_w(f)<=w1+w2, so v_w(f)>w1+w2 for some w => NO mate. Catches x^2 (w~(2,1)), x^2+y^3 (w~(3,2)), x^3; finds none for coordinates. = the JC analogue of S106/THM-2067's v(Pi)=1!=0; generalizes codex THM-2045's grading obstruction to every weight.
- **(2) S106 DvdK object REALIZED = the resonant faces:** xy/homog dodge (1), sitting ON a resonant w-face (v_w=w1+w2); there '1 in im(D_face)' for the w-quasi-homogeneous face collapses (C*-weight action) to a 1-VARIABLE CONSTANT-TERM condition = the DvdK/S101/S106 object; the sweep over weighted faces = boxeph S225 descent-termination; S106 orbit-product / boxeph-S231 monomial-certificate apply per face.
- **(3) RESIDUAL PINNED:** x+x^2y (fiber C*, no mate) DODGES every weight and is no single face => the obstruction is the GLOBAL coker(D_f)=H^1(fiber), not any one valuation = the JC(2) core; the weight/face obstructions are LOCAL shadows of the global fiber~=C.
- **PAYOFF (honest):** the S106/S107 route and the classical fiber-~=C route are the SAME obstruction coker(D_f) -- this DE-SPECULATES S107 by grounding it in established theory, and pinpoints where valuations bite (high/resonant strata) vs the C*-fiber residual (global, needs descent to terminate). NEXT tool-matched sub-target: prove local=>global for Keller components with a SINGLE resonant face, combining the (2) DvdK-face nonvanishing (S106/S231) with a (S225) coprime-interval Frobenius descent bound.
- **HONEST:** mostly a synthesis of known objects (Hamiltonian reformulation, Brieskorn module, Kaliman fiber theorem, THM-2045 grading); value = realizing S106/S107 as VERIFIED obstructions + unifying with fiber theory + locating the residual. JC(2) remains open. HYP-8950.
## death-star-2026-07-22-S108 -- asymptotic unit-distance problem: a positive count with a cancelling spectral kernel (exploration, NOT a new bound)
## boxeph-2026-07-22-S233 -- THM-2067 abstract contradiction assembled + t-adic closing over C(t); gap refined to THM-1550 alone (kernel-pure, HYP-8942)

**Owner:** long session, keep going until GMC(2) formalization complete; pull from agents.

**NEW kernel-pure ([propext,Classical.choice,Quot.sound]):**
- GMC2OrbitProduct.orbit_product_contradiction: transitive G-action + equivariant f + G-fixed subset product p + valuation v, v(prod_Omega f)=0 but v(p)!=0 => False. The full abstract THM-2067 contradiction (packages prod_pow_card_group_eq + valuation_zero_of_prod_fixed).
- GMC2RatFuncClosing.monomial_pow_ne_const: a*t^N != const in F(t) for a!=0, N>=1 (pullback F[t]<->F(t) injective + natDegree). The concrete t-adic closing.

**REFINEMENT of S232:** valuation is t-adic on C(t), NOT the splitting field -- both Pi=c*t and C_Phi=(-1)^d r0/rd are in C(t), so the closing is an elementary degree argument inside C(t). Gap shrinks to exactly THM-1550. Confirmed Mathlib wrapper hooks exist: Polynomial.Gal.galAction_isPretransitive, Gal.smul_def, MulDistribMulAction Phi.Gal SplittingField, Fintype instances -- so my abstract core DOES plug into Mathlib's Galois machinery.

**REMAINING INTERFACE:** (A) irreducibility of X^M-tR over C(t) [death-star; Gauss], (B) THM-1550 Pi=c*t [death-star's unramified-Hensel = THE gap; gives hfix + v(Pi)=1], (C) Vieta v(C_Phi)=0, (D) equivariance/wrapper, (E) Check A CT(Lambda^m)=[u^Mm]R^m. (A)+(B) hard/gap (death-star owns); (C)-(E) mine gap-free. Messaged death-star.

**Honest:** NOT complete -- completion gated on THM-1550 (death-star's Hensel piece). This session assembled the abstract contradiction + proved the C(t) closing kernel-pure, refined the gap, and specified the interface. 3 pushes. Artifacts: reflection thm2067-abstract-contradiction-...-boxeph-S233.md, HYP-8942, GMC2OrbitProduct.lean (+orbit_product_contradiction), GMC2RatFuncClosing.lean.

**Honest after correction:** the proposed LRC reduction is withdrawn; the
completed kernel-pure single-character DvdK1 leaf survives. Artifacts:
reflection doubling-homeomorphism-plus-mirror-parity-...-boxeph-S227.md,
HYP-8920, script (+.out), Lean GMC2DvdKTwoCharge.lean.
## codex-2026-07-21 -- HYP-8905 bridge audit and MISTAKE-229

The exact binary homogeneous Hessian calculation from S103 survives and lands
inside THM-2063. The claimed `NC2 -> GMC(2) -> JC(2)` chain does not: the
general symmetric reduction from JC(2) lands in four variables, and no
Gaussian-to-Laplacian predicate map was supplied. S225's VC(4), planar
leading-form/Jelonek, and Lame-for-polygons programs remain separate routes
sharing a descent heuristic, not equivalent formulations. The dimension-three
rank/cycle comparison is heuristic, not a classification of collisions.

## boxeph-2026-07-21-S225 -- working JC(2): the obstruction is a descent termination; coprime intervals are its tool (HYP-8905)

**Owner:** long session working to prove the planar Jacobian Conjecture; pull past threads creatively.
## death-star-2026-07-21-S103 -- Planar JC IS NC2 one-sidedness: 2D nilpotent Hessian => one-sided; the JC-true/false boundary = the GMC2 unique-vs-coincident-cycle threshold. HYP-8910.

**Owner directive:** work to prove the planar Jacobian conjecture, pull in past threads creatively.

- **CONTEXT:** target = JC(2) (Alpoge THM-1300 killed JC dim>=3); THM-1830 puts my NC2/GMC2 upstream (NC2=>GMC2=>...=>JC(2)); Zhao VC (Hess P nilpotent <=> Delta^m(P^m)=0) = the SAME moment-vanishing as GMC2's E[P^m]=0.
- **PROVED/VERIFIED (planar_jc_..._S103.py):** 2D symmetric Hess P nilpotent <=> harmonic (trace 0) AND det 0; for harmonic P=A z^d+B zbar^d, det(Hess)=-4 d^2(d-1)^2 A B |z|^{2(d-2)} => det=0 <=> A B=0 => P ONE-SIDED (= verbatim NC2 conclusion). One-sided P is harmonic (Zhao VC trivial) + ONE-FIBER-LINEAR (F2-iF1=-iz => codex THM-2063 tame). So the SYMMETRIC planar JC / 2D Zhao-VC case IS NC2 one-sidedness, PROVED.
- **THE BOUNDARY (the unification):** 2x2 nilpotent has RANK<=1 = ONE isotropic direction (one-sided); dim>=3 reaches rank>=2 = MULTIPLE isotropic dirs = RESONANCE = Alpoge counterexample (JC false). rank<=1 => parallel gradients => functional dependence => one-fiber pencil (codex). SAME threshold as GMC2 S101 (unique vs coincident cycle) + boxeph S217 entropy (rigid=zero-entropy=one-sided). Planar = the last one-direction dimension.
- **HONEST:** proves the symmetric case (=NC2 one-sided); NOT full JC(2) (non-symmetric rank-1 parallel-gradient = codex THM-2063 open crux). A unification placing planar JC in my GMC2/NC2/resonance framework + a bridge to codex's pencil.
- **ADOPTED codex MISTAKE-227:** my S102 strict measure detects M(S)>1/14 only (vanishes on the tight AP = measure-zero boundary packet, not a counterexample); the integer-kernel=GMC2 unification SURVIVES as a strict-BULK observable, boundary handled by THM-2058 exact packets. reflection planar-jc-is-nc2-one-sidedness-...-S103. HYP-8910.

## codex-2026-07-21 -- HYP-8879 strict-measure boundary correction

**HONEST:** JC(2) OPEN (JC(n>=3) FALSE, Keller THM-1300). Worked it via 3 reduction PROGRAMS (none complete): (A) de Bondt+Zhao => VC(4) (JC(2)<=>a dim-4 Laplacian moment nullcone, doubling n->2n); (B) klein-S329 Euler-Zariski cover-degree-3 bootstrap (cuspidal Jelonek curve, ramification-parabola escape, NO CF); (C) mac-mini-S137 golden-corner/Lame (subtractive Euclid on Newton slopes).
- **MISTAKE-227:** S102's integral with `g=1[||x||>1/14]` detects only
  `M(S)>1/14`. It vanishes on the tight AP `{1,...,13}`, whose weak safe set is
  the six-point packet `U_14/14`, so zero measure is not an LRC counterexample.
- **Surviving bridge:** after Fejer convergence control, the integer-kernel
  expansion is a strict-bulk observable. THM-2058's exact denominator packets
  handle equality; finite Sidon/AP ratios prove no AP-core reduction.
- **THM-2060 proved with prior-art boundary:** the clock-independent order is
  `q=a/gcd(a,w)`, and every tail histogram bin has the sharp floor
  `q-ceil(q/7)`. The qualitative one-tail dodge was already THM-760/761/765;
  the new gain is exact CRT support, multiplicity, and clock composition.
- **THM-2061 proved reduction:** the only primitive imprimitive-eleven-core
  residue is the dyadic two-odd-tail seam. Strict failure is an open folded-
  diamond cover of the closed core-safe set; it forces divisor completeness
  through `14`, tails below `12 max(C)`, measure at most `4/63`, and has no
  survivor for normalized cores in `{1,...,19}`.
- **THM-2062 proved atlas sieve:** deletion determinantal indices turn
  hereditary primitivity on every saturated two-anchor interval into an exact
  squarefree CRT wheel with at most two bad projective directions per prime;
  rank-one deletions become affine `+-1` terminals plus a 1D coprime wheel.
- **THM-2064 incoming synthesis:** the independent common-clock capacity
  theorem proves the full multi-tail union bound and the same unique dyadic
  two-tail seam; THM-2060 is now explicitly its sharp histogram specialization.
- **THM-2065 proved circuit-ray collapse:** THM-2051's bounded support-three-
  to-five relation pulls back to either a persistent coefficient-row circuit
  or one primitive projective parameter. Thus every circuit-free two-anchor
  template has only finitely many strict-null rows, filtered exactly by the
  THM-2062 wheel. Persistent height-`2^20` marked circuits are the residual.

**VERIFIED (jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.py):** 2D symmetric case EASY -- nilpotent Hessian <=> P prop (x+iy)^d (harmonic), Delta^m(P^m)=0, invertible. Lame worst-case = Fibonacci (longest Euclid chain <200 = (144,89)). [Restricted JC(2) proved elsewhere: THM-1345 equivariant, THM-1370 elliptic, THM-1365 poly-Galois, geom-degree<=2, THM-2063.]

**CREATIVE UNIFICATION:** the JC(2) obstruction has THREE equivalent forms -- (A) VC(4) both-signs radial nonvanishing; (B) leading-form descent {P_A,Q_B}=0 propagating to a coordinate (THM-1345 s5); (C) Lame-for-polygons CF bound -- and ALL THREE are DESCENT/RETURN-TERMINATION problems = exactly what my coprime-interval/numerical-semigroup/Frobenius engine (S223 DvdK, S224 Wall A) handles; Lame-Fibonacci = the effective bound. codex THM-2045 (smooth factorized R has NO planar Jacobian mate, an exponent-semigroup/Newton-edge obstruction) already uses this engine.

**BOUNDED:** VC(4) GMC-like (E=L o CT) but GMC(2)NOT=>JC(2) (doubling to rank>=2, S205); JC<->LRC 'shared n=12' WITHDRAWN (S137's 12 = Fibonacci proxy). POSITIVE PRIOR: Keller collision min = dim 3, so n=2 is BELOW threshold -> expect JC(2) TRUE + provable by low-rank/coprime-interval means.

**Corrections adopted (mining):** CF thread = mac-mini-S137 (not klein-S329); reify-ladder = deathstar-S75 (not defunct THM-1750); apolarity/Fischer != THM-1685/1710/1735 (those are TNC). NOT a proof of JC(2) -- an assembled route map + verified pieces + the unified descent-termination obstruction. Artifacts: reflection working-jc2-the-obstruction-is-a-descent-termination-...-boxeph-S225.md, HYP-8905, script (+.out).
- **PROVED:** for `S=aC union {w}`, safe phases on any clock `N` are the CRT
  fiber product of a core packet modulo `N` and a tail packet modulo
  `Na/gcd(w,Na)`. Reducing both to their common gcd turns existence into an
  exact histogram dot product and counts every safe `k/(Na)` grid phase.
- **Reach:** unlike THM-2057's automatic nonzero-residue argument, the formula
  works for `N>14`. The exact replay checked `53760` direct identities,
  `2903040` grid indices, `44761` small missing-clock specializations, and
  found `34854` larger-clock certificates.
- **Typed bulk/obstruction split:** the histogram dot product is exactly its
  positive zero mode plus nontrivial finite Fourier channels. The integer
  Cauchy test alone proves `14195` of the `14978` positive audit rows. This is
  the rigorous replacement for the unsupported modular-cusp language.
- **Finish route:** THM-2058's proved primitive phase-order counts populate the
  core histogram, and its longitudinal interval populates the tail histogram.
  A zero dot product exports disjoint residue
  supports to the signed Euler/deletion layer.
- **Assumption challenge:** the lossless carrier is a bipartite CRT
  compatibility graph. It is not a runner tournament or a modular cusp;
  orienting its symmetric ties destroys the theorem.
**Owner:** leverage the recent ideas to make progress toward LRC. Uses the RIGOROUS tools (not the cusp metaphor -- codex MISTAKE-226 accepted).

## codex-2026-07-21 -- concrete GMC(2) residue and conditional NC2 capstone checked

- `GMC2NormalizedResidue` now proves the complete normalized three-case
  calculation: non-dilated channels vanish by their multinomial, dilated
  off-face channels vanish by the factorial gap, and face channels reindex to
  the exact Frobenius power of the undilated face constant term.
- `GMC2SupportFaceBridge` proves exact geometric-face/support-face reindexing,
  identifies the specialized lifted seed with that constant term, and derives
  the global scaled floor, exact face height, and strict off-face gap from a
  concrete reference channel.
- `GMC2NC2` checks the whole post-specialization contradiction and exports
  `nc2_of_dvdK1_of_heightWitnessSupplier` plus the GMC(2) endpoint. The old
  one-`sorry` `GMC2NC2Capstone` WIP is replaced by a sorry-free compatibility
  surface. All audited declarations use only Lean's standard axiom trio.
- Direct checks of the concrete residue, support bridge, conditional capstone,
  and compatibility surface pass at the default heartbeat. The aggregator
  check stopped before elaboration at its pre-existing missing
  `GMC2GoodReduction.olean`; that dependency chain was not built under the
  no-big-build constraint.
- **Exact remaining boundary:** the reference-channel extractor and height
  theorem are separately green, but their direct existential wrapper into
  `HeightWitnessSupplier` deterministically exhausts elaboration (also at
  800k heartbeats) without a type or mathematical error. DvdK remains a
  separate published external premise. A final default-budget redesign
  compressed the reference to mass and balance: the compact extractor, base
  obligations, and base-contradiction adapter all elaborated, but the final
  `nc2_of_dvdK1` compositor still timed out at `whnf`; the experiment was
  reverted. No repository-wide build was run.
- **Live-main connection audited:** S223/HYP-8895 recasts the positive-
  coefficient one-variable DvdK return set as a numerical semigroup. That is
  a useful future formalization route, but it does not replace `DvdK1` here:
  the mixed-sign cancellation case still rests on S222's unfinished saddle
  argument. S103/HYP-8910 independently finds the same one-sided charge shape
  in a symmetric planar-Jacobian subcase, but explicitly supplies no transfer
  theorem to Gaussian moments and therefore does not change this Lean spine.
**THE MOVE:** assemble the structural theorems into an EXACT rational covering-min (upgrade over S206's float grid):
- THM-2047 s2 (PROVED): every maximizer t*=a/q has q|v_i+v_j (q<=2max) => M(S)=max over pair-sum vertices a/q, exact rational.
- S212/HYP-8845: covering => chi(G_delta) EVEN + mirror-symmetric => scan a/q in (0,1/2] (HALVING, verified).
- S223: candidate a/q coprime (three-distance/CF).

**VERIFIED:** M(deep well {1..12,182}) = 14/183 EXACTLY at t*=14/183, q=183=182+1=Phi_6(14) (pair-sum vertex, coprime CF [0;13,14]); 14/183>1/14 => LRC holds rigorously. SHARPENED disproof search (exact M): deep well 14/183, AP12+364 28/365, non-AP {1..11,13,168} 14/173, 2*AP 7/92 (non-primitive) -- ALL >= 1/14. No disproof.
## boxeph-2026-07-21-S221 -- the cusp frame is a repo-wide difficulty-LOCATOR (HYP-8885)
## death-star-2026-07-21-S100 -- CONFINE the GMC(2) DvdK dependency: elementary except in the resonant-signed corner (a verified sharpening, not a full bypass). HYP-8877.

**REDUCTION (the progress):** Wall A <=> every PRIMITIVE covering 13-set has some pair-sum vertex a/q in (0,1/2] with min_v||v a/q|| >= 1/14 -- the exact-arithmetic (residues mod q) form of the n=12 AP-core rigidity (S214 rank-11 vertex). Mirror halves the domain, pair-sum finitizes the vertices, coprime/CF names the target (q=Phi_6(14)).

**Honest:** rigorous covering-min tool + halving + finite exact-arithmetic reduction of Wall A + disproof-free confirmation of the tested class; NOT a proof of Wall A (the AP-core rigidity -- 'every primitive covering core has a lonely pair-sum vertex' -- is still the open crux). Converges with death-star-S101 (DvdK-free) + my S222/S223: both GMC and LRC reduce to exact residue/coprime-interval combinatorics. Artifacts: reflection leveraging-the-toolkit-an-exact-rational-covering-min-and-a-sharpened-wall-a-boxeph-S224.md, HYP-8900, script (+.out).
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
**Owner:** go back through the repo, apply the cusp frame to under-attended problems, show its power.
**Owner directive:** continue, merge in scaled cores and clocks.
## codex-2026-07-21-LRC-primitive-packets -- exact one-dimensional deck/fan carrier

- **THM-2058 PROVED:** every transverse safe packet splits uniquely by reduced
  phase order, with divisor summation, Mobius inversion, labelled unit
  transport, an exact Ehrhart/Beatty law, and a primitive discrepancy bound.
  For fixed bad `N` and THM-2055 hull owner, positivity and determinant failure
  leave one interval of coprime longitudinal coordinates, minus collision walls.
- **Load-bearing sidecar:** THM-2055 hull deletion applies only to determinant
  ownership. On the one-tail plane, hull representative `r=13` leaves strict
  `S_38` unresolved (`D_51=1/17`), while non-hull `r=10` recovers its exact
  `1/12` exit (`D_48=1/12`). Primitive packets are phase-order, not THM-2041
  frequency projectors, and they are not CRT-multiplicative.
- **Orbit-product transfer:** the TNC monodromy norm becomes an exact
  unit-stabilizer identity for primitive packets. At `N=27` the stabilizer is
  `{1,-1}`, giving nine packet images with incidence one; the same calculation
  proves that an unlabelled packet loses signed orientation and creates no seed.
- **Exact cusp replacement:** phase height gives a proved bulk/boundary/null
  trichotomy. Strict templates occur primitively on all sufficiently large
  prime grids; tight templates have finite support on level-14 pair-sum clocks;
  subthreshold templates have no packets. This is the lawful pointwise carrier,
  unlike the retracted modular-form attachment in MISTAKE-226.
- **Unique-channel no-go:** at `N=29`, the template
  `(1,2,3,4,5,6,8,9,10,11,12,13,14)` has unique antipodal maximizers `+/-4`
  but height only `2/29`. HYP-8878's uniqueness principle therefore schedules
  singleton cells but cannot replace the LRC height coordinate.
- **Referee repair:** the theorem now scopes transverse coefficients to nonzero
  integers, cites the arity-free THM-2047 pair-sum theorem and the exact
  `THM-1065-doubling-family-mod-six-characterization.md` source for the
  Goddyn--Wong boundary, and tests the zero-measure `S_24` period-14 packet.
- **Open residual:** the fixed-star enumeration is exact and one-dimensional;
  uniform THM-2052 atlas compression and pair-sum/Euler/relative-Fejer
  discharge of surviving rows remain open.

## death-star-2026-07-21-S98 -- NC2 capstone: skeleton of `DvdK1 → NC2` typechecks (architecture validated); full completion plan worked out + reference-channel friction resolved. HYP-8805.

**REDUCTION:** DvdK -> standard analytic combinatorics (dominant-saddle nonvanishing) + THM-1840 periodicity + S208 confluent cusp -- all effective, none DvdK's machinery. Makes GMC(2)'s angular/Eisenstein floor (S221) DvdK-free + effective (yields the open effective-DvdK bound m0).
## boxeph-2026-07-21-S221 -- the cusp frame is a repo-wide difficulty-LOCATOR (HYP-8885)

**THE FRAME as a diagnostic:** object = EISENSTEIN (computable floor/main term/local) + CUSP (hidden obstruction = genus = deep arithmetic entropy S218). Difficulty is always the CUSP; the frame localizes it + predicts the first-hard-case = first positive cusp dim.

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
**HONEST:** split is weight-2 X0(14) NOT L(S); VALUE modular / EXISTENCE+constant COMBINATORIAL (floor constant NOT L(14a,1)/L(sym^2)/period -- in the descent); MISTAKE-087 (Phi_6 non-extremal); no weight-1 dihedral construction. Artifacts: reflection the-lrc14-obstruction-is-the-first-cusp-form-scaled-cores-clocks-are-cusps-boxeph-S220.md; HYP-8880; scripts S220 + S219 (+.out).
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
**Honest:** THM-2053 + kps-S17 standing; Heegner/Paley facts classical+verified; my contribution = the synthesis/leverage-framing (residual = binary form, Heegner rigidity, rank-or-Euler = isotropic/anisotropic). A target, not a closure. Ties THM-2053+THM-2052+kps-S17+S212+S214+S215. Artifacts: reflection the-anisotropic-determinant-gate-...-boxeph-S216.md, HYP-8865, script (+.out).
**INCOMING SYNTHESIS:** THM-2052's new two-anchor refinement says the
rank-eleven branch is a finite list of one-projective-parameter stars. THM-2053
now cuts every star to the explicit finite primitive disk `||d||<91L`.
HYP-8855 correctly identifies the AP as the achiral rank-eleven boundary
vertex, but its tournament orientation remains a lens: finite discharge must
retain the signed coefficient heights, pair-sum rulers, and endpoint owners.

**ASSUMPTION CHALLENGED:** a twelfth bounded relation is not the only route to
finiteness, and HYP-4346's existence of some escaping direction has the wrong
quantifier for a specified row. The geodesic estimate is uniform in the target
direction. The remaining problem is finite but not yet executed: reduce bases,
compress the two-anchor atlas, and prove HYP-2896-style resonance fans or exact
Euler certificates throughout the disks. LRC(14) remains open.

**ARTIFACTS:** THM-2053, HYP-8846, updated HYP-8841/frontier/reflection.

**LATE PULL:** THM-2054 reserves relative Fejer decorrelation along a character
line. It is the right analytic partner for the tangent disks: discharge whole
off-resonance cells, then spend exact pair-sum/Euler work only on resonant
lattice points. HYP-8860 organizes the useful primes but remains a diagnostic
modulus lens; it does not preserve signed star coefficients or endpoint owners.

## boxeph-2026-07-21-S214 -- the rank-11 AP-core is the achiral vertex where codex's rank-or-Euler frontier meets (HYP-8855)

**Owner:** understand primes 7,5,11 as well as 2,3, through 'a set of pairwise relations is a tournament'.

**UNIFYING LAW (verified, understanding_primes_via_paley_tournaments_boxeph_S215.py):** prime p IS its Paley object (i->j iff j-i is a QR mod p). p=3 mod4 (3,7,11) -> Paley TOURNAMENT: self-converse, vertex-transitive, doubly-regular, char_A=(x-(p-1)/2)(x^2+x+(p+1)/4)^((p-1)/2), quad disc -p, roots (-1+-i sqrt p)/2 = quadratic Gauss sum i sqrt p. p=1 mod4 (5,13) -> self-complementary Paley GRAPH, real spectrum (-1+-sqrt p)/2 = Gauss sum sqrt p. p=2 = the reversal INVOLUTION (Phi_2=x+1, chirality, S210-S213). Verified: Paley-3=(x-1)(x^2+x+1)=Eisenstein; Paley-7=(x-3)(x^2+x+2)^3 (THM-1830); Paley-11=(x-5)(x^2+x+3)^5; |g_p|^2=p (i sqrt p for 3,7,11; real sqrt p for 5,13).

**PERIODIC TABLE for LRC(14)=2*7:**
- 2 = the INVOLUTION (chirality/reversal iota:t->1-t, S212/S213; Phi_2 the antipode).
- 3 = the ATOM: Paley-3 = the 3-cycle, char x^2+x+1 = Eisenstein; the argmax Phi_6(14)=183, t*=14/183.
- 5 = the GOLDEN FOIL (=1 mod4): self-comp GRAPH, real sqrt5 = Q(sqrt5) = Fibonacci = the LRC loosest/foil (S206); also Bonferroni depth.
- 7 = the APEX (14=2*7): Paley-7 tournament, i sqrt7 = my S212 Euler-branch index; Phi_7/Phi_14 hardness; F_7[C_14]=F_7[X]/(X+-1)^7 (THM-2043, verified x^14-1==(x-1)^7(x+1)^7 mod 7); cap field Q(cos2pi/7) cubic disc 49.
- 11 = the RANK (=3 mod4): Paley-11 (i sqrt11), the rank-11 AP-core/relation code (S214); SCARCE (only multiple<=14 is 11) => forced rigid speed.
1. TOURNAMENTS EVEN: A000568=SC+2*(chiral pairs). THM-587 (PROVED): P_n(1)=A000568, P_n(-1)=SC=self-converse=antipodal Euler/Lefschetz number; SC=2,2,8,12,88,176 (n=3..8) all EVEN => A000568 EVEN for n>=3. Enumeration n<=6: (total2/SC2/0pairs),(4/2/1),(12/8/2),(56/12/22), count=SC+2pairs holds. THM-1830 BLUE self-comp atom=1 iff n odd = local single-fixed-point face.
2. TOOTHPICK A139250: sim=OEIS exactly; D4 mirror sym; =(#axis)+2*(#pairs); #axis ODD (central seed) => A139250 ODD all n>=1; first differences have exactly ONE odd term (seed); A139250(2^k)=(2^(2k+1)+1)/3=1,3,11,43,171 Jacobsthal.
3. LRC (S212): chi(G_delta) == #iota-fixed lonely pts (mod 2).
- **HONEST FRAME:** GMC2Reduction.lean already had the easy direction + `NC2At`/`ChargeOneSided`; the whole remaining problem is `nc2 : ∀P, NC2At P` (= THM-2022). Full sorry-free completion is BLOCKED on §2 (number-field descent, HEAVY) + §3 (Duistermaat-van der Kallen THM-1630, NOT in Mathlib = a citation) -- genuinely multi-session, NOT done this session. What I did: formalize the self-contained arithmetic engine kernel-pure.
- **FORMALIZED kernel-pure** ([propext, Classical.choice, Quot.sound], no sorry/native_decide): (1) **architecture `gmc2_of_nc2`** -- GMC(2) ⟸ NC2, sorry-free (the whole problem rests on the single theorem NC2); (2) **§1 `wick_expansion`** -- E(P^m) = Σ multinomial·∏coeff^k·wt(radial), the exact channel sum M_m, via Mathlib `Finset.sum_pow_eq_sum_piAntidiag` + E-linearity (E_add/E_monomial/E_sum) + prod_monomial; (3) **§5 `sum_natCast_mul_pow_char`** -- (Σ w_s g_s)^p = Σ w_s g_s^p in char p (Frobenius fixes natCast weights ⟹ face survives as Q̄^p); (4) **§4 `multinomial_dilate_modEq`** -- multinomial Lucas (Mathlib had only binomial), assembled via multinomial_insert + Choose.choose_mul_mul_modEq_choose_nat.
- **FIXED codex's `GMC2FrobeniusFace.lean`** (§3-4 face geometry) -- it FAILED to build (linarith nonlinear gap: lambda·charge i vs lambda·charge j; substitute the charge equality first). File not in aggregate build, so break was invisible to `lake build`.
- **Mathlib API survey (verified, recorded in reflection):** Lucas/Kummer/Frobenius/multinomial-theorem/Zariski all located; multinomial-Lucas + DvdK were the only NOT-IN-MATHLIB gaps (multinomial-Lucas now assembled; DvdK stays a citation).
- **REMAINS for complete nc2:** §2 descent (Zariski `finite_of_finite_type_of_isJacobsonRing` + number-field residue Frobenius, HEAVY but Mathlib-stocked); §3 DvdK (cite); §4 no-carry channel-survival wrapper (MODERATE); final contrapositive assembly. reflection formalizing-thm-2022-...-S94. HYP-8805.
## death-star-2026-07-21-S93 -- PR packaging of three-term/Hermite formalization (HYP-8805)

**SYNTHESIS:** p mod 4 decides tournament(i sqrt p)-vs-graph(sqrt p), tight(Eisenstein 2/3/7)-vs-slack(golden 5); the Gauss sum IS the Paley spectrum, fixing each prime's role. Ties S206(foil)+S212(i sqrt7)+S213(chirality)+S214(rank-11)+opus-S434(Paley symmetric-intransitive pole, opposite the transitive AP nullcone vertex).

**Honest:** all rows are verified/classical facts; the contribution is the SYNTHESIS -- one Paley construction giving each small prime a tournament personality via p mod 4, mapped onto LRC(14)'s tight/slack/apex/rank structure. Artifacts: reflection each-prime-is-its-paley-tournament-...-boxeph-S215.md, HYP-8860, script (+.out).
**Honest:** THM-587, the toothpick OEIS/oddness/Jacobsthal facts, and S212 are each proved/verified; the contribution is the UNIFICATION (one reversal, one Lefschetz parity count==#fixed) + the explicit toothpick<->chirality dictionary. Artifacts: reflection chirality-toothpicks-and-why-tournament-counts-are-even-one-lefschetz-parity-boxeph-S213.md, HYP-8850, script (+.out).
**Owner:** push the DC(2)/planar-JC thread to its next decisive target and transfer the proved/disproved mechanisms toward LRC(14), while repeatedly integrating incoming work.

**DC(2) result (THM-2049, PROVED local/formal statement; not DC(2)):** in the exact Ore algebra `Q[x,q][ell;delta]`, `beta(sum a_k ell^k)=min_k(v_x(a_k)-2k)` is multiplicative and commutators raise beta by two. The associated bracket is `{ell,q}_0=2`. For the Weyl boundary symbols, the simultaneous grade-`g` correction map is `(A,B)->(8/3)(u-2)A+(2u^2-10u+9)B/9`; it is surjective because the two `u` polynomials are coprime. Thus the grade-six residual is exact. An exact ladder advances grades `6,...,13` to `14`; a formal beta-adic `[S,T]=1` lift exists. The open gates are polynomial termination and the coupled `D` relations. This corrects HYP-8802/8803's earlier suggestion that the first invariant grade might carry the obstruction.

**LRC no-go (THM-2050, PROVED):** AP13 and `AP13` with `12->26` have identical full local phase-height function germs on `|h|<1/728` at every unit point `a/14`, yet `M=1/14` and `M=1/12`. Local top data, even as a full germ, cannot determine global loneliness.

**Incoming synthesis:** THM-2047 supplies the lossless signed phase-height/Euler carrier; THM-2048 supplies the fiber-quantization pruning tax; HYP-8840 identifies GMC's constant-term/volume leverage and its zero-volume ceiling. Later pulls supplied THM-2048's genuine Cover14 gain and strengthened THM-2051: after paying pair covariance exactly, the no-support-`3..5`-relation branch at height `2^21` has positive safe volume. THM-2052 then proves every hypothetical counterexample already has eleven independent bounded three-support relations and lies in a finite two-dimensional rational atlas. The exact transfer is `volume/tax -> strict branch`, `Euler signed wall word -> tight branch`, and a labelled Noetherian rank-or-Euler rule inside the relation branch as the missing glue. No literal algebra map between GMC/DC and LRC is asserted.

**Exact termination-sidecar audit (HYP-8841):** pair-sum maxima, threshold interval/point topology, complete first exits, and every peel tax were computed on AP/GW, `12->26`, `12->36`, `12->96`, `12->84`, P10+K33, and the incoming Cover14 tax-gain row. It exactly reproduces the latter's peel-`93` excess `2413467317/235670635200`. The tax fires on deep/covering controls but misses the smallest hostile/K33 controls and is not a scalar termination height. THM-2047 proves the strict search is complete by `q<=2 max(S)`. With THM-2051 now proved, the remainder lies in the bounded small-relation branch, but every control already has a height-one support-three circuit (the first seven share `1+2=3`). Raw circuit existence ranks last in the enlarged carrier tournament. The next decisive Wall-A clause is therefore an active-owner circuit/Euler endpoint-survival lemma when neither a tax violation nor a positive pair-sum margin occurs.

**Past-work pull into the decisive clause:** HYP-2108 already gives the exact endpoint-cover functional. If a core-safe component has midpoint `m_i` and length `l_i`, it is swallowed by the open danger arcs of peel `w` iff `||wm_i||+(w/2)l_i<1/14`. Hence the active-owner target is precisely `P_w(C)=max_i(||wm_i||+(w/2)l_i-1/14)>=0`. With THM-2052, the next theorem becomes a rank-or-Euler alternative: some peel has `P_w>=0`, or active owner data supplies a twelfth independent bounded relation and reaches the finite maximal-minor terminal. HYP-3117/HYP-3120 identify endpoint incidence as the missing proof-circuit input; HYP-8845 halves the covering case because a survivor has a mirror partner and `chi>=2`.

**Artifacts:** THM-2049, THM-2050, HYP-8841, the updated exact Ore script/output, `lrc14_termination_sidecar_codex_20260721.py/.out`, and reflection `from-Ore-boundary-acyclicity-to-LRC14-Euler-termination-codex-20260721.md`.
**Owner:** look for other topological advances the repo has made; come up with creative LRC arguments combining and extending them.

**PULL (repo topological toolkit, credited):** THM-2047 (codex) chi(G_delta)=#components, LRC(14)<=>chi(G_{1/14})>0 [PROVED]; HYP-3015 (codex-S179) {G_delta}=superlevel filtration of f_S, M(S)=top death, persistence barcode; opus lonely-set Euler-char certificate + kps cohomological_three_distance Alexander duality (b0(lonely)=b0(cover), arcs alternate); HYP-3025 arc-Cech nerve + Betti-defect sidecar, HYP-3101 normal-fan Cech barcode component bound (open); kps-S19 LRC Lefschetz = free iota:t->1-t + Gauss sum i*sqrt(7) (ordinary Lefschetz blind); THM-587 metagraph reversal Lefschetz (tournament side, PROVED). P1-P5 re-verify these.

**NEW (P6, verified):** f_S(1-t)=f_S(t) => G_delta is iota-INVARIANT (iota = my S210 involution). iota's fixed points {0,1/2}: f_S(0)=0 always, f_S(1/2)=0 iff some speed EVEN. Every COVERING set has an even speed => both fixed points dangerous => iota acts FREELY on G_delta => chi(G_delta) EVEN. So codex's chi>0 sharpens to:
  LRC(14) for covering S  <=>  chi(G_{1/14})>=2  <=>  a MIRROR PAIR {t*,1-t*} of lonely windows survives.
Verified: deep well {1..12,182} chi=24 (12 pairs); tight (1,2,3)@1/4={1/4,3/4} chi=2; all-odd (1,3,5,7) = iota-FIXED exception (1/2 lonely, chi=1 ODD = Borsuk-Ulam fixed point, the classical all-odd-lonely case).

**LEVERAGE:** (1) equivariant HALVING of Wall A -- find one lonely window in [0,1/2], mirror automatic; (2) parity obstruction -- chi even => never 1; a disproof needs chi=0 = every mirror pair killed simultaneously (iota-symmetric covering); (3) kps-S19's Lambda(iota)=0 blind BECAUSE free; the odd-equivariant index = Gauss sum i*sqrt(7) is the Borsuk-Ulam obstruction on G_delta/iota. The equivariant chi (even) + odd index (i*sqrt7) are the two halves of the Z/2-equivariant Euler class of the good set. Topological form of THM-1820 mirror pairs (B3) + S210.

**Honest:** P1-P5 re-verify existing fleet toolkit (credited, not claimed); the equivariant even-chi mirror-parity sharpening (P6) is new and verified. Forcing chi>=2 for every 13-speed covering core = LRC itself, OPEN; reduces Wall A to 'G_{1/14}(C)/iota on [0,1/2] nonempty'. Artifacts: reflection the-good-sets-reversal-symmetry-...-boxeph-S212.md, HYP-8845, script (+.out).
**Honest:** a verified DECOMPOSITION of the remaining work + a proven ceiling on the GMC/volume half, NOT a proof or new implication. Corrects+adopts codex THM-2047; retracts the S209 toric-complement and S210 'torus needs saddles' overclaims. Artifacts: reflection where-gmc2-reaches-lrc14-...-boxeph-S211.md, HYP-8840, script (+.out).
> **CURRENT-TRUTH WARNING:** this is a concise provenance log, not a truth
> source. Later audits can correct an earlier session on the same page. Start
> with `START-HERE.md`, `CURRENT-FRONTIER.md`, `ACTIVE-GUARDRAILS.md`, and the
> linked theorem or mistake entry.
> **CURRENT-TRUTH WARNING (2026-07-21):** This file is a concise provenance ledger, not a proof authority. Read `CURRENT-FRONTIER.md`, `01-canon/ACTIVE-GUARDRAILS.md`, and the cited theorem or mistake file before reusing a claim. Historical session narratives live in git and may contain claims later corrected below.

# Recent research sessions

This ledger records only changes that affect how a new mathematical session should
start. It deliberately omits chronological play-by-play. The proof files are the
authority; this file says why the current routing looks the way it does.

## 2026-07-21 — audited LRC(14) terminal stack

The live theorem chain was reread from its proofs and computations rather than
copied from earlier summaries.

- `THM-2051` gives a finite-geometry reduction: a counterexample has a genuine
  support-3, -4, or -5 integer relation of height at most `2^20`.
- `THM-2052` turns rank 11 into finitely many two-anchor stars and rank 12 into a
  finite box. These are reductions, not a completed enumeration.
- `THM-2053` constructs the exact transverse residue deck. Adjacent normalized
  columns are independent of Ungar; an exposed pair with sum divisible by `N`
  gives `M(v) >= 1/13 - R/(2N)`. Its determinant gate and 26 tangent disks are
  sufficient certificates. Gate failure means “uncertified,” not “bad.”
- `THM-2054` proves an abstract relative Fejer inequality. Its `H=2^19`
  calculation is numerically compatible with the desired scale, but the
  LRC-specific plateau and routing hypotheses remain open.
- `THM-2055` supplies the precise normal-fan reduction, including basis, deletion,
  tie, non-hull, and radial-sector scope.
- `THM-2056` identifies the determinant defect after Kelvin inversion with an
  exact Farey-neighbor defect and supplies an existential finite acute unimodular
  fan. It is a certificate geometry theorem, not an LRC closure and not a
  Heegner phenomenon.
- `THM-2057` proves a general missing-clock tax and closes two full families:
  `{a,2a,...,11a,13a,w}` through the `12a`/`14a`/`84a` leaves, and
  `{a,2a,...,12a,w}` through the `13a`/`14a`/`182a` leaves. These are two
  structured planes, not the general 14-runner case.
- `THM-2059` gives the exact generalized-CRT join of core-safe and tail-safe
  packets for any one-tail clock. Its histogram dot product counts compatible
  classes; the `Na/lcm` lift counts safe grid phases. Zero overlap rejects only
  that clock grid.
- `THM-2058` is an empty reservation. It contains no theorem, proof, script, or
  result and must never be cited as progress.

Therefore **LRC(14) remains open**. The most compressed current obligation is to
route every residual two-anchor star either to the transverse-deck determinant
certificate, to the relative-Fejer hypotheses, or to a genuinely new certificate.
`HYP-8871` records one open typed deck/Farey/clock state graph for doing so.
`HYP-8900` adds an exact per-row pair-sum evaluator, but does not make the
unbounded covering family finite or prove its restated Wall-A condition.

## 2026-07-21 — corrections to recent synthesis sessions

The recent sequence, entropy, and modular explorations produced useful objects,
but several interpretations outran their evidence. The correction ledger now
routes them as follows.

- `MISTAKE-227`: the arithmetic-progression computation used index `11!`, not the
  intended quantity. Do not reuse its numerical conclusion.
- `MISTAKE-228`: Paley tournaments do not assign LRC runner roles by themselves.
- `MISTAKE-229`: the observed determinant/Farey geometry does not justify a
  `-7`, Heegner, or complex-multiplication interpretation.
- `MISTAKE-230`: inverse classes are inversion orbits; rational primes do not
  support the claimed universal “arithmetic entropy.”
- `MISTAKE-231`: entropy is observable-relative. The surviving exact data include
  distinct fibers and characteristic polynomials, not a universal scalar law.
- `MISTAKE-232`: there is no Paley tournament at composite order 14;
  `lcm(12a,14a)=84a`; and the proposed Frobenius language was inapplicable. The
  exact `S_7`/`A_13` matrix facts remain separate facts.
- `MISTAKE-233`: scaling an LRC speed set does not create a `Gamma_0` level;
  modular cusps are not missing clocks; the weight-two level-14 newform has
  rational coefficient field, not a forced `Q(sqrt(-7))` period field; and
  `X_0(11)` already has genus one. No modular bridge to LRC(14) was proved.
- `MISTAKE-234`: a return semigroup records support reachability, not
  mixed-sign noncancellation. The aperiodic polynomial
  `z-z^(-1)+z^2-z^(-2)` has every return length from two onward but zero
  constant term at every odd power.
- `MISTAKE-235`: LRC and GMC have differently weighted kernel-indexed sums,
  not one noncancellation object. S102's four truncated low-dimensional rows
  prove no Fourier-tail bound or AP-core reduction; the tight non-AP row
  `{1,...,11,13,24}` already defeats the proposed premise.

The S99, S102, and S217–S223 artifacts are retained as quarantined historical
experiments. Their corrected exact computations may seed new work, but their
refuted interpretations may not be promoted without a new bridge.

## 2026-07-21 — NC2 formalization boundary

The imported module `TournamentH7.GMC2Formalization` is sorry-free at the proved
nodes and formalizes the root-imported part of the GMC2/NC2 argument. The separate
`TournamentH7.GMC2NC2Capstone` module is a typechecking capstone skeleton with one
remaining `sorry` and is not imported by the proved root. The cited external
theorem has not thereby been formalized. Startup text must preserve all three
distinctions.

### Citation-minimization routes landed later that day

- `HYP-8878` isolates an elementary DvdK-free case: a unique minimum-mass
  balanced channel makes the first relevant constant term one nonzero monomial.
  The mechanism is proved; `98/116` (`84%`) is only its finite size-3/4 support
  census in `[-4,4]`, not an asymptotic density.
- `HYP-8877` is the coarser edge/positive-coefficient precursor. Both results
  can shrink the cited branch in a future formalization but do not eliminate
  the coincident-minimum-channel cancellation locus.
- `HYP-8890` proposes a saddle/Watson replacement for DvdK. The positive-real
  saddle case and periodic examples are useful controls; the general complex
  dominant-saddle and noncancellation theorem is still missing.
- `HYP-8895` exactly packages support return lengths as a semigroup. By
  MISTAKE-234 it is citation-free only for positive coefficients or after a
  separate noncancellation certificate such as HYP-8878.
- `HYP-8885` reuses “Eisenstein/cusp” as a difficulty-locator. Its exact
  tournament and kernel computations survive, but MISTAKE-231/233 forbid
  treating cusp, genus, or entropy as a common object or `f14` as LRC data.

## Durable adjacent results worth connecting

- `THM-2047` and `THM-2050` delimit what the current carrier/no-go machinery can
  and cannot certify. A failed carrier is information about the missing invariant.
- `THM-2044`, `THM-2045`, and `THM-2049` provide useful tournament and sequence
  structure only within their stated scopes; none is an automatic LRC bridge.
- Sequence profiles should be treated as families of observables—fibers, moments,
  recurrences, Dirichlet series, transforms, and singularity data—not as one
  privileged entropy number.
- Tournament analysis is most useful when the vertex set, pairwise observable,
  switch/gauge, preserved predicate, lost information, and tie Hamiltonian path
  are all declared. Runners are only one possible vertex set.

## Process change made durable

Every exploratory mathematical session now carries three simultaneous lanes:

1. an **Anchor** tied to the main live obstruction;
2. a **Niche** thread selected because it is underexplored or structurally odd;
3. a **Wildcard** prompted by a genuine mathematical compulsion.

For each new object, ask what it preserves, what it forgets, why the observed
claim is true or false, and which other active object it can compose with. Record
mechanisms and failed bridges, not just outcomes. Turn certificate failure into
an address with sidecar data; promote successful meta-patterns into
`META-PATTERNS.md`; and keep current-truth summaries short enough that future
agents read the proof rather than inherit folklore.
