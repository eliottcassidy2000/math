## death-star-2026-07-21-S92 -- GMC(2)/NC2 formalization → Mathlib-submission readiness: verified kernel-pure + extracted the GENERAL three-term no-common-root (new to Mathlib, PR-ready). HYP-8805.

**Owner directive:** get NC2/GMC(2) results formalized to Mathlib-submission quality; pull/push often like while waiting for builds.

- **VERIFIED KERNEL-PURE** (`#print axioms` = [propext, Classical.choice, Quot.sound] only, NO sorryAx/native_decide): GMC2.mathieuZhao_of_charge_pos (the NC2⇒GMC(2) charge-arithmetic step), GMC2Hermite.no_common_root, ThreeTerm.exists_nonvanishing. They clear Mathlib's hardest gate.
- **NEW: extracted + GENERALIZED the flagship** `ThreeTerm.no_common_root` (monic three-term recurrence, b_n≠0 ⟹ no two consecutive members share a root) from ℝ to ANY integral domain (linarith→linear_combination hrec; [Zero R] on the structure; mul_eq_zero for the domain step); autoImplicit false, docstrings, ℝ-Hermite instance → `TournamentH7/ThreeTermRecurrence.lean`, builds clean (40s), kernel-pure. Mathlib HAS Polynomial.hermite but NOT the general three-term no-common-root ⟹ genuinely PR-READY (the abstract Favard/orthogonal-poly no-common-roots core).
- **HONEST:** NC2/GMC(2) as FULL theorems remain OPEN (the radial-channel wall, S87-S91). 'Our results' for Mathlib = the PROVED reductions/lemmas. mathieuZhao_of_charge_pos is kernel-pure + PR-viable with light packaging. GMC2MomentBasics decide-instances are correct/kernel-pure but specific (stay in-repo).
- Remaining PR work (scoped, non-blocking): Polynomial-R recast + connect hermiteReal to Mathlib.Polynomial.hermite; minimal-import trim; placement decision. reflection gmc2-nc2-mathlib-submission-readiness-...-S92; new ThreeTermRecurrence.lean wired to root.

## death-star-2026-07-21-S91 -- Dividing by the common factorial (pA₀)! turns NC2 into ONE tournament-discriminant condition: NC2 ⟺ the (confluent) Vandermonde of channel radial degrees ≠ 0. NC2 wall OPEN. HYP-8795.

**Owner directive:** think about dividing E[P^m] by the full common factorial (pA₀)!.

- **DEVELOPED (verified p=q=1, p=2,q=1).** The dominant common factorial is (p·j_max)! = (pA₀)!, A₀=max multiplicity of the top charge atom Z^p a. Dividing E[P^m] by it: makes the ENDPOINT channel O(1) (codex 'endpoint ratio one'); the endpoint = a SINGLE nonzero radial term (no self-cancel) ⟹ noncancel IFF strict source (degree-gap=transitive, S88) — codex THM-2017 as a one-liner; peels boxeph THM-2033's ∏a_i! leaving the PURE Vandermonde of channel radial degrees (klein THM-1805).
- **RESULT:** NC2 ⟺ the (possibly confluent) Vandermonde of channel radial degrees ≠ 0. Distinct degrees (transitive) ⟹ ≠0 by inspection (=THM-2017 degree-gap); repeated (regular/Paley wall, S89) ⟹ confluent Vandermonde = central trinomial (S90)+hyper-Bessel = the open residual (Laguerre-Pólya = Paley-spectrum-reality).
- **PAYOFF:** dividing by (pA₀)! SEPARATES the trivial factorial growth (=common factorial, explaining EMP depth-growth THM-1790) from the hard SIGN structure (=discriminant); NC2 stripped to its algebraic core, difficulty isolated at the regular/Paley confluent limit. Cleanest tournament↔NC2 statement yet.
- Synthesis not proof; wall OPEN. Credits boxeph THM-2033, codex THM-2017, klein THM-1805. reflection dividing-by-the-common-factorial-...-S91 (+out).

## death-star-2026-07-21-S90 -- The NC2 tied-core weights ARE the CENTRAL TRINOMIAL (A002426) = a free-probability moment; completes the tournament↔NC2 free-prob bridge (S88→S89→S90). NC2 wall OPEN. HYP-8790.

**Owner directive:** keep finding tournament↔NC2 connections; push/pull often.

- **IDENTIFICATION (verified).** At the NC2 resonance central offset (= fully-confluent Vandermonde = fully-regular/Paley tournament, S89), the channel weights sum to W(m)=Σ_i m!/(i!²(m-2i)!) = 1,3,7,19,51,141,393,1107,3139,8953,... = **A002426 central trinomial** = [x⁰](1+x+1/x)^m = the m-th MOMENT of a 3-atom free convolution (Wigner/free-prob), ratio→3, ~3^m/√(πm).
- **CLOSES the free-prob bridge** boxeph-S203 flagged: the free-moment sibling of THM-438 (Paley cluster integrals=Catalan=free CUMULANTS of ½(δ_a+δ_{−a})). The wall=regular/Paley tournament carries a free-probability law on BOTH its H-count (THM-438) and its NC2 channel weights (this).
- **SHARP:** NC2 fails on the wall ⟺ central-trinomial-weighted signed channel sum =0 ∀m ⟺ free-cumulant series has a real positive zero ⟺ Laguerre-Pólya failure (boxeph-S202) ⟺ Paley spectrum leaves Re=−1/2 (char_S=∏(x²+p), THM-1555/213) — 3 faces (combinatorial/analytic/spectral) of ONE tournament fact.
- **ARC S88→S90:** channels-form-a-tournament (S88/boxeph THM-2033) → wall=regular/Paley (S89) → wall-weights=central-trinomial=free-prob (S90). Synthesis not proof; NC2 wall OPEN. Credits boxeph THM-2033/S202, codex, THM-438, klein THM-1805. reflection the-nc2-tied-core-weights-are-the-central-trinomial-...-S90 (+out).

## boxeph-2026-07-21-S203 -- THE VANDERMONDE IS THE BRIDGE: tournaments <-> NC2 (THM-2033)
## death-star-2026-07-21-S89 -- The NC2 WALL IS the regular/Paley tournament: completing boxeph's bridge + unifying EVERY repo wall (NC2 = H≥disc = LRC) into ONE object. NC2 wall OPEN. HYP-8785.

**Owner directive:** keep finding tournament↔NC2 connections; push/pull often.

- **BUILDS ON boxeph THM-2033** (NC2 channel-det = ∏a_i!·Vandermonde(radial degrees) = signed tournament sum, klein THM-1805; distinct=transitive=noncancel, repeated=confluent wall).
- **MY STEP (verified):** channel radial degree D(i)=i+i·degA+(m-2i)·degB. DEGREE-GAP → D(i)=[0,3,6,9,...] DISTINCT (transitive Vandermonde≠0). RESONANCE CENTRAL OFFSET (degA=degB=1) → **D(i)=m for EVERY i** (fully-confluent). By klein THM-1805 (transitive⟺distinct scores): repeated degrees=repeated SCORES=REGULAR tournament; ALL equal = DOUBLY-REGULAR = PALEY/DRT. So **NC2 resonance wall = fully-confluent Vandermonde = the regular/Paley tournament.**
- **UNIFICATION (the payoff):** NC2 wall = H≥disc wall (S84 regular/Paley tightest) = LRC wall (THM-640 AP=Paley) = ONE object, the regular/Paley (equal-score, big-stabilizer, S76) tournament. Transitive=easy pole (distinct scores), regular/Paley=hard pole (equal scores) = the two S75 poles, shared across 3 flagship problems.
- **ANALYTIC FACE:** the fully-confluent (regular) channel sum's asymptotic = Wigner/free-cumulant (THM-438, H(Paley)~e·avg); codex hyper-Bessel + boxeph Laguerre-Pólya boundary = the Paley char_S=∏(x²+p) real-rooted spectrum = Re=−1/2 critical line = quasirandomness (S85). NC2 noncancel on the wall = confluent Paley Vandermonde/Wronskian≠0 = real-rootedness = reality of the Paley spectrum (tournament-spectral, THM-1555/213).
- Unification not proof; NC2 wall OPEN. Credits boxeph THM-2033/S202, codex, klein THM-1805. reflection the-nc2-wall-is-the-regular-paley-tournament-...-S89; script nc2_confluent_vandermonde_is_regular_S89 (+out).

## death-star-2026-07-21-S88 -- The CHANNEL-TOURNAMENT LENS: NC2 is a tournament-nullcone on its radial channels; the regular-channel core is the wall; explains why domination (MISTAKE-202) was refuted. NC2 OPEN. HYP-8772.
## codex-2026-07-21-support-dirichlet-atlas — THM-2000/2005: reciprocal sequences become support profiles

- **Archived paper and corrected object.** Read the complete Downey--Ong--Sellers preprint.  Its polygonal table is one digamma-clock axis of the master figurate array; odd polygons follow without a new argument, the square row is a digamma-to-trigamma confluence, and centered polygonals supply a second clock.  The literal harmonic-subset object is the support profile `D_A(z)`, with multiplicities isolated as a collision tax; one scalar mass is only a valuation.
- **Exact analytic extension.** THM-2000 now has stored, byte-identical normal/`-O` evidence and fixed hashes.  THM-2005 proves the profile/abscissa law, integer Abel--Dini lift and Bertrand recursion, Egyptian conservation hyperplane, ordinary/centered polygonal clocks, maximum-`c3` and Forcade profiles, fibbinary/Moser Mahler laws and block tails, primitive-residue/LRC measure, Sylvester remainder, and the alternate-vertex/tournament audit.
- **Tournament census closed to a rigorous interval.** Burnside's all-odd-cycle formula is reconstructed exactly through A000568(20).  Orbit growth bounds the remaining deduplicated reciprocal tail by `<3.11e-44`, placing the support mass between `1.853534132290116317333715265823800054971816577029176...` and `1.853534132290116317333715265823800054971816608078302...`.  This also corrects the old non-A000568 terminal value and cleanly separates the indexed collision tax `2`.
- **Concurrent tournament integration.** After rebasing, THM-2016's reducibility ceiling gave a new exact interpretation: with `tau_c(n)=M_3(n-1)/M_3(n)`, `product_(j=4)^N tau_c(j)=1/M_3(N)`.  Its normalized defects have denominator word `(6,5,8,7,...)`, exactly the adjacent-pair reversal of `{5,6,7,...}`, so their profile is `zeta(s)-H_4^(s)` with abscissa one.  Sorting the same hazards changes the prefix-product mass from `75/4-24 log 2` to `2`; the exact parity-shuffle tax is `67/4-24 log 2`.
- **Formalization and historical repair.** The sorry-free `SupportHarmonicFigurate.lean` module now includes ordinary-polygonal factorization/splitting, both parity branches of the maximum-`c3` denominator algebra, and both condensation ratio/defect identities; its isolated build passed all 8,475 jobs with no `sorryAx`.  THM-1985/1990, both legacy scripts/outputs, MISTAKE-209, the hypotheses/results indexes, and both reflections now retract the false `c3` maximum, indexed/support conflation, impossible `H=2` neighbor, linear-growth iff, scalar-fingerprint/transcendence claims, and conjectural H-spectrum divergence while preserving the exact simplex core.  MISTAKE-217 additionally repairs THM-2016's false largest-SCC shortcut: the ceiling follows from SCC additivity plus discrete convexity/superadditivity of `M_3`, not by discarding smaller SCC contributions; the referee checks 6,634 exact finite rows of that repaired core.  The formerly empty THM-2016 result is now regenerated and byte-identical under normal/optimized Python.
- **Honest frontier.** This closes a broad reciprocal-sequence packet, not LRC-14 itself.  The primitive-residue profile exactly recovers THM-819's LRC measure, but the live LRC-14 terminal-word/large-modulus mathematics remains separate.  Next profile targets are normalized tails for the other census rows, centered/automatic Lean recurrences, and growth laws for actual LRC witness/modulus supports.

## codex-2026-07-21-gmc2-degree-gap — THM-2017 closes the strict degree-gap three-weight region; HYP-8766 isolates resonance

**Owner:** find connections between tournaments and NC2; long session, push/pull often.

**THE BRIDGE (verified).** NC2 noncancellation = det[(a_i+k)!] (THM-1815) = ∏a_i!·Vandermonde(radial degrees) = Σ_T sgn(T) x^{score(T)} = the SIGNED TOURNAMENT SUM (THM-1805/my THM-1925; transitive tournaments survive). Verified exactly: det=∏a!·Vand over 7 degree sets; Vandermonde=tournament-sign-sum (n=3:30,n=4:1440). So the object deciding NC2 IS the tournament sign-sum over the channel degrees.

**TWO REGIMES.** DISTINCT degrees => Vandermonde≠0 => TRANSITIVE channel (death-star HYP-8772) => noncancellation (codex THM-2017 degree-gap). REPEATED degrees = the resonance WALL => ordinary Vandermonde VANISHES but the CONFLUENT Vandermonde (derivative/Wronskian row) SURVIVES nonzero (verified) = codex's 1/m hyper-Bessel correction / my Laguerre-Pólya boundary ODE θ²Φ=ξΦ (HYP-8775). Respects MISTAKE-212 (tied channels = CONFLUENCE not intransitivity; degree preorder ties, derivative order breaks it).

**codex's λ IS the node-spacing.** Channel factorial-degree D(k)=dm+λk (λ=e-rd), so the transitivity Vandermonde = λ^{C(nch,2)}·(int): |λ| = node spacing, r-|λ| = confluence order. |λ|>=r+1 separated (THM-2017), |λ|=r boundary (my L-P), 0<|λ|<r band (HYP-8766), λ=0 FULLY CONFLUENT = central resonance (codex HYP-8771) = my τ=1 regular core. codex's regime map = the confluence-order stratification of ONE discriminant.

**META (fermionic/bosonic).** The tournament sign-sum = Vandermonde = a DETERMINANT (fermionic, signed) -- nonvanishing is free (distinct nodes; the sign-involution collapses to the transitive core). NC2's E[P^m] = a PERMANENT (bosonic, all +, no sign) -- doesn't collapse, hence hard. THM-1815's miracle: the bosonic noncancellation reduces to a fermionic Vandermonde on the channel degrees -- NC2 borrows the fermion's rigidity. The wall (regular/tied core) = my invariant-resistant hot center (THM-2016): maximally symmetric, confluent, no cheap invariant separates it -- WHY domination failed (MISTAKE-202: regular channels have no source).

**Unifies FOUR threads as ONE object (confluence of the tournament sign-sum):** THM-1815 transitivity Vandermonde / THM-1805-1925 tournament sign-sum / death-star HYP-8772 channel lens / codex THM-2017 hyper-Bessel + my HYP-8775 L-P boundary. Plus my continuum (THM-1979/2013/2016) reads the same axis as cyclic temperature.

**Push/pull:** checkpointed twice; codex reserved THM-2023 (proving the L-P claim I raised) -- convergent. NEW connection for the fleet.
> **Post-push upgrade (THM-2023):** HYP-8775's general `Phi_(p0,q0)`
> Laguerre--Polya claim is proved by Gauss multiplication plus the
> Baricz--Singh positive-parameter `0F_Q` zero theorem. The `rd-e=r` boundary is
> NC2-clear off the negative real `xi` axis for all primitive charges. This does
> not cover the opposite `Psi_r` boundary or remove the negative-root exceptions.

**Watson machinery MAP (Explore sweep).** The cross-shell descent splits by degree gap λ=e-rd into 3 regimes (codex THM-2017): (a) |λ|>=r+1 degree-gap-dominant = PROVED (mixed-factorial + dominated convergence); (b) |λ|=r sharp boundary = hyper-Bessel Phi_{(p0,q0)}(ξ), clear iff nonzero; (c) 0<=|λ|<r interior = full-entropy saddle = OPEN (codex HYP-8766/8771 central resonance). Single-shell radial closed by my THM-1565 Radial Lemma (Watson-Nevanlinna) + klein THM-1665 (per-component Watson) + EMP; the symmetric-top charge-0 projection IS a modified Bessel I_0 (already canon THM-1835).

**Honest scope:** SYNTHESIS/identification (verified computations + THM-1805/1815), NOT a proof of NC2. The residual is ONE discriminant at distinct nodes (opus THM-1710 multinomial-ratio, shared with TNC) + its confluent limit at the tied core (hyper-Bessel/L-P). Backlog lead: the free-probability/semicircle bridge (death-star S88: channel weights=free cumulants, tied core=semicircle=regular tournament=my τ=1; Catalan/Wigner = THM-438 H(Paley)~e·avg).

**Next:** develop the free-prob/semicircle bridge (regular tournament spectrum <-> NC2 tied-core entropy saddle); connect the confluent Vandermonde to the Hermite/Wronskian radial closure (THM-1615). Artifacts: THM-2033, HYP-8780, reflection the-vandermonde-is-the-bridge-tournaments-and-nc2-boxeph-S203.md, script confluent_transitivity_vandermonde_boxeph_S203.py (+.out).

**Coordination:** this COMPLEMENTS codex THM-2017/2018 (I resolve the boundary zero-loci to a classical L-P question; codex owns the interior central resonance HYP-8771). Updated GMC2-FINISH-MAP with the regime map + the boundary addendum. @codex: the boundary is a Pólya-Schur multiplier-sequence problem, not an open zero-locus -- and it's unconditional for real-positive data.

**Honest scope:** did NOT prove GMC(2) or the general L-P claim; I mapped the Watson machinery, corrected my own refuted dominance claim, and reduced codex's boundary residual to a named classical (Laguerre-Pólya/Pólya-Schur) problem with strong numerical evidence + rigorous I_0 base. Artifacts: HYP-8775, HYP-8770 correction, reflection watson-estimates-for-gmc2-and-the-laguerre-polya-boundary-boxeph-S202.md, script hyperbessel_boundary_zeros_boxeph_S202.py (+.out), finish-map addendum.

**Next:** prove Phi_{(p0,q0)} in L-P via Pólya-Schur (or find the first complex zero = the inner resonance band boundary); connect the multiplier-sequence closure to codex's central resonance.

**Two-front repo map (2 Explore agents) => the precise state.** GMC(2) = [broad proved skeleton] + [ONE residual].
- PROVED: DvdEZ/NC2 => GMC(2) (Lean, no sorry, mathieuZhao_of_charge_pos); polar bridge E=L∘CT_u makes the ANGULAR layer the Duistermaat–van der Kallen theorem (THM-1630/1645) -- the gap is PURELY RADIAL (Laplace determinacy, ker L != 0). Sound closed strata include sign-coherent, two-charge all-degree, pure radial/EMP, span-2, constant-endpoint arbitrary-middle (THM-2014), strict three-weight degree gaps plus generic sharp boundary (THM-2017), single-straddle, the single-character pair base, bounded span<=4, and the constant span-6 certificate. The arbitrary-radial `{−1,0,1}` claim and pair-radical supply are open as corrected above.
- THE ONE RESIDUAL BRIDGE: radial-channel noncancellation for every two-sided P. Its named, overlapping forms are (1) symmetric-top Watson dominance (HYP-8770), (2) three-weight resonance asymptotics (HYP-8766), (3) asymmetric-top/bottom-up descent, and (4) multilevel pair-radical elimination after first-return cancellation (HYP-8765). No span-uniform finite bound is available because detection depth grows with radial degree.

**DELIVERABLE 1: the GMC(2) FINISH-MAP** (00-navigation/GMC2-FINISH-MAP-2026-07-21.md) -- assembles the skeleton and radial bridge. Corrected reading: pair-straddle radical membership is the desired combinatorial form, not a consequence already supplied by THM-1770; the cross-shell Laplace descent is the analytic form.

**DELIVERABLE 2: HYP-8770 -- I OWN the symmetric-top Watson dominance** (my THM-1565 Radial-Lemma territory). Precise statement: E[P^m]=Σ_V a_V V!, unique top a_max=C(m,m/2)(αβ)^{m/2}; a_max≠0 => E[P^m]≠0 (m≫0) => αβ=0. It is a genuine Borel/Watson determinacy -- the V-level factorial grading is only POLYNOMIAL while level coeffs are EXPONENTIAL, so the crude gap bound fails; closed only in the coefficient-dominant sub-case. Built an exact no-sympy moment engine (charge-balanced-tuple formula) that VALIDATES death-star's span-6 certificate and finds no two-sided nullcone member across 5 strata (constant+radial).

**Push/pull:** finish-map pushed (checkpoint); broadcasting the division of labour at close. Corrected split: symmetric-top Watson => HYP-8770; three-weight resonance => HYP-8766; multilevel pair radicals after first-return cancellation => HYP-8765; stratum certificate bank => death-star+codex.

**Honest scope:** I did NOT prove GMC(2). I assembled the finish-map, validated the moment machinery, and pinned + took ownership of the symmetric-top sub-residual with the coefficient-dominant partial. After the correction, Watson dominance is one major radial subproblem alongside resonance and multilevel cancellation, not a completed reduction of every other radial case. No canon overridden; HYP-8770.

**Next:** prove the symmetric-top Watson dominance for two shells (extend THM-1565), then general; the coefficient-dominant case is already clean. Artifacts: GMC2-FINISH-MAP-2026-07-21.md, HYP-8770, gmc2_symmetric_top_dominance_boxeph_S201.py (+.out).

## codex-2026-07-21-gmc2-proportional-central-resonance -- THM-2021 / HYP-8771

**Owner directive:** finish NC2 and mine seemingly unrelated repo work for ideas.

- **FULL CLOSURE INCOMING / INDEPENDENT HOSTILE AUDIT: PASS.** THM-2022 proves
  arbitrary-support NC2 by algebraic descent plus Frobenius amplification of
  the lowest balanced Wick face. I checked the exact-support localization,
  rational LP dual, no-projection-collision lemma, DvdK input, good-prime
  selection, common-factorial normalization, Kummer no-carry iff componentwise
  `p`-dilation, strict off-face quotient gap, multinomial Lucas congruence, residue-field
  Frobenius, neutral case, and NC2=>GMC(2) finish. No load-bearing gap found.
  Crucially the proof retains the whole minimum face as `Q^p`; it does not
  repeat MISTAKE-211's atomwise separation.
- **Unrelated-work challenge (MISTAKE-214).** The incoming Vandermonde identity
  is useful, but HYP-8785/8790 confused node values with score exponents when it
  declared the tied NC2 wall regular/Paley. The `m=2` tied core has two equal
  nodes and no regular two-vertex tournament. I retained the confluent
  Vandermonde and central-trinomial identities, downgraded Paley/free-probability
  iff claims to analogy, and linked the actual tied-face certificate `Q^p`.

> **Post-incoming correction.** The completed THM-2018 proves the full
> all-charge proportional hypersurface unconditionally. Its root-of-unity EGF
> argument only needs toral nonzeros arbitrarily far out, which must intersect
> EMP's cofinite radial nonzero tail. THM-2021's Legendre analysis is a sharper
> zero-geometry refinement, and HYP-8771 is no longer an NC2 blocker.

- **Incoming synthesis/correction (MISTAKE-213).** Concurrent THM-2018 first reserved the symmetric proportional identity `h=kappa*b^2`; its completed theorem proves the all-charge factorization and the sufficient property that the toral factor is nonzero arbitrarily far out. The initial claim that no-consecutive-zero is insufficient was backward: it already supplies unbounded toral nonzeros, which meet EMP's cofinite radial tail. Concurrent THM-2019 (affine-height supports) is transverse, while THM-2020's finite-place channel separation supplies a route to the stronger zero-profile question.
- **THM-2021 (proved self-contained part, corrected by MISTAKE-213).** For `P=Z^p a(s)+b(s)+Zbar^q c(s)`, if `h=s^(pq/g)a^(q/g)c^(p/g)=kappa*b^r`, then `E[P^m]=A_m^(p0,q0)(kappa)L(b^m)` exactly. At `p0=q0=1`, `A_m` is a Legendre transform, every zero parameter is negative real, and its three-term recurrence forbids consecutive zero levels. That already closes the full symmetric proportional slice: EMP makes the radial factor nonzero eventually, while the toral factor is nonzero arbitrarily far out. THM-2018 supplies the all-charge version.
- **Announced external refinement, honestly labeled.** A February 2026 IAS announcement for Mangoubi--Kadets--Weller Weiser states that a fixed nonzero point is a zero of only finitely many Legendre polynomials. No public proof/preprint was located. This sharpens the symmetric zero profile but is not required for NC2 after THM-2018.
- **HYP-8771.** The higher-charge finite-zero-recurrence statement is now a stronger toral-sequence target, not the proportional-slice NC2 gate. Exact gcd scans found zero shared-root events across 1,379 level pairs for `(1,1),(1,2),(1,3),(2,3)`. Tournament vertices are moment levels; gcd-degree flips the chronological edge. All four tournaments are transitive (zero flips/cycles, singleton SCCs, one Hamiltonian path). This is evidence, not proof.
- **Artifacts.** THM-2021, HYP-8771, reflection `nc2-central-resonance-is-legendre-finite-zero-recurrence-...`, exact script/output `gmc2_proportional_legendre_finite_recurrence_thm2021.py/.out`; HYP-8766 and the finish map updated. THM-2022 now closes full NC2 beyond these slices.
## codex-2026-07-21-NC2-followup -- THM-2022 proves full NC2 and GMC(2) by Frobenius amplification of the lowest balanced Wick face

- **RESULT:** For arbitrary finite support and arbitrary complex coefficients, the Gaussian moment nullcone consists exactly of the two strict one-sided charge loci. Hence GMC(2) follows by the already-proved charge reduction.
- **MECHANISM:** Descend a fixed-support torus null point to algebraic coefficients. The lowest balanced face has a nonzero DvdK constant term `Q` at some level `m0`, with common Wick height `A0`. At a good prime `p`, normalize `M_(p*m0)` by `(p*A0)!`. Kummer/Lucas remove every channel outside the `p`-dilated face; Frobenius makes the entire tied residue exactly `Q^p != 0`. No first-return atom, unique channel, coefficient sign, radial-degree gap, or asymptotic dominance is assumed.
- **UNRELATED-WORK TRANSFER:** This is the exact realization of the repo's “retain a side channel under Frobenius twist” and tropical carry-wall ideas. The preserved side channel is the whole face sum, so MISTAKE-211's cross-atom cancellation is absorbed rather than denied. Dimension two is essential because balance leaves one monotone Wick factorial.
- **AUDIT:** Three independent hostile reconstructions plus the parallel GMC2 task found no gap. Exact verifier covers a cancelled earlier face layer, two colliding minimum channels, 126 off-face channels, a rational intercept, and neutral `A0=0`; stored transcript matches. Commits `fab3b36b7`, `055b68a94`, `c50d86f5f`, `0f0f8e6d4`, `b7dba73a5` (with later integration commits possible).

## codex-2026-07-21-NC2-transfer -- exact-period Frobenius packets and the honest LRC boundary

- **PROVED (THM-2041):** In the cyclic group algebra of period `q` and
  characteristic `p` coprime to `q`, every `p`-stable integral twist satisfies
  `ct_q(Theta Lambda^(pm))=ct_q(Theta Lambda^m)^p`. Primitive-support and
  Ramanujan exact-period kernels are stable because multiplication by `p`
  permutes the units modulo `q`. This includes composite period `q=14` and
  propagates any already-nonzero exact-period packet through all `p`-power
  levels.
- **LRC SYNTHESIS (HYP-8800):** THM-2022's four services separate into
  descent, face exposure, off-face filtration, and whole-face Frobenius
  preservation. LRC finite checking and packet reductions resemble the first
  two, and THM-2041 supplies the fourth. The missing content is a nonzero
  integer base safe-count packet, an off-packet filtration, and a lift that
  retains safe inequalities plus endpoint/boundary ownership.
- **THREE NO-GOS:** finite phase-function algebras have zero divisors and no
  DvdK base theorem; LRC danger-count factorial moments terminate at order 13,
  so prime amplification exits the alphabet; and a nonzero signed Ramanujan
  trace need not imply any safe phase. The recommended attack amplifies
  labelled period packets, not intersection order.
- **OTHER FRONTS:** the same whole-initial-layer principle can preserve unit
  root orbits in Hensel/CRT work, toral character packets, resultant-unit
  layers, and finite cyclic counts. Each still needs a problem-specific base
  nonvanishing theorem and filtration. It does not automatically control
  Paley tournament scores or higher-dimensional products of Wick factorials.
- **CORRECTION (MISTAKE-215 / THM-2040):** the incoming global
  de-factorialization/Paley equivalence was false. A common factorial exists on
  the exact even symmetric monomial wall, where the quotient is a central
  trinomial polynomial, but THM-2022 uses only a prime-local minimum-valuation
  normalization. Scalar moments are not the special moment-matrix
  determinants.
- **COMPUTATION:** 71,484 Ramanujan dilation identities and 1,740 cyclic
  twisted-moment identities passed, with a failing non-invariant-twist control;
  the corrected de-factorialization audit passed 72 exact monomial-wall cases.
- **INDEPENDENT PACKET AUDIT:** the merged finite-abelian formulation passed
  1,479 period/order/prime triples and 4,437 tests each of exact-period
  projector commutation, nonzero-layer preservation, and Ramanujan-energy
  dilation. Five bad-prime examples verify the sharp nilpotent boundary.
- **SURPRISING EXACT TRANSFERS:** THM-2022 now covers every positive rational
  Gamma radial shape by the same lowest-face argument, while THM-346 gains the
  congruence `A_M^p=A_M mod p` for odd-prime cube-walk transport between
  arbitrary buckets. A three-factor Wick example shows why a scalar face does
  not automatically extend the proof to higher Gaussian dimension.
- **LRC TARGETS AND REPAIR:** the strongest next seed routes are THM-671's
  resolved-modulus `B5` supply and a familywise Fejer/Toeplitz theorem. Two
  deeper sidecars are now explicit: cubing is a single six-cycle on `U_14`,
  and characteristic 7 retains parity times Hasse-jet depth. Older claims that
  2 ramifies in `Q(zeta_14)` were corrected; the actual obstruction is
  nonunit/non-etale group-algebra behavior.
- **EXTERNAL CHECK:** current primary-source frontier was cross-checked against
  Sungkawichai--Trakulthongchai (LRC through 13 runners),
  Malikiosis--Santos--Schymura (finite checking), Giri--Kravitz and
  Jain--Kravitz (spectra), and Bedert (Riesz products).
- **PART 1 — CRUX further reduced (S83's (i) H(reg)≥Szele avg → quasirandomness).** Binding case = DOUBLY-REGULAR (Paley, max disc): smallest H/(n·disc) among regulars (Paley-7: 3.38 vs rot 25; Paley-11: 35.6 vs 8457). Paley is QUASIRANDOM (K-spectrum ±i√n, ratio →0) ⟹ H=(1+o(1))·avg (quasirandom Ham-path counting lemma); measured H/avg≈2.0-2.4 (bounded). n·disc/avg=n(n+1)^{(n-1)/2}/n!→0 super-exp (0.71 at n=7, 0.069 at n=11), so large n is loose (needs only a tiny fraction of avg); small n direct. Min-strong (Busch) route EXCLUDED — doubly-regular disc ~(√n/2)^n is too big — so the crux genuinely needs regular=quasirandom=near-average. The crux is now standard pseudorandomness, not an eigenvalue-product mystery.
- **PART 2 — GMC(2) positivity.** The S81 Pell identity sharpens: for REAL-coeff P, E[sym(P)²]−E[alt(P)²]=E[P·P̃]=E[|P|²]≥0 (Bargmann norm) ⟹ E[sym²]≥E[alt²] — a RIGOROUS proof of klein THM-1810's bosonic≥fermionic at the squared-moment level. HONEST: orthogonal to the nullcone (one-sided P=Z has E[sym(P^m)²]=m!/2>0 despite E[P^m]=0); a Bargmann-PD handle on the TORAL side (S67/S77 toolkit), the open RADIAL gap unaffected — does NOT close GMC(2).
- No open problem closed; a genuine crux reduction + a small rigorous positivity. Credits klein THM-1950/1810. reflection crux-reduced-to-quasirandomness-...-S84; script crux_reduction_and_gmc2_positivity_S84 (+out). LRC(≤13) not re-audited.
## opus-2026-07-20-S448 - Folded the sequence-sweep findings into THM-1985: the LRC deep-well = harmonic/triangular (THM-805 bridge)

Folded the S447 background sequence-sweep (28-sequence master table + the repo's existing harmonic/
gamma/zeta appearances) into THM-1985 s4. THE KEY SYNTHESIS (verified exact):
- THM-805: the LRC deep-well base measure = a HARMONIC number over a TRIANGULAR number:
  m({1..k}; 1/(k+2)) = H_k / C(k+2,2). At k=12: H_12/C(14,2) = 6617/194040 (exact, = mac-mini
  THM-736's |G'({1..12})|). So the divergent-edge object (harmonic H_k) and the convergent figurate
  size (triangular C(k+2,2), sum 1/·=2 of THM-1985 s1) combine into the single most important LRC
  extremal measure -- the reciprocal-sum program meets the LRC flagship exactly.
- The LONELINESS FLOOR spectrum IS the harmonic series: M_floor(n)=1/n => sum = sum 1/n (diverges).
- THM-1926 tournament zeta = 1/det(I-uA) = Euler product over primitive cycles (cycles = primes).
- Staircase 3-cycles k(k-1) telescope to 1. Deep-well denominators Phi_6(n)=n^2-n+1 (cyclotomic),
  sum 1/Phi6=1.798 converges.
- REFINEMENT (agent): only pure binomials C(n,k) give rational k/(k-1); the 1+2^d / 2^{n-2}+1 /
  Cayley-Dickson 2^k+1 families converge to IRRATIONAL Erdos-Borwein-type constants. A002088, A001764
  absent from repo.

THE LOOP CLOSES: the figurate 2 (triangular=tournament size), the harmonic H_k (resistance/divergent
edge), and the zeta_T Euler product are three faces of one structure; the LRC deep-well = H_k/(triangular).

Files: THM-1985 s4 + depends_on (added THM-805/1926); harmonic_triangular_bridge_opus_S448.py(out).
Cites THM-805/1926/1370/1920/1930/1970/1975, kps THM-1980.

## opus-2026-07-20-S447 - Reciprocal sums = the harmonic-scale face of the poly/#P tower (THM-1985)

Owner: the reciprocal of an integer sequence is a subset of the harmonic numbers; study sum 1/a_n
for as many repo sequences as possible; figurate reciprocals (triangular=2), Abel-Dini, Bertrand;
1+1/2+..+1/5 > 2 already while sum 1/triangular = 2.

THM-1985: a sequence's GROWTH (its poly-tower position) = its reciprocal sum's CONVERGENCE. THREE
STRATA. (1) FIGURATE invariant-SIZES = char_S coefficients (THM-1920): sum_{n>=k} 1/C(n,k)=k/(k-1)
(exact k=2..6). Tournament sizes -> RATIONALS: arc=C(n,2)=triangular => sum 1/arc = 2 (the Downey-
Ong-Sellers triangular identity realized on the tournament -- the char_S subleading series sums to
exactly 2); c3max=C(n,3)->3/2; var-max=2C(n,3)->3/4. Degree-k invariant -> reciprocal-sum k/(k-1).
(2) COUNTING seqs (super-exp) -> fast transcendentals: sum1/A000568=2.8535, sum1/A038375=2.6293,
sum1/A051337=2.198, sum1/A002854=1.062; Cayley-Dickson sum1/(2^k+1)=0.7645 (Erdos-Borwein cousin),
H=1+2^(n-2) SC-neighbor -> 1.2645. (3) H-VALUE spectrum (odds minus {7,21}, THM-1370) ~linear =>
sum 1/H-value DIVERGES (harmonic-slow) = H's VALUE SET sits at the convergence/divergence boundary =
the reciprocal-sum face of THM-1970's formula/#P edge. ABEL-DINI: no series at the exact boundary =
kps THM-1980's 'Redei parity is the last formula'. BERTRAND boundary = sum 1/(n ln n).

THE PICTURE: rational k/(k-1) [figurate invariant sizes, poly, deep convergence] | transcendental
[counting sequences, the census] | DIVERGES [H-value spectrum, #P, the edge]. The reciprocal sum is
the harmonic-scale coordinate that recovers the poly/#P tower.

Concurrent: a background agent is sweeping the full repo sequence list (30+ OEIS A-numbers, growth
rates, existing harmonic/gamma appearances like THM-805 resistance=harmonic-number) -- findings fold
into THM-1985 next session.

OPEN: identify the counting constants (2.85, 2.63 -- e/pi/new?); the H-value density c in c*ln x;
Bertrand-scale repo sequences.

Files: THM-1985; HYP-8745; reciprocal_sums_of_repo_sequences_opus_S447.py (+out). Namespace clean
(1985/8745). Cites THM-1920/1930/1970/1975/1370, kps THM-1980, Downey-Ong-Sellers.

## kind-pasteur-2026-07-21-S128c144 - The figurate reciprocal ladder and the harmonic edge: our sequences as sub-series of the harmonic series (THM-1990)

Owner: reciprocal of an integer sequence = a subset of the harmonic numbers; 1+1/2+..+1/5 > 2 while
sum 1/T_n = 2; study our sequences' reciprocal sums extensively and extend. Looked back at every
integer sequence in the corpus through the reciprocal-sum lens.

**PROVEN CORE (THM-1990) -- the figurate ladder.** Telescoping (exact rational, verified):
sum_{n=k}^N 1/C(n,k) = k/(k-1)*(1 - 1/C(N,k-1)) => sum_{n>=k} 1/C(n,k) = k/(k-1) for k>=2, from
1/C(n,k) = k/(k-1)[1/C(n-1,k-1) - 1/C(n,k-1)]. THE LADDER by dimension:
  k=1 vertices n = HARMONIC = DIVERGES (the edge); k=2 arcs/triangular = 2 (the owner's value);
  k=3 tetrahedral = 3/2; k=4 = 4/3; ... -> 1.
Arc count of K_n = C(n,2), so the corpus' central staircase sequence has reciprocal signature EXACTLY
2, and its telescoping 1/C(n,2)=2(1/(n-1)-1/n) IS the arc->vertex (dim2->dim1) reduction.

**THE HARMONIC EDGE = p=1 = dimension 1, UNIFYING THREE LENSES.** The ladder diverges only at k=1
(the ground set). This is the SAME p=1 boundary as THM-1980 (H's formula content collapses to one
bit = 'H at p=1') and THM-1870 (cycle counts turn #P at the Hamiltonian length). Three independent
lenses -- reciprocal convergence, 2-adic depth, cycle length -- place the marginal case at the same
p=1 corner. The project's dimensional ladder n->C(n,2)->C(n,3) IS the figurate reciprocal ladder
crossing that edge.

**RECIPROCAL SIGNATURE table (VERIFIED to named constants):** arcs 2, tetrahedral 3/2,
var(lambda^2)=2C(n,3) -> 3/4, squares pi^2/6, factorial e, central binomial 4/3+2pi sqrt3/27 (EXACT
18 digits), Catalan 2+4pi sqrt3/27, Fibonacci=reciprocal-Fibonacci const, Mersenne 2^n-1=Erdos-Borwein,
2^n=1. 2-ADIC/THETA: labeled tournaments sum 1/2^{C(n,2)}=1.6416325607... (partial theta q=1/2);
switching classes 2^{C(n-1,2)} = 1 + that (the +1 = extra n=1 term, PROVEN exact); SIGNED pentagonal
= Euler prod(1-2^-n)=(1/2;1/2)_inf (pentagonal number thm, THM-488 hub). Census fingerprints:
A000568=3.8535, A002854=3.0618, score=3.9325, tangent=1.5663, secant=2.2171.

**CONVERGENCE DICHOTOMY:** sum 1/a diverges iff a grows at most linearly; EVERY combinatorial repo
sequence (degree>=2) converges, only the linear ground-set ones (vertices, odds, H-spectrum) sit on
the harmonic edge -- a clean 'ground set vs structure' separation.

**EXTEND:** sigma(a)=sum 1/a_n as a sequence fingerprint invariant; telescoping = the a-monoid
(THM-1885, a:n->n+1) Mode-A action on figurate dimensions; a transcendence gradient (poly->rational,
exp->irrational analytic, lacunary->theta). NEXT: identify the theta constant 1.6416...; signed
reciprocal sums Sum(-1)^n/a_n for tournaments/even-graphs; Sum 1/H(T) as bridge to THM-1980.

Reframing + a proven ladder + a verified signature table; classical constants (figurate, Basel,
Erdos-Borwein, pentagonal thm) unified as one harmonic-subset classification of the corpus.
**Files:** THM-1990; reflection the-harmonic-boundary-p-equals-1-recurs-kps-S128c144; HYP-8750;
script reciprocal_sums_of_our_sequences_kps_S128c144 (+out). Namespace: THM-1990/HYP-8750 (clean).
## death-star-2026-07-21-S83 -- H≥disc: the REGULAR SUB-BASE reduced to ONE average (H(reg)≥Szele avg), proved n≥7 modulo that crux; the Pfaffian injection IS the even/odd (Ω/E_n) duality. HYP-8698.

**Owner directive:** work the Pfaffian injection and the regular sub-base (S82 handles toward klein THM-1950's open base H(C)≥max(1,s)disc(C) for HYP-8636).

- **REGULAR SUB-BASE H(reg)≥n·disc(reg): PROVED for n≥7 modulo one crux.** Chain H(reg) ≥(i) n!/2^{n-1} ≥(ii) n(n+1)^{(n-1)/2}/2^{n-1} ≥(iii) n·disc(reg).
  - **(iii) PROVED (AM-GM):** disc(reg)=∏(1+μ_j²)/2^{n-1}, Σμ_j²=C(n,2) fixed ⇒ ∏ maximized at equal μ_j²=n (doubly-regular=Paley) ⇒ disc(reg)≤(n+1)^{(n-1)/2}/2^{n-1}. Tight Paley-3,7,11.
  - **(ii) PROVED (elementary):** (n-1)!≥(n+1)^{(n-1)/2} — fails n=3,5, holds n=7 (720≥512), ratio increasing ⇒ all n≥7.
  - **(i) THE CRUX (conjecture, strongly evidenced):** every regular tournament has ≥ the Szele average n!/2^{n-1} Hamiltonian paths. Exhaustive n=3 (H=3≥1.5), n=5 (all 24 regular H=15≥7.5); samples n=7 (min 171≥79), n=9 (min 3243≥1418), huge margins. n=3,5 direct.
  - So the "regular is the wall" crux (S75/S76) is now a SINGLE tractable Ham-path statement (plausibly Moon/Alon/Busch), far easier than the eigenvalue-product original. Doubly-regular (Paley) tightest.
- **PFAFFIAN INJECTION:** aggregate 2^{n-1}H≥Σ_{S even}Pf(K[S])²=det(I+K) confirmed with room (slack 112,416 at n=5,6). STRUCTURAL READING: disc=Σ Pf² counts EVEN cycle-covers (cycle-space); H=I(Ω,2) counts via ODD cycles (OCF); so H≥disc = "the ODD (OCF) count dominates the EVEN (Pfaffian) count" — the even/odd, cut/cycle, E_n/Ω duality (S80). Per-subset injection open (Pf(K[S])²≤H(T[S])H(T\S) not clean; the right compatibility is subtler).
- Reduces HYP-8636's open crux to (i); no new theorem. Credits klein THM-1950. reflection the-regular-sub-base-...-S83; script h_ge_disc_regular_subbase_S83 (+out). GMC(2)/LRC(14) untouched.

## opus-2026-07-20-S446 - The path-cover polynomial is the refined compositional invariant; the formula/#P edge is REAL (THM-1975)

## codex-2026-07-21 -- THM-2000 support-harmonic transform and full figurate mass surface

Interpreted an integer sequence literally as a harmonic **support**, separating
it from indexed multiplicity by the collision tax.  Proved the Abel--Stieltjes
counting identity, logarithmic-block and partial-sum Dini laws, full Bertrand
boundary, regular-variation tail, heat/Mellin transform dictionary, and the
random needle-width threshold.  The master figurate array now has both a beta
integral and a finite partial-fraction form at every point: off resonance it is
rational-argument digamma/cyclotomic-log; exactly when `s-2>=2` divides `d`, one
double pole gives `Q+Q*pi^2` (hence transcendental).  Added exact polygonal,
Faulhaber p=1..5, balance, quarter-square, Gauss-theta, run-filtration, and
ratio/sum/product laws; THM-1127's affine toothpick phases force divergent
support via an arithmetic progression.  Corrected THM-1985/1990, the raw
tournament Dirichlet/EGF reversal, finite-prefix/offset errors, the open global
H-spectrum, and the false all-n self-line extrapolation.  Referee normal/`-O`
outputs are byte-identical and `RESULT=PASS`; a new sorry-free Lean module
certifies the algebraic core.  Files: THM-2000, MISTAKE-209/210/211,
`support_harmonic_abel_dini_figurate_surface_thm2000.py/.out`,
`SupportHarmonicFigurate.lean`, and reflection
`the-sequence-is-its-logarithmic-occupancy-not-its-index-codex-20260721.md`.

## boxeph-2026-07-21-S198 -- THM-1979 TOURNAMENT SPACE IS A SPECTRUM (single point -> continuum)

**Owner:** understand tournament space on n vertices as a spectrum from a single point (transitive) to a continuum; the maximally-spread score classes house the different structure.

**THE FRAME (verified n<=7).** Tournament space fibers over the SCORE SEQUENCE (Landau polytope; counts 2,4,9,22,59 for n=3..7). The spectral coordinate is score spread sigma^2=Var(scores) in [0,(n^2-1)/12]. Cyclicity is its EXACT affine image: c3 = n(n^2-1)/24 - (n/2)sigma^2 (= the classical c3=C(n,3)-sum C(s_i,2) restated). So score-spread and cyclicity are ONE axis, opposite directions.

- SINGLE POINT = TRANSITIVE (sigma^2 max, scores 0..n-1): fiber=1, c3=0, char_A=x^n (nullcone vertex), zeta=1 (my THM-1926, no periodic structure), reducible. The ordered structureless pole.
- CONTINUUM = REGULAR/near-regular (sigma^2->0): fibers SWELL -- max fiber 1,1,3,12,47 (n=3..7), -> inf with n; at low sigma^2 every class is STRONG and mostly MODULAR-PRIME; c3 max = n(n^2-1)/24. This is where the different structure lives (circulant/Paley thread THM-1955, |R|-independence mac-mini THM-1966 first at n=7, the whole irreducible interior THM-1978).
- MONOTONE LAW: fiber size, strong-frac, modprime-frac all run OPPOSITE to sigma^2. High-spread fibers are singleton reducible chains (strong=0); the regular center is all-strong. The structurally-richest score class is NEAR (not exactly at) the center: n=7 peak fiber 47 at sigma^2=4/7, score (2,2,3,3,3,4,4), all 47 strong, 29 modular-prime.
- LIMIT: the tournamenton spectrum -- transitive W=1_{x>y} (single ordered point) to quasirandom W=1/2 (positive-entropy continuum); degree function d(x): x -> 1/2; sigma^2 = integral (d-1/2)^2.

**UNIFIES the recent arc:** reduction principles = statements about the TOP (transitive/high-sigma^2, reducible); hardness = the BOTTOM (regular/low-sigma^2, irreducible). n=7 wall = first crack of the point opening into the continuum (THM-1830 atom). The clean theorems live at the point; the mathematics lives in the continuum.

**Housekeeping:** opus-S445 first-pushed THM-1970 (and THM-1975) so I renumbered my S197 n>=7-regime theorem THM-1970 -> THM-1978 (off the round-number grid to stop the collision churn). This session: THM-1979, HYP-8732. Integrated mac-mini THM-1966 (|R| new invariant from n=7). No overrides.

**Honest scope:** the fibration, the affine c3-sigma^2 law, the fiber/strong/modprime monotonicity, and the max-fiber sequence 1,1,3,12,47 are VERIFIED exact n<=7. The affine law is a restatement of a classical identity; the contribution is the SPECTRUM framing (single point <-> continuum via score spread) that unifies reducibility, cyclicity, structural diversity, and the n=7 breaks on one axis.

**Next:** (1) the max-fiber growth 1,1,3,12,47 -- asymptotics (does it track the # of tournaments with the balanced score seq, ~ 2^{C(n,2)}/(Landau count)?); (2) the exact location of the structural-diversity peak vs the regular center (off-center at n=7); (3) the quasirandom-neighborhood entropy as the continuum's "dimension". Artifacts: THM-1979, HYP-8732, reflection tournament-space-as-a-spectrum-single-point-to-continuum-boxeph-S198.md, script tournament_space_spectrum_boxeph_S198.py (+.out).
