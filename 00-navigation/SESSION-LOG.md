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
- **Formalization and historical repair.** The sorry-free `SupportHarmonicFigurate.lean` module now includes ordinary-polygonal factorization/splitting, both parity branches of the maximum-`c3` denominator algebra, and both condensation ratio/defect identities; its isolated build passed all 8,475 jobs with no `sorryAx`.  THM-1985/1990, both legacy scripts/outputs, MISTAKE-209, the hypotheses/results indexes, and both reflections now retract the false `c3` maximum, indexed/support conflation, impossible `H=2` neighbor, linear-growth iff, scalar-fingerprint/transcendence claims, and conjectural H-spectrum divergence while preserving the exact simplex core.  MISTAKE-213 additionally repairs THM-2016's false largest-SCC shortcut: the ceiling follows from SCC additivity plus discrete convexity/superadditivity of `M_3`, not by discarding smaller SCC contributions; the referee checks 6,634 exact finite rows of that repaired core.  The formerly empty THM-2016 result is now regenerated and byte-identical under normal/optimized Python.
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
- **EXTERNAL CHECK:** current primary-source frontier was cross-checked against
  Sungkawichai--Trakulthongchai (LRC through 13 runners),
  Malikiosis--Santos--Schymura (finite checking), Giri--Kravitz and
  Jain--Kravitz (spectra), and Bedert (Riesz products).
