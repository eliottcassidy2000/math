> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## death-star-2026-07-21-S93 -- Mathlib-PR packaging of the three-term no-common-root: Polynomial-R recast + minimal imports + the Mathlib-MISSING three-term Hermite recurrence proved, with "consecutive Hermite share no root" as the flagship application. All kernel-pure. HYP-8805.

**Owner directive:** work the natural next continuations of the S92 formalization toward finalized completion; pull/push often.

- **Absorbed the fleet correction:** NC2/GMC(2) are PROVED on paper by THM-2022 (codex; Frobenius amplification of the lowest balanced Wick face, certificate `Q^p`, DvdK constant-term THM-1630 as citation input). My S89-S91 Vandermonde/Paley synthesis is superseded (MISTAKE-214/215). **Honest vindication:** THM-2022 §4 IS the owner's S91 "divide by (pA0)!" directive -- that instinct was the crux; my error was the residue (it is `Q^p` via Kummer/Lucas, not a Vandermonde).
- **`ThreeTermRecurrence.lean` upgraded to PR-quality:** (1) **Polynomial-R recast** -- added `ThreeTerm.poly : ℕ → R[X]` (`noncomputable`), the `eval_poly` bridge, and `no_common_root_poly` in `Polynomial.IsRoot` form; (2) **minimal imports** (`import Mathlib` → 4 modules; build **8475→1202 jobs**); (3) **`hb` weakened to `∀ n, b (n+1) ≠ 0`** (only b at ≥1 occurs; Hermite has b0=0) -- strict generalization.
- **NEW `HermiteThreeTerm.lean` (the compelling application):** proved **`derivative_hermite_succ`** (`(hermite (n+1))' = (n+1)•hermite n`) and **`hermite_recurrence`** (three-term `hermite(n+2)=X·hermite(n+1)-(n+1)·hermite n`) -- **both Mathlib-MISSING** (Mathlib ships only the derivative-form `hermite_succ`); then **`hermite_no_common_root`** (no two consecutive Hermite polys share a root) by exhibiting `Polynomial.hermite` as the `ThreeTerm ℤ` instance `hermiteZ`. `module` tactic closes the ladder identity (H_n' terms cancel).
- **All 8 results verified kernel-pure** (`[propext, Classical.choice, Quot.sound]`, no sorry/native_decide), warning-free. Two PR-ready files, wired to root. Lean pitfall saved (Module-ℤ smul ≠ SMulZeroClass smul: use `map_smul`/`show…from`, not bare `rw`; `module` needs `mul_smul_comm` first).
- **Remaining (non-blocking):** `mathieuZhao_of_charge_pos` minimal-import trim; formalizing THM-2022 itself (DvdK not in Mathlib) is the real multi-session target. reflection gmc2-nc2-mathlib-submission-readiness-...-S92 (updated). NC2/GMC(2) FULL = THM-2022 (proved on paper), not this Lean payload.

## boxeph-2026-07-21-S205 -- JC<->LRC = one n=12 AP-rigidity; comprehensive view; Keller counterexample verified; red-team suite (HYP-8810)
## boxeph-2026-07-21-S207 -- cake, bagel, Moser and Fibonacci are ONE Pascal triangle (HYP-8820)

**Owner:** relate the repo's polygonal/polyhedral (figurate) work to Fibonacci and the cake & bagel cutting sequences.

**MINED:** the repo's figurate framework = opus-S317 (Vandermonde-truncation law: polygonal vs polyhedral; polygonal = first two Vandermonde layers of Pascal; polyhedral row-sum=2^n, shallow-diagonal=Fibonacci; polygonal row-sum=A000127 Moser circle), klein-S313 (the (r,g) shadow lattice, g-bonacci kernels 1−x−x^{g+1}, missing-region DEFICIT-1), mac-mini-S137 (the Hurwitz golden-corner principle: JC₂'s golden Fibonacci-degree corner + LRC's penultimate-convergent extremality + the g-bonacci kernels = one shape).

**SYNTHESIS (verified exact):** everything is ONE Pascal triangle read three ways.
- Full row sums = 2^n. Shallow-DIAGONAL (skip) sums = FIBONACCI. Truncated-row sums = the figurate CUTTING sequences.
- lazy caterer A000124 = C(n,0..2) (2D disk); CAKE A000125 = C(n,0..3) (3D ball); Moser A000127 = C(n,0)+C(n,2)+C(n,4) (polygonal row-sum); BAGEL (solid torus) = C(n,3)+n(n+1) = 1,2,6,13,24,40,62 (3 cuts->13).
- THE SURPRISE: bagel − cake = T_n − 1 (triangular minus one) = the DEFICIT-1 = klein-S313's g-bonacci-kernel missing-region boundary effect. The torus's topological hole IS the g-bonacci kernel's off-by-one. Genuine bridge between the cutting geometry and the Fibonacci-kernel side.
- g-bonacci kernels 1/(1−x−x^{g+1}) (verified): g=1=Fibonacci exactly; g=2,3 = the shadow-lattice family. The generating-function bridge between the row (cutting) and diagonal (Fibonacci) readings.

**So cake/bagel/Moser (rows) and Fibonacci (diagonals) are two projections of one Pascal/figurate triangle** — the same golden/figurate scaffold on which JC₂ (golden-degree corner) and LRC(14) (anti-golden Eisenstein extremal, the penultimate-convergent it forbids) sit (mac-mini-S137).

**Honest:** synthesis + verified figurate identities (cake/bagel as Pascal truncations, bagel−cake=T_n−1=deficit-1, Fibonacci skip-sum, g-bonacci kernels), tying opus-S317 + klein-S313 + mac-mini-S137 + my S206 LRC-Fibonacci-foil into one picture. (My polygonal-skip sub-computation had an indexing bug; cite opus-S317's verified version.) Artifacts: reflection cake-bagel-and-fibonacci-are-one-pascal-triangle-boxeph-S207.md, HYP-8820, script cake_bagel_figurate_fibonacci_boxeph_S207.py (+.out).
## boxeph-2026-07-21-S206 -- what an LRC(14) disproof must be, and why Fibonacci was proposed as a foil (HYP-8815)

> **Post-session audit:** MISTAKE-221 retracts the advertised near-AP,
> anti-golden, and full-autocorrelation “characterization.” The denominator
> scan gives exact rational lower bounds, not exact maxima in general.

**Owner:** mine connections to twelve and Fibonacci, consider disproof
constructions for LRC(14), and refine their necessary structure.

**Rigorous survivor.** A counterexample can be divided by its gcd, must be
Cover14 with `M<1/14`, and cannot have a dilated-AP maximum-deletion core by
THM-1017. THM-730 therefore gives that twelve-set a strict additive-triple
deficit. These are necessary conditions; the resummation from triple deficit to
loneliness remains open.

**Exploratory frame.** The deep well `{1,...,12,182}` has the Eisenstein
identity `183=Phi_6(14)` and maximizing time `14/183=[0;13,14]`. Fibonacci/
golden continued fractions were tested as an opposite pole. This motivates an
anti-golden hostile-control program, but does not prove that a counterexample
is near-AP, anti-golden, or controlled by the maximum runner alone.

**Computation.** The finite candidate scan found a rational witness above
`1/14` for every displayed AP, near-AP, Fibonacci-flavored, and sampled covering
packet. It recovers the canonical deep-well value and the imprimitive
`2*AP` witness `7/92`. After MISTAKE-221 the script reports these as
denominator-truncated lower bounds; values above threshold exclude only the
listed packets.

**Artifacts:** HYP-8815; MISTAKE-221;
`04-computation/lrc14_disproof_search_boxeph_S206.py` and matching output;
`07-reflections/what-an-lrc14-disproof-must-be-and-why-fibonacci-is-the-foil-boxeph-S206.md`.

## boxeph-2026-07-21-S205 -- JC/LRC AP-rigidity analogy and counterexample audit (HYP-8810)

> **Post-session scope:** HYP-8810 is a wildcard analogy, not a proved common
> reduction. THM-1017 is one-way on the LRC side. MISTAKE-205 keeps public
> provenance for the verified Keller map unsettled.

**Owner:** relate the Jacobian Conjecture to LRC, challenge assumptions, and
survey open problems while pulling concurrent work.

**Exact verification.** The owner-supplied three-dimensional Keller map has
constant Jacobian determinant `-2` and an exact collision, so it is a verified
counterexample object in repo canon (THM-1300/1315). Discovery and publication
credit are recorded separately from this calculation.

**Connection proposed.** AP/transitive/one-sided extremals appear as cold
vertices in several quotient pictures. The session proposed comparing the open
planar-JC reduction with the LRC(14) AP-extraction supplier. The shared language
is productive only after an exact map, preserved predicate, and loss ledger are
given; current canon does not identify the two residual problems.

**Artifacts:** HYP-8810 and
`07-reflections/jacobian-and-lonely-runner-two-nullcones-that-diverge-boxeph-S205.md`.
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
certifies the algebraic core.  Files: THM-2000, MISTAKE-209/210/213,
`support_harmonic_abel_dini_figurate_surface_thm2000.py/.out`,
`SupportHarmonicFigurate.lean`, and reflection
`the-sequence-is-its-logarithmic-occupancy-not-its-index-codex-20260721.md`.

Continuation: promoted THM-2005's full support-Dirichlet profile, integer
Abel--Dini representatives, fibbinary/Moser automatic dimensions and tail
bounds, collision-safe Egyptian continuum, Sylvester remainder, primitive
residue fibers, and Forcade divisor profiles.  The exact new tournament
criterion is
`sigma_-1(Div(C(2^p,2)))>=2`, with equality exactly for a Mersenne prime.
Kakeya achievement sets expose topology hidden by equal mass: `k`-simplex
reciprocals form exactly `2^(k-2)` interval components, while powers of `k`
give a Cantor set for `k>=3`; both triangular and binary atoms fill `[0,2]`.
Corrected max-`c3` to its parity splice and mass `75/4-24log2`, then used its
discrete convex increments to repair THM-2016's invalid dropped-SCC-summands
proof (MISTAKE-216).  A score-sidecar audit also narrowed the local-profile
claim: profiles add beyond score at `c3=12`, whereas score already explains
the old `c3=11` gain.

## boxeph-2026-07-21-S198 -- THM-1979 TOURNAMENT SPACE IS A SPECTRUM (single point -> continuum)

**Owner:** understand tournament space on n vertices as a spectrum from a single point (transitive) to a continuum; the maximally-spread score classes house the different structure.

**THE FRAME (corrected by MISTAKE-218; verified n<=7).** Tournament space fibers over the SCORE SEQUENCE (Landau counts 2,4,9,22,59 for n=3..7). Score variance lies in `[epsilon_n,(n^2-1)/12]`, with `epsilon_n=0` odd and `1/4` even. Cyclicity is its exact affine image, and `tau=(sigma_tr^2-sigma^2)/(sigma_tr^2-epsilon_n)`. This is a shell coordinate, not a complete structural axis.

- SINGLE POINT = TRANSITIVE (sigma^2 max, scores 0..n-1): fiber=1, c3=0, char_A=x^n (nullcone vertex), zeta=1 (my THM-1926, no periodic structure), reducible. The ordered structureless pole.
- BALANCED EDGE = REGULAR (odd) / near-regular (even), at variance `epsilon_n`; c3 max is `n(n^2-1)/24` odd and `n(n^2-4)/24` even, and every maximum-cyclic class is strong.
- NON-MONOTONE FIBERS: max fiber sizes are 1,1,3,12,47, but at n=7 size 47 occurs at variances 4/7 and 10/7 while the regular fiber has size 3. At n=6, variance 5/4 carries both size-3 all-strong and size-1 reducible fibers. The old monotone-richness claim is false.
- LIMIT: `W -> d_W` is a projection. Variance is `integral (d_W-1/2)^2`; many structured regular tournamentons share `d_W=1/2` with the quasirandom point.

**UNIFIES the recent arc:** reduction principles = statements about the TOP (transitive/high-sigma^2, reducible); hardness = the BOTTOM (regular/low-sigma^2, irreducible). n=7 wall = first crack of the point opening into the continuum (THM-1830 atom). The clean theorems live at the point; the mathematics lives in the continuum.

**Housekeeping:** opus-S445 first-pushed THM-1970 (and THM-1975) so I renumbered my S197 n>=7-regime theorem THM-1970 -> THM-1978 (off the round-number grid to stop the collision churn). This session: THM-1979, HYP-8732. Integrated mac-mini THM-1966 (|R| new invariant from n=7). No overrides.

**Honest scope:** the fibration, affine c3 law, parity seam, counterexamples to structural monotonicity, and max-fiber sequence 1,1,3,12,47 are exact through n=7. The durable contribution is the projection/fiber framing: score variance locates cyclicity, while the hard structure is exactly what that projection forgets.

**Next:** (1) max-fiber asymptotics without assuming the balanced score sequence is extremal; (2) the shell/fiber entropy profile; (3) invariants of the regular-degree tournamenton fiber beyond the degree map. Artifacts: THM-1979, HYP-8732, reflection tournament-space-as-a-spectrum-single-point-to-continuum-boxeph-S198.md, script tournament_space_spectrum_boxeph_S198.py (+.out).

## codex-2026-07-21-DC2-JC2 -- rank-two Poisson reconstruction, quantum anomaly, and planar no-mate theorem

- **THM-2044 / PC(2) FALSE:** independently reconstructed the four-polynomial
  rank-two Poisson counterexample from the supplied `R=x(2-3xq)`. It is a
  one-variable symplectic suspension of THM-1300 after the shear `q=y+xz/3`.
  All six bracket identities pass exactly; expanded term counts are
  `(2,35,246,78)`; three explicit rational points form one exact fibre.
- **THM-2045 / PLANAR POSITIVE RESULT:** for every `ab!=0`,
  `R=x(a-bxq)` has no polynomial planar Jacobian mate. Laurent `x`-sectors
  isolate the only constant-producing coefficient, and its ODE forces a
  constant not divisible by the required `s=xq` factor. The suspension cannot
  be de-stabilized to a planar counterexample through its first coordinate.
- **HYP-8802 / DC(2) GATE:** naive Weyl-symmetric quantization has a nonzero
  cubic Moyal anomaly in five of six relations, with term counts
  `42,42,0,3,165,273`. The `(D,R)` anomaly is central over `R`; a finite two-step
  correction gives an exact 332-term symbol with `M(Dq,R)=1`. Simultaneous
  terminating correction is open and is the precise direct-`A_2` target.
- **CORRECTION:** rank-two Poisson means four commutative generators; it does
  not automatically mean a Weyl `A_2` endomorphism or a planar Keller map.
  Cotangent/Hamiltonian-dual quantization safely lands in conventional `A_4`.
- **EXTERNAL CHECK:** BKK stable equivalence, Adjamagbo--van den Essen's
  equivalence framework, Lee--Li's 2024 planar Newton-polygon program, and the
  2024 Joseph-normal-form Dixmier sieve were checked from primary sources. No
  public copy of the supplied exact abstract/formula was located.
