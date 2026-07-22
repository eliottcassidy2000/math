## death-star-2026-07-21-S96 -- NC2 formalization: §4 channel-survival COMPLETE + §1 balanced-channel form + contrapositive entry; the entire self-contained arithmetic engine of THM-2022 (§1/§4/§5) is now kernel-pure (16 theorems). HYP-8805.

**Owner directive:** work the nc2 formalization, pull from other agents often.

- **§4 COMPLETE (all 3 residue-layer cases of THM-2022, kernel-pure):** (1) `multinomial_dvd_of_exists_not_dvd` — non-p-dilated channel ⟹ p∣multinomial (no-carry), via NEW `dvd_choose_of_dvd` (p∣n,¬p∣k ⟹ p∣C(n,k), from the absorption identity k·C(n,k)=n·C(n-1,k-1)) + the `multinomial_insert` recursion — MUCH cleaner than the Legendre digit-sum route I'd feared; (2) `factorial_dilate_dvd` — off-face dilated channel ⟹ p·(p·A0)! ∣ (p·A')! (factorial ratio killed by p), via `ascFactorial`; (3) `multinomial_dilate_modEq` (on-face survive) was already done. So the surviving residue layer = exactly the on-face dilated channels.
- **§1 completed:** `wick_expansion_balanced` — E(P^m) = the balanced-channel sum M_m (only charge-0 radial channels contribute; drops wt=0 terms). Plus `not_chargeOneSided_iff` — ¬one-sided ⟺ 0∈conv(charges), the contrapositive entry point.
- **STATE:** §1 (Wick + balanced + charge_radial), §4 (all 3 cases), §5 (Frobenius + face_sum_frobenius=Q̄^p + face_sum_ne_zero), architecture gmc2_of_nc2, contrapositive entry — ALL kernel-pure. The self-contained arithmetic engine of THM-2022 is DONE.
- **REMAINS:** §2 number-field descent (heavy, Mathlib-stocked); §3 DvdK (cite, ~person-months to formalize — see S95 roadmap); the number-field assembly (E(P^{pm0})/(pA0)! ≡ Q̄^p mod 𝔭) which inherently needs §2 (the normalization + mod-𝔭 reduction cross char 0→p). memory nc2-gmc2-lean-formalization-state updated. HYP-8805.

## boxeph-2026-07-21-S214 -- the rank-11 AP-core is the achiral vertex where codex's rank-or-Euler frontier meets (HYP-8855)

**Owner:** work incoming LRC progress (pull often); explore rank-11 / '11 private-coordinate relations' through 'relations = a tournament'.
> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## boxeph-2026-07-21-S213 -- chirality + toothpick A139250 + tournament evenness = one Lefschetz parity (HYP-8850)

## codex-2026-07-21-LRC-triangle-pair -- THM-2051 cutoff sharpened to `2^20`

- **PROVED refinement:** THM-965's exact pair covariance has only seven
  reduced ratios below `-1/220`. Packing pair covariances over all speed
  triangles gives the universal triangle floor `-1297/70070` and the global
  thirteen-speed pair-sum floor `-1297/2695`.
- **CLOSURE:** the resulting pair reserve is `1668/18865`. At Fejer height
  `2^20`, the support-`3..5` whole-product tail is below
  `221806728/2517633377`, leaving exact margin
  `43815012/138469835735>0`. Thus every hard row has a genuine support-`3..5`
  relation already at height `2^20`, matching THM-2052's first relation-code
  scale.
- **REFEREE:** a new frozen script checks the seven exceptional ratios, all
  `182` oriented exceptional triangles, the edge/triangle double count, and
  the final rational margin in normal and optimized Python. The older
  `2^21` referee remains intact as a reproducibility checkpoint.

## codex-2026-07-21-LRC-two-anchor -- THM-2052 explicit star-atlas refinement

**PULLED (incoming):** codex THM-2052 (PROVED) relation-rank descent -- pigeonhole harvests dim W_{91^6,3}(v)>=11 support-<=3 bounded relations (RANK 11); 12th pins <v> (finite terminal). Rank-or-Euler frontier HYP-8841 (OPEN): each peel => rank 11->12 OR Euler survivor chi>0 (my HYP-8845 = the Euler branch). Also pulled: rank-two geodesic terminal (dd899d4f2), pointed plane-transport HYP-8846, Fejer-BV THM-2051, active-owner circuit localization.

**THREE rank objects (the collision is the substance):** A=harvested code W (11->12, THM-2052 PROVED), B=ambient {k.v=0} (rank 12), C=AP-core L(AP) (rank 11, mac-mini S25).

**VERIFIED (object C, tournament lens; rank11_relation_lattice_is_the_transitive_tournament_boxeph_S214.py):** the 12 AP speeds {1..12} = transitive tournament T_12. L(AP) rank-11; Gram TRIDIAGONAL = path/Jacobi (three-term backbone); 11 basis d_k = adjacent-pair covering relations = transitive Hasse spine = the '11 private-coordinate relations'; minimal vectors = additive triples (30, kissing 60) = additive energy (S211); score 0..11=AP, char_A=x^12 = reify-ladder nullcone vertex (THM-1750); transitivity Vandermonde = braid A_11 (rank 11). PALINDROMIC v_i+v_{13-i}=13 => SELF-CONVERSE=ACHIRAL fixed point of reversal (S213), feeds q|v_i+v_j (THM-2047).

**SYNTHESIS:** the AP is the achiral iota-fixed point where BOTH frontier branches meet -- rank descent bottoms at the rank-11 AP vertex; Euler chi maximal at the AP (deep well chi=24, S212). Wall A = 'the rank-11 achiral transitive vertex is the unique optimizer' = reify-ladder + chirality.

**HONEST CAVEAT (THM-2052 + MISTAKE-224):** the tournament ORIENTATION is a structural/diagnostic LENS, NOT the descent carrier -- a binary orientation discards the signed coefficients+heights (THM-2052 = signed coding theory); surviving antisymmetric content into the proof = pair-sum law (THM-2047) + t->-t mirror-pairing (S212/S213/HYP-8845). Lens for locating the vertex; coding for the descent; chi for the Euler branch -- complementary. Not a new proof step. Artifacts: reflection the-rank11-ap-core-...-boxeph-S214.md, HYP-8855, script (+.out).
1. TOURNAMENTS EVEN: A000568=SC+2*(chiral pairs). THM-587 (PROVED): P_n(1)=A000568, P_n(-1)=SC=self-converse=antipodal Euler/Lefschetz number; SC=2,2,8,12,88,176 (n=3..8) all EVEN => A000568 EVEN for n>=3. Enumeration n<=6: (total2/SC2/0pairs),(4/2/1),(12/8/2),(56/12/22), count=SC+2pairs holds. THM-1830 BLUE self-comp atom=1 iff n odd = local single-fixed-point face.
2. TOOTHPICK A139250: sim=OEIS exactly; D4 mirror sym; =(#axis)+2*(#pairs); #axis ODD (central seed) => A139250 ODD all n>=1; first differences have exactly ONE odd term (seed); A139250(2^k)=(2^(2k+1)+1)/3=1,3,11,43,171 Jacobsthal.
3. LRC (S212): chi(G_delta) == #iota-fixed lonely pts (mod 2).

**TOOTHPICK<->CHIRALITY DICTIONARY:** central seed toothpick (axis-fixed) <-> achiral/self-converse classes; mirror toothpick PAIR <-> chiral pair {T,T^op}; D4 mirror <-> R=converse=complement=antipode; A139250 parity=#axis(mod2) <-> A000568 parity=SC=P_n(-1)(mod2); doubling A139250(2^k) <-> level-graded signed-cycle-index recursion. The toothpick fractal = the PICTURE of 'one achiral seed generating chiral mirror-pairs'. Unifies with S212: (P(+1),P(-1)) = Burnside count / Lefschetz fixed set = the chirality Euler class (S212 even chi + odd Gauss sum i*sqrt7).

**Honest:** THM-587, the toothpick OEIS/oddness/Jacobsthal facts, and S212 are each proved/verified; the contribution is the UNIFICATION (one reversal, one Lefschetz parity count==#fixed) + the explicit toothpick<->chirality dictionary. Artifacts: reflection chirality-toothpicks-and-why-tournament-counts-are-even-one-lefschetz-parity-boxeph-S213.md, HYP-8850, script (+.out).
**Owner:** push the DC(2)/planar-JC thread to its next decisive target and transfer the proved/disproved mechanisms toward LRC(14), while repeatedly integrating incoming work.

**DC(2) result (THM-2049, PROVED local/formal statement; not DC(2)):** in the exact Ore algebra `Q[x,q][ell;delta]`, `beta(sum a_k ell^k)=min_k(v_x(a_k)-2k)` is multiplicative and commutators raise beta by two. The associated bracket is `{ell,q}_0=2`. For the Weyl boundary symbols, the simultaneous grade-`g` correction map is `(A,B)->(8/3)(u-2)A+(2u^2-10u+9)B/9`; it is surjective because the two `u` polynomials are coprime. Thus the grade-six residual is exact. An exact ladder advances grades `6,...,13` to `14`; a formal beta-adic `[S,T]=1` lift exists. The open gates are polynomial termination and the coupled `D` relations. This corrects HYP-8802/8803's earlier suggestion that the first invariant grade might carry the obstruction.

**LRC no-go (THM-2050, PROVED):** AP13 and `AP13` with `12->26` have identical full local phase-height function germs on `|h|<1/728` at every unit point `a/14`, yet `M=1/14` and `M=1/12`. Local top data, even as a full germ, cannot determine global loneliness.

**Incoming synthesis:** THM-2047 supplies the lossless signed phase-height/Euler carrier; THM-2048 supplies the fiber-quantization pruning tax; HYP-8840 identifies GMC's constant-term/volume leverage and its zero-volume ceiling. Later pulls supplied THM-2048's genuine Cover14 gain and strengthened THM-2051: triangle-packing the exact pair covariances proves that the no-support-`3..5`-relation branch already at height `2^20` has positive safe volume. THM-2052 then proves every hypothetical counterexample already has eleven independent bounded three-support relations and lies in a finite two-anchor, one-projective-parameter atlas. The exact transfer is `volume/tax -> strict branch`, `Euler signed wall word -> tight branch`, and a labelled Noetherian rank-or-Euler rule inside the relation branch as the missing glue. No literal algebra map between GMC/DC and LRC is asserted.

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
## codex-2026-07-21 -- NC2/GMC(2) formalization: descent and face geometry closed; normalized residue assembly isolated

**Owner directive:** complete NC2/GMC(2) formalizations without big builds; work open mathematics instead of launching a broad build.

- Pulled/fetched the live mainline repeatedly and integrated the concurrent
  `wick_expansion`, multinomial-Lucas, `face_sum_frobenius`, non-cancellation,
  and `charge_radial` APIs into the modular library. No default-target or
  repository-wide build was run; checks were isolated modules or narrow
  dependency targets only.
- Moved `GMC2Reduction` into the library namespace while retaining its historical compatibility entry point. Added the negative-charge power bound and Mathieu--Zhao branch, moment-nullity for both strict orientations, `ChargeOneSided`, `NC2At`, and the exact `mathieuZhao_of_nc2At` wrapper.
- Added `GMC2FrobeniusFace`: balanced multiplicities on one exposed face have common Wick height; using a strict off-face term raises height; integral heights satisfy `A0+1<=A`; exact exponent support projects injectively to charge on the face.
- Added `GMC2FrobeniusResidue`: multivariate Frobenius proves any non-`p`-dilated multinomial channel of total mass divisible by `p` vanishes mod `p`; componentwise dilation satisfies multinomial Lucas; weighted face sums become `Q^p`; every strict normalized factorial quotient is divisible by `p`. This removes `p>m0` from the formal isolation lemma, though the paper retains it for the elementary carry narration.
- Added `GMC2WickChannels`: `E` additivity/finite sums, exact monomial evaluation and monomial-product arithmetic, the full `piAntidiag P.support m` expansion, and equation (1) with an explicit balance test and Wick factorial.
- Closed the former algebraic-descent gap end to end. The rational route builds
  `Q[z_i,z_i^-1] subset C`, descends every moment relation to a number field,
  and preserves a selected nonzero face seed. The integral route builds
  `Z[z_i,z_i^-1,Q(z)^-1]`, proves `Z` Jacobson in the pinned Mathlib, and
  specializes directly to a finite field of an output prime characteristic.
  Its universal property preserves every integral zero relation, so the
  p-dependent normalized moment can be selected after p is known.
- Formalized rational lowest-face existence and its exact Finset package,
  charge/coordinate dictionaries, the weighted tilted-height identity,
  balanced floor/equality, scaled natural height floor, strict off-face gap,
  and extraction of a concrete balanced reference channel from the DvdK seed.
- Added rational and integral universal moment/constant-term polynomials,
  exact support transport, integral-vs-rational seed cast, and the normalized
  moment polynomial with pre-reduction factorial cancellation.
- Added seed-preserving number-field and direct finite-field packages,
  `GMC2GoodReduction` as an independently checked longer route, and complete
  channel extension/dilation: mass, charge, bidegree, height, coefficient
  powers, injectivity, exact antidiagonal image, and sum reindexing.
- Expanded `GMC2Formalization` from five modules to the full checked spine. No
  `sorry`, `native_decide`, or custom axiom declarations; printed audits use
  only standard Mathlib foundations. DvdK is an explicit proposition passed
  to downstream theorems.
- Hostile audit PASS on THM-2022 and the new kernels. Precision repair: choose `pfrak` in `O_K`, use local units, and cancel `(p*A0)!` termwise in `O_(K,pfrak)` before reduction; never invert its zero residue.
- Honest remaining internal work is now the concrete instantiation of the
  abstract three-case residue sum with normalized Wick channels, followed by
  the final composition `DvdK1 -> NC2 -> GMC(2)`. Formalizing the published
  one-variable DvdK theorem itself remains a separate major project. Reflection:
  `the-prime-is-an-output-not-an-input-gmc2-integral-descent-codex-20260721`.
