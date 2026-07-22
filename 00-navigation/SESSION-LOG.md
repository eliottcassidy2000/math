## death-star-2026-07-22-S108 -- asymptotic unit-distance problem: a positive count with a cancelling spectral kernel (exploration, NOT a new bound)

**Owner directive:** work the unit distance problem; creative new routes; mine unrelated threads.

- **OPENED the asymptotic Erdos unit-distance problem in the repo** (max u(n): upper O(n^4/3) Spencer-Szemeredi-Trotter, lower n^{1+c/loglog n} grid; 40-year gap). Fresh: the repo had only small-n exact (THM-431 u(21)=57, THM-432 Moser-lattice) + lattice-density (THM-412); ZERO on Szemeredi-Trotter / Guth-Katz / sums-of-two-squares. This is their n->inf continuation (THM-412's unit-distance layer IS the r2 count).
- **VERIFIED ANCHORS** (unit_distance_asymptotic_probe_deathstar_S108.py): grid count at squared-dist k ~ (1/2)r2(k)n; the r2-maximizer over k<=N is ALWAYS a product of primes =1 mod4 (25,65,325,5525=5^2.13.17), max r2 ~ N^{(log2+o1)/loglog N} = the n^{c/loglog n} growth; U/n runs 1.9->9.8 for n up to 1e4; upper bound = point-vs-unit-circle incidence, K_{2,3}-free (max common unit-neighbor=2 verified) => KST O(n^3/2) refined to ST O(n^4/3).
- **CROSS-THREAD PLACEMENT (the point of the session):** (a) extremal distance = product of primes =1 mod4 = smooth-number/clock object, same SHAPE as LRC divisor-complete cores (analogy, not transfer). (b) KEY -- unit distances are on the CANCELLATION side of the repo's dichotomy: U(P) is a positive edge-count, but the SPECTRAL bound (the only route to n^{1+o1} for structured sets) pairs |P^|^2 against the SIGN-INDEFINITE Bessel kernel sigma^=J_0, so the positive-definite/Hankel manoeuvre that unlocked GMC does NOT transfer (records WHY); files unit distances with LRC covering (signed sinc, boxeph S228/S229) and DvdK coincident-cycle (S101/S106) -- a positive count with a cancelling spectral certificate. (c) upper-bound obstruction = unexploited translation-CONGRUENCE (all circles are translates of one; ST is generic; Guth-Katz cracked distinct distances via SE(2) but not unit distances for exactly reason (b)).
- **ROUTE-SHAPE:** the sharp bound wants a cancellation-surviving arithmetic/topological certificate (the LRC chi-criterion THM-2047 / THM-2067 orbit-product template), NOT a positivity argument -- the structure-vs-randomness step for translation-covariant incidence.
- **HONEST:** anchors classical (r2 lower, ST upper), verified in-repo + connected to THM-412/431/432; (b)/(c) are a placement + route, not a theorem; the n^4/3 vs n^{1+o1} gap remains open. HYP-8945.
## boxeph-2026-07-22-S232 -- THM-2067 orbit-product core in Lean; GMC(2) now hinges solely on DvdK1 (kernel-pure, HYP-8941)

**Owner:** long session to fully formalize GMC(2); push/pull often.

**STATE UPDATE (crucial):** death-star-S105 discharged HeightWitnessSupplier kernel-pure => clean endpoints GMC2NC2.nc2_of_dvdK1 : DvdK1 -> NC2 and gmc2_of_dvdK1, NO height hypothesis. **GMC(2) <= NC2 <= DvdK1 ALONE.** Uniform DvdK1 = codex THM-2067 (Galois orbit-product); death-star-S106 map = 4 Mathlib-ready pieces + 1 valued-field gap (THM-1550 = unramified Hensel).

**FORMALIZED (GMC2OrbitProduct.lean, kernel-pure [propext,Classical.choice,Quot.sound], = S106 §5 target 1, gap-INDEPENDENT abstract core, 4 thms):**
- prod_smul_eq_prod_pow_card_stabilizer: transitive action => prod_g f(g.x) = (prod_a f a)^|Stab x| (Fintype.prod_fiberwise + fiber<->stabilizer bijection).
- card_stabilizer_eq_card_stabilizer: stabilizer order constant on transitive action (conjugation bijection).
- prod_pow_card_group_eq: the ORBIT-PRODUCT EQUATION -- equivariant f + G-fixed subset product p => p^|G| = (prod_Omega f)^(|S|*|Stab|).
- valuation_zero_of_prod_fixed: the CONTRADICTION ENGINE -- v(prod_Omega f)=0 => v(p)=0 (THM-2067 gets contradiction from v(p)=v(ct)=1).

**Gap-independence verified:** Mathlib has Polynomial.Gal.galAction_isPretransitive + MulDistribMulAction p.Gal SplittingField (transitivity/action free), but instantiation needs a VALUATION ON THE SPLITTING FIELD of C(t) = the ramified extension = exactly the S106 gap (-> unramified Hensel via X=sZ, t=s^M). So the abstract core is the clean gap-independent boundary.

**Remaining to full GMC(2)** (coordinating with death-star): Galois wrapper (instantiate via galAction_isPretransitive), Check A (CT(Lambda^m)=[u^Mm]R^m), irreducibility of X^M-tR/C(t), Vieta, THM-1550/Hensel (the gap).

## codex-2026-07-22 -- rank-one wheel boundary and residue-incidence transfer

**THE 12%:** unique-channel bypass at ANY mass (S230) reclassifies 4 supports free -> 87.9% (102/116); residual = exactly 14/116 = 12.1%, the coincident-channel/symmetric supports (involution u->-1/u, f(-1/u)=-f(u), THM-2070, pairs compositions => card>=2 at every mass, no single-power certificate).

**THE TOOL -- MONOMIAL CERTIFICATE:** CT(f^m)=constantTermRelation(q,m) homogeneous deg m in coeffs. Support DvdK-free elementarily IFF ideal <CT(f^m)> contains a monomial mu=prod x^e (unit on the coeff torus): mu=sum g_m CT(f^m) => some CT(f^m)!=0 on torus, same finite mass set for all torus points. = V(I) cap torus empty = DvdK1 (true, THM-2067) but as an explicit finite certificate, replacing Galois by exact Q-linear-algebra (graded per degree, no numpy/sympy). RESULT: ALL 14 residual supports certified at degree<=6 => every straddling support size 3-4 in +-4 is DvdK-free with NO Galois (88% unique-channel + 12% monomial-cert). Paradigm {-2,-1,1,2}: 12 b^2c^2 = (3ad+9bc)CT(f^2) - CT(f^4), two-mass {2,4}.

**FORMALIZED (kernel-pure [propext,Classical.choice,Quot.sound]):**
- GMC2DvdKMonomialCertificate.dvdk1_of_monomialCertificate: the ENGINE -- identity prod X^e = sum_{m in M} g_m constantTermRelation(q,m) + c_i!=0 => exists m>=1, CT(f^m)!=0 (any field Q-algebra). No positivity/uniqueness/Galois. Unique-channel = the M={m0} special case; this generalizes to the coincident-channel stratum.
- GMC2DvdKResidualExample.dvdk1_neg2_neg1_1_2: DvdK1 fully discharged for the symmetric paradigm {-2,-1,1,2}. CT(f^2), CT(f^4) evaluated by decide-on-piAntidiag-filter; certificate via linear_combination; two-mass disjunction.

## codex-2026-07-22 -- live-pull correction and the true depth-five exit

- MISTAKE-231 retracts the first, colliding THM-2080 terminal-resonance gate:
  it reversed the implication `G_Q subset E_h`, which covers `E_h^c`, not
  `E_h`. Its exact mixed-radius fold formula survives inside the valid THM-2080.
- The first-reserved THM-2080 is valid: for odd `h`, every 1/14 danger comb
  overlaps the 1/7 guard comb by at least `1/42`, with equality only at
  `q=6h`. Hunter's star makes six distinct danger combs insufficient to cover
  the guard complement. Hence terminal size is at least seven and the dyadic
  tower has depth at most four.
- For THM-2081's rank-seven deficit, the THM-2080 fold formula shows that
  negative mod-14 correlation pays exactly for multiplicity excess in a cover
  of `E_h^c`; this is the live scalar boundary, not a small-ratio gate.

## boxeph-2026-07-22-S229 -- kernel-pure Lean: the unique-channel DvdK-free criterion (any support) + the cancellation/inclusion-exclusion dictionary (HYP-8930)
## boxeph-2026-07-22-S230 -- bypassing the GMC2 DvdK dependency for the unique-channel class (kernel-pure Lean, HYP-8931)

**Owner:** creatively bypass the GMC2 dependency on DvdK, or find an easier formalization.

**THE LOCALIZATION:** codex's spine consumes DvdK1 in exactly ONE place -- GMC2FaceSeed.exists_nonzero_lowest_face_seed, only to produce "∃m>=1, CT(lowest_face^m)!=0". Everything around it (slope lambda, level delta, exact face F, straddling, charge-injectivity, coeff-nonzero) is DvdK-free Newton-polygon geometry (GMC2.exists_rational_lowest_face_finset). So the whole GMC2 DvdK dependency is ONE seed implication.

**THE BYPASS (kernel-pure, [propext,Classical.choice,Quot.sound]):**
- dvdk1_of_uniqueChannel (GMC2DvdKUniqueChannel.lean): the exact DvdK1 conclusion (∃m>=1, CT(f^m)!=0), NO premise, whenever some size has a unique balanced composition. Discharges the interface input for the unique-channel class.
- exists_nonzero_lowest_face_seed_of_uniqueChannel (new file GMC2DvdKUniqueChannelBypass.lean): a DROP-IN replacement for codex's exists_nonzero_lowest_face_seed -- identical conclusion, DvdK1 premise replaced by LowestFaceUniqueChannel P; reuses codex's geometry lemma verbatim, swaps only the final DvdK call for S229 ct_ne_zero_of_unique_balanced.
NET: every P whose lowest face has a unique channel needs NO DvdK axiom for NC2, only HeightWitnessSupplier. Covers death-star-S101/HYP-8878's 84%.

**HONEST BOUND:** residual = coincident-channel stratum (card>=2, symmetric/resonant 16%). The involution u->-1/u (f(-1/u)=-f(u), THM-2070) pairs compositions => even multiplicity at every mass => never unique (e.g. {-2,-1,1,2}); irreducible DvdK = codex THM-2067 (Galois orbit-product). BLOCKED all elementary erasure routes: face-simplification (THM-2070 any Laurent poly is a GMC lowest face), saddle (S222 retracted), char-p (harder: multinomials vanish mod p, Frobenius gives CT=0), genericity (feasibility!=cancellation).

**Next step (proposed to codex):** parameterize the descent by the seed lemma (take exists_nonzero_lowest_face_seed's conclusion as input) so both DvdK1 and my unique-channel seed drive it => DvdK-axiom-free NC2 for the 84%.

### Relation-free all-height lane

- THM-2083 proves a uniform short-relation alternative. If every triple
  `(h,q_i,q_j)` has relation height tending to infinity, character convergence
  sends mixed overlaps to `2/49`, restricted pair weights to `5/343`, the
  scalar deficit to zero, and the maximum tree to `30/343`. Such packets have
  a large positive relative-Hunter margin.
- Consequently there is an absolute `H_7` such that every rank-seven terminal
  containment satisfies `a h+bq_i+cq_j=0` for some nonzero coefficients of
  height at most `H_7`. This is uniform and rigorous but currently
  ineffective; the next task is a Fejer constant chase or template-by-template
  endpoint/CRT discharge.
- A 500-row structured stress test with maxima up to `1200` found no relative-
  tree failures; generic margins clustered near the predicted `30/343`.
  This is evidence only and is not used in THM-2083.

### Effective Selberg relation gate

- THM-2085 replaces THM-2083's ineffective constant by `H_7=57`. Vaaler's
  degree-57 interval pair is assembled into the signed box minorant
  `prod U-sum_r(U_r-L_r)prod_(s!=r)U_s`; this is pointwise below the box
  indicator even though the minorant may be negative.
- If every guard/two-speed triple has relation height greater than `57`, finite
  Fourier exactness gives `I_i>=5363/164836` and every outside-guard edge
  `w_ij>=655135/66923416`. Hence
  `tau-Delta>=6435/8365427>0`, contradicting containment by THM-2081.
- This does not conflict with THM-537's analytic-minorant wall: no nonnegative
  polynomial supported inside an arc is asserted. The signed coordinate-
  labelled correction is the entire repair. The exact rational referee passes
  normally and under `python -O`; the same certificate is nonpositive through
  degree `56`.
- Incoming THM-2082 sharpens the next split. Translated prime grids retain
  projective residue incidence, whereas the rank-one code wheel retains only
  deletion primitivity. Its frozen unbounded family rules out extracting a
  height cutoff from scalar code/cogirth, fold, and divisor sidecars alone.

### Short-relation cut, not one relation

- THM-2087 observes that THM-2085 needs only a relation-free spanning tree.
  Under containment the graph `ij <=> lambda(h,q_i,q_j)>57` is therefore
  disconnected. The complementary short-relation graph contains a spanning
  `K_(a,7-a)` and at least six indexed height-57 relations.
- If any speed has a two-term guard relation, it is `q/h=r/s` with coprime
  `r,s<=57`; terminal guard oddness makes `s` odd. Otherwise all cross-relation
  speed coefficients are nonzero. Eliminating one opposite-component anchor
  puts every speed in `A_q h+B_q x+C_q q=0`, `C_q!=0`, with coefficient height
  at most `2*57^2=6498`.
- This is the guard-localized two-anchor star missing from the global relation
  alternative. It is still unbounded and does not prove LRC(14), but it exposes
  a local coefficient-row interface for THM-2062/2065 and the residue-incidence
  interface for THM-2082. The prior dyadic guards and original seam tails still
  need to be reconstructed in the same two parameters before the thirteen-row
  circuit-ray theorem can be applied.

### Cut-matrix rank terminal

- THM-2088 treats the complete short-relation cut as a coefficient-row
  matroid. Because both speed coefficients on every cross edge are nonzero,
  leaf elimination proves every cut forest is independent. A spanning tree
  supplies rank six; the positive terminal tuple caps rank at seven.
- At rank seven, seven original height-57 support-three rows recover the
  primitive `(h,Q)` by maximal minors. Hadamard gives the explicit finite
  bound `floor(sqrt(3^7*57^14))=91421508108581` on every coordinate.
- At rank six, the six tree rows define a two-parameter template and every
  chord is persistent on it. The cut types `1+6`, `2+5`, `3+4` have
  respectively `0`, `4`, `6` chord tests. Thus the all-height rank-seven
  residual is now: bounded guard ratio, finite rank-seven terminal, or
  persistent rank-six cut template.

### Persistent means flat affine holonomy

- THM-2089 imports THM-1226's forgotten relation-cycle holonomy into the new
  cut. Writing each edge as `q_j=alpha_e q_i+beta_e h`, a fundamental chord is
  persistent exactly when `product alpha=1` and the transported offset sum is
  zero. Clearing denominators gives the two exact integers `D_C-N_C` and
  `R_C`; a nonflat simple cycle gives a guard/base relation of height at most
  `6*57^6=205778683494` and is already THM-2088's finite rank-seven branch.
- In the flat branch, a tree gauge gives `q_i=u_i(z+v_i h)` with rational
  gauge height at most `6*57^6`. Integrality is the explicit congruence lattice
  `D_i | N_i z+R_i h`, of index at most `57^36`; positivity is one interval in
  `z/h`, and distinctness removes at most 21 rational walls unless a pair is
  identically equal, in which case the template is empty.
- The exact referee checks 1,500 signed clearing identities, all 28 squares of
  flat `K_(2,5)`/`K_(3,4)` gauges, and a perturbed connection with six detected
  nonflat squares under both normal and optimized Python.

### Global splice or frozen terminal

- THM-2090 lifts a cut spanning tree to six independent support-three rows on
  the original thirteen speeds: `a h+bq_i+cq_j=0` becomes
  `2a(16h)+b(32q_i)+c(32q_j)=0`, of height at most `114`.
- These six rows sit inside THM-2052's rank-at-least-eleven height-`91^6`
  relation code. Rank twelve is finite. At rank eleven, restrict its
  two-dimensional kernel to the last-guard/terminal block. Restriction rank
  two makes the last guard and some terminal speed independent anchors, so
  THM-2052's triple pigeonhole gives a height-`91^6` star for all thirteen
  original speeds.
- Restriction rank one is sharper: after undoing fixed dyadic scaling, every
  admissible terminal vector is a rational multiple of one primitive `(h,Q)`.
  Integrality makes the multiplier integral, terminal primitivity makes it
  `+-1`, and positivity makes it `+1`. Hence `(h,Q)` is literally frozen and
  the remaining integer points form an affine lattice line moving only the
  three earlier guards and two original tails.
- The six local height-114 circuits already satisfy THM-2065's persistent-
  circuit alternative. Asking Fejer for one unspecified extra relation is now
  provably vacuous; the finish must use circuit multiplicity/location, owner
  addresses, translated-grid incidence, or exact relative-Hunter phase data.

### Frozen terminal lines are finite

- THM-2092 reuses THM-2077's all-depth needle invoice at `r=4`:
  `max(S)<=32B*max(1,8/6)=128B/3`, where `B=max(Q)`.
- THM-2090's restriction-rank-one branch literally fixes `(h,Q)`, hence fixes
  `B`; its five-coordinate affine outer line therefore meets a finite absolute
  box. Across the finite relation-template bank this whole branch is finite.
- THM-2088's uniform terminal bound `B<=91421508108581` now lifts to the
  explicit original-row bound `max(S)<=3900651012632789`. These rows remain
  unenumerated; “finite” is a reduction, not a safety verdict.
- The only potentially unbounded no-pair branch left by this pipeline is
  THM-2090's rank-eleven, restriction-rank-two star anchored at the last guard
  and one terminal speed. It retains the six height-114 cut circuits,
  THM-2086's modular/nonlacunary filters, and the dyadic owner address.

### Dyadic global-star finite terminal

- THM-2093 identifies the original depth-four speed valuations as an exact
  nested evaluation-word support flag of weights `2,3,4,5,6` modulo
  `2,4,8,16,32`. At the binary level this is THM-2069's weight-two deletion
  codeword; the higher ring layers retain the guard staircase that the binary
  matroid forgets.
- If an outer coordinate is `2^t u`, `0<=t<=3`, its primitive star relation
  has private coefficient divisible by
  `2^(4-t) gcd(h,q)/gcd(gcd(h,q),u)`. In particular the two tails and three
  earlier guards force private factors `16,16,8,4,2`, with normalized relation
  `A h+2B q+c u=0` and `A==c mod 2`.
- Rooting the six persistent terminal cut rows at the chosen terminal anchor
  gives a prime-by-prime path invoice. Terminal primitivity forces
  `v_p(gcd(h,q))<=6 floor(log_p 57)`, hence
  `gcd(h,q)|lcm(1,...,57)^6`; no prime greater than `57` can be shared.
- Clearing the eleven private star denominators gives an integral two-
  coordinate presentation with coefficient radius at most
  `sqrt(2)*(91^6)^12`. THM-2053 then bounds its primitive anchor direction,
  and the whole original row satisfies
  `max(S)<2912*lcm(1,...,57)^6*(91^6)^13`.
- This removes the last potentially unbounded no-pair branch. The separate
  THM-2087 bounded guard-ratio branch and all finite-box discharge remain.
  The exact referee checks 1,000 valuation flags, 253,680 small staircase
  relations, sharp path budgets at all 16 primes through 57, and 500
  denominator-clearing stars under normal and optimized Python.

## death-star-2026-07-22-S105 -- GMC2 formalization CAPSTONE CLOSED: the whnf wall is SOLVED; HeightWitnessSupplier discharged kernel-pure; clean DvdK1 -> NC2 and DvdK1 -> GMC(2) now compile.

**Owner directive:** work creatively at whnf-tuning, or how it can be bypassed.

- **SOLVED the S104 blocker.** The direct existential wrapper into `HeightWitnessSupplier` -- which codex's own aggregator docstring recorded as "exceeds Lean's elaboration budget" (>6.4M heartbeats, 32x default) -- now compiles at DEFAULT 200k heartbeats, kernel-pure.
- **ROOT CAUSE (isolated by binary search + LOW-heartbeat fast-fail):** the `whnf` explosion is in APPLYING `GMC2FaceReferenceChannel.exists_reference_channel_of_nonzero_face_seed`. The elaborator repeatedly reduces `P.coeff` (a `Finsupp` lookup) to weak-head-normal-form while unifying the supplied coefficient `fun s : F => P.coeff s` against the extractor's parameter. (`hface_tilted` alone always compiled; the fault was purely the coefficient unification.)
- **FIX (the bypass):** seal the coefficient behind an opaque local def -- `set c : F -> C := fun s => P.coeff s with hc` at the proof start -- removing every `P.coeff` occurrence `whnf` could unfold, so unification succeeds structurally. One line. No `maxHeartbeats` bump, no `irreducible`, no axiom, no `sorry`. The downstream obligation/package lemmas are coefficient-free, so the opaque `c` never needs unfolding.
- **DELIVERED** (`TournamentH7/GMC2HeightWitness.lean`, now imported by the `GMC2Formalization` aggregator): `heightWitnessSupplier_holds : HeightWitnessSupplier`, and the clean endpoints `nc2_of_dvdK1 : DvdK1 -> NC2` and `gmc2_of_dvdK1`. All three `#print axioms` = [propext, Classical.choice, Quot.sound]. Full `lake build` of the aggregator succeeds (8509 jobs). Corrected the aggregator docstring that documented the wall as unresolved.
- **NET:** the GMC(2) descent endpoints now depend on ONLY the one published analytic input (`GMC2DvdKInterface.DvdK1`) -- no `HeightWitnessSupplier` hypothesis. Complements boxeph's S226-S229 (kernel-pure DvdK1 for two-charge, positive-coefficient, and unique-channel/arbitrary support -- S229 mechanizes my S101/HYP-8878 unique-cycle criterion), which removes the OTHER hypothesis on the 84% DvdK-free stratum; the residual is the card>=2 general-complex DvdK1 = codex THM-2067. GENERAL LESSON: an unexplained `whnf`/isDefEq timeout applying a lemma to an argument built from `P.coeff`/`Finsupp`/projections is a defeq-unfolding blowup -- seal the subterm opaque with `set`/`generalize`, don't raise heartbeats. HYP: none new (closing the existing capstone).
## death-star-2026-07-22-S104 -- GMC2 formalization: pinpointed + wrote the last capstone discharge (HeightWitnessSupplier); structurally correct + statements axiom-checked, but the proof hits a pathological whnf wall (>6.4M heartbeats). One perf-fix from clean DvdK1 -> NC2.
**Owner:** aim earnestly at formalizing DvdK; make it simpler / circumvent it; spill over to LRC.

**Honest:** the DvdK-free (unique-channel) side is now kernel-pure in Lean for arbitrary support, subsuming S226 and complementing S228; the coincident-cycle (card>=2) stratum remains the THM-2067 Galois frontier. The synthesis is a reading of proved theorems (THM-1820/1810/406/515/671), not a new theorem. Artifacts: reflection the-unique-channel-dvdk-in-lean-...-boxeph-S229.md, HYP-8930, Lean GMC2DvdKUniqueChannel.lean (5 theorems).
**Honest scope:** the DvdK dependency is localized to one seed implication and discharged kernel-pure for the unique-channel 84%; NOT a full bypass (coincident-channel 16% = THM-2067). NC2-level wiring is a one-line codex-owned change. Artifacts: reflection bypassing-the-gmc2-dvdk-dependency-...-boxeph-S230.md, HYP-8931, Lean GMC2DvdKUniqueChannelBypass.lean + dvdk1_of_uniqueChannel.
**Honest:** general engine + 1 concrete instance in Lean; other 13 certificates script-verified (same decide+linear_combination recipe, mechanical). NOT a uniform degree bound (= effective DvdK, Sturmfels/ESV open) -- but any explicit support is a finite, Galois-free certificate computation. Artifacts: reflection eliminating-dvdk-for-the-residual-12-percent-...-boxeph-S231.md, HYP-8932, dvdk_monomial_certificate_residual_boxeph_S231.py, GMC2DvdKMonomialCertificate.lean, GMC2DvdKResidualExample.lean.
**Honest:** NOT full GMC(2); the kernel-pure abstract orbit-product core of THM-2067 (the sole remaining GMC(2) input) + contradiction engine, + a verified map of the remaining instantiation to the one S106 gap. 2 mid-session checkpoints pushed. Artifacts: reflection the-thm2067-orbit-product-core-...-boxeph-S232.md, HYP-8941, Lean GMC2OrbitProduct.lean.

