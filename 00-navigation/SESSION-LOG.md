> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## death-star-2026-07-22-S107 -- JC(2): the resonance dictionary + a manufactured-valuation orbit-product route (exploration, NOT a proof)

**Owner directive:** work the 2D Jacobian Conjecture; try creative new routes; explore unrelated threads for ideas.

- **VERIFIED ANCHORS** (exact bivariate poly arithmetic, jc2_resonance_dictionary_leadingform_deathstar_S107.py): det(D(I+H))=1+div(H)+jac(H1,H2); **homogeneous 2D Keller maps are EXACTLY shears** H=c(bx-ay)^d(a,b) (=> homogeneous 2D JC trivial; all difficulty is non-homogeneous, why BCW cubic reduction must raise dimension); leading forms of a Jacobian pair are **powers of a common binary form h**; the **descent reduces unless neither degree divides the other** (Abhyankar-Moh stuck stratum = boxeph S225 coprime/Lame); **places at infinity = distinct roots of h**, so single-root h=l^k => 1 place => Abhyankar-Moh => TAME, and a **JC(2) counterexample requires multi-root h (>=2 places)**.
- **THE RESONANCE DICTIONARY** (cross-thread synthesis): the tame-vs-hard split of THREE repo problems is ONE single-root/multi-root phenomenon -- DvdK unique vs coincident cycle (S101/S106); Hessian rank-1 one-sided (x+iy)^d vs rank>=2 (THM-1300 dim-3 ctrex, S103); JC(2) single-root h (tame) vs multi-root (open). Not a formal chain (MISTAKE-229; GMC=/>JC, S205) -- a structural claim that the same TOOL should hit all three resonant residuals.
- **NEW ROUTE (sketch, crux open):** transfer S106's manufactured-valuation Galois orbit-product (THM-2067) to JC(2)'s places at infinity. (a) Lefschetz reduces a counterexample to Q-bar coefficients [valid]; (b) roots of h lie in a number field so Gal permutes the places at infinity in orbits [valid -- the arithmetic supplies the orbit THM-2067 needs, which S205's LRC packet lacked]; (c) OPEN CRUX: exhibit a product of local invariants over a Galois orbit of places, rational yet with a valuation incompatible with jac=1 (the analogue of Pi(t)=ct vs the Vieta norm). Framing: boxeph S146 Euler ledger = topological reciprocity; this route = its arithmetic (valuation/norm, Weil-reciprocity-flavored) refinement; it sharpens boxeph S225's descent-FEASIBILITY engine to a CANCELLATION/valuation one.
- **Parallel tool (boxeph S231, same day):** F invertible <=> a coordinate certificate x=P(f,g) exists (Nullstellensatz), the EFFECTIVE face; boxeph's DvdK monomial certificate is the same coin. Dictionary gains a unified "certificate exists <=> tame; effective degree bound = open part" column.
- **HONEST:** anchors classical (Abhyankar-Moh / van der Kulk), verified in-repo + assembled; the dictionary is a synthesis not an implication; the route's crux is open. JC(2) remains open. HYP-8940.
## boxeph-2026-07-22-S231 -- eliminate DvdK for the residual 12% of straddling supports via a monomial certificate (kernel-pure Lean, HYP-8932)

**Owner:** get rid of DvdK for the remaining ~12% of straddling supports.
## codex-2026-07-22 -- the missing relative-Hunter channel and two all-height branches

- THM-2086 gives the exact Fourier split
  `w_pq=(5/7)rho_pq-(epsilon_p+epsilon_q)/7+R_h(p,q)`. The genuine term is an
  absolutely convergent sum over nonzero relations `ap+bq+ch=0`; the other
  channels are the global pair bulk and the two mixed-fold axes. Summed over a
  tree this is an identity, not an asymptotic.
- On THM-2081's sharp packet, the three contributions are
  `100421/1177176`, `-16117/4512508`, and `-2833331/203062860`, summing to the
  independently atomized positive margin `561797/8288280`.
- When `7|h`, every cross edge between a 7-nonmultiple and a 7-multiple has
  `epsilon_p=0`, `rho_pq=1/49`, and `R_h(p,q)=0`. Divisor completeness and
  hereditary primitivity force a spanning `K_(L,H)` tree; its margin is at
  least `5/294`. The whole apex-divisible guard branch is closed at all heights.
- When `7 not|h` and five speeds are divisible by seven, the high-high genuine
  channel vanishes. THM-1234 plus uniform `K_5` tree averaging gives restricted
  weight `88/1911`; paying the two low mixed folds leaves `23/1911`. Thus the
  live modular profiles have only one through four 7-divisible speeds.
- The unrelated high-frequency/BV route gives a second all-height result:
  `|w(B,q;h)-(1/7)mu(D_q cap C_h)|<=(q+h)/(3B)`. The exact overlap spectrum
  then closes `sum_(q!=B)q+6h<(17/1078)B`.
- The rank-seven residual is now simultaneously `7 not|h`, nonlacunary, and,
  by THM-2083, supported on a bounded guard/two-speed relation template. This
  restores exactly the incidence erased by THM-2082's scalar code wheel and
  mirrors GMC's support-versus-genuine-channel split without solving DvdK1.

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

