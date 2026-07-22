> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## codex-2026-07-22 -- rank-one wheel boundary and residue-incidence transfer

- THM-2082 specializes THM-2069 exactly: a scalar row has code enumerator
  `1+(p-1)z^w_p`, and hereditary primitivity is `w_p>=2` for every prime.
  Thus the one-deletion wheel collapses to modulus one and cannot manufacture
  a terminal height variable.
- A translated-prime-grid lemma is the positive replacement: carrier ratios
  are uniformly safe, each noncarrier consumes at most `ceil(p/7)` residues,
  and the guard consumes at most `ceil(2p/7)`. A strict leftover is an actual
  point of `G_Q minus E_h`; failure is an explicit carrier-ratio branch.
- The unbounded rows `{1,...,s-1,360360*32^j}`, `s=7,...,10`, freeze their full
  prime code profiles while passing hereditary/divisor, quarter, relative
  height, and scalar-fold filters. Every row nevertheless escapes at `3/31`.
  Hence those scalar sidecars alone cannot prove a height cutoff.
- Two divisor-complete hereditary rows have the same full-support `p=17` code
  `1+16z^8` but zero versus twelve safe grid residues. This is the exact
  support-incidence loss familiar from tournament score fibers, the
  `[72,36,16]` cocircuit-design gate, and GMC's multiple-channel cancellation:
  the live LRC carrier is labelled projective residue incidence, concretely
  THM-2081's guard-restricted event graph.

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

