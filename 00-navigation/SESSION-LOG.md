> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

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

**Owner:** work the next Lean DvdK target; mine past results on cancellation and inclusion-exclusion.

**FORMALIZED (GMC2DvdKUniqueChannel.lean, kernel-pure [propext,Classical.choice,Quot.sound], 5 thms):**
- ct_ne_zero_of_unique_balanced: a UNIQUE balanced composition of size m => CT(f^m) = one multinomial term => != 0 for EVERY complex c (any support, any #charges). = death-star-S101/HYP-8878 unique-primitive-cycle criterion, mechanized, coefficient-independent (the DvdK-free 84%).
- two_balanced_of_ct_zero: contrapositive dichotomy -- CT=0 + one channel => a SECOND channel; cancellation needs >=2 coincident cycles.
- ct_ne_zero_of_card_eq_one + two_le_card_balanced_of_ct_zero: cardinality form -- |balancedSet q m|=1 => CT!=0 (84%), CT=0 => card>=2 (hard 16%). Coincident-cycle count = a Lean Finset.card.
- two_charge_via_unique: S226 pair theorem = the Fin 2 instance (subsumed).
DvdK-free zone in Lean now = every support with a unique minimal channel; residual card>=2 = codex THM-2067 Galois.

**MINE -- cancellation and inclusion-exclusion are ONE object:** THM-1820 (bridge: same moment-nullcone, discriminant=|alphabet|; LRC bounded {0..13} => Bonferroni IE terminates at k=13 => finite B5 cert; GMC unbounded => depth>=d+1). THM-406 M1b (loneliness p0=Sum_j(-1)^j S_j = IE alternating overlap-volume sum) + M2 ({p0=0} all-orders, Bonferroni one-sided = Vitali wall). THM-515 (theta/sinc lattice sum alternates, THM-504 cross-level). THM-1810 (bosonic/fermionic permanent/determinant = the shared sign-involution; THM-2070 dihedral witness). THM-671 (B5=2052/7^5~+0.1221 first positive Bonferroni truncation).

**HONEST SELF-CORRECTION of S228:** LRC covering is NOT "all-cancellation, no positive regime". There IS a positive regime -- odd-Bonferroni B5 (THM-671) + codex guard-capacity (THM-2076 union bound + THM-2080 overlap floor mu(D_q cap E_h)>=1/42) -- handling low-resonance cores (the 84% analogue), which provably cannot reach the all-orders floor on resonant/AP cores (THM-406 M2; S_5 explodes 156-172x iid, THM-686), EXACTLY as positive-coeff/unique-channel DvdK cannot reach the general complex case (THM-2067). Same positive-vs-cancellation split; only the alphabet (finite vs unbounded) differs -- why LRC cancellation terminates and GMC's does not.
### Rank-seven continuation

- THM-2081 restricts the seven danger combs to the guard complement and runs
  Hunter on that relative event graph. If `I_q=measure(D_q intersect E_h)`
  and `tau_h` is the maximum restricted spanning-tree weight, containment
  forces `tau_h<=2/7-sum I_q`; a strict reverse inequality gives actual safe
  mass outside the guard.
- Exact rational atomization checks all `131` hereditary divisor-complete
  rank-seven cores through height `24` and `4,120` allowed odd guards. The
  mixed-overlap scalar leaves `1,322`; the relative tree leaves zero. Its
  worst positive margin is `561797/8288280` at
  `Q=(1,9,10,11,13,14,24), h=23`.
- The all-height inequality is open. The new carrier is projective and
  scalable--three-frequency outside-guard atoms plus a weighted graphic
  matroid--but the common-gcd/endpoint phase classification has not yet been
  proved.

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

