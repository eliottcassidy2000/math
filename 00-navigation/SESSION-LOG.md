> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## death-star-2026-07-22-S104 -- GMC2 formalization: pinpointed + wrote the last capstone discharge (HeightWitnessSupplier); structurally correct + statements axiom-checked, but the proof hits a pathological whnf wall (>6.4M heartbeats). One perf-fix from clean DvdK1 -> NC2.

**Owner directive:** finish the GMC2 formalization.

- **STATE (codex's spine):** `GMC2NC2.nc2_of_dvdK1_of_heightWitnessSupplier : DvdK1 -> HeightWitnessSupplier -> NC2` is DONE sorry-free; the one remaining input is `HeightWitnessSupplier` (produce A0 + NormalizedHeightPackage from the face + nonzero seed) -- exactly the reference-channel/height-floor friction I flagged in S98.
- **I WROTE THE DISCHARGE** (`heightWitnessSupplier_holds`): wire `exists_reference_channel_of_nonzero_face_seed` (hface_tilted via `GMC2FaceDictionary.tiltedHeight_eq`) -> `normalized_height_obligations_of_face_reference` (hlower,hface,hrefBalanced,hrefMass,hrefHeight) -> `normalized_height_package_of_base` -> `⟨A0, ·⟩`. Types all verified; the derived `nc2_of_dvdK1 : DvdK1 -> NC2` + `gmc2_of_dvdK1` STATEMENTS axiom-checked correct.
- **BLOCKER:** the discharge PROOF hits a pathological `whnf` elaboration timeout -- fails at >6.4M heartbeats (32x default), while codex's spine compiles at default 200k. Almost certainly WHY codex left HeightWitnessSupplier as a hypothesis. So the capstone is one PERFORMANCE fix from clean DvdK1 -> NC2 (only DvdK1 hypothesis).
- **HONEST:** structural discharge identified + written, statements correct, but not compiling (whnf wall). Removed my GMC2HeightWitness.lean to keep the build clean; documented the exact composition + candidate perf fixes (prove NormalizedHeightPackage fields directly instead of composing the intermediate structure; irreducible/instance guards) in memory + a letter to codex. HYP: none new (documenting the existing capstone).
## codex-2026-07-21 -- THM-2069/2074 code-wheel and density-one LRC transfer

- **Higher deletion wheel:** THM-2069 proves that failure of some `k`-deletion
  gcd modulo `p` is exactly evaluation-code weight at most `k`. Cogirth is the
  first failure radius; full-rank deletion indices give exact CRT and primitive
  density products; rank deficiency is bad at every prime. The exact referee
  passes normally and under `python -O`.
- **Unrelated payoff:** the Paley-`e8` application closes THM-211. Its fourteen
  triangles are two distinct cyclic Fano planes, not two orientations of one
  Fano block. HYP-2430 is now the precise radius-16 first-cocircuit realization
  gate for a hypothetical `[72,36,16]` code.
- **Generic LRC theorem:** THM-2074 turns THM-2051's finite relation trap into
  density-one strict LRC(14): at most `R B^12` exceptions, with exact
  `R=25173854387233097811887443361297472`, versus
  `B^13/(13! zeta(13))+O(B^12)` increasing primitive rows. Any prescribed
  prime-power tower also contains whole certified speed congruence packets.
- **Honest residual after pulling THM-2073--2079:** LRC(14) remains open on the
  structured relation arrangement. MISTAKE-230 forbids descending the empty
  full-row safe set after deleting the two odd tails; terminal arguments must
  retain THM-2077/2079 owner/address sidecars. Hyperplane density and code-wheel
  primitivity do not supply phase height or those owners.

## boxeph-2026-07-21-S227 -- doubling homeomorphism + mirror-parity (LRC reduction); full two-charge DvdK in Lean (HYP-8920)

**Owner:** complete GMC(2) formalization; LRC math -- doubling as a continuous bijection + the unique safe-child condition.

**LRC (mapping the owner's prompts to codex's frontier):** 'making doubling a continuous bijection' = THM-2075 (doubling D:t->2t is a HOMEOMORPHISM between consecutive safe sets, chi/#comp/endpoint INVARIANT per address-sheet); 'the unique hild-safe condition' = THM-2073 (unique safe-child law, dyadic tower). VERIFIED (doubling_homeomorphism_meets_mirror_parity_boxeph_S227.py): phi_{2C}(t)=phi_C(2t) => G_{2C}=D^{-1}G_C; D 2-to-1 on S^1 but homeomorphism on each binary-address half (deck tau:t->t+1/2, distinct from my mirror iota); chi even for dyadic-seam covering sets.
- **Minimum bounded owner bank:** THM-2068 turns the THM-2066 census into an
  exact set-cover problem. Inside clocks `15..34`, seven clocks
  `{25,26,27,28,32,33,34}` are necessary and sufficient for all `59,880`
  primitive divisor-complete eleven-cores through maximum `24`; all banks of
  at most six undominated clocks were exhausted and every chosen clock has a
  private core.
- **Uniform structural descent:** after pulling THM-2072's fixed-bank no-go,
  THM-2073 transfers THM-775's forgotten safe-child mechanism to the strict
  `1/14` seam. Every imprimitive deletion has gcd two, the first four lifts
  are partitioned `2+1+1`, and descent iterates through divisor-complete
  quotient cores (including the new denominator-`26` shell) to a hereditary
  terminal. THM-2076's Haar-capacity lemma forces terminal size at least six,
  sharpening depth to at most five. Exact referees pass normally and under
  `python -O`. THM-2075 then proves that doubling is a homeomorphism along the
  tower: component/Euler counts persist, lengths and measure halve exactly,
  each component carries one constant binary address, and every endpoint has
  an inherited terminal-core owner. THM-2078 uses the terminal guard-height
  bound and an exact denominator-`8192` bitset audit to close every nontrivial
  tower with terminal maximum at most `24`: all `4,484,931` cores of ranks
  `6,...,10` were filtered, `30,594` were hereditary/divisor-complete, and no
  allowed guard survived even the necessary rational grid. LRC(14) remains
  open on the unbounded hereditary terminal lane and its address assignment.
- **Hostile parity correction:** MISTAKE-230 retracts HYP-8920's claimed
  `chi(G_S)=0 -> chi(G_terminal)=0` descent. THM-2075 starts at `G_C`, not at
  the empty `G_S`; the original tails kill both outer children. Mirror
  evenness of the nonempty terminal survives, but it does not discharge tail
  coverage. The S227 two-charge DvdK Lean work is independent and unaffected.
- **Correct mirror repair:** THM-2079 proves the safe-child section is
  reversal-equivariant. Mirror terminal components have bitwise-complementary
  addresses `a` and `2^r-1-a`; inherited endpoint owners match, while each
  original odd tail flips its owner parity. This halves the address search
  but explains why mirror parity alone cannot contradict the outer cover.

**SHARPENED REDUCTION (my contribution = S212 mirror-parity + codex THM-2073/2075):** a dyadic-seam disproof (S=2C ∪ {x,y}) has chi(G_{1/14})=0; chi doubling-invariant => terminal core also chi=0; chi even => needs chi=0 exactly. So Wall A (dyadic-seam case) <=> no hereditarily-primitive TERMINAL core has chi(G_{1/14})=0. Honest: full doubling 2-to-1 (chi doubles) vs per-sheet homeo (chi fixed); disproof chi=0 preserved either way.
## boxeph-2026-07-21-S226 -- kernel-pure Lean proof of the two-charge DvdK seed (HYP-8915)
**Owner:** work on completing the formalization of GMC(2).

**GMC(2) LEAN:** extended GMC2DvdKTwoCharge.lean with exists_nonzero_ct_pair' (SWAP orientation, index 0 = -n charge). Both pair orientations + dvdk1_pair KERNEL-PURE ([propext, Classical.choice, Quot.sound]) -- the full single-character DvdK1, DvdK-premise-free (any injective straddling Fin 2->Z up to relabeling). General DvdK1 (>=3 charges) = codex THM-2067 (Galois orbit-product) + height package = remaining boundary.

**Honest after correction:** the proposed LRC reduction is withdrawn; the
completed kernel-pure single-character DvdK1 leaf survives. Artifacts:
reflection doubling-homeomorphism-plus-mirror-parity-...-boxeph-S227.md,
HYP-8920, script (+.out), Lean GMC2DvdKTwoCharge.lean.
