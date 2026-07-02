/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S42)
-/
import TournamentH7

/-!
# The global axiom audit (the submission-state certificate)

One `#print axioms` per pillar theorem of the LRC(14) corpus.  Reading guide:
* `[propext, Classical.choice, Quot.sound]` — the mathlib-standard trio: fully kernel-checked.
* `+ Lean.ofReduceBool` — uses `native_decide` (compiled evaluation trusted); every such use
  is a finite rational computation independently mirrored in Python in `04-computation/`.
* `sorryAx` appearing ANYWHERE below is a release blocker — the audit exists to prove it absent.

The endgame theorem `lrc14_endgame` takes its two remaining obligations as PARAMETERS, so its
axiom profile certifies that everything BENEATH the final surface is proof, not promise.
-/

-- the endgame surface
#print axioms TournamentH7.LRC14Assembly.lrc14_endgame

-- the mathlib-track core
#print axioms LonelyRunner.conjecture_one
#print axioms LonelyRunner.conjecture_two
#print axioms LonelyRunner.isLonelyAt_of_forall_not_dvd
#print axioms LonelyRunner.exists_residue_one_of_tight
#print axioms LonelyRunner.isLonelyAt_of_unit_residue_missed

-- the certificate engine (modules 6 core + sharp form) and its packs
#print axioms LonelyRunner.WitnessCert.cert_lonely_tail
#print axioms LonelyRunner.WitnessCert.cert_ladder
#print axioms LonelyRunner.WitnessCert.cert_ladder'
#print axioms LonelyRunner.WitnessCert.depth3_pack_row
#print axioms LonelyRunner.WitnessCert.depth3_perLevel_row
#print axioms LonelyRunner.WitnessCert.certAP_all
#print axioms LonelyRunner.WitnessCert.cert3_all
#print axioms LonelyRunner.WitnessCert.cert4_all

-- the pattern layer (module 2, both directions) and module 3 (both frames)
#print axioms LonelyRunner.RatIntervals.length_inter_le_left
#print axioms LonelyRunner.RatIntervals.overlapQ_1_7
#print axioms LonelyRunner.RatIntervals.originNest_attained
#print axioms LRCCommensuration.seven_commensuration
#print axioms LRCCommensuration.seven_mul_volume_inter_periodic
