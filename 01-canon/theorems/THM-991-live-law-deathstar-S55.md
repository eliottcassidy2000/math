# THM-991 — The live law of the canonical family (death-star-2026-07-17-S55)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCLiveLaw.lean`,
standard trio ×7). Source: HYP-7280. With THM-987's deep count, the canonical
census is now closed-form on both ends.

## The law

**liveCount((1,…,13), q) = 6 if 14 ∣ q, else 0** (recon: count AND exact set
{o·(q/14) : o a unit mod 14} verified for every q ∈ [2, 3000)).

Lean delivers the census-critical directions:
1. `resonant_live` + `liveCount_ge_six_of_dvd`: at q = 14m the six unit
   multiples o·m are live — the mod-scaling `(x·m) % (14m) = (x % 14)·m`
   reduces every band check to `14 ∤ (i+1)·o`.
2. `live_gap`: a live multiplier's 14 residues {k·p mod q : k ≤ 13} are
   pairwise ≥ ⌈q/14⌉ apart AND co-apart — both edges of the safe band through
   the DIFFERENCE SPEED (the canonical family is difference-closed: that is
   the property the rigidity consumes).
3. `live_forces_dvd` + `liveCount_eq_zero_of_not_dvd`: **off resonance the
   live set is EMPTY** — the block injection (pairwise gaps ≥ c pack 14
   residues into M/c + 1 blocks, so the spread is ≥ 13c; the wrap gap adds c)
   closes q ≥ 14⌈q/14⌉, forcing 14 ∣ q. No sorting machinery: the block map
   x ↦ x/c is injective outright, and the same-quotient cancellation keeps
   every step linear for omega.

The exact-6 upper bound at resonance (the m-net/δ-chain extraction) is the
named remaining piece; recon confirms it.

## Significance

The tight family's census is now: deep(q) exact (THM-987), live(q) = 6·[14∣q]
(this theorem, ≥/0 directions). Together with `canonical_lonely` (THM-987)
the q = 14 story is fully explained IN-KERNEL: six live equality witnesses,
one deep half-point, race won. The rigidity argument (difference-closure ⟹
gap packing ⟹ resonance) is the discrete Steinhaus/three-distance shadow and
should generalize to any difference-closed prefix family — the next
structural question alongside the generic Wagner-circle program.
