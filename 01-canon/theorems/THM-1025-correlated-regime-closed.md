---
id: THM-1025
title: THE CORRELATED REGIME IS CLOSED BY THE SAWTOOTH — the complementary floor pair: the INDEPENDENCE floor (THM-1012) is sharp when b ≫ a, the SAWTOOTH floor ρ ≥ 1/49 − g²/(4ab) (from THM-965's fold bound, g = gcd) is sharp when ab is large after gcd reduction, which comparable blocks automatically satisfy; with dilation invariance ρ(ga,gb) = ρ(a,b) every pair reduces to a primitive one, and only 15 primitive pairs have a′b′ ≤ 12 (all with ρ > 0, minimum 1/84) — so the pair floor is now available in EVERY regime, and on comparable 7-blocks the combined floor sum is positive in 400/400 (range 0.106–0.150) exactly where THM-1012 failed on all 400
status: verified exactly — sawtooth floor 4000/4000, dilation invariance 300/300, the 15-pair finite table computed, comparable-block floor sums 400/400 positive and each verified ≤ the exact Hunter floor; LEAN: the coprime case is ALREADY kernel-pure (LRCFloorTable.overlap_floor_rat, S347); remaining is the general-g version and dilation invariance
source: opus-2026-07-17-S357 (owner: attack the correlated regime with the sawtooth or beat structure)
depends_on: [THM-965 (the folded identity — the deviation IS the two folds), THM-1021 (the negative this answers), THM-1012 (the complementary independence floor), LRCFloorTable (the coprime case, already in the kernel), LRCSharpWallBound (the consumer)]
scripts: 04-computation/correlated_regime_opus_S357.py -> 05-knowledge/results/correlated_regime_opus_S357.out
---

# THM-1025 — the correlated regime, closed by the sawtooth

THM-1021 showed the dense core resists the INDEPENDENCE floor: its defect
2λ(a/b) → 2λ as a → b, swamping the leading 4λ² sevenfold. The answer is
that the sawtooth supplies the *complementary* floor.

**The sawtooth floor.** THM-965 gives
μ(D_a ∩ D_b) = [4ab + fold_M(S mod M) − fold_M(D mod M)]/(196ab) with
M = 14g, fold_M(r) = r(M−r) ∈ [0, M²/4]. Bounding the folds:

> **ρ(a,b) ≥ 1/49 − g²/(4ab)** (verified 4000/4000).

This is sharp exactly when **ab is large** — the comparable regime —
whereas THM-1012 is sharp exactly when **b ≫ a**. The two floors are
complementary, and between them they cover the whole ratio range.

**Reduction to primitive.** ρ is dilation-invariant, ρ(ga,gb) = ρ(a,b)
(verified 300/300), so every pair reduces to a coprime one, where the
floor reads ρ ≥ 1/49 − 1/(4a′b′). This is positive iff a′b′ > 49/4, i.e.
**a′b′ ≥ 13**. Only **15** primitive pairs fall below that threshold; all
have ρ > 0, with minimum ρ(1,12) = 1/84. The finite table closes them.

**The dense core regime.** Combining (sawtooth above threshold, table
below), the floor sum over the six consecutive pairs of a comparable
7-block is **positive in 400/400**, range 0.106–0.150 — and each was
checked to lie below the exact Hunter floor. THM-1012's floor sum was
≤ 0 on all 400. The regime THM-1021 identified as resistant is closed.

**Why this had to be the sawtooth.** The independence heuristic
ρ ≈ (2λ)² approximates the correlation away, so it can only be sharp when
the scales genuinely decouple. The folded identity does not approximate:
it *is* the correlation, written as two folds. Comparable combs are
correlated by construction, so only the exact form can see them — which is
precisely what THM-1021 predicted would be needed.

**Lean status — and a note against myself.** The coprime case of this
floor is ALREADY kernel-pure: `LRCFloorTable.overlap_floor_rat` (S347)
states exactly `1/49 − 1/(4ab) ≤ muNum a b/(14ab)`. I proved it three
sessions ago while thinking of it as a table entry for THM-964's
*separated* regime, and did not recognise that it is the tool that closes
the *correlated* one. What remains in Lean is the general-g form (fold
≤ 49g²) and dilation invariance — both short, given `muNum_folded`.
