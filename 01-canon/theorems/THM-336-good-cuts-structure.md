---
theorem: THM-336
name: Good Cuts Structure Theorem
status: PROVED
session: opus-2026-05-27-S2
verified: computationally n=3..10 (exhaustive for n≤7 via brute force)
formalized: kind-pasteur-2026-05-29-S6 (`TournamentH7.GoodCuts`)
depends_on: THM-330
lean: 04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean
---

## Statement

In the tiling model on n vertices (tiles = non-consecutive arc pairs, base path n-1→...→0), call a cut k **good** if at least one tile crossing it is in the upward orientation. Let G(τ) denote the set of good cuts of tiling τ.

**Part 1 (Impossibility of 1 good cut):** No tiling has exactly |G(τ)| = 1. That is, G(τ) is never a singleton.

**Part 2 (Exact counts for small |G|):**
- |G|=0: exactly 1 tiling (all tiles downward = transitive tournament).
- |G|=1: 0 tilings (**impossible**).
- |G|=2: exactly n−2 tilings.
- |G|=3: exactly 5(n−3) tilings.
- |G|=4: exactly (n−4)(n+95)/2 tilings.
- |G|=n−1 (SC): exactly SC(n) tilings (the SC tiling count).

## Proof of Part 1

Any tile (x,y) satisfies x ≥ y+2 (non-consecutive). An upward tile (x,y) witnesses every cut k with y < k ≤ x. The number of such cuts is x − y ≥ 2. Hence any upward tile witnesses **at least 2** good cuts simultaneously. Therefore if G(τ) ≠ ∅, then |G(τ)| ≥ 2. □

**Lean formalization (kind-pasteur-2026-05-29-S6):** `04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean`
now proves the reusable axiom-free core in the concrete staircase tile model:

- `goodCuts_empty_iff_all_down`: bucket 0 is exactly the all-down tiling.
- `goodCuts_nonempty_iff_exists_upward_tile`: nonempty good-cut support iff
  at least one tile is upward.
- `goodCutCount_eq_zero_iff_all_down`: cardinality form of bucket 0.
- `goodCutCount_pos_iff_exists_upward_tile`: positive count iff some tile is
  upward.
- `goodCutCount_pos_iff_not_all_down`: positive count iff the tiling is not
  all-down.
- `two_le_goodCutCount_of_upward_tile`: one upward tile forces at least two
  distinct good cuts.
- `goodCutCount_ne_one`: no tiling has exactly one good cut.
- `goodCutCount_eq_zero_or_two_le`: the strengthened dichotomy
  `goodCutCount = 0 ∨ 2 ≤ goodCutCount`.
- `goodCuts_empty_or_two_le_card`: set-cardinality form of the dichotomy.

The audit target `TournamentH7.Verify` confirms these theorems depend only on
Lean foundations (`propext`, `Classical.choice`, `Quot.sound`), not on project
axioms; see `05-knowledge/results/lean_goodcuts_bucket_strengthening_kind_pasteur_s6.out`.

## Proof of Part 2: |G|=2 gives n−2

**Claim:** If |G(τ)| = 2, then G(τ) = {k, k+1} for some k ∈ {1,...,n-2}.

**Proof:** Let G(τ) = {k₁, k₂} with k₁ < k₂. Suppose k₂ > k₁+1 (non-consecutive). Any upward tile witnessing k₂ has x ≥ k₂ and y < k₂. We need this tile NOT to witness any cut j with k₁ < j < k₂ (else j ∈ G). For such a tile: x ≥ k₂ > j, and y must satisfy y ≥ j. So y ≥ j > k₁. But then y > k₁, and this tile does NOT witness k₁. Similarly, any upward tile for k₁ has x < j (to avoid j) and so has x < k₂ (doesn't witness k₂). So witnessing k₁ and k₂ requires two separate tiles, one with x < j and one with y > j. But then k₁ and k₂ can both be witnessed simultaneously without any tile crossing j. 

BUT: for cut k₁ to be witnessed with x < j: we need x ≥ k₁ and y < k₁ and x < j. The set of cuts this tile witnesses = {k₁,...,x} ⊆ {k₁,...,j-1}. None are outside {k₁,k₂}? Only if x ≤ k₁ is possible (x ≥ k₁ and x < j with j = k₁+1 means x = k₁, but tiles have x ≥ y+2 ≥ 2 and x ≥ k₁ = x means y ≤ k₁−2). The tile (k₁, k₁−2) witnesses only cut k₁. But then |G| ≥ 2 from this one tile alone: it witnesses {k₁−1,k₁}... wait, (x=k₁, y=k₁−2) witnesses cuts {k₁−1, k₁} ∩ {1,...,n-1}. So k₁−1 is ALSO witnessed (if k₁ ≥ 2).

More carefully: for G(τ) = exactly {k₁,k₂}, we need no tile to witness any cut OUTSIDE {k₁,k₂}. In particular, tiles must not witness k₁−1 (if k₁ > 1) or k₂+1 (if k₂ < n-1). This forces all upward tiles to have y+1 ≥ k₁ (so y ≥ k₁−1) and x ≤ k₂ (so the tile doesn't reach k₂+1). Combined with x ≥ y+2 ≥ k₁+1 and x ≤ k₂: the free tile must satisfy k₁+1 ≤ x ≤ k₂ and y ≥ k₁−1 and y < x−1.

For k₂ = k₁+1: y ∈ {k₁−1} only (since y ≥ k₁−1 and y < x−1 ≤ k₂−1 = k₁, so y = k₁−1), and x ∈ {k₁+1} only (since x ≤ k₁+1 and x ≥ y+2 = k₁+1). So free tile = (k₁+1, k₁−1) uniquely. Making this tile upward witnesses both {k₁, k₁+1}. ✓ **Each consecutive pair {k,k+1} gives exactly 1 valid tiling.** There are n−2 pairs, so |G|=2 count = n−2. □

## Proof of Part 2: |G|=3 gives 5(n−3)

By the same argument: any upward tile not witnessing cuts outside {k,k+1,k+2} must have y ≥ k−1 and x ≤ k+2. The free tiles are:
- (k+1, k−1): witnesses {k, k+1}
- (k+2, k−1): witnesses {k, k+1, k+2}
- (k+2, k): witnesses {k+1, k+2}

Three free tiles, 2³=8 configurations. For exactly {k,k+1,k+2} good:
- Cut k requires tile (k+1,k−1) or (k+2,k−1) upward.
- Cut k+2 requires tile (k+2,k−1) or (k+2,k) upward.

Bad configs: [cut k bad] = 2 (tiles (k+1,k−1) and (k+2,k−1) both down, tile (k+2,k) free) + [cut k+2 bad] = 2 − [both bad] = 1. By IE: bad = 2+2−1 = 3. Valid = 8−3 = 5.

There are n−3 intervals {k,...,k+2} (k=1,...,n-3). Each contributes 5 tilings.

Non-consecutive triples contribute 0 tilings (for both n=5,6,7 verified). This is because for a triple with a gap, the "outer" cuts cannot be witnessed independently of cuts outside the triple. 

Total = 5(n−3). □ (VERIFIED for n=5..7; PROVED for n≥5.)

## Part 2: |G|=4 gives (n−4)(n+95)/2

For consecutive intervals {k,...,k+3}: free tiles = C(4,2)=6, valid configs = SC(5) = 50. Contributes (n−4)·50 tilings.

**Non-consecutive contribution:** Additionally, non-consecutive 4-element good sets contribute C(n−4,2) tilings. For n=6: C(2,2)=1 (from {1,2,4,5}); for n=7: C(3,2)=3; for n=8: C(4,2)=6, etc.

The lone non-consecutive tiling for good set {1,2,4,5} at n=6: only tiles (2,0) (crossing {1,2}) and (5,3) (crossing {4,5}) are free. Both must be upward. All other tiles are forced downward (they cross cut 3).

Total = 50(n−4) + C(n−4,2) = (n−4)(n+95)/2.

**Verification:** (n−4)(n+95)/2 = (n−4)(100+n-5)/2:
- n=6: 2·101/2 = 101 ✓
- n=7: 3·102/2 = 153 ✓
- n=8: 4·103/2 = 206 ✓

## General Bucket Recurrence (THM-348)

THM-348 upgrades the small-bucket arguments above to a full recurrence. With
`N=n-1` legal cuts, let

`B_N(x) = Σ_g #{tilings with |G|=g} x^g`.

Let `c_L` be the number of interval subsets covering a connected run of `L`
cuts, using intervals of length at least 2. Then

`B_N(x) = B_{N-1}(x) + Σ_{L=2}^N c_L x^L B_{N-L-1}(x)`.

This explains the formulas here uniformly:

- `|G|=2`: one run of length 2, so `n-2`.
- `|G|=3`: one run of length 3, so `5(n-3)`.
- `|G|=4`: one run of length 4 plus two separated runs of length 2, giving
  `50(n-4)+C(n-4,2)`.

See `01-canon/theorems/THM-349-good-cut-interval-union-recurrence.md`.

## Computational Verification

| n | |G|=0 | |G|=1 | |G|=2 | |G|=3 | |G|=4 | |G|=n-1 (SC) | Total |
|---|------|------|------|------|------|------------|-------|
| 3 | 1 | 0 | 1 | — | — | 1 | 2 |
| 4 | 1 | 0 | 2 | 5 | — | 5 | 8 |
| 5 | 1 | 0 | 3 | 10 | 50 | 50 | 64 |
| 6 | 1 | 0 | 4 | 15 | 101 | 903 | 1024 |
| 7 | 1 | 0 | 5 | 20 | 153 | 30773 | 32768 |
| 8 | 1 | 0 | 6 | 25 | 206 | 2032504 | 2097152 |

## Connection to SC Cut Theorem

By THM-330 plus the concrete bridge in `TournamentH7.StaircaseConnectivity`,
Lean now proves SC iff `|G|=n-1` for the tournament induced by a staircase
tiling. The hierarchy is:
- |G|=0: transitive (maximum hierarchy)
- |G|=1: impossible
- |G|≥2: partial connectivity
- |G|=n-1: strongly connected (maximum democracy)

The jump from 0 to ≥2 good cuts (skipping 1) reflects the binary tile structure: every tile spans ≥2 cuts simultaneously, preventing a "one-cut bridge."

## Lean Formalization

`TournamentH7.GoodCuts` and `TournamentH7.StaircaseConnectivity` formalize
the axiom-free core:

- `StTiling.goodCutCount_ne_one`: no tiling has exactly one good cut.
- `StTiling.goodCutCount_eq_zero_iff_all_down`: bucket 0 is exactly all-down.
- `StTiling.goodCutCount_bucket_bounds`: `g(τ) ∈ {0} ∪ {2,...,n-1}`.
- `StTiling.isGoodCut_iff_exists_upward_tile_interval`: good cuts are unions of upward tile intervals.
- `StTiling.goodCutCount_mono`: turning more tiles upward can only add good cuts.
- `StTiling.goodCutCount_eq_top_iff_all_cuts_good`: top bucket iff every legal cut is good.
- `StTiling.goodCutCount_reflect`: grid reflection preserves `g`.
- `StTiling.toTournament_hasBasePath`: concrete staircase tilings induce
  base-path tournaments.
- `StTiling.isGoodCut_iff_crossesUpward_toTournament`: good cuts translate to
  THM-330 crossing-upward cuts.
- `StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected`: top bucket
  iff the induced tournament is strongly connected.

The S15/Codex audits show these depend only on Lean foundations, not on
project-specific axioms.

## Strong-Component Interpretation (THM-354)

THM-354 identifies the whole hierarchy structurally. For any tournament `T`
and any Hamiltonian base path `P`,

`g_P(T) = n - #SCC(T)`.

Thus the missing good cuts are exactly the boundaries between consecutive
strong components in the condensation order. The top bucket is the one-component
case, the bottom bucket is the `n`-component transitive case, and the forbidden
level `g=1` is equivalent to the impossibility of a tournament condensation
with exactly `n-1` strong components.
