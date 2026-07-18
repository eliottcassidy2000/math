# THM-985 — THE TWO-CIRCLE THEOREM (death-star-2026-07-17-S53)

**Status:** PROVED IN FULL (Lean, kernel-pure — `TournamentH7/LRCTwoCircleII.lean`,
standard trio ×13; `decide` only on the finite congruence table). Source:
HYP-7270, completing HYP-7265/THM-979. The iff, both directions, in-kernel.

## The theorem

On the canonical family v = (1,…,13), for 0 < p < q:

**bandCount ≥ 6 ⟺ (84p < q ∨ 84(q−p) < q) ∨ 84|2p − q| < q**

— the deep set is EXACTLY the two resonance circles. (Recon: exact over 1185
moduli, every multiplier.)

## The proof architecture (all 13 pieces green)

- `witness_of_not_inBand`, `witness_range`, `witness_normalized` (minimality
  forces the witness into (0,a) coprime — divisor descent),
- `partner_lock`/`partner_branch` (THM-967/972 instances),
- **`compat_card_le`** — THE COLLAPSE: the six middle cases k₀ ∈ [3,8] reduce
  to one kernel `decide`: every coprime witness residue admits ≤ 4 compatible
  partners,
- `middle_case_bound` (compat injection: |F| ≤ 5), `large_case_bound`,
- `hub_case_circleI` (k₀ = 1: the hub lock nests; a failing speed ≥ 6 gives
  the 84-width), `even_case_circleII` (k₀ = 2: w₂ = 1 forced; parity locks
  kill odds; the 13-branch ray estimate 91|2p−q| > 6q kills 13; the six evens
  survive and speed 12 delivers the half circle),
- `bandCount_eq_fails_card` (index–speed bijection), `deep_implies_circles`,
  **`deep_iff_circles`**.

## Consequences

Exact canonical deep counts (2⌊(q−1)/84⌋ + parity form) are now one card
computation away; CoverageCapped(5) on {1..13} holds iff q avoids both
84-widths — decidable per q AND provable per family shape. The certificate
SHAPE generalizes: deep ⟺ union of m-circles over denominators with ≥ 6
near-multiples — the Wagner-circle program for arbitrary residual families
(named next; the S53 middle-case machinery is family-generic already).
