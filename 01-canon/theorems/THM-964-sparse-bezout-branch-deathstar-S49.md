# THM-964 — The sparse Bézout-branch decomposition (death-star-2026-07-17-S49)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCSparseBranch.lean`,
standard trio ×6). Source: HYP-7245. Consumes THM-963. Closes the structural
layer of the 29 sparse pairs of {1,…,13}.

## Statement (Lean, ×6)

1. `witness_unique`: a failing runner's band witness is unique — so the Bézout
   residue `k = j′w_a − i′w_b` is a well-defined function of the jointly-failing
   multiplier, and the branches `k ∈ {−1,0,+1}` (THM-963) PARTITION the set.
2. `branch_zero_iff` / `branch_zero_count`: the k = 0 branch is exactly the
   GCD-speed narrow band — for EVERY reduced pair, locked or sparse (no
   i′+j′ ≤ 13 needed once k = 0 is imposed): count `2⌊(q−1)/(14·max(i′,j′))⌋`.
3. `branch_one_iff`: the k = 1 branch, with Bézout `j′X₀ − i′Y₀ = 1`, is
   exactly the two-band condition on `Z = g·p − t·q`.
4. `branch_no_qmultiple`: the branch interval never contains a multiple of q
   (it would force `j′X₀ − i′Y₀ = 0`) — no p ≡ 0 boundary correction.
5. `branch_mirror_card`: `p ↦ q−p` negates k (witness `w ↦ g·v′ − w`), so
   N⁺ = N⁻ by explicit bijection.

## The complete pair ledger of the canonical family (recon-exact)

With A = 14X₀−1, B = 14Y₀+1, C = 14i′, D = 14j′:
- **49 locked pairs** (i′+j′ ≤ 13): N = 2⌊(q−1)/(14·max)⌋ (THM-963, Lean).
- **29 sparse pairs** (14 ≤ i′+j′ ≤ 25): N = 2⌊(q−1)/(14·max)⌋ + 2·N⁺ with
  **N⁺ = ⌊(Bq−1)/D⌋ − ⌊Aq/C⌋** — verified 232/232 across all 29 pairs
  (recon (f)); the Lean interval-count of this floor form is the named next
  step (structure fully proved; the count is a lattice-interval evaluation).
- Limit law: N/(q−1) → 1/(7max) + 2(i′+j′−14)⁺/(14i′j′); reproduces boxeph
  LEM-044's μ exactly (pair (8,9): 77/3528 ✓; (7,8) → 1/49, the 7∣k law ✓).

## The honest B5 scope finding (recon B, `b5limit_recon_deathstar_S49.out`)

Assembled S2 limit on {1..13} = **462821/210210 ≈ 2.2017** — confirmed against
the census (2.1971 at q = 20011, converging ✓; equilibrium 78/49 ≈ 1.5918).
BUT **B5/(q−1) ≈ −9.6 < 0 and liveCount ≈ 0 at generic prime q**: the
canonical family is maximally resonant — its lonely instants live only at
resonant moduli (q = 14: p = 1 is live — the equality witness). Therefore the
B5-funnel live floor MUST carry a dissociation hypothesis; compressed
near-tight families belong to direct-witness/citation routes, not the funnel.
This properly scopes the nucleus: the a-priori live floor is a statement about
DISSOCIATED strata.
