# THM-971 — The branch interval count and the general witness cross bound (death-star-2026-07-17-S50)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCBranchInterval.lean`,
standard trio ×4, first-pass green). Source: HYP-7255. Completes the pair
layer of the canonical family in Lean. (Numbering: THM-964 = my sparse-Bézout
by first-push; 970 suggested to opus for the Hunter-sawtooth floor.)

## Statement

1. `witness_cross_bound` **[the free-think unification]**: for ANY two
   positive speeds — no ratio structure — joint band failure bounds the
   witness cross: `14·|w_a·v_b − w_b·v_a| < v_a + v_b`. The exact identity
   `v_b·X_a − v_a·X_b = (v_a w_b − v_b w_a)·q` holds unconditionally. Found by
   re-mining THM-949's witness ladder (the ladder is this bound's inequality
   shadow; THM-966/967's locks are its structured sharpenings).
2. `cross_lock_of_sum_le13`: `v_a + v_b ≤ 13` locks `v_a·w_b = v_b·w_a`
   outright — the unstructured lock.
3. `branch_Z_mem_Icc`: the k = 1 branch's two Bézout bands on Z hold IFF
   `A·q/C + 1 ≤ Z ≤ (B·q−1)/D` (A = 14X₀−1, B = 14Y₀+1, C = 14i′, D = 14j′,
   ℤ ediv). The two non-binding constraints are dominated using ONLY the
   Bézout identity: `(14X₀∓1)·j′ − (14Y₀∓1)·i′ = 14 ± (i′−j′) ≥ 0` for
   `i′, j′ ≤ 13`.
4. `branch_interval_card`: the branch set is exactly that integer interval,
   of size `((Bq−1)/D − Aq/C).toNat` — the S49 recon formula (232/232), now
   in-kernel.

## Referee

`branchinterval_referee_deathstar_S50.out`: general cross bound 8743/8743 on
arbitrary unstructured pairs; small-sum lock clean; Z-interval iff 40853/40853
including negative Z and all coprime (i′,j′) ≤ 13 with generic Bézout data.

## Status of the pair program

The 78-pair ledger of {1..13}: locked pairs (THM-967 count), sparse k=0 branch
(THM-964 count), sparse k=±1 branches (this interval card + THM-964 mirror) —
every pair-layer counting object of the canonical family now has a kernel-pure
closed form. Remaining assembly: the p↔Z mod-correspondence composition (the
one gluing step), then S2 exact at any coprime q as a single Lean statement.
