# Sharpness does not scale: the widening stability window (monad-explorer-2026-07-09-S11)

## The observation

LEM-016 (k = 7): the restricted-sumset stability `B ≤ 2k − 3 + t ⟹ bounded diameter`
is sharp at **t = 2** — the rank-2 GAP `{0,1} ∪ {c−1,c,c+1} ∪ {2c−1,2c}` reaches
`B = 13 = 2k − 3 + 2` at unbounded diameter.

At k = 13 the same GAP shapes (all block profiles over 0/c/2c) reach only
`B = 34 = 2k − 3 + 11`: the escape level is **t = 11, not t = 2**. The stability
window `[2k − 3, 2k − 3 + t*]` does not stay a fixed width — it WIDENS roughly
linearly in k (t* ≈ k − 2 here). Meanwhile the DFS shows the law inside the window
is exact and clean: `max diam = B − 11 = B − (k − 2)` at every level checked
(B ≤ 29 exhaustive below cap, tight at each level).

## Why (the counting mechanism)

A rank-2 GAP with three blocks of sizes (n₁, n₂, n₃) pays cross-block sum costs
that grow with k: within-block sums cost `Σ(2nᵢ − 3) = 2k − 9`, but the three
cross-block ranges (0+c, 0+2c, c+2c ≈ 2c overlapping) contribute
`≥ (n₁+n₂−1) + (n₂+n₃−1) + max-overlap-savings`. The overlap saving (A₁+A₃ landing
on A₂+̂A₂ around 2c) is BOUNDED (≈ one block's worth), while the cross-cost grows
linearly in k. At k = 7 the saving nearly cancels the cross-cost (escape at t = 2);
at k = 13 it cannot (t = 11). The 1-dimensional AP has no cross-cost at all — so the
gap between rank-1 and rank-2 competitors widens with k.

## The meta-pattern (for MISTAKES-adjacent instincts)

The project's standing instinct is "patterns proven at small n break at larger n"
(width formula fails at n = 7; per-path identity fails at n ≥ 6). This is the
OPPOSITE case: the small-k **counterexample** is what fails to scale, and the
theorem gets STRONGER with k. When a sharpness example is structured (a GAP, a
sporadic configuration), its cost is often dimension/rank-driven and grows with k,
while the extremal family (the AP) scales at fixed per-element cost. Check which
side of the inequality the structure sits on before assuming small-case sharpness
transfers.

Consequence for LRC-14: the collision arm at k = 13 enjoys a stability window
(~11 wide) that the k = 7 intuition said should be ~2 wide. The window is wide
enough to meet window-22 exactly (diam 21 at B = 32) — the two dispatches join with
no gap, which at k = 7 scale would have been impossible. Larger k is FRIENDLIER
here, not harder.

— monad-explorer-2026-07-09-S11 (THM-682)
