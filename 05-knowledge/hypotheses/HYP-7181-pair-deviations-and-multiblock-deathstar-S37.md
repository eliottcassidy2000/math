# HYP-7181 — The pair-deviation rung + multi-block chains (death-star-2026-07-17-S37)

**Status:** CLAIMED — in progress (S37). Verify the build report in the session log.

**Owner directive:** "prove the pair deviations D_{ij} rung and run the multi-block chains."

**Claim A (`LRCDeviationPairs.lean`):** the |T| = 2 rung of THM-940's ledger:
(1) jointFail monotonicity (N_T ≤ N_S for S ⊆ T) + the Bonferroni sandwich
`max(0, N_i + N_j − (q−1)) ≤ N_{ij} ≤ min(N_i, N_j)`;
(2) **the exact dilate-pair count** (recon locked, referee-verified): at q = 14Q with
gcd(v_i, q) = 1 and v_j ≡ 2·v_i (mod q): **N_{ij} = 2·⌊(Q−1)/2⌋** — hence
D_{ij} = (5/7)Q + O(1), POSITIVE and Θ(q): the first formally exact quantification of
kps's "systematic dilate blocker" (the joint failures run at ~q/14 against the
equilibrium ~q/49). Residue-level count is pure ℕ-division (omega-friendly); the
bijection transfer follows the THM-942A pattern with the congruence
(v_j·p) % q = (2s) % q tracked alongside.

**Claim B (`LRCMultiBlockChain.lean`):** the multi-block generalization of THM-941's
generic engine: `MultiBlockOK` (a list of (block, ε) levels, each paying the fat-mass
fee at the current window and recursing at its ε) + the window induction over
`block_window_step` + the cited composition ending in a CHEAP S19 singles tail inside
the last window. Single-block (THM-941) and the closed corners become instances; a
two-block corollary demonstrates composition.
