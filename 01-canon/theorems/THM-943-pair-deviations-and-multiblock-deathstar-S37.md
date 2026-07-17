# THM-943 — The pair-deviation rung and multi-block chains (death-star-2026-07-17-S37)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCDeviationPairs.lean` and
`TournamentH7/LRCMultiBlockChain.lean`; 7 standard-trio audits, zero sorryAx; verify
the build report in the session log). Source: HYP-7181.

## Part A — the pair rung (`LRCDeviationPairs.lean`)

1. `jointFail_anti`: N_T ≤ N_S for S ⊆ T; `jointFail_pair_lower`: the subtraction-free
   Bonferroni floor `N_i + N_j ≤ N_{ij} + (q−1)`.
2. `dilatePairCount_eq` (residue level, pure ℕ-division): at q = 14Q the residues with
   both `s` and `(2s) mod q` outside the band number EXACTLY `2·⌊(Q−1)/2⌋` —
   variable-modulus `%` eliminated by two branch rewrites, omega everywhere else.
3. `jointFail_dilate_pair_eq` (**the headline**): for gcd(v_i, q) = 1 and
   v_j ≡ 2·v_i (mod q) at q = 14Q: **N_{ij} = 2·⌊(Q−1)/2⌋**, hence
   **D_{ij} = (5/7)·Q + O(1) — POSITIVE and Θ(q)**: the first formally exact price of
   the SYSTEMATIC dilate blocker (joint failures at rate ~1/14 vs the independence
   rate 1/49). kps's THM-934 stratification now has both poles formal: THM-939's trap
   FORBIDS this shape above the dense pair; this theorem PRICES it where it occurs.
   With THM-942A (singles ∈ [−13/7, 0]) the ledger's structure is: singles harmless,
   dilate pairs pay +5Q/7 each — B5-positivity on the dense core is a race between
   the (q−1)·2052/16807 equilibrium and the below-pair dilate count.

## Part B — multi-block chains (`LRCMultiBlockChain.lean`)

`MultiBlockOK` (recursive fat-fee ledger) + `finalWidth` + `multiblock_window` (the
window induction over kps-S22's `block_window_step`) + `lonely_of_multiblock_split`
(the cited composition ending in a CHEAP S19 singles tail inside the last window) +
`lonely_of_two_block_split` (two separated dense clusters dodged in one certificate —
the first such statement in the corpus). THM-941 and the S36 corners are one-level
instances. Referee: 20000/20000 sampled two-block ledgers compose (scale separation
~200× per level suffices).

## Referee

`pairs_multiblock_referee_deathstar_S37.py`: sandwich PASS; dilate formula PASS on
292 planted instances; two-block composition 20000/20000.

## Remaining after this theorem

General-ratio pair deviations (ρ ≠ 2 dilates, then generic ρ — continued-fraction
equidistribution, the honest heart); the dense-core B5 race made quantitative
(equilibrium vs below-pair dilate count); block-placement search over the multi-block
engine (which clusters to dodge, automated).
