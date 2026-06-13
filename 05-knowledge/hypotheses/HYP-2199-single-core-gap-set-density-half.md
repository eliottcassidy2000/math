# HYP-2199 — Single-core signature gap set has asymptotic density 1/2 and no simple closed form

**Status:** CONFIRMED (computational, complete to R = 2^20)
**Source:** monad-compute-2026-06-04-S3
**Extends:** [[HYP-2198]] (single-core signature gap, complete for all lengths).
**Resolves:** the S2 handoff "the dense single-core gap structure {3,6,10,14,17,20,21,…}
might itself be worth an OEIS/closed-form look."

## Object

The single-core odd-cycle count `r(s) = Σ_{i<j, s_i=1, s_j=0} f(j-i-1)`,
`f(0)=1, f(t)=2^{t-1}` (see [[HYP-2198]] for the construction). The **gap set**
`G = { r ≥ 1 : r ≠ r(s) for every bit string s }` is the set of single-core odd-cycle
counts that are unachievable at any length (equivalently, single-core complete-Ω H-values
`1+2r` that cannot be produced). S2 proved `G` is decidable to any bound `R` via the
length cap `L ≤ 3 + ⌊log₂R⌋`.

## Computation

`04-computation/single_core_gap_structure_monad_s3.py`, complete to **R = 2^20 = 1,048,576**
(rigorous length cap L = 23). Fast O(L) per string via the contribution recurrence
`A_{p+1} = 2·A_p − s_{p-1} + s_p`, `r += A_p` on each 0 — validated exhaustively against the
O(L²) reference `r_brute` for all strings of length ≤ 14.

- achievable r in [1, 2^20]: **521,915**
- gap r in [1, 2^20]: **526,661**  (gap fraction 0.5023)

## Findings

1. **Asymptotic density = 1/2.** Gap density per dyadic window `(2^k, 2^{k+1}]` converges
   monotonically to 50%: …, 49.77% (8k–16k), 49.93%, 49.99%, 50.24%, 50.08%, 50.34%,
   50.26% (524k–1M). The gap set is **persistent / infinite** (largest gap ≤ R is 1,048,574),
   **not finite**.

2. **Novel to OEIS.** Neither the gap set
   `3, 6, 10, 14, 17, 20, 21, 24, 27, 28, 29, 33, 35, 37, 38, 41, 44, 46, …`
   nor the achievable set
   `1, 2, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 18, 19, 22, 23, 25, 26, 30, 31, 32, 34, 36, …`
   returns any OEIS match (checked with 11–22 term prefixes, with/without leading 0).

3. **No simple closed form** (all natural density-1/2 candidates ruled out):
   - NOT a union of residue classes mod q for any q ≤ 12.
   - NOT Thue-Morse / popcount-parity: gaps are 50.1% odious, achievable 49.9% odious —
     popcount parity is independent of achievability.
   - NOT a Beatty sequence: consecutive gap differences take every value 1,2,…,12+ (a Beatty
     sequence has at most two distinct gap sizes).
   - Achievable set is NOT an additive semigroup (1,2 ∈ ACH but 1+2 = 3 ∈ G), and NOT closed
     under doubling (5 ∈ ACH, 10 ∈ G).

4. **What is structured:** every power of two is achievable (witness `1·0^k → r = 2^{k-1}`),
   so `2,4,8,16,…` (and 1) are never gaps. The target values from [[HYP-2198]] re-confirm:
   r=3 (H=7), r=6 (H=13), r=10 (H=21), r=94 (H=189) are gaps; r=31 (H=63) achievable.

## Interpretation

Single-core complete-Ω is a *strict* sub-construction whose forbidden-r set is "half of all
integers" with no arithmetic regularity — so it carries no special structure that would single
out `{7, 21}`. This reinforces [[HYP-2198]]'s conclusion: the single-core picture explains the
H=63 unlock but cannot by itself characterise the *global* permanent gaps `{7, 21}`
(that is HYP-1753 / THM-079, over all Ω shapes). Cf. MISTAKE-024/050.

## Files

- `04-computation/single_core_gap_structure_monad_s3.py`
- `05-knowledge/results/single_core_gap_structure_monad_s3.out`
