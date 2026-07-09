# THM-666 — The pair-sum ruler theorem: every local maximum of the reach function lives at t = p/(v_i+v_j)

**Status:** CLAIMED (stub, proof + verification landing this session)
**Source:** mac-mini-2026-07-09-S65
**Depends on:** none (elementary). Precedents: THM-420 (empirical pair-sum witnesses on the shell-free residual), mac-mini-S64 / klein-S206 (used as an exact-computation device, unproved).

## Claim (to be filled with proof this session)

Let `S = {v_1 < … < v_k} ⊂ Z_{>0}`, `m(t) = min_i ‖v_i t‖`, `M(S) = max_t m(t)`.

1. Every local maximum of `m` occurs at a rational `t = p/(v_i+v_j)` for some `1 ≤ i ≤ j ≤ k`
   (the case `i = j` giving single-runner peaks `p/(2v_i)`, `p` odd).
2. Hence `M(S)` is attained at a **pair-sum ruler** rational with denominator `≤ 2·Vmax`,
   and `M(S) ∈ Q` with denominator dividing some `v_i + v_j`.
3. Realization consequence: pair-(i,j) center events have spacing `1/(v_i+v_j)`; a component of
   the wide-gap set on which teeth `i,j` bound a `>1/7` arc, of length `≥ 2/(v_i+v_j)`, contains
   a properly-oriented event and hence a lonely time.

Verification script: `04-computation/lrc14_pair_sum_ruler_macmini_S65.py` (landing).
