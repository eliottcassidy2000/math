# THM-950 — The pinned-p counting and the census criterion (death-star-2026-07-17-S42)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCDeepCount.lean`; verify the
build report in the session log). Source: HYP-7196. The composition THM-949 named:
CoverageCapped-style control on explicit strata by counting.

## Statement

1. `window_unique_p`: for q ≤ 7·v, two bad multipliers sharing a witness coincide —
   each witness value admits AT MOST ONE bad p (14·v·|Δp| < 2q ≤ 14v ⟹ Δp = 0).
2. `deepSet_card_le` (**the bottom injection**): for any runner set S containing a
   windowed runner i, `#{p : S ⊆ bad(p)} ≤ v_i` — deep sets inject into witness
   ranges. With THM-949's tower determinism the whole ladder-deep set is checkable
   per bottom witness: **CoverageCapped-style statements on explicit strata are
   finite checks.**
3. `bled_ge_neg792` / `bled_eq_low` (decide): the uncapped pointwise ledger is 1 at
   c = 0, vanishes for 1 ≤ c ≤ 5, and never drops below −792 = −C(12,5).
4. `B5_ge_live_sub_deep` (**the census criterion**): UNCONDITIONALLY,
       **B5 ≥ liveCount − 792·#{p : bandCount ≥ 6}**;
   `B5_pos_of_live_beats_deep`: live > 792·deep ⟹ B5 > 0. On recon-typical strata
   the deep count is ZERO (4000 planted ladder-7 families at q ≤ 14·v_bot: 99.5%
   zero, max 3 — only at resonant q), so the criterion is simply liveCount > 0.

## The closing shape

The race is now: count live multipliers (a census), bound deep multipliers (the
bottom injection + tower determinism per stratum — finite), win by 792:1. No caps,
no moments, no analysis in the criterion itself.

## Referee

`deepcount_referee_deathstar_S42.out`: window uniqueness PASS (30k runner/q pairs,
zero witness collisions); census floor PASS (60 full families vs direct B5). Recon
distribution: {0: 3979, 1: 17, 2: 3, 3: 1} deep p's over 4000 ladder-7 plantings.
