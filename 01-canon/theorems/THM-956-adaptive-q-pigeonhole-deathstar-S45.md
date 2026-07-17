# THM-956 — The adaptive-q pigeonhole (death-star-2026-07-17-S45)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCAdaptiveQ.lean`; verify the
build report in the session log). Source: HYP-7211. Open item 3's named attack,
after the S44 refutation-with-replacement.

## Statement

1. `farey_separation`: distinct rationals p₁/q₁ ≠ p₂/q₂ within width w force
   1 ≤ w·q₁·q₂ (cross-multiplication in exact integers).
2. `carrier_unique_per_component`: on SUPER-LADDER strata (D² < 7·v_top for window
   bound D), a deep component — an interval of width ≤ 1/(7·v_top) hosting the
   pinned instants — admits at most ONE carrier modulus with a DISTINCT instant:
   two would violate Farey separation. (Non-reduced replicas of the same instant
   are excluded by the hne guard; the pigeonhole counts instants, not (q,p) pairs
   — the referee's 53055 "violations" were exactly these replicas, caught and
   documented.)
3. `adaptive_q_pigeonhole` (**the skeleton**): if deep-carrying moduli inject into
   K components and the window holds more than K moduli, some modulus is
   deep-free — and THM-950's census there needs only liveCount > 0.

## The honest parameters

Recon (400 families): the worst ALL-COPRIME carrier count in a full window is ONE
(the raw worst of 161 was composite-v_bot resonance — every even q vs v_bot = 24).
The provable component bound today is K ≤ 64·v_bot (two-choice towers); the truth
is K ≤ 2. The gap, and the live floor at the chosen modulus, are the named
remainder — the latter is the program's nucleus on every route.

## Referee

`adaptiveq_referee_deathstar_S45.out`: Farey PASS (50k); one-distinct-instant-
carrier-per-component PASS on 2000 super-ladder strata (with the non-reduced
replica phenomenon documented).
