# THM-949 — The q-window lemma and the witness-ladder rigidity (death-star-2026-07-17-S41)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCQWindow.lean`; verify the
build report in the session log). Source: HYP-7190. Discharges both opens THM-947
named: the q-window and the constraint-system rigidity on trapped strata.

## Statement

1. `failWitness_pos_of_window` (**the q-window lemma**): for q ≤ 14·v_i, every band
   failure at p ≥ 1 has witness n ≥ 1 — the near-zero catastrophe cannot occur
   inside the window. `failWitness_le`: n ≤ v_i for p < q. **Witnesses live in
   [1, v] inside the window** — THM-947's constraints become near-proportionality
   of positive integers.
2. `witness_ladder` (**the rigidity engine**): simultaneously bad runners with
   3·v_i ≤ v_j ≤ 13·v_i and n_i ≥ 1 force **n_j ≥ 3·n_i** — the witness vector
   inherits the ladder. The 13 is sharp for the argument (slack
   (v_i+v_j)/(14·v_i) < 1 ⟺ v_j < 13·v_i; at 13 equality still lands
   n_j > 3n_i − 1) — **and 13 is the residual's own compression constant: the
   rigidity window IS the residual's ratio window.**
3. `witness_chain_geometric`: k+1 simultaneously bad ladder runners force
   n_top ≥ 3^k·n_bottom (List.Chain' induction); seven force n_top ≥ 729.
4. `witness_pin`: witness n ≥ N forces (14N−1)·q < 14·v·p — **ladder 7-overlaps
   are pinned far from zero** (p/q > 10205/(14·v_top)). Deep coverage requires a
   geometric witness vector the family's own magnitudes must accommodate inside
   p < q: the constraint-system rigidity in exact integers, no analysis.

## Referee

`qwindow_referee_deathstar_S41.out`: window/ceiling/ladder/pin PASS (40k samples
each); the 7-chain geometric bound confirmed on 18 REAL overlap instances found by
search (n_top ≥ 729·n_bot live, not just planted).

## Remaining after this theorem

The cap-supply composition: window + rigidity + trap ⟹ CoverageCapped(6) on
explicit (q, v)-strata via counting the pinned p's (each 7-overlap p confined to a
width-(q/(7·v_top)) window near n_top·q/v_top with n_top ≥ 729); mixed-ratio chains
through the dense-core certificate's per-step windows; the fleet composition
(kps exhaustion × codex lacunary × this rigidity).
