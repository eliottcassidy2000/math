---
id: THM-697
title: The first-window witness — the packed top-block supply gap closed by a route SIMPLER than the test points: the multiplier w = 7j+6 puts the cluster phase at EXACTLY 6/7 (V·w mod 7V = 6V), so cluster residues are 6V − e·w with ONE inequality per speed (2ew < 11V), while small speeds ride the OPEN FIRST WINDOW W_P = (1/(14·pmin), 13/(14·pmax)) — nondegenerate for EVERY P with pmax + 1 ≤ 13·pmin, in particular ALL top-blocks [b,13] — no test point, no missed class; leftmost-j constructive (bracketed by 112·pmin); COMPLEMENTARITY: the degenerate-window P's (containing {1,13}) are exactly non-top-blocks covered by THM-693/696's test-point route — the two witnesses jointly cover every small part with packed clusters. FORMALIZED SORRY-FREE, GREEN ON FIRST COMPILE (LRCFirstWindowWitness.lean); demo: the k = 8 top-block {9..13} with the CONSECUTIVE APEX cluster strictly live at (70000, 559), and its certificate-supply data produced through the fattening bridge
status: PROVED AND FORMALIZED (cluster_band, firstWindow_strictlyLive, firstWindow_j_exists, demoTopBlock_strictlyLive, demoTopBlock_supply — all kernel-pure [propext, Classical.choice, Quot.sound], green on first compile, root-wired, project green). With THM-696: the certificate-supply citation's remaining content shrinks to [spread clusters at top-blocks (e ≳ 10·pmin — censused all-positive at S239/S240)] ∪ [non-taxonomy realization shapes (the coverage audit)] ∪ [bounded V (THM-686 windows/banks)].
source: klein-2026-07-11-S249 (HYP-5970; continuing S248's remainder (a) with a reroute to the simpler uniform route)
depends_on:
  - kps LRCStrictRuler (StrictlyLive), LRCTestDataSupply (the fattening bridge consumed by the demo)
related:
  - THM-691(B) (the packed sliver this makes constructive at finite V), THM-693/696 (the complementary test-point witnesses)
---

# THM-697 — the first-window witness

## The construction

At t = w/(7V) with w = 7j+6: the top phase is frac(V·t) = frac(w/7) = 6/7
EXACTLY, so every cluster speed V−e has residue ((V−e)·w) mod 7V = 6V − e·w
(one ring identity + the emod window), strictly in-band iff 2·e·w < 11·V.
Small speeds p need only the no-wrap band 7V < 14·p·w < 91V, which holds
throughout the OPEN first window W_P = (1/(14·pmin), 13/(14·pmax)) — an
interval of length ≥ (13·pmin − pmax)/(14·pmin·pmax), nondegenerate whenever
pmax + 1 ≤ 13·pmin: for top-blocks [b,13] the room is 13(b−1)/(14·13b) —
enormous. The leftmost w ≡ 6 (mod 7) past the window's left edge is
constructed by ediv + mod-7 rounding, bracketed within 112·pmin of the edge,
keeping e·w ≤ ~e·V/(2·pmin): the packed range e ≤ 10·pmin.

## Why this closes S248's gap (a)

The test-point route (THM-693/696) requires a modulus q* ∈ [8,13] outside P
— exactly what top-blocks [b,13] deny. The first window is the OPPOSITE
resource: it is degenerate exactly at ratio-13 P's (pmin = 1 with pmax = 13)
— which always admit q* (they are never top-blocks). THE TWO WITNESSES ARE
COMPLEMENTARY: every small part with a packed cluster is covered by one of
them, and both feed the fattening bridge (THM-696) into mac-mini's
THM527ACertificateSupply.

## The supply's remaining content after this theorem

[spread clusters at top-blocks, e ≳ 10·pmin — the S239/S240-censused sliver,
per-class decidable] ∪ [realization shapes outside the two-/multi-scale/ray
taxonomy — the shape-coverage audit, the critical open mapping] ∪
[bounded V — THM-686 windows and the banks]. Plus the other two citations:
the ≤7-arcs pigeonhole and THM-661's interior.

## Formalization & files

`04-computation/lean/TournamentH7/TournamentH7/LRCFirstWindowWitness.lean`
(kernel-pure, green on first compile, root-wired). Demo cross-check:
`05-knowledge/results/lrc14_first_window_witness_demo_klein_S249.out` — the
top-block apex family {9,…,13} ∪ {V−7,…,V} at V = 10000, all thirteen
residues strictly in-band at (70000, 559), cluster phase exactly 6/7.
