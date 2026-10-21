---
id: HYP-9346
title: "After THM-4731: (1) REFUTED before publication: the odd-n rigidity of self-line endpoints (true at n = 5, 7, 9) fails at n = 11, where eight endpoints have |Aut| = 3; (2) the line metagraph and the blue graph on the self-converse classes are connected for every n (HYP-9345 (b), re-observed here) and all non-pure-blue classes lie in one black component (true n = 4..8); (3) the self-line sequences |Q_n| = 2, 8, 16, 88, 412, 3408, 22376, ..., |Q'_n| = 2, 8, 40, 176, 852, 14864, 103264, ... and the merged black self-loops 0, 4, 24, 114, 618, 8998, 62666, ... have no closed form in sight (at n = 11 the black self-lines 116096 fall below SC(11)/2 = 139984); the blue self-lines at even n are 1, 2, 4, 24 (n = 4, 6, 8, 10; HYP-9345 (c)) and n = 12 is the next test; (4) self-lines concentrate on near-regular classes: the regular multiset carries none at odd n (leg law), and the complement-symmetric multiset nearest to regular with s_0 != s_(n-1) carries the largest fraction of self-line classes (4/5, 6/15, 43/85 at n = 6, 7, 8) — conjecture: for every n"
status: >
  (1) REFUTED at n = 11 (FINITE-EXACT: two independent implementations; kept as the record of a three-point
  mirage). (2) NUMERICAL (n = 4..8 exact; the first two statements are HYP-9345 (b)).
  (3) data, no conjecture; the n = 12 self-line run (12! systems, compiled) is the next test. (4) NUMERICAL (n = 6..8).
session: opus-2026-10-20-S46
source: 05-knowledge/results/tiling_line_metagraph_self_lines_twins_depth_20261020.md
scripts:
  - 04-computation/experiments/self_lines_twisted_burnside_20261020.py, self_lines_structure_20261020.py, tiling_line_metagraph_census_20261020.cpp, tiling_line_metagraph_analysis_20261020.py
related:
  - THM-4731, THM-643, THM-648, THM-649, THM-793, THM-849, THM-2454, HYP-9345 (mac-mini, reserved, concurrent)
---

# HYP-9346 — what the self-line engine and the n = 8 census leave open

## Statements

1. **(Odd-n rigidity; REFUTED.)** At `n = 5, 7, 9` every endpoint `t in Q_n` has `|Aut(T_t)| = 1` (and is black, THM-648), so `|Q_n| = P_kappa(n) = sum_g 2^(dim S_g)` there; at `n = 11` eight of the `232192` endpoints have `|Aut| = 3` (two Klein orbits with scores `(3^3, 4^2, 5, 6^2, 7^3)` and `(3^4, 4, 5, 6, 7^4)`, found by the session's engine and by the audit's independent solver, VF2-verified), so `|Q_11| = 232192 < P_kappa(11) = 232208`. Even `n` has endpoints with `|Aut| = 3` (n = 6, 8, 10) and `5` (n = 10). The permutations realising self-lines at `n = 5, 7` have cycle types `(3,1,1), (5)` and `(4,3), (5,2), (6,1)` (`self_lines_structure_20261020.py`).

2. **(Connectivity; NUMERICAL.)** For every `n`: the line metagraph (all lines, classes as nodes) is connected and the blue line graph on the self-converse classes is connected (HYP-9345 (b), from THM-4729's census; re-observed here from an independent census), and all non-pure-blue classes lie in one black component (the pure-blue classes are black-isolated by THM-643 T2, so the black graph has `1 + pureblue(n)` components). True at `n = 4..8` (components `1 / 1 / 2, 4, 3, 5, 4`). *Update (opus-2026-10-20-S47, THM-4732): true at `n = 9` (191,536 classes; black components `7 = 1 + 6`), and at `n = 8, 9, 10` every closed blue line between non-pure-blue classes has a literal three-black-move detour (HYP-9347).*

3. **(The self-line sequences; data.)** `|Q_n| = 2, 8, 16, 88, 412, 3408, 22376, 232192` (n = 4..11), `|Q'_n| = 2, 8, 40, 176, 852, 14864, 103264, 1009360`, converse-merged black self-loops `0, 4, 24, 114, 618, 8998, 62666, 620472`; ratios `|Q_n|/SC(n) = 1, 1, 4/3, 1, 2.34, 1.24, 2.55, 0.83`. Blue self-lines at even `n`: `1, 2, 4, 24` (n = 4, 6, 8, 10; THM-649 killed `2^(n/2-2)` at `n = 10`); `n = 12` is the next test, reachable by a compiled version of the solver (`12!` systems). No closed form is proposed.

4. **(Near-regular concentration; NUMERICAL.)** Self-line classes have complement-symmetric score multisets (THM-4731 (4)); the regular multiset carries none at odd `n` (the leg law forces `s_0 != s_(n-1)`), and among the realised multisets the one nearest to regular with `s_0 != s_(n-1)` carries the largest fraction of self-line classes (`4/5` at n = 6 for `(2^3, 3^3)`, `6/15` at n = 7 for `(2, 3^5, 4)`, `43/85` at n = 8 for `(3^4, 4^4)`). Conjecture: this holds for every `n`, and that fraction does not tend to `0`.

## Cheapest decisive tests

- Statement 1 is settled negatively; the compiled solver at `n = 12, 13` remains the test for statement 3 (blue self-lines at `n = 12`, HYP-9345 (c3)).
- Statement 2 at `n = 9` needs the class-level census (`2^28` tilings with canonical forms); the blue part is cheaper: the `2^16` blue tilings at `n = 9` and their classes (the pure-blue engines of THM-2444/2454 already touch them).
- Statement 3: `n = 12` for the blue self-line sequence.
- Statement 4 at `n = 9` from the same blue/black census.

## Why

The owner asked for comprehensive patterns and forced arrangements of the merged metagraph as `n` grows. THM-4731 settles the parity-level structure of the black matching (twins, `2B`, handshake), computes the self-lines past the enumeration wall and refutes the last outstanding self-line law (its own rigidity guess died at `n = 11` the same day); what remains are the global shape statements above, which the parities do not force.
