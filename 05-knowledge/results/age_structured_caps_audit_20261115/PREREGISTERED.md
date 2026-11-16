# Pre-registered predictions (klein-2026-11-15-age)

Written to `prediction_caps.txt` at 2026-11-15 21:43 MST, before `cone_test.py` was run and before the length-9 and
length-10 cap files were read. (The maximiser-structure predictions P-A/P-B of the sibling packet
`chanpak_span3_growth_two_apex_seam_20261115` are reproduced in `prediction_maximisers.txt`; both confirmed there.)

| # | prediction | verdict |
|---|---|---|
| P1 | at history length `K = 0` the coordinate-cap cone reproduces the earlier failure (hundreds of violating pairs, excess of order tens) | **confirmed**: 598 violations in round 1 (138 features), largest excess 256 |
| P2 | violations decrease as `K` grows from 0 to 3 but stay positive at every `K ≤ 3` | **refuted in its first clause, confirmed in its second**: 598 (K = 0), 5,328 (K = 1), 42,806 (K = 2), 341,961 (K = 3) violations in round 1 — they grow with the number of keys; none vanishes |
| P3 | value iteration does not converge within ten rounds for `K ≤ 2` | **confirmed**: K = 0 stalls at 550 violations on dead directions (infinite relative caps); K = 1 diverges (largest `Z` cap 38 → 49,325); K = 2 diverges in three rounds (largest `Z` cap 38 → 8,375) |
| P4 | the per-state maxima over words of length ≤ 9 equal those over length ≤ 7 for every state, and the largest is 14 | **first clause refuted, second confirmed**: states 81, 101, 114 grew (13/9 → 25/17, 8/5 → 44/27, 7/5 → 25/17); largest 14 (states 30, 99) |
| P5 | the deletion bound `κ(s)` equals the observed maximum for most two-fragment states and exceeds it for some three- and four-fragment states | **refuted**: κ is sharp only for the 16 trivial states (14 complete, two singleton pairs/triples); it exceeds the observed maximum for 102 states, and is below it for none |
| P6 (written 2026-11-15 22:20 MST, before the length-11 run `hratio5_11_*` was read) | at length 11 at least one of the three growing caps (states 81, 101, 114) moves again, and state 81 reaches 49/33 (the next term of (3·2^k + 1)/(2^(k+1) + 1)) | PENDING (run launched at close-out; verdict to be recorded by the next session) |

Computed after the pre-registration and not predicted: the per-letter maxima to length 9 (1, 9, 17, 19, 33, 33, 39, 31,
with `Z` one less), the value `max κ = 30`, the length-10 continuation of the growing caps, and the Zsigmondy reading of
the covering wheel.
