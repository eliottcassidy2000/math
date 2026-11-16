# Pre-registered predictions (klein-2026-11-16-hex)

Written to `prediction_hex.txt` at 2026-11-16 10:37 MST, before any run of this session.

| # | prediction | verdict |
|---|---|---|
| P-A | every rooted base (order ≤ 6 two added vertices; order 7 one added vertex) whose three digits share a common `(X′, g′)` has `H = 5` | **confirmed** (SINK_DIGITS.md: 3 + 2 bases, all of count 5) |
| P-B | at order 6, bases with `D ≠ 0` have a `3H + 2` extension at least twice as often as bases with `D = 0` | **refuted**: 83 % against 80 % (SINK_DIGITS.md) |
| P-C | with three added vertices on orders ≤ 5 the common-continuation bases still have `H ≤ 5` | **confirmed** (orders 4, 5: counts 5 only) |
| P-D | Collatz tree tournament, even-first order: `H(T_N) = 1` for `N ≤ 8`, `H(T_9) ≥ 3`; at `N = 20` even-first and odd-first differ, both odd, both > 100 | **confirmed**: `H(T_9) = 5`; 87,545,175 against 105,295,407 (COLLATZ_TOURNAMENT.md) |
| P-E | inflating the first 10 nodes into transitive 2-blocks multiplies `H` by a non-power of two | **confirmed**: 11 → 6,069 (factor 551.7) |

Not predicted: the generalised binary odometer search (PORT_MOVES.md), run after these.
