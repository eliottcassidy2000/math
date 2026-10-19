---
id: HYP-9335
title: "The Farey / mediant / Mandelbrot dictionary across Collatz, Hadwiger-Nelson, the planar Jacobian conjecture, Goldbach and twin primes: (a) Collatz cycle shapes are Stern-Brocot nodes of log_3 2 and the sporadic ones sit at semiconvergents; (b) the Goldbach / twin-prime singular series truncated by Stern-Brocot depth does NOT converge (REFUTED guess: the circle-method cancellation is per denominator, Ramanujan sums, not per mediant depth); (c) the Farey graph is 3-coloured by the parity class of (a, b) and the rational unit-distance graph is bipartite (Woodall), so the Farey skeleton contributes nothing to chi(R^2) >= 5; (d) the Bass-Serre tree of Aut(C^2) = A *_(A cap E) E is the Jacobian analogue of the Farey tree (dual to the Farey tessellation for PSL2(Z)), and HYP-9070's Euclidean depth is the Stern-Brocot depth of the degree ratio; (e) codex's parabolic clocks at the cusp c = 1/4 are the Mandelbrot side of the Collatz odd anchor -1"
status: >
  OPEN as a dictionary; each row typed. (a) OPEN: FINITE-EXACT placements in THM-4715 (3) and the CITED Legendre
  mechanism of THM-4633; the conjecture is that every sporadic cycle of 3x+d with minimum m > (d/(3 ln 2)) K^2 is at a
  convergent, and those with smaller minimum at semiconvergents or off the path. (b) REFUTED as a convergence statement by
  singular_series_farey.py: for N = 2 the depth-16 truncation gives 0.885 against 1.320 (denominators <= 30 give 1.313);
  for N = 210: 4.205 at depth 16 against 4.225 (denominators <= 100: 4.238). (c) KNOWN (Woodall 1973; the three cusps of
  Gamma(2)). (d) KNOWN structure (Jung-van der Kulk; Serre); the identification with HYP-9070's search order is a
  reading, not a theorem. (e) a reading of codex's parabolic_scaling_session_20261019 through THM-4714 (1b,c).
session: mac-mini-2026-10-19-farey
source: 05-knowledge/results/farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md
related:
  - THM-4714, THM-4715, HYP-9334, THM-4633, THM-998 (the LRC deep set is a Farey major-arc dissection), THM-501/503 (the LRC singular series), HYP-9070
  - 04-computation/experiments/farey_mandelbrot_20261019/singular_series_farey.py (+ .out)
---

# HYP-9335 — the Farey / mediant / Mandelbrot dictionary

| conjecture | Farey / mediant object | Mandelbrot object | status |
|---|---|---|---|
| Collatz | parity number `ρ` (dyadic) ↔ Farey address `?^(−1)ρ` (continued fraction of run lengths); cycle shapes `(K,X)` on the Stern–Brocot path of `log_3 2`; concatenation = mediant of slopes | external angles; cardioid bulbs = balanced cycles; cusp `1/4` = anchors `0, −1`; tip `−2` = `−1/3`; `1/a`-bulbs = anchors `c_a` | THM-4714/4715 PROVED; HYP-9334 OPEN; (a) OPEN |
| Goldbach, twin primes | circle method = Farey dissection; singular series = sum over Farey fractions of Ramanujan sums | bulb sizes `~ q^(−2)` like Ford circles (numerology) | (b) REFUTED: mediant-depth truncation diverges from the series; the conductor filtration is the right one |
| Hadwiger–Nelson | rational unit vectors = Farey parameters (THM-4057, THM-3333); the Farey graph is 3-coloured by `(a mod 2, b mod 2)` | none found | (c) KNOWN: `Q^2` is bipartite (Woodall); irrational directions are needed for `χ ≥ 5` |
| planar Jacobian | Jung–van der Kulk amalgam, Bass–Serre tree (the Farey tree is PSL2(Z)'s); Euclidean depth of `(deg f, deg g)` | parabolic cotangent germs and `π/√ε` clocks (codex) are the same cusp as THM-4714's anchor | (d), (e) readings |

## Open questions worth a session each

1. Prove HYP-9334 (an identity for the Christoffel carry modulo `2^q − 3^p`).
2. (a): quantify, with THM-4633's budget, which Stern–Brocot nodes can host a sporadic cycle of `3x + d` with minimum `m`; the census to `d ≤ 49` suggests only convergents and first semiconvergents.
3. (b): which ordering of Farey fractions preserves the Ramanujan-sum cancellation while refining the denominator order (e.g. by Ford-circle size `1/(2q^2)` = by denominator; the mediant tree is incompatible).
