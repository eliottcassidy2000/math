# Hadwiger–Nelson type-M experiments (mac-mini-2026-10-14-wilson)

Exploratory and not canon. Goal: small 5-chromatic unit-distance graphs. Current records:

* **Overall:** 509 vertices (Parts, arXiv:2010.12665).
* **Moser-spindle-free:** 1,441 vertices (Heule, Geombinatorics 31(2), 2021).

Every graph here has exact coordinates in a real multiquadratic field, with `x, y ∈ Q(√p : p ∈ HN_PRIMES)`. Every unit distance is decided by exact arithmetic. Colourability is decided by CaDiCaL 1.9.5 through pysat. UNSAT answers are solver verdicts, not DRAT-checked certificates.

## Toolkit

| file | role |
|---|---|
| `hn_udg.py` | Field arithmetic (`HN_PRIMES`, default `3,5,11`), Minkowski sums, rotations, exact unit-distance graphs, a selector-based colouring oracle, core shrinking and vertex-critical deletion |
| `type_m_probe.py` | Parts' type M, `G = L727 ∪ w·S361`, for `w ∈ {ρ, w_7/3, η²}` |
| `rotation_probe.py` | `w_t` for a squared reference radius `t`: `cos = 1 − 1/(2t)`, `sin = √(4t−1)/(2t)` |
| `rotation_scan.py` | 25 base pairs `A ∪ w_t B` for each `t` |
| `forced_radii.py`, `forced_radii_scan.py` | Radii `t` at which every 4-colouring repeats the centre's colour |
| `mono_sets.py` | Every realizable set of centre-coloured points on a reference orbit |
| `minimize.py` | Randomized core-guided deletion until vertex-critical; writes `best_*.json` with exact coordinates |

## Findings so far (2026-10-14)

1. **Reproduction.**
   * Parts' pieces have exactly his sizes: `V25, V31, V37, V49`, `L727 = V49 ⊕ V37`, and `S361 = V31 ⊕ V25`. S361 reaches 361 only with his second definition, `V31 = H¹ ∪ i√3·H^{−1,1}`.
   * `L727 ∪ ρ·S361`, with `ρ = (7 + i√15)/8`, has 1,087 vertices and 5,892 edges. It is not 4-colourable (CaDiCaL, 13 s).
   * One greedy deletion pass gave 599 vertices, not yet vertex-critical. The multi-seed run is in `minimize_L727_S361_rho.log`.
2. **Forced radii** (`forced_radii_scan.out`).
   * `t = 4` (radius 2) is forced in nearly every base graph with V31, V37 or V49 summands.
   * `t = 4/9` and `t = 20/9` are forced only in `V49 ⊕ V37` and `V49 ⊕ V49`. Their rotation fields are `Q(√−7)` and `Q(√−71)`.
   * `t = 7/3` is forced nowhere. Its rotation `(11 + 5i√3)/14` would stay in the Moser field, but `L727 ∪ w_7/3 S361` is 4-colourable.
3. **Mono-set rigidity** (`mono_sets_*.out`). At `t = 4`, L727 admits only 19 realizable mono-sets: six of size 6 (three antipodal pairs each), twelve of size 12, and the full orbit. S361 admits 19,916, the smallest being 6 antipodal pairs at the non-hexagonal angles. No single rotation joins every pair of L727's minimal mono-sets by radius-2 unit chords, so the type-M contradiction also uses the other cross edges.
4. **Rotation scan at `t = 4`** (`rotation_scan_1.out`). `A ∪ ρB` is 5-chromatic exactly when one piece is rigid (`V49 ⊕ V37` or `V49 ⊕ V49`) and the other forces radius 2. Two pieces that merely force radius 2 (`V37 ⊕ V37` with `V31 ⊕ V25`) give a 4-colourable graph.
