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
| `smus_probe.py` | Smallest non-k-colourable induced subgraph via pysat's OptUx (implicit hitting sets) |
| `spindle_free_build.py` | Exact rebuild of the Voronov–Neopryatnaya–Dergachev spindle-free graph (series 2, case 10) in `Q(i, √2, √3, √5)`, plus its 4-colourability check |
| `minimize_npz.py` | Core-guided minimisation of a large saved graph (phase 1: assumption cores; phase 2: deletion until vertex-critical) |
| `verify_graph.py` | Independent re-check of a candidate: all-pairs exact unit distances, two solvers, vertex-criticality |

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
5. **Other rotations** (`rotation_scan_2.out`).
   * 225 graphs (`t = 4/9, 11/9, 4/3, 5/3, 20/9, 7/3, 3, 5, 7`, 25 base pairs each) are all 4-colourable.
   * Among these hinges only Parts' `ρ` works, even where both pieces force the radius (`t = 4/9`).
6. **Exact smallest MUS.** OptUx gives V25's smallest non-3-colourable subgraph (9 vertices) instantly, but it did not finish `H1 ⊕ H1` within 10 minutes. It is not usable at the 1,087-vertex scale.
7. **Synthesis.** See `05-knowledge/results/hn_type_m_hinge_rotations_and_forced_radii_20261014.md` for the synthesis and the record-attempt plan.
8. **Spindle-free rebuild** (`spindle_free_build.log`).
   * `|M1| = 73`, `|M2| = 865`, `|M3| = 32,257`, `|W| = 64,513`.
   * `|E| = 542,472` exact unit edges, i.e. the paper's 542,352 common edges plus the 120 new edges of case 10.
   * The rebuild therefore matches arXiv:2106.11824 exactly. Common denominator 96.
   * The graph itself (an 8 MB npz) is kept out of git. `spindle_free_build.py` regenerates it in about a minute.
   * **Not 4-colourable** (CaDiCaL 1.9.5, 1,974 s). This independently confirms the paper's 5-chromatic, spindle-free `G_64513` (case 10) on an exact rebuild.
9. **Vertex-critical baselines** (`minimize_L727_S361_rho.log`, `crit_*.json`). Four seeds give vertex-critical subgraphs of `L727 ∪ ρS361` with 614, 611, 617 and 618 vertices, against Parts' 509. Their union has 784 vertices; `smus_union4.log` holds a bounded OptUx run on it.
10. **Spindle-free pieces** (`spindle_free_pieces.log`, `spindle_free_forced.log`, `spindle_free_mixed.log`).
    * `ψ = φ_0^(−1) ρ` is Parts' hinge.
    * `M_3` and `M_2 ⊕ M_1^clock` force radius 2; `M_1 ⊕ M_1` does not.
    * All four symmetric reductions tested are 4-colourable.
