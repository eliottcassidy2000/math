# Which tournaments get designated special: from the bypass `H_n` to two 3-cycles glued three ways

**opus S15, 2026-10-01.** This continues the shaved-tournament note
[`shaved_tournaments_unavoidable_cores_20261001.md`](shaved_tournaments_unavoidable_cores_20261001.md).

- Script: `04-computation/experiments/tournament_designations_20261001.py` (+ `.out`).

**Status.**
- FINITE-EXACT: every isomorphism class with `n ≤ 8`, and every object listed.
- Two-block cycles for all `n`: Benhocine–Wojda (1983).
- Oriented Hamiltonian paths for all `n`: Grünbaum (1971), Rosenfeld, Havet–Thomassé (2000).
- No new theorem for general `n`. The results here are the small-`n` census and its structure.
- Independent audit: OWED.

## 0. The question

> the 4-vertex object [is] in every tournament except C3 and one 5-vertex tournament ... do the same analysis on
> part of that one 5-vertex tournament that is made special. see if it designates a 6-vertex tournament as special.

**The mechanism.** A spanning oriented graph `D` *designates* the `n`-tournaments that avoid it. The owner's
`H_4` (`A→B→C→D` plus `A→D`) is the oriented 4-cycle with a source `A`, a sink `D`, and directed paths of lengths
1 and 3. Its family `H_n` (Grünbaum's bypass) designates exactly `C3` and `T*_5 = C3[1,C3,1]` (vertex `a` beats a
3-cycle, the 3-cycle beats `b`, and `b → a`).

**The special part of `T*_5`.** It is its inner 3-cycle, the vertex of `C3` that was blown up. The analysis below
shows what the next objects designate at 6 vertices.

## 1. The same family, one step on: two-block cycles

`D(n; p, n−p)`: a source and a sink joined by directed paths with `p` and `n − p` arcs. Here `p = 1` is `H_n`.

| n | (p, n−p) | designated (avoiding) tournaments |
|---|---|---|
| 3 | (1, 2) | `C3` |
| 4 | (2, 2) | `TT2[1,C3]` (source + 3-cycle) and `TT2[C3,1]` (3-cycle + sink) |
| 5 | (1, 4) | `C3[1,C3,1] = T*_5` |
| 6 | (2, 4) | **`TT2[C3,C3]`: two 3-cycles, one beating the other completely** (`|Aut| = 9`, scores `1,1,1,4,4,4`) |
| 7, 8 | any | none |

**Answer: yes.** The same analysis designates exactly one 6-vertex tournament, `C3 ⇒ C3`. That is the special part
of `T*_5` (its inner 3-cycle) doubled and placed in series.

Every tournament designated by a two-block cycle is a composition of 3-cycles. The reason is local. A 3-cycle
gives each of its vertices in- and out-degree exactly 1 inside it, so it can host neither the source's two
out-arcs nor the sink's two in-arcs. (§1 of the check gives the case analysis for `C3 ⇒ C3`.) The general-`n`
statement, that there are no exceptions beyond these, is Benhocine–Wojda's theorem.

## 2. The Paley line: nearly alternating objects

All oriented Hamiltonian paths and all non-directed oriented Hamiltonian cycles, `n ≤ 8`:
- **Paths.** Only antidirected paths fail, and only in `C3`, `R5` (the regular 5-tournament) and `P7` (Paley). This
  is the classical Grünbaum–Rosenfeld–Havet–Thomassé exception list, reproduced.
- **Cycles.**
  - Two-block words designate the 3-cycle compositions of §1.
  - Four-block and antidirected words designate the "Paley line":
    - `00101` is avoided only by `R5` (`n = 5`);
    - `001011` and `001101` are avoided only by **`G_par`** (`n = 6`);
    - `0010101` is avoided only by `P7` (`n = 7`).
  - Antidirected cycles are avoided by 3, 13 and 19 classes at `n = 4, 6, 8`.
  - Nothing else is avoidable for `n = 7, 8`.
- **Path + one or two extra arcs** (`n ≤ 7`):
  - `H_6 + (0→3)`, the bypass plus one more arc from the first vertex, is avoided only by **`G_par`**.
  - `H_7 + (0→5)` is avoided only by `P7`.

## 3. Two 3-cycles, glued three ways

Every special 6-tournament found is two copies of the 3-cycle (the special part of `T*_5`) glued together.

| gluing of `A = C3`, `B = C3` | cross arcs | what it is | designated by |
|---|---|---|---|
| in series | all `A → B` | `TT2[C3,C3]`, `|Aut| = 9` | `D(6;2,4)`, the two-block family |
| matching, **parallel** rotation | `a_i → b_i`, the rest `B → A` | **`G_par`**: scores `2,2,2,3,3,3`, `|Aut| = 3`, `H = 45`, `hc = 5` | `H_6 + (0→3)`; cycles `001011`, `001101` |
| matching, **antiparallel** rotation | same | **`P7` minus a vertex**: same scores, `|Aut|`, `H`, `hc` | not singled out by any object here |

**`G_par` and `P7` minus a vertex are twins.** They share every invariant in the table, including the maximum
`H = 45` at `n = 6`, yet they differ:
- In `P7` minus a vertex every arc lies on an odd number of Hamiltonian paths (15 of 15; mac-mini's THM-4524). In
  `G_par` only 9 of 15 do.
- Only `P7` minus a vertex extends to the Paley heptagon.

A first version of this analysis identified the designated tournament as `P7` minus a vertex using invariants
alone. The exact class computation corrected that before anything was written up.

## 4. The pattern

**Two families.**
- **Rigid compositions of 3-cycles** (`C3`, `TT2[1,C3]`, `T*_5`, `TT2[C3,C3]`) are designated by objects with few
  blocks.
- **Symmetric "near-Paley" tournaments** (`R5`, `G_par`, `P7`) are designated by objects that alternate.

**At 6 vertices the two families meet.** Every designated 6-tournament is two 3-cycles glued: in series, in
parallel, or antiparallel.

**At 7 vertices only `P7` survives** (for these objects). This matches its role as the obstruction for the largest
unavoidable subgraphs (THM-4526).

At 8 vertices only the antidirected cycle is avoidable. Rosenfeld conjectured that from 9 vertices on, every
non-directed oriented Hamiltonian cycle is unavoidable. Thomason proved this for very large `n`, and Havet (JCTB 80
(2000)) for `n ≥ 68`. For paths, Havet–Thomassé settled all `n`.

## 5. Reproduction

```bash
python 04-computation/experiments/tournament_designations_20261001.py
```

This takes about 15 minutes.
