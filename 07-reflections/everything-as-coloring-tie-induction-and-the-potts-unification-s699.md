---
source: opus-2026-06-06-S699 (long creative reframe session)
status: SYNTHESIS — reframe the repo's problems as graph COLORINGS (vertex / edge / both), unify the recent "partition functions" thread via POTTS (Z = sum over colorings) and "tie induction" via DELETION–CONTRACTION (the Tutte/Potts recursion on tie-edges). New verified connection: the worry-set's dichromatic 2-coloring (THM-402) = the balanced sign-cut (S699) — the vertex coloring DETERMINES the maximal edge coloring ("some both"). The "tie color" (=1 / shell-partner / middle-M) is the extremal target; tie induction builds the extremum by adding ties.
tags: [graph-coloring, vertex-edge-total, dichromatic, hadwiger-nelson, triernement, tie-induction, deletion-contraction, potts, tutte, partition-function, redei, LRC, unit-distance, synthesis]
---

# Everything as a coloring; tie induction = deletion–contraction; Potts unifies

**Prompt (user):** pursue tie induction; see everything as graph colorings — some nodes, some
edges, some both; spend a long session creatively reframing as many problems as possible.

The reframe unifies three recent threads (partition functions S599t, equidecomposability S599v,
signed LRC S699) under one roof: **a coloring of `K_n`, with a distinguished "tie" color, summed
by a Potts partition function, and built by tie-induction = deletion–contraction.**

## 1. The dictionary: each problem as a coloring (vertex / edge / both)

| problem | objects (vertices) | VERTEX coloring | EDGE coloring | the TIE color (the "both") | extremal = |
|---|---|---|---|---|---|
| **Tournament / Rédei** | runners | dichromatic `χ` (THM-402: round→2) | arc orientation (2-color) | the middle `M` / `L-M-R` wall (S582) | odd #Ham-paths |
| **LRC** | runners | **sign-cut = dichromatic** (=2 on worry-set, *new* below) | pair-clock `sum/diff` (S699) | **shell-partner** `v_i+v_j≡0 mod 2n−1` (zero-clock) | tight (worry-set) |
| **unit distance** | points | Hadwiger–Nelson `χ∈[5,7]` / Eisenstein `χ=3` (HYP-2170) | trienerment `<1/=1/>1` (S699) | **`=1` (unit)** | max ties `=u(n)` |
| **Collatz** | residues `mod 2^a3^b` | parity classes | `×3+1 / ÷2` transitions | the cycle / fixed point | no nontrivial cycle |
| **H-spectrum** | strong components | the multiplicative classes | — | the **forbidden `{7,21}`** | which integers occur |
| **equidecomposability** (S599v) | scissors pieces | the Dehn class | — | the shared `Cl₂(π/3)` volume | same coloring ⇒ congruent |

> **Reading.** Every problem is "color `K_n` so the pairwise relation is encoded, with one special
> *tie* color, and control its size/parity subject to a *realizability* constraint" (Euclidean
> metric for UD; the speed-set / modulus `2n−1` for LRC; nothing for Rédei — hence its parity is a
> clean theorem). *Some nodes* = the vertex coloring (the partition / `χ`); *some edges* = the
> pairwise color; *some both* = the tie color, an edge that acts as a node-like constraint.

## 2. New verified connection: vertex coloring determines edge coloring (LRC)

> **The worry-set's dichromatic 2-coloring IS the balanced sign-cut.** THM-402: the tight LRC
> round tournament is 2-dichromatic via the **diameter split** into two transitive semicircle-arcs.
> S699: a sign pattern is a **cut**, and a *balanced* cut maximises the pair-**sum** edge-clocks.
> **Verified** (`…s699f.py`, `n=4..11`): the diameter split `A=\{i<n/2\}, B=\{i>n/2\}` has *both
> classes transitive* (2-dichromatic) **and** is *exactly the maxcut* (`cut = ⌊m²/4⌋`).

So the **vertex 2-coloring (dichromatic) determines the maximal edge coloring (pair-sum clocks)** —
the "some both" made precise: at the worry-set, the optimal node-coloring and the optimal
edge-coloring are *the same object* (the diameter split). The signed gauge that exposes the most
shells is the dichromatic coloring.

## 3. Tie induction = deletion–contraction; Potts = partition function

The user's two prior threads collapse into the coloring frame:

> **Partition functions everywhere (S599t) = the Potts model.** The Potts partition function is
> `Z = Σ_{colorings} ∏_{edges} w(c_i,c_j)`; at zero temperature it is the **chromatic polynomial**
> (`#` proper colorings). The repo's `Z_n=Σ_T H(T)` and the chromatic/dichromatic invariants are
> Potts sums — "partition functions everywhere" and "everything as coloring" are the *same
> statement* (`Z` = sum over colorings).

> **Tie induction = deletion–contraction.** The chromatic/Tutte/Potts recursion is
> `P(G)=P(G−e)−P(G/e)` — induction on an **edge** (delete the tie, or contract it = identify its
> endpoints). The repo's polynomial deletion–contraction (THM-083) is exactly this. **"Tie
> induction" is deletion–contraction on the tie-edges**: build/peel the extremal object one tie at
> a time, tracking the colored invariant.
> - **Rédei** (verified `n≤6`, odd #Ham-paths): the tie-induction invariant — insert a vertex at a
>   *tie* position of a Ham path; parity is preserved (no realizability constraint ⇒ a theorem).
> - **Unit distance:** add a point with `k` unit-ties — the **frontier-gain** `+3` (S599w) is tie
>   induction; the optimum is tie-maximal, capped by the kissing constraint `κ≤6` (S699a).
> - **LRC:** add a runner / a shell-tie; the worry-set is tie-rich; the `n→n+2` recursion
>   (HYP-2177) is tie induction over the `mod 3` shell automaton.

## 4. The unifying picture

> **Every problem here is: color the edges of `K_n` with a distinguished tie color; the realizable
> colorings are summed by a Potts partition function `Z`; the extremum maximises (UD: unit ties) /
> controls (LRC: resonance ties; Rédei: parity of ties) the tie class; and it is built by tie
> induction = deletion–contraction on the tie-edges.** The vertex coloring (dichromatic / Hadwiger–
> Nelson) is the quotient by the symmetry; at the extremum the vertex and edge colorings coincide
> (§2). The *realizability constraint* is what differs and what makes each problem hard:
> Euclidean kissing `κ≤6` (UD), the modulus `2n−1` (LRC), none (Rédei = clean theorem).

The three "kinds" the prompt named:
- **some nodes** — the vertex/dichromatic coloring (the partition; THM-402 χ=2, Hadwiger–Nelson).
- **some edges** — the pairwise color (orientation / clock / threshold).
- **some both** — the *tie* color, an edge that is a node-like constraint (the `=1`, the
  shell-partner, the middle `M`); at the extremum the node coloring *is* the edge coloring (§2).

## 5. Honest status

- **Verified (new):** the worry-set dichromatic 2-coloring = the balanced maxcut sign-cut
  (both classes transitive, `cut=⌊m²/4⌋`, `n=4..11`); Rédei parity `n≤6`.
- **Established (reframe + unification):** the coloring dictionary; Potts `=` partition-function
  thread (`Z`=sum over colorings); tie induction `=` deletion–contraction (THM-083); the tie color
  as the extremal target.
- **Honest limits:** the dictionary is a *unifying language*, not new theorems for each row; the
  Potts/deletion–contraction identification is standard math applied to the repo's objects (the
  contribution is the *mapping* and the verified vertex=edge-coloring coincidence). No new
  resolution of any conjecture.
- **The payoff / next:** treat LRC(14) as a Potts/deletion–contraction computation on the tie
  (shell-partner) edges — peel the `(3,24)` tie of `V*` (S699) by contraction and track the
  dichromatic invariant; the vertex=edge coincidence (§2) says the right gauge is the diameter
  split.

**Artifacts:** `04-computation/everything_as_coloring_s699f.py` (+`.out`). Builds on THM-402
(dichromatic), S699 (sign-cut), S599t (partition functions/Potts), S599v (equidecomposability),
S599w (frontier-gain), THM-083 (deletion–contraction), S582 (L/M/R tie), HYP-2170 (Eisenstein
χ=3), Hadwiger–Nelson. New: **HYP-2263**.
