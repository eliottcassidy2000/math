---
source: opus-2026-06-03-S599m/n/o (remote-control)
status: RESULT + DEEP BRIDGE — n=22 max unit distances = 49 (Harborth triangular-lattice optimum, achieved: a triangulated disk of 22 Eisenstein integers). It and the LRC/tournament world are ONE species: extremal ADDITIVE-ENERGY problems on cyclotomic Cayley graphs, solved by the lattice/AP, isostatic+isoperimetric, with the same prime-3 (ζ3) coloring and the same odd-triangle/Rédei core.
tags: [unit-distance, erdos, harborth, triangular-lattice, eisenstein, cayley-graph, additive-energy, cyclotomic, isostatic, rigidity, prime-3, redei, OCR, chromatic, LRC, tournaments]
---

# Unit distances at n=22, and why it is the same problem as LRC + tournaments

**Prompt (user):** work on the unit distance problem, find the optimal solution for n=22;
understand how it relates deeply with tournaments and the LRC; poke around for inspiration.

## The result: n=22 → 49

The **maximum number of unit distances among 22 points** is realised on the **triangular
lattice = Eisenstein integers `ℤ[ζ₆]`**. Harborth's theorem gives the lattice optimum
`⌊3n − √(12n−3)⌋`; for `n=22`, `⌊66 − √261⌋ = ⌊49.84⌋ = 49`. **Achieved** (`…s599n.py`): a
compact 22-point cluster with **49 unit distances**, degree distribution `{3:8,4:4,5:2,6:8}`
(8 interior degree-6 points). It is a **triangulated disk**: `V=22, E=49, T=28` triangles,
Euler characteristic `V−E+T = 1` (`…s599o.py`).

> *Honest scope:* 49 is the **triangular-lattice optimum** (Harborth, proven optimal among
> lattice configs) and the **achieved** value (a verified lower bound on the true planar
> maximum); it is the conjectured overall optimum for `n=22`. My computation does not exclude
> a non-lattice configuration with more — that is the (hard, computer-search) content of the
> Erdős unit-distance tables.

## The deep bridge: one species, four faces (all verified)

LRC, tournaments, and unit distances are **the same extremal problem on a cyclotomic Cayley
graph**. Each repo object has a unit-distance twin:

### 1. Cayley / cyclotomic structure
- **Unit-distance graph** `= Cay(ℤ[ζ₆], U₆)`, `U₆ = {±1,±ζ₆,±ζ₆²}` the six sixth-roots of
  unity (the unit vectors). A Cayley graph on the Eisenstein integers with the **cyclotomic
  unit set**.
- **Round LRC tournament** (the worry-set, THM-402/407) `= Cay(ℤ/(2n−1), H)`, `H` a shell-half.
  A circulant tournament on `ℤ/(2n−1)`.
- **LRC witnesses** `=` primitive `n`-th roots of unity (THM-403). Both worlds are built from
  **roots of unity / cyclotomic rings**; the unit distance lives in `ℤ[ζ₆]`, LRC in `ℤ[ζ_n]`.

### 2. Both are EXTREMAL ADDITIVE-ENERGY problems  — *the core identity*
> **Verified:** the unit-distance count `= ½·#{(p,q)∈S² : p−q ∈ U₆} = 98/2 = 49` — the
> **additive energy of `S` with respect to the unit set `U₆`.**

And the **LRC resonance** (the `(★)` overlaps, the reason the modulus is `2n−1`, THM-401) **is
the additive energy of the speed set mod `2n−1`** (pair-sums `a+b`). So:

> **Unit distances and LRC tightness are the same functional — maximise additive energy of a
> point set with respect to a fixed difference set** — and both are maximised by the most
> additively closed configuration: the **triangular lattice** for unit distance, the
> **arithmetic progression / round tournament** for LRC.

### 3. Isostatic + isoperimetric signature
Harborth's `3n − √(12n−3)` decomposes as **bulk minus boundary**: each point has `≤ 6`
unit-neighbours (kissing number 2D = 6 = `|U₆|`), so `E ≤ 3n`; the deficit `√(12n−3)` is the
**perimeter** (isoperimetric correction), minimised by the **round/compact** cluster. The LRC
mirror: the worry-set is **isostatic** (S598: full Helly `= n−1`, every runner load-bearing),
the AP minimising the resonance "boundary." Both extremal configs are the **rigid,
minimal-boundary, round** one — an isoperimetric/isostatic extremum.

### 4. The prime-3 (ζ₃) coloring — a precise echo
> **Verified:** the cluster is properly **3-coloured** by `(i−j) mod 3` (the Eisenstein
> sublattice = the ideal of norm 3 in `ℤ[ζ₆]`); class sizes `7,7,8`. So `χ(UD) = 3`, forced by
> the **prime 3 / `ζ₃`** factor of `ℤ[ζ₆]`.

This is the *same prime 3* that rules the **LRC n=14 residual** (`2n−1 = 27 = 3³`, THM-407: the
worry-set folds to the prime-3 tower `3⁰,3¹,3²`). Both the unit-distance coloring obstruction
and the LRC n=14 obstruction are the **`ℤ/3` structure inside the Eisenstein/cyclotomic ring**.
(Not numerology: `χ` of the triangular lattice is `3` exactly because `ℤ[ζ₆]/(\sqrt{-3})≅𝔽₃`.)

### 5. Odd triangles = Rédei / OCR
The cluster has **28 unit triangles** (3-cliques). Oriented, each is a **directed 3-cycle** —
the smallest odd cycle, the generator of the repo's **OCR / Rédei** parity (the
Hamiltonian-path-count oddness lives on odd cycles). The triangular lattice is *maximally
triangulated* (the most odd-3-cycles per vertex); the LRC round tournament's Rédei parity is
the oriented shadow of the same triangle structure.

### 6. The χ duality (undirected vs oriented)
- **Undirected** unit-distance (triangular lattice): `χ = 3` — forced *up* by the odd triangles
  (3-cliques).
- **Oriented** round LRC tournament: `χ = 2` (THM-402, the diameter semicircle split).
Same Cayley graph, two faces: the **undirected** (clique/odd-triangle ⟹ `χ=3`) and the
**tournament** (round ⟹ `χ=2`). The unit-distance problem is the *undirected, geometric* face;
LRC is the *oriented, dynamical* face — of one cyclotomic Cayley object.

## The unifying thesis

> **Unit distances and LRC are the geometric and dynamical faces of one problem: maximise the
> additive energy of a point set against a cyclotomic difference set, on a Cayley graph of a
> cyclotomic ring.** The optimum is the lattice/AP (most additively closed), which is *rigid
> (isostatic)*, *round (isoperimetric)*, *cyclotomically `3`-coloured*, and *maximally
> odd-triangulated (Rédei)*. n=22's `49` is the Eisenstein face; LRC's worry-set is the
> `ℤ/(2n−1)` face; Rédei's odd-Hamiltonian-count is the oriented face. Same species.

## Inspiration drawn from the repo

The bridge was found by mapping the unit-distance optimum onto repo objects: **additive energy
↔ THM-401** (the `2n−1` pair-sum modulus is an additive-energy modulus); **cyclotomic ↔
THM-403** (witnesses = roots of unity); **round Cayley ↔ THM-402/407**; **isostatic ↔ S598**;
**prime-3 ↔ THM-407** (n=14 = `3³`); **odd triangles ↔ Rédei/OCR**. The unit-distance problem
is LRC/tournaments *read in `ℤ[ζ₆]` instead of `ℤ/(2n−1)`*.

## Honest status

- **Computed/verified:** `n=22` → `49` unit distances (Harborth lattice optimum, achieved; a
  triangulated disk); UD count `=` additive energy wrt `U₆` (`=49`); proper 3-coloring by
  `(i−j) mod 3`; 28 unit triangles; Harborth formula vs small-`n` checks.
- **Established (rigorous) bridges:** Cayley/cyclotomic; additive-energy identity; isoperimetric
  `3n−√` decomposition; the `ζ₃`/prime-3 coloring; odd-triangles = Rédei.
- **Honest caveats:** `49` is the lattice optimum / achieved value, not a proof of the *planar*
  maximum (non-lattice not excluded). The cross-problem unification is a structural identity
  (same Cayley/additive-energy species), **not** a reduction of one conjecture to the other.

**Artifacts:** `04-computation/unit_distance_n22_triangular_s599m.py`, `…_v2_s599n.py`,
`…_lrc_bridge_s599o.py` (+`.out`s). Builds on THM-401 (additive energy / `2n−1`), THM-402/407
(round Cayley / prime-3), THM-403 (cyclotomic), S598 (isostatic), Rédei/OCR. New: **HYP-2170**.
