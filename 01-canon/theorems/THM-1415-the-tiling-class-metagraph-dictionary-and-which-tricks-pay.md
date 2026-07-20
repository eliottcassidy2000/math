---
id: THM-1415
title: "THE TILING / ISO-CLASS / MERGED-METAGRAPH DICTIONARY, EXACTLY — with the mutual-computation tricks audited for what actually pays. (1) The fibration is complete at all three levels, with a one-line double-count proving the global checksum Σ_c H(c)/|Aut(c)| = 2^{C(n−1,2)}. (2) THE HAM-PATH TRICK: fibre(c) ↔ HamPaths(T_c)/Aut(T_c), so every fibre is enumerable from ONE representative and the 2^m cube never has to be touched. (3) TILINGS ARE THE SWITCHING CLASSES: the tiling cube is arc space modulo cut(K_n), the base path is a spanning-tree transversal, and every switching class contains exactly 2^{n−1} tournaments — a base-path-INDEPENDENT description of the tiling set. (4) THE BASE-PATH-INDEPENDENT SUBGROUP IS DEAD IN BOTH DIRECTIONS: ∩_P Γ_P = 0 and ⟨∪_P Γ_P⟩ = the whole arc space, at n = 4,5,6,7 — so there is neither a common star transformation nor a common star invariant, and kind-pasteur's proposed repair cannot exist. (5) THE d=1 WIGGLY QUOTIENT IS EXACTLY G_n (5/30/290 edges at n=4/5/6), because of the AVOIDABLE-ARC LEMMA: whenever every Hamiltonian path of T uses arc a, flipping a is a self-loop (verified n ≤ 7, zero counterexamples). (6) TRICK AUDIT: the Aut-orbit refinement decays 33%→7% across n=4..7 and is asymptotically worthless because almost all tournaments are rigid; the class-BFS is the trick that pays."
status: >
  MIXED, stated per part.
  (1) PROVED (double count) + VERIFIED exhaustively n = 4,5,6.  The orbit-stabiliser
      half is CANON (CLAUDE.md, LEM-003) and was verified per class by boxeph-S157 the
      same day -- it is NOT claimed here; §1 adds the merged/sigma-orbit refinement and
      the blue-count identity.
  (2) PROVED (orbit-stabiliser + LEM-003 freeness) + VERIFIED against exhaustive n=4,5,6
      and run at n=7.
  (3) PROVED (spanning-tree transversal for the cut space) + VERIFIED n = 4,5,6 with
      every switching class of size exactly 2^(n-1).
  (4) VERIFIED n = 4,5,6,7 (exact GF(2) linear algebra over ALL n!/2 Hamiltonian paths).
      A decisive negative for a proposed method, not a theorem about tournaments.
  (5) PROVED at n = 4,5,6 by direct comparison; the avoidable-arc lemma that would give
      it for ALL n is VERIFIED n <= 7 (zero counterexamples) and PROVED for the
      transitive class only.  Stated as a lemma-conditional general claim, not as proved.
  (6) MEASURED, n = 4..7.
  Nothing here advances a named open problem; it is infrastructure plus one decisive
  negative.
source: klein-2026-07-20-S336 (owner: figure out how tiling sets map to iso class nodes in the merged metagraph exactly; use the edges and nodes and tilings to all compute each other efficiently, look for tricks; consider even more creative ideas than a descending star-type invariant coming from a base-path-independent subgroup -- the natural candidate being the intersection of Gamma over all spanning paths)
depends_on:
  - THM-549   # sigma = the complement in tiling coordinates; Fix(sigma); h(n)
  - THM-1382  # star flips span cut(H)
  - THM-1405  # rank Gamma = n-1 exactly (kind-pasteur / mac-mini)
  - THM-1410  # klein-S335: sigma permutes, phi translates
prior_art:
  - "CLAUDE.md + LEM-003: tilings * |Aut| = H (orbit-stabiliser on the tiling fibration)"
  - "boxeph-2026-07-20-S157: same owner directive; verified fibre*|Aut| = H and sum fibres = 2^m per class, and computed HALF = 40/544, QUARTER = 21/276 at n=5/6. §1 overlaps their work and credits it."
  - "death-star: G_7/Z_2 with V = 272 -- reproduced exactly here, validating the build."
related:
  - THM-1400  # kind-pasteur: base-path-relativity of cut-space constructions
script: 04-computation/tiling_class_dictionary_klein_S336.py (+ .out)
---

# THM-1415 — the dictionary, and which tricks pay

## 0. What is mine and what is not

The orbit-stabiliser fibre law `|fibre(c)| · |Aut(c)| = H(c)` is **canon** (CLAUDE.md, LEM-003), and
**boxeph-S157 verified it per class on the same owner directive, hours before this session.** It is
not claimed here. §1 records the dictionary in full because the rest depends on it, and adds only
the merged/σ-orbit refinement and the blue-count identity. §§2–6 are the new content.

## 1. The dictionary, all three levels

Write `c` for an iso class, `H(c)` for its Hamiltonian-path count, `B(c)` for its number of
**blue** (grid-symmetric, σ-fixed) tilings, and `h = ⌊(n−1)²/4⌋` (THM-549).

```text
   2^m tilings   --->   A000568(n) iso classes   --->   (A000568 + SC)/2 merged nodes
        m = C(n−1,2)              fibre = H/|Aut|            fibre = (2−[SC])·H/|Aut|
```

| statement | status |
|---|---|
| `\|fibre(c)\| = H(c)/\|Aut(c)\|` | canon (LEM-003); re-verified |
| `Σ_c H(c)/\|Aut(c)\| = 2^{C(n−1,2)}` | **proved below**; verified n=4,5,6 |
| merged fibre `= (2−[c SC])·H/\|Aut\|` | verified n=4,5,6 |
| σ-orbits over a merged node `= H/\|Aut\|` (non-SC), `= (H/\|Aut\| + B)/2` (SC) | verified n=4,5,6 |
| `B(c) = 0` for every non-SC class | verified n=4,5,6 |
| `Σ_{SC classes} B(c) = 2^h` | verified n=4,5,6 |

**Proof of the checksum.** Count pairs `(T, Q)` with `T` a labelled tournament on `[n]` and `Q` a
Hamiltonian path of `T`. Summing over `T` gives `Σ_T hp(T)`. Summing over `Q` instead: each of the
`n!` vertex sequences is a Hamiltonian path of exactly `2^{C(n,2)−(n−1)} = 2^{C(n−1,2)}` tournaments
(the path pins `n−1` arcs, the rest are free). So `Σ_T hp(T) = n! · 2^{C(n−1,2)}`. Grouping the
left side by iso class, each class contributes `(n!/|Aut|)·H`, and dividing by `n!` gives
`Σ_c H(c)/|Aut(c)| = 2^{C(n−1,2)}`. ∎

The last two rows say the **blue tilings live exactly over the SC classes** and total `2^h` — which
is THM-549's `|Fix(σ)|`, now *distributed* across the SC spine. So the half tiling is precisely the
fibration over the merged metagraph, and `B(c)` is the exact SC correction term.

## 2. The Ham-path trick — fibres without touching the cube

**Proposition.** `fibre(c) ↔ HamPaths(T_c) / Aut(T_c)`, for any single representative `T_c`.

*Proof.* A tiling in `fibre(c)` is a labelled copy of `c` in which `n→n−1→…→1` is a Hamiltonian
path. Choosing such a copy is the same as choosing a Hamiltonian path of `T_c` to *become* that
sequence, and two choices give the same tiling iff they differ by an automorphism. `Aut` acts freely
on Hamiltonian paths (LEM-003), so the count is `H/|Aut|`. ∎

This is the operational form: **the fibre is not just counted but enumerated, from one
representative.** Everything below is downstream of it — the whole dictionary can be built from a
list of class representatives, never touching `2^m`.

| `n` | classes | SC | merged | `G_n` edges | canonicalisations | time | naive `2^{C(n,2)}` |
|---|---|---|---|---|---|---|---|
| 4 | 4 | 2 | 3 | 5 | 24 | <0.1 s | 64 |
| 5 | 12 | 8 | 10 | 30 | 120 | <0.1 s | 1,024 |
| 6 | 56 | 12 | 34 | 290 | 840 | <0.1 s | 32,768 |
| 7 | 456 | 88 | **272** | 4,086 | 9,576 | 0.6 s | 2,097,152 |

The `n=7` merged count **272** reproduces death-star's independently built `G_7/Z_2` exactly, which
validates the build. Checksum `Σ H/|Aut| = 32,768 = 2^{15}` holds at `n=7`.

## 3. Tilings *are* the switching classes

**Theorem.** The tiling cube is arc space modulo the cut space: `GF(2)^{E(K_n)} / cut(K_n)`. The base
path, being a spanning tree, meets every coset exactly once, so tilings are a **transversal for
switching** — and every switching class contains exactly `2^{n−1}` labelled tournaments.

*Proof.* `cut(K_n)` acts simply transitively on the orientations of any spanning tree; the base path
is one. So each coset has a unique member whose path arcs all point `k → k−1`, which is a tiling.
Counting: `2^{C(n,2)}/2^{n−1} = 2^{C(n−1,2)} = 2^m`. ∎

Verified `n = 4,5,6`: `8/64/1024` distinct switch-canonical tournaments, every fibre of size exactly
`2^{n−1}`.

**Why this matters.** It gives a description of the tiling *set* that does not mention the base path:
*tilings = switching classes of tournaments.* Equivalently, every tournament can be switched, in
exactly one way, to make a prescribed Hamiltonian path descend. The base path enters only in *naming*
the coordinates — which is exactly the base-path-relativity kind-pasteur's THM-1400 diagnosed, now
localised: **the set is canonical, the coordinates are not.**

## 4. The base-path-independent subgroup — dead in both directions

kind-pasteur (THM-1400) proposed that a descending star-type invariant must come from a
base-path-independent subgroup, "e.g. the intersection of `Γ` over all spanning paths", and
boxeph-S157 handed this on unrun as "the canonical invariant candidate". **It does not exist.**

Embed every `Γ_P = ⟨star_{K_n∖P}(v)⟩` in the common ambient `GF(2)^{E(K_n)}` and range over all
`n!/2` Hamiltonian paths:

| `n` | Ham paths | `C(n,2)` | `dim Γ_P` | `dim ⟨∪_P Γ_P⟩` | `dim ∩_P Γ_P` |
|---|---|---|---|---|---|
| 4 | 12 | 6 | 3 | **6** | **0** |
| 5 | 60 | 10 | 4 | **10** | **0** |
| 6 | 360 | 15 | 5 | **15** | **0** |
| 7 | 2,520 | 21 | 6 | **21** | **0** |

Both directions fail, and they are the two things one could have wanted:

- `∩_P Γ_P = 0` — there is **no transformation** that is a star flip for every presentation.
- `⟨∪_P Γ_P⟩` is everything — hence `∩_P (invariants of Γ_P) = 0`: there is **no invariant** shared
  by all presentations either.

So the repair is not merely hard, it is impossible: the star construction has nothing
base-path-independent in it at all. Combined with §3 this is a clean split — *the tiling set is
canonical, the star group on it is not*.

## 5. The `d=1` wiggly quotient is exactly `G_n`

The tiling model's single-tile flips only see the `C(n−1,2)` **non-path** arcs, so *a priori* the
`d=1` wiggly quotient could be a proper subgraph of the single-arc-flip metagraph `G_n`. It is not:

| `n` | `G_n` edges | `d=1` wiggly quotient edges | equal? |
|---|---|---|---|
| 4 | 5 | 5 | **yes** |
| 5 | 30 | 30 | **yes** |
| 6 | 290 | 290 | **yes** |

**The reason — the avoidable-arc lemma.** To realise a `G_n` edge `{[T],[T+a]}` inside the tiling
model one needs a Hamiltonian path of `T` avoiding `a` (then relabel it to be the base path, and `a`
becomes a tile). So the equality follows from:

> **Lemma (avoidable arc).** If every Hamiltonian path of `T` uses the arc `a`, then `T + a ≅ T`.

| `n` | forced `(class, arc)` pairs | of those, flip is a self-loop | counterexamples |
|---|---|---|---|
| 4 | 3 | 3 | **0** |
| 5 | 6 | 6 | **0** |
| 6 | 13 | 13 | **0** |
| 7 | 36 | 36 | **0** |

Proved for the transitive class: its unique Hamiltonian path forces all `n−1` consecutive arcs, and
flipping a consecutive arc of a linear order gives another linear order — a self-loop. The general
case is **verified, not proved**, and is the honest gap: with it, `d=1` wiggly quotient `= G_n` for
all `n`.

## 6. Trick audit — what actually pays

Stated as a negative-included audit rather than a list of wins:

- **PAYS: the class-BFS.** Building `G_n` by BFS over iso classes costs `|classes| · C(n,2)`
  canonicalisations — 9,576 at `n=7` against 2,097,152 tournaments naively, and against `2^{21}`
  tilings for the cube route. This is the dominant saving.
- **PAYS: the Ham-path fibre formula (§2).** It removes the `2^m` cube from the pipeline entirely;
  at `n=8` the class route needs `6880 × 28 = 192,640` canonicalisations against `2^{21}` tilings
  and `2^{28}` tournaments.
- **PAYS: the checksum `Σ_c H/|Aut| = 2^{C(n−1,2)}`** as a cheap total-correctness test on any
  class census — it caught nothing here, which is the point of running it.
- **DOES NOT PAY: the Aut-orbit refinement.** Flipping one arc per `Aut(T)`-orbit instead of per arc
  saves `33.3%, 26.7%, 16.2%, 6.9%` at `n = 4,5,6,7`. It **decays**, because the fraction of rigid
  classes rises (`2/4, 7/12, 41/56, 399/456` = 50%, 58%, 73%, 88%) and almost all tournaments are
  rigid asymptotically. Recorded so nobody implements it expecting a win at large `n`.

## 7. Scope

Infrastructure and one decisive negative (§4). No `H`-value, LRC, or extremal claim. §5's general
form waits on the avoidable-arc lemma.

*Files: `04-computation/tiling_class_dictionary_klein_S336.py` (+ `.out` in `05-knowledge/results/`).*
