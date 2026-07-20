---
id: THM-1410
title: "THE HALF FOLD PERMUTES, THE QUARTER FOLD TRANSLATES — the two generators of the quarter tiling sit on opposite sides of the star (cut/cycle) split. σ is the tile permutation induced by the vertex involution ρ(v)=n+1−v, a graph automorphism of H = K_n∖P, so σ(star v) = star(n+1−v) and σ NORMALISES the star group: the half fold is INTERNAL to the star algebra, acting on the holonomy word by a permutation. φ (the antipodal, flip-all-tiles map) lies in cut(H) iff H is bipartite iff n ≤ 4; for n ≥ 5 the triangle {1,3,5} ⊆ H obstructs, and φ acts on the holonomy word by TRANSLATION by the nonzero functional ε(C) = |C| mod 2, hence FREELY on the 2^{m−n+1} star orbits. So the quarter fold costs EXACTLY ONE holonomy bit and the surviving conserved quantities are exactly the EVEN-cycle parities, of dimension m−n. Quantitatively the defect of φ from the star group is m − maxcut(H) = C(n−1,2) − ⌊n²/4⌋ + 1 = ⌊(n−3)²/4⌋, carried by an explicit tile set: the tiles internal to the two halves of the base path. That number equals the THM-549 half-tiling size h(n−2), but the Mode-B reading of that coincidence is REFUTED in §4c — both sides are quarter-squares for the elementary two-triangles-make-a-half-square reason, and the tile sets are not in bijection. And the triangle's three sides read: the two LEGS are star_H(1) and star_H(n), swapped by σ; the HYPOTENUSE is Fix(σ), the ρ-mirror matching."
status: >
  PROVED.  Every part is elementary and complete: (1) ρ ∈ Aut(H) is immediate;
  (2) the bipartiteness criterion plus the explicit triangle; (3) linearity of
  size-parity on the cycle space; (4) maxcut(K_n∖P) = ⌊n²/4⌋ − 1 with matching
  interval bipartition and a one-line connectivity upper bound; (5) direct.
  INDEPENDENTLY VERIFIED n = 3…11 (all five parts, plus the translation law on
  2,000 random tilings per n, zero violations) and the defect closed form to n = 15.
  n = 3 is degenerate (m = 1 < n−1) and is excluded from the rank statements.
  This is STRUCTURE, not a bound: it says nothing about H-values, LRC, or metagraph metrics.
source: klein-2026-07-20-S335 (owner: consider half and quarter tilings, how they relate to the map-graph work, and whether the new insights sharpen them — build on the repo's existing results, do not duplicate)
depends_on:
  - THM-549   # the half fold: sigma, Fix(sigma) = {x+y=n+1}, h(n) = floor((n-1)^2/4)
  - THM-1382  # klein-S333: star flips span cut(H), invariants are cycle(H)
  - THM-1405  # kind-pasteur-S128c109: rank Gamma = n-1 exactly; the holonomy quotient
related:
  - THM-550   # codex: the half-tiling parity recurrences
  - THM-790   # the LEG LAW — the score-level shadow of part (5)
  - THM-584   # complement = antipode on the arc hypercube
  - MISTAKE-033  # sigma (complement tournament) vs phi (complement tiling) — the distinction this theorem is built on
script: 04-computation/halfquarter_starsplit_klein_S335.py (+ .out)
---

# THM-1410 — the half fold permutes, the quarter fold translates

## 0. What was already known, and what this adds

The repo has two folds of the tiling cube `F = GF(2)^m`, `m = C(n−1,2)`:

- the **half fold** by the grid reflection `σ(x,y) = (n+1−y, n+1−x)` — THM-549/550. `Fix(σ)` on
  tiles is the antidiagonal `{x+y = n+1}`, of size `⌊(n−1)/2⌋`, and the half-region has size
  `h(n) = ⌊(n−1)²/4⌋` (A002620). THM-549 §2: **`σ` *is* the tournament complement** in tiling
  coordinates (`σ(tiling) = ρ(T^op)`, `ρ(i) = n+1−i`).
- the **quarter fold** by the Klein 4-group `⟨σ, φ⟩` — kind-pasteur S15/S16 — where `φ` is the
  **antipodal** map, flip *all* tiles. MISTAKE-033 exists precisely to keep `φ ≠ T^op` apart.

Separately, THM-1382 (klein-S333) and THM-1405 (kind-pasteur-S128c109) established the **star
split**: with `H = K_n ∖ P` (so `E(H)` = the tiles), the vertex-star flips span `cut(H)`, of rank
exactly `n−1` for `n ≥ 4`, and the invariants are `cycle(H)`, of dimension `m−n+1` — mac-mini's
*holonomy word*.

**Nobody has asked where `σ` and `φ` sit relative to that split.** They sit on opposite sides, and
that is the whole content below. It is the sharpening the owner asked for: the two generators of
the quarter fold are not two symmetries of the same kind.

Throughout, `hol(t) ∈ cycle(H)*` denotes the holonomy word, `hol(t)(C) = ⟨t, C⟩`; two tilings share
a star orbit iff they share a holonomy word.

## 1. The half fold is internal: σ normalises the star group

**Proposition.** `σ` is the tile permutation induced by the vertex involution `ρ(v) = n+1−v`, and
`ρ ∈ Aut(H)`. Hence `σ(star(v)) = star(n+1−v)`, so `σ(cut(H)) = cut(H)` and `σ(cycle(H)) = cycle(H)`.

*Proof.* Reading a tile as an edge `{x,y}` of `H`, `σ{x,y} = {n+1−x, n+1−y} = ρ{x,y}`. The base
path `P` has edge set `{{k,k−1}}`, and `ρ{k,k−1} = {n+1−k, n+2−k}`, again consecutive, so `ρ(P) = P`
and therefore `ρ(H) = H`. A graph automorphism carries `δ_H(v)` to `δ_H(ρ(v))`, which is the star
statement; the stars span `cut(H)`, so `cut(H)` is `σ`-stable. `σ` is a coordinate permutation,
hence orthogonal for the standard `GF(2)` form, so it also preserves `cycle(H) = cut(H)^⊥`. ∎

So `σ` descends to the star-orbit space and acts there **linearly, by a permutation of the holonomy
word**. The half fold never leaves the star algebra.

## 2. The quarter fold is transverse: φ is not a star flip

**Proposition.** `φ` = translation by the all-ones vector `1 = E(H)`. Then

```text
1 ∈ cut(H)   ⟺   H is bipartite   ⟺   n ≤ 4.
```

*Proof.* `1 = δ_H(S)` says every edge of `H` crosses `S`, i.e. `H` is bipartite with parts `S`,
`V∖S`. For `n ≥ 5` the vertices `1, 3, 5` are pairwise non-consecutive, so `{1,3}, {3,5}, {1,5}`
are all tiles: an odd cycle. For `n = 4`, `H = {41, 31, 42}` is the path `3–1–4–2`; for `n = 3`,
`H` is a single edge. ∎

## 3. What φ does to the holonomy word: translation by ε, and freeness

Define `ε : cycle(H) → GF(2)`, `ε(C) = |C| mod 2`. It is **linear**, since
`|C₁ Δ C₂| = |C₁| + |C₂| − 2|C₁ ∩ C₂|`.

**Theorem.** `hol(t + 1) = hol(t) + ε`, and `ε ≠ 0` iff `H` is non-bipartite, i.e. iff `n ≥ 5`.
Hence for `n ≥ 5`:

1. `φ` acts on the `2^{m−n+1}` star orbits as translation by a **nonzero** vector, therefore
   **freely** — it pairs orbits, fixing none, and the orbit count is even.
2. The quarter fold costs **exactly one holonomy bit**.
3. The conserved quantities surviving the quarter fold are exactly the **even-cycle parities**
   `ker ε`, of dimension `(m−n+1) − 1 = m − n`.

*Proof.* `hol(t+1)(C) = ⟨t,C⟩ + ⟨1,C⟩ = hol(t)(C) + |C| mod 2`. And `ε = 0` iff every cycle-space
element has even size iff `H` has no odd cycle. Translation by a nonzero vector of an `F₂`-affine
space is a fixed-point-free involution. A functional `C` gives a `φ`-invariant conserved quantity
iff `⟨1,C⟩ = 0`, i.e. `C ∈ ker ε`, a hyperplane. ∎

**Verification** (`n = 3…11`): the translation law held on 2,000 random tilings at each `n` with
zero violations; `1 ∈ cut(H)` agreed with bipartiteness at every `n`; and

| `n` | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|
| `m` | 1 | 3 | 6 | 10 | 15 | 21 | 28 | 36 | 45 |
| holonomy dim `m−n+1` | 0 | 0 | 2 | 5 | 9 | 14 | 20 | 27 | 35 |
| **quarter** dim `m−n` | 0 | 0 | **1** | **4** | **8** | **13** | **19** | **26** | **34** |
| defect `⌊(n−3)²/4⌋` | 0 | 0 | 1 | 2 | 4 | 6 | 9 | 12 | 16 |

## 4. How far φ is from being a star flip — and the number that comes out

The natural quantitative form of §2 is the Hamming distance from `φ` to the star group:

```text
defect(n) := min_{g ∈ cut(H)} wt(1 + g) = m − maxcut(H).
```

**Lemma.** `maxcut(K_n ∖ P) = ⌊n²/4⌋ − 1` for `n ≥ 3`.

*Proof.* For `S` neither empty nor everything,
`|δ_H(S)| = |S||S^c| − |δ_P(S)| ≤ ⌊n²/4⌋ − 1`, because `P` is a spanning connected subgraph and so
has at least one edge crossing any nontrivial bipartition. Taking `S = {1,…,⌊n/2⌋}`, an interval of
the path, gives `|δ_P(S)| = 1` and `|S||S^c| = ⌊n²/4⌋`, attaining it. ∎

**Corollary.** `defect(n) = C(n−1,2) − ⌊n²/4⌋ + 1 = ⌊(n−3)²/4⌋`.

*(Parity check: `n = 2k` gives `(k−1)(k−2)`; `n = 2k+1` gives `(k−1)²`; both equal `⌊(n−3)²/4⌋`.)*

So **the antipodal map is a star flip plus exactly `⌊(n−3)²/4⌋` extra tiles, and no fewer.**

### 4b. The defect *set*: the tiles that do not cross the middle of the path

The number has an exact carrier. The maximising bipartition is the balanced **interval**
`S = {1,…,k}`, `k = ⌊n/2⌋`, so a minimum-weight representative of the coset `φ + cut(H)` is
supported on the non-crossing tiles:

```text
supp = E(H[S]) ⊔ E(H[S^c]),        |supp| = C(k−1,2) + C(n−k−1,2).
```

> **The obstruction to `φ` being a star flip is carried entirely by the tiles internal to the two
> halves of the base path** — the tiles that do not cross its midpoint. Those two tile sets are
> themselves the staircases of the two induced sub-tournaments.

Verified `n = 3…13` against brute-force max-cut: the interval bipartition is optimal at every `n`,
and `|supp|` matches `defect(n)` exactly.

### 4c. Why `defect(n) = h(n−2)` — and why that is *less* than it looks

THM-549 §3 gives the half-tiling size `h(n) = ⌊(n−1)²/4⌋`, so numerically

```text
defect(n) = ⌊(n−3)²/4⌋ = h(n−2).
```

It is tempting to read this as a **Mode-B** fact (`n → n−2` is exactly "both legs removed", and §5
identifies the legs as the two endpoint stars). **That reading is wrong, and I am recording the
deflation rather than the temptation.** §4b gives the honest explanation: with `a = k−1`,
`b = n−k−1` and `j = n−3` we have `a + b = j+1` and `|a−b| ≤ 1`, so

```text
C(a,2) + C(b,2) = ⌊j²/4⌋
```

is the classical **two-triangles-make-a-half-square** identity. Both sides of `defect(n) = h(n−2)`
are quarter-squares for the *same elementary reason* — "half of a triangle" and "two complementary
triangles" are both `A002620` — and the two tile sets are **not** naturally in bijection: the defect
set is two disjoint sub-staircases of the `n`-staircase, while `h(n−2)` is a folded half of the
single `(n−2)`-staircase. There is no Mode-B descent here, only `A002620` appearing twice. Logged
in the hypothesis index as raised-and-resolved in the same session.

## 5. The triangle's three sides, in star language

| side | tiles | star reading |
|---|---|---|
| vertical leg | `{(x,1)}`, `n−2` tiles | `star_H(1)` — a star-flip generator |
| horizontal leg | `{(n,y)}`, `n−2` tiles | `star_H(n)` — a star-flip generator |
| hypotenuse | `{x+y = n+1}`, `⌊(n−1)/2⌋` tiles | `Fix(σ)`, the `ρ`-mirror matching `{v, n+1−v}` |

and `σ` **swaps the two legs**, since `ρ(1) = n` (verified `n = 3…11`).

> **Two of the triangle's sides are cut generators; the third is the fixed set of the fold.**

This is the star-algebra reading of "everything is the triangle", and it is consistent with — and
explains the shape of — **THM-790's LEG LAW**, which says the `d=m` line (that is, `φ`) moves the
centered scores by `Δx = 8(e₁ − e_n)`, the difference of the two *leg* out-degrees. The legs are
exactly the two stars at the endpoints of the base path, so the leg law is the score-level shadow
of the fact that `φ`'s relation to the star group is controlled at the path endpoints. *(Consistency
observation, not a new derivation of THM-790.)*

## 6. The one-line summary, and what it corrects in intuition

> The quarter tiling is **not** "the half tiling folded again". `σ` and `φ` act on the invariant
> content by operations of different types — a **permutation** and a **translation** — and only the
> translation destroys information.

That asymmetry is invisible in the Klein-4 presentation `⟨σ, φ⟩`, where the two generators look
interchangeable, and it is the reason opus-S20's finding that *"the fold is a single terminal step,
no quarter-tiling"* has a companion on the other generator: the `σ`-fold cannot be iterated because
`σ` is already an involution with a large fixed set, while the `φ`-fold cannot be iterated because
it is **free** — for a completely different reason. Two different obstructions, one Klein group.

## 7. Scope

Structure only. Nothing here touches `H`-values, iso classes, LRC, or metagraph metrics — and in
particular THM-1382 §"Both open questions ANSWERED" and THM-1405 §I still apply: **holonomy is a
tiling invariant, not a tournament invariant.** Everything above therefore lives at the tiling
level, below iso classes, and inherits that ceiling. What it buys is that the two folds the repo
already uses are now placed exactly within the star algebra, with the cost of each fold named.

*Files: `04-computation/halfquarter_starsplit_klein_S335.py` (+ `.out` in `05-knowledge/results/`).*
