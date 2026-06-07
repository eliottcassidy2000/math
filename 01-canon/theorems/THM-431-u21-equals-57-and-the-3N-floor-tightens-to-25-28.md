---
id: THM-431
name: u(21)=57 (proven) and the 3N-beating floor tightens to N* in [25,28]
status: VERIFIED (external proof cited for u(n); the N*-interval synthesis is an
        elementary rigorous consequence proved here)
date: 2026-06-07
session: monad-explorer-2026-06-07-S710
depends_on:
  - THM-421   # unit-distance 3N common-neighbour floor: N* in [17,32]
  - THM-412   # density quantization of the 2D unit-distance layer
  - HYP-2267  # triangular sqrt(7) Eisenstein construction line
  - HYP-2285  # exact value of N* in [17,32] (this sharpens it)
external:
  - "Alexeev, Mixon, Parshall, 'The Erdos unit distance problem for small point
     sets', arXiv:2412.11914 (2024)."   # exact u(n) for n<=21, upper bds 22..30
  - "Schade (1993); Engel ('Moser lattice'), lower bounds n<=30 (via AMP Table 1)."
---

# THM-431: u(21) = 57, and N* (first N beating 3N) lies in [25, 28]

Let `u(N)` = the maximum number of **unit distances** among `N` points in the plane
(the Erdos unit-distance maximum; OEIS A186705). This resolves the dispatched
n=21 campaign and sharpens THM-421.

## (A) u(21) = 57  [VERIFIED — external proof]

Alexeev–Mixon–Parshall (arXiv:2412.11914, Theorem 1) prove **u(n) exactly for all
`n ≤ 21`** by matching the Schade lower bounds with new computer-assisted upper
bounds (forbidden-subgraph `F`-free enumeration + a custom unit-distance embedder
that sidesteps cylindrical algebraic decomposition), and they **fully enumerate**
the densest graphs. In particular

```
        u(21) = 57.
```

Before AMP the best bounds were `57 ≤ u(21) ≤ 68`. The exact small-N table:

```
 n   16  17  18  19  20  21
u(n) 41  43  46  50  54  57
```

### The extremal n=21 graph is a Cartesian product (clean structure)

AMP's first n=21 extremal embedding is the **generic-angle Minkowski sum of a unit
equilateral triangle `T = K₃` and a unit wheel `W = W₇`** (hub + `C₆`, i.e. the
Eisenstein rosette/flower: 6 spokes + 6 rim = 12 edges). At a generic relative
angle the Minkowski sum realises the **graph Cartesian product `K₃ □ W₇`** as a
unit-distance graph, so

```
   u(21) = e(K₃ □ W₇) = e(T)·n(W) + n(T)·e(W) = 3·7 + 3·12 = 21 + 36 = 57.
```

(21 distinct vertices `= 3·7`; no edge collisions at generic angle — the standard
Erdos product construction.) **Note vs canon:** the S630 reflection's reading
`57 = 20 + 37` (Hamiltonian-spine + centered-hex bulk) is *not* the structure of
the proven extremal graph; the correct split is the product split `57 = 21 + 36`
(`e(T)n(W) + n(T)e(W)`). Most AMP extremal graphs `n ≤ 21` live in Engel's "Moser
lattice"; the lone exception is one `n = 16` graph.

## (B) The 3N-beating floor: N* ∈ [25, 28]  [PROVED here]

THM-421 defined `N*` = the smallest `N` with `u(N) > 3N` (average degree `> κ = 6`)
and proved `N* ∈ [17, 32]` (combinatorial cherry/common-neighbour floor 17;
`√7` Eisenstein construction at 32). AMP's table tightens **both** ends.

AMP Theorem 1(b) gives proven bounds for `22 ≤ n ≤ 30`:

```
 n     22  23  24  25  26  27  28  29  30
u≥     60  64  68  72  76  81  85  89  93     (Schade/Engel constructions)
u≤     61  66  72  78  84  90  96 103 110     (AMP upper bounds)
3n     66  69  72  75  78  81  84  87  90
```

**Floor `N* ≥ 25`.** For every `n ≤ 24`, the proven *upper* bound satisfies
`u(n) ≤ 3n`:
 - `n ≤ 21`: exact, and `u(n) < 3n` (largest is `u(21)=57 < 63`);
 - `n = 22`: `u(22) ≤ 61 < 66`;
 - `n = 23`: `u(23) ≤ 66 < 69`;
 - `n = 24`: `u(24) ≤ 72 = 3·24` (so not *strictly* greater).
Hence **no set of `≤ 24` points beats `3N`**: `N* ≥ 25`.

**Ceiling `N* ≤ 28`.** The Engel/Schade lower bound `u(28) ≥ 85` is *realizable*
(an explicit 28-point unit-distance graph with 85 edges), and `85 > 84 = 3·28`.
Hence **28 points beat `3N`**: `N* ≤ 28`.

```
        ==>   N* ∈ [25, 28].          (was THM-421: [17, 32])
```

This shrinks the THM-421 / HYP-2285 gap from 16 candidate values to **4**.

## (C) The √7 Eisenstein lane is the WRONG family for N*

Independent exact-integer search in the `√7` Eisenstein graph
(`unit_distance_3n_floor_sharpen_s710.py`, all counts via the 12-vector `Q=7`
shell) finds the densest `k`-subsets first exceed `3k` only at **`k = 39`**
(disk, `U = 119 > 117` — exact recount confirmed; the repo's THM-421 anneal does
better, `k = 32`, `U = 97`, but that is still inside the `√7` family). Engel's
**Moser-lattice** constructions beat `3N` already at `n = 28`. So the `√7` disk
family that THM-421(B) used to *upper-bound* `N*` is **not** where the true
first-beating happens — the Moser lattice is strictly better at the small-`N`
frontier. THM-421's "`N* ≤ 32` via `√7`" is correct but loose; AMP/Engel give
`N* ≤ 28`.

## Open (handed to HYP-2298)

Pin `N*` exactly in `[25,28]`. The sharpest target is **`n = 27`**, where the best
known construction *ties* (`u(27) ≥ 81 = 3·27`) but it is unknown whether
`u(27) > 81`. Either a construction beating `3N` at `n ∈ {25,26,27}` (lowers the
ceiling) or an upper bound `u(n) ≤ 3n` for `n ∈ {25,26,27}` (raises the floor)
would settle it. AMP's deficits to `3n` are: `25:[−3,+3]`, `26:[−2,+6]`,
`27:[0,+9]`, so the lower bounds run `−3, −2, 0` — **closing to a tie exactly at
27** before the construction breaks through at 28.

## Relation to canon

Consistent with and **strictly sharpens** THM-421 (no contradiction: `[25,28] ⊂
[17,32]`; no court case needed). The combinatorial floor 17 (cherry/Shrikhande,
THM-421(A)) is unaffected — it is the *design-theoretic* floor; `[25,28]` is the
*Euclidean-realizable* floor, so the THM-421 "cost of realizability" gap `17→32`
is now the narrower `17→25..28`. Credit: AMP (arXiv:2412.11914), Engel, Schade.
