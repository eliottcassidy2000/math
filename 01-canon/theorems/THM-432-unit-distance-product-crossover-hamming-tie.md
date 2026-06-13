---
id: THM-432
name: the Erdos-product 3N-crossover; the n=27 unit-distance tie is the Hamming graph H(3,3)
status: PROVED (product formula + criterion, classical) + VERIFIED exact (Q(sqrt3),
        explicit algebraic angles); uses cited proven u(a) values (AMP 2412.11914)
date: 2026-06-07
session: monad-explorer-2026-06-07-S711
depends_on:
  - THM-431   # u(21)=57; N* in [25,28]; the tie at n=27
  - THM-421   # common-neighbour 3N floor (N* in [17,32])
  - THM-412   # density quantization of the unit-distance layer
  - HYP-2170  # UD graph = cyclotomic Cayley graph; "3" = kissing/2 = Eisenstein rosette
external:
  - "Alexeev, Mixon, Parshall, arXiv:2412.11914 (2024) — exact u(n), n<=21."
  - "Erdos (1946) — the Minkowski-sum / Cartesian-product unit-distance construction."
verified_by: 04-computation/unit_distance_product_crossover_monad_s711.py
---

# THM-432: The Erdos product and the n=27 tie = the Hamming graph H(3,3)

Let `u(N)` = the Erdos unit-distance maximum (OEIS A186705) and `N*` = the smallest
`N` with `u(N) > 3N` (average degree `> kissing number κ = 6`). THM-431 proved
`N* ∈ [25,28]` and that the best known construction **ties** (`u(27) ≥ 81 = 3·27`),
leaving open *why* the deficit closes to a clean tie exactly at `n = 27 = 3³`. This
theorem answers that, and shows the first crossing is **not** a product.

## (A) The product formula and the 3N criterion [PROVED — classical]

For finite planar unit-distance graphs `G, H` placed at a **generic relative angle**,
the Minkowski sum `G ⊕ H = {g + h}` realises the graph **Cartesian product `G □ H`**
as a planar unit-distance graph (Erdos 1946): the `n(G)n(H)` sum-points are distinct,
and the unit distances are exactly the "axis-parallel" ones — `g~g'` in `G` (fixed
`h`) and `h~h'` in `H` (fixed `g`) — with no diagonal coincidences at generic angle.
Hence

```
   n(G □ H) = n(G)·n(H),     e(G □ H) = e(G)·n(H) + n(G)·e(H).
```

Writing the **edge density** `ρ(G) = e(G)/n(G) = ½·avgdeg(G)`, the product density is
`e(G□H)/n(G□H) = ρ(G) + ρ(H)`. Since "beats 3N" means density `> 3` (avg degree `> κ = 6`):

```
   G □ H beats 3N   ⟺   ρ(G) + ρ(H) > 3   ⟺   avgdeg(G) + avgdeg(H) > κ = 6.
                  (TIE iff ρ(G)+ρ(H) = 3, i.e. the avg degrees sum to exactly κ.)
```

The kissing cap `κ` reappears as the additive budget for the two factors' densities.

## (B) The n = 27 tie IS the Hamming graph H(3,3) = K₃^□3 [VERIFIED exact]

`u(9) = 18` is realised by `K₃ □ K₃` (two unit triangles): `n = 9`, `e = 3·3+3·3 = 18`,
**4-regular**. Iterating, the `n = 27` tie is the 3-fold product of unit triangles

```
   K₃ □ K₃ □ K₃  =  K₃^□3  =  Hamming graph H(3,3)  =  the 3×3×3 rook's graph,
   n = 3³ = 27,   e = 81,   6-REGULAR,   so  e = 81 = 3·27  EXACTLY.
```

Each vertex has exactly `2` neighbours per coordinate × `3` coordinates = degree
`6 = κ`. So the construction hits the kissing threshold **on every vertex
simultaneously** — it ties `3N`, and **cannot beat it**, precisely because a product
of triangles is forced 6-regular (no degree can exceed `κ`). This is the structural
reason the THM-431 deficit closes to a clean `+0` at `n = 27 = 3³`: the "`3³`" is
literally `K₃^□3`, and the "tie" is its 6-regularity.

**Exact verification** (`unit_distance_product_crossover_monad_s711.py`, all
coordinates in `Q(√3)`, rotations by exact Pythagorean unit complex numbers
`(3/5,4/5),(5/13,12/13),…` so a squared distance is unit iff it equals exactly `1`):

| product | N | U | 3N | regularity |
|---|---|---|---|---|
| `K₃□K₃` (= u(9) max) | 9 | 18 | 27 | 4-regular |
| `W₇□K₃` (= AMP's proven u(21) extremal `K₃□W₇`) | 21 | **57** | 63 | deg 5..8 |
| `K₃□K₃□K₃ = H(3,3)` | 27 | **81** | 81 | **6-regular**, `U = 3N` |
| `W₇□W₇` (a product that beats 3N) | 49 | 168 | 147 | deg 6..12 |

The `W₇□K₃ = 57` row independently **reproduces AMP's proven `u(21) = 57` extremal
graph** as an exact `Q(√3)` point set, validating the product machinery against a
known optimum.

## (C) N* is NOT a product — it is a rigid 2D blob [VERIFIED]

Census over all factorisations `N = a·b` using the **proven** `u(a)` values (AMP, exact
`a ≤ 21`; best factors = `u`-maximisers, since `e(G□H)` is increasing in each `e`):

```
   smallest product TIE  (U = 3N):   N = 27  (= K₃^□3 = H(3,3))   and   N = 30 (3×10)
   smallest product BEAT (U > 3N):   N = 32  (= K₂ □ G₁₆,  U = 16 + 2·41 = 98 > 96)
```

Since the true `N* ∈ [25,28]` and `28 < 32`, **no generic Cartesian product attains
the first crossing**: the optimal `N*`-construction is necessarily a non-product,
irregular **blob** (some vertex of degree `> κ`), consistent with AMP/Engel's extremal
graphs living in the "Moser lattice" rather than as products. The product is the
*decomposable / maximally-symmetric* rung of the construction ladder; it saturates
`3N` (tie) but the average degree of a product is dragged below `κ` by each factor's
low-degree boundary, so it cannot cross.

**Tightness at 27.** The best product per `n` (n=22..30) gives deficits to `3n` of
`−9,·,−6,−5,−5,+0,−1,·,+0` (· = prime `n`, no product), vs AMP's `−6,−5,−4,−3,−2,
+0,+1,+2,+3`. The product is **tight with the global best exactly at `n = 27`**
(both `81`) and loses by only `1–3` elsewhere — the small gap is the edge the
irregular blob buys by concentrating degree past `κ` on interior vertices.

## Consequence for OPEN-Q-057 / N*

Strong structural evidence (not a proof) that `u(27) = 81` (hence `N* = 28`): the most
efficient *structured* 27-point construction, `H(3,3)`, is forced 6-regular and ties;
beating `81` would require an irregular 27-point blob, and AMP's *upper* bound at 27 is
only `90`. The honest open question is whether such a degree-concentrating 27-blob
exists. To **lower** the ceiling one would seek an 82-edge 27-point unit-distance graph
(necessarily non-product); to **raise** the floor, an upper bound `u(27) ≤ 81`.

## Relation to canon

Sharpens the THM-431 interpretation of the `n=27` tie (now identified as `H(3,3)`).
Consistent with THM-421/THM-412 (the `κ = 6` budget governs both the floor `C(κ,2)+2`
and the threshold `(κ/2)N`; here it is the *additive density budget* `ρ(G)+ρ(H) ≤ 3`).
No contradiction, no court case. Credit: AMP (arXiv:2412.11914), Erdos (product).
The reflection
`07-reflections/symmetry-saturates-irregularity-violates-the-hamming-tie-s711.md`
develops the cross-domain pattern (unit-distance ↔ LRC worry-set).
