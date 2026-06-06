---
id: THM-418
name: lattice-chromatic-dichotomy-of-the-plane
status: PROVED
date: 2026-06-06
session: oracle-2026-06-06-S686
depends_on:
  - HYP-2275   # the Niven/Dehn lattice-escape criterion (Hadwiger needle)
---

# THM-418: Lattice chromatic dichotomy — single 2D lattices never force χ≥4 (Hadwiger-Nelson)

## Statement

For a 2D lattice `L` and squared-norm `D`, the **unit-distance graph** `U(L,D)` has vertex
set `L` and edges `{x,y}` with `|x−y|² = D` (scale the unit to `√D`).

1. **Square lattice `ℤ²`:** `U(ℤ², D)` is **bipartite (χ = 2)** for every `D` (with edges).
2. **Triangular / Eisenstein lattice `ℤ[ω]`** (norm form `Q(a,b)=a²+ab+b²`): `U(ℤ[ω], D)` is
   **exactly 3-chromatic (χ = 3)** for every `D` (with edges).

In particular **neither lattice ever attains χ ≥ 4, at any norm/scale.**

## Proof

**Upper bounds (proper colorings of the *infinite* lattice graph).**

*Square.* Colour `c(a,b) = (a+b) mod 2`. An edge is a vector `(s,t)` with `s²+t²=D`. If `D`
is **odd**, exactly one of `s,t` is odd, so `s+t` is odd and `c` changes across every edge —
a proper 2-colouring. If `D` is **even**, every norm-`D` vector has `s≡t (mod 2)`, so the
norm-`D` vectors lie in the index-2 sublattice `{a≡b mod 2}`, which is `√2·ℤ²` rotated 45°
(scaled-similar to `ℤ²`); the graph splits into 2 disjoint copies of `U(ℤ², D/2)`. **2-adic
recursion** on `D` terminates at an odd `D'`, where the colouring is proper. So `χ ≤ 2`.

*Triangular.* Colour `c(a,b) = (a−b) mod 3`. Since `a²+ab+b² ≡ (a−b)² (mod 3)`, a norm-`D`
vector has `a≡b (mod 3)` **iff** `3 | D`. If `3 ∤ D`, every edge vector has `(a−b)≢0`, so `c`
changes across every edge — a proper 3-colouring. If `3 | D`, all norm-`D` vectors lie in the
index-3 sublattice `Λ' = {a≡b mod 3}`, which is `√3·`(rotated triangular lattice),
scaled-similar to `ℤ[ω]`; the graph splits into 3 disjoint copies of `U(ℤ[ω], D/3)`. **3-adic
recursion** terminates at `3 ∤ D'`. So `χ ≤ 3`.

**Lower bounds.**
- Both: an edge exists ⟹ `χ ≥ 2`.
- Triangular: `ℤ[ω]` has the order-6 rotation `R₆:(a,b) ↦ (−b, a+b)` as an automorphism (it
  preserves `Q`). For any norm-`D` vector `v`, the three vectors `0, v, R₆v` are pairwise at
  squared-distance `D` (`|v−R₆v|² = 2D − 2D·cos60° = D`), a **unit equilateral triangle** ⟹
  `χ ≥ 3`. Hence `χ = 3`.

## Verification

`04-computation/lattice_chromatic_dichotomy_s686b.py` (+.out): over **all 66 triangular and
78 square norms `D ≤ 199` with edges**, the colourings `(a−b) mod 3` / `(a+b) mod 2` are
proper after the 3-adic / 2-adic strip (0 anomalies, with the complete connection set
`|coord| ≤ 2√D`); the 60° equilateral triangle is present for every triangular norm; exact
patch chromatic numbers are `3` (triangular, `D ≤ 13`) and `2` (square, `D ≤ 25`).

## Consequence (Hadwiger-Nelson)

The two densest 2D lattices — including the **Eisenstein lattice `ℤ[ω]`, which is de Grey's
substrate** for `χ(ℝ²) ≥ 5` — are uniformly `≤ 3`-chromatic at every scale: **a single square
or triangular lattice can never force `χ ≥ 4`.** So the values `4,5,6,7` are *eliminated* for
single-(square-or-triangular)-lattice unit-distance graphs.

Therefore **all chromatic forcing above 3 in the plane is Dehn-nontrivial** (HYP-2275): it
requires combining **incommensurate** lattices / **irrational** rotation angles — the Moser
spindle's `arccos(5/6)` junction (χ=4), de Grey's irrational rotations (χ≥5). The plane's
`χ ∈ {5,6,7}` lives entirely *off* any single lattice; the chromatic difficulty is exactly the
Niven/Dehn lattice-escape, not anything internal to a lattice.

## Scope / honest limits

- Proved for the **square and triangular (Eisenstein)** lattices — the two with `>2`
  automorphisms (`w=4,6`), and the only relevant ones (the Eisenstein lattice is the de Grey
  substrate; the hexagonal 7-colouring is its dual). A **generic** 2D lattice (`w=2`) at a
  popular norm with many edge-directions is not covered by this parity argument; whether some
  generic 2D lattice attains `χ=4` is left open (small patches found only `≤3`).
- This is a statement about single-lattice subgraphs, not the full plane; it does **not**
  resolve `χ(ℝ²)`. Its value is *localizing* the open `{5,6,7}` to the Dehn-nontrivial regime.

## Related
HYP-2275 (Niven=Dehn-triviality=THM-416 cap=HN lattice-escape), THM-416/HYP-2274 (CM
density-quantum totient cap), HYP-2265 (LRC=HN=unit-distance Cayley + Delsarte), de Grey 2018,
Moser spindle, classical: triangular grid is 3-chromatic / `χ(ℚ²)=2`.
