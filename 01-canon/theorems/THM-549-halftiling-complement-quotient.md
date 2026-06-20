---
id: THM-549
title: The half-tiling is the tiling-model fundamental domain of the complement involution — the y=x reflection σ(x,y)=(n+1-y,n+1-x) of a tiling equals φ(T^op) (reverse all arcs, relabel i→n+1-i), with ⌊(n−1)²/4⌋ cells, a ⌊(n−1)/2⌋-cell fixed diagonal (the SC spine), and even/odd size recursions A+B−C / A+B−C+D−E−F+G
status: PROVED (the symmetry σ(tiling)=φ(T^op) verified EXHAUSTIVELY n=3..6 over all 2^m tilings; sizes/fixed-count/recursions verified n=2..12). Connects the user's S4 half-tiling framework to the complement quotient G_n/Z_2 and the one/two/three-far Newton recursion (codex HYP-2680/2681).
source: mac-mini-2026-06-20-S4
depends_on: []
prior_art:
  - THM-280   # reflection = complement (pre-existing; THM-549 gives the tiling-coordinate σ=φ(T^op) form + exhaustive proof)
  - THM-442   # full-tiling third-difference recursion (the half-tiling parity recursions refine this)
related:
  - THM-548   # one/two/three-far boundary-value decomposition (the coverage incarnation of the 7-term recursion)
  - HYP-2680  # codex Phi_s Stirling multi-far hierarchy
  - HYP-2681  # codex cube-root / Eisenstein A+B+C-D-E-F+G = pair-tax shadow
external: Rédei tournaments; complement (converse) involution; A002620 (quarter-squares).
---

# THM-549 — The half-tiling = complement-quotient fundamental domain

## Setup
Tiling model: tournament on `{1..n}`, fixed base path `n→…→1` (arc `k→k-1`); **tiles** = arcs
`(x,y)`, `x>y`, `x-y≥2` (`m=C(n-1,2)`); a **tiling** orients each tile. Reflection
`σ(x,y)=(n+1-y, n+1-x)` (the grid map).

## Proved
1. **σ is an involution** on tiles; **fixed tiles** `= {x+y=n+1}`, count `⌊(n−1)/2⌋`.
2. **`σ(tiling) = φ(T^op)`** where `T^op` reverses all arcs (complement) and `φ(i)=n+1-i`
   relabels to restore the base path. VERIFIED exhaustively `n=3..6` (all `2^m` tilings).
   ⟹ the y=x reflection IS the complement involution in tiling coordinates.
3. **Half-region size** `= (m + ⌊(n−1)/2⌋)/2 = ⌊(n−1)²/4⌋` (A002620, quarter-squares) — the
   Burnside orbit count of σ. Full tiling is triangular `C(n-1,2)`; the half is square-ish.
4. The fixed diagonal `{x+y=n+1}` = the **self-complementary spine** (tiles invariant under
   complement) — the tiling-geometric realization of the metagraph SC spine (CLAUDE.md).
5. **Size recursions:** EVEN `h(n)=2h(n-1)-h(n-2)` (= A+B−C); ODD
   `h(n)=2h(n-1)-2h(n-3)+h(n-4)` (= A+B−C+D−E−F+G, the −C+D equal-size but geometrically
   distinct: C an edge-overlap, D a corner). Both satisfied by `⌊(n−1)²/4⌋`.

## Verified extensions (mac-mini-S4)
6. **Square/pronic shape (the even/odd "different shape"):** odd `n=2k+1` ⟹ half-tiling `= k²`
   (perfect SQUARE); even `n=2k` ⟹ `= k(k−1)` (PRONIC rectangle). The square (odd) admits the
   3-corner decomposition; the pronic (even) does not — exactly the user's "even tournaments do
   not produce half-tilings with 3 corners." (VERIFIED n=2..14.)
7. **Recurrence order:** the order-4 recurrence `h(n)=2h(n−1)−2h(n−3)+h(n−4)` (char poly
   `(x−1)³(x+1)`) holds for ALL n; EVEN n *additionally* satisfies the order-2 `h(n)=2h(n−1)−h(n−2)`.
   So even is degenerate (extra recurrence, pronic), odd is generic (order-4, square).
8. **Complement-invariance ⟹ 2× computation (engineering):** `c3(T)=c3(T^op)`, `HP(T)=HP(T^op)`
   for ALL tournaments (VERIFIED n=4,5,6, 100%); scores complement to `(n−1)−s`. Hence any
   complement-invariant invariant (cyclic-triangle count, Ham-path count, H, OCF, Walsh order-2)
   is computable by summing over the `⌊(n−1)²/4⌋`-cell half-region instead of the full
   `C(n−1,2)` — a literal 2× saving, the diagonal SC spine handled once. **The OCF (odd-cycle
   counts of ALL lengths) is complement-invariant too** (VERIFIED 100% n=4,5,6), so the CENTRAL
   Claim-A / μ computation also halves on the half-tiling — not just c3/HP.

## Significance
The half-tiling factors out the SAME Z_2 (complement = reverse-all-arcs) as the merged metagraph
`G_n/Z_2`, but in tiling coordinates where the involution is a literal mirror over `y=x` and the
SC spine is the literal diagonal. The 7-term odd recursion (3 corners + 3 edges + 1 center)
mirrors the one/two/three-far Newton packets (3 one-far + 3 two-far + 1 three-far, THM-548 /
codex HYP-2680/2681), with the 3-fold structure carried by the cube roots of unity (Eisenstein).
