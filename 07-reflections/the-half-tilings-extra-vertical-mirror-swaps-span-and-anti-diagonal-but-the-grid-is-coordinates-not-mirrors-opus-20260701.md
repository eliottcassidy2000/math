# The half-tiling's extra vertical mirror swaps span and anti-diagonal — but the "grid" is the coordinate lattice, not a family of mirror lines

*opus-2026-07-01-S21. The owner asked: are there meaningful operations on the half-tiling model that aren't
proper folds; does it have perpendicular parallel lines of symmetry along its diagonals; what grid do they form?
Computed the exact symmetry group and characterized the extra operation.*

## The two single reflections (both regions are D₁)
Computed the full isometry symmetry group of the tile point-set (n=4..9):
- **Full staircase**: symmetry `= ⟨R⟩`, a single reflection across the **anti-diagonal `x+y=n+1`** — this is
  `R = the complement fold` (S18), and it **fixes the span** `d=x−y` (reflects the anti-diagonal `s=x+y`).
- **Half-region**: symmetry `= ⟨V⟩`, a single reflection across the **vertical `x=(n+3)/2`**.

**Why the half gains V:** in the half-region, row `y` runs `x ∈ [y+2, n+1−y]` (staircase edge to anti-diagonal),
whose midpoint `((y+2)+(n+1−y))/2 = (n+3)/2` is **independent of `y`**. So every row is centered on the same
vertical line, and `V : x ↦ (n+3)−x` is a reflection — one the *full* staircase lacks (there the row-center
`(y+2+n)/2` drifts with `y`). This is the concrete "extra symmetry" of the folded region.

## What V is — and isn't
`V` is a genuine **involution of the half-tiling cube** (`V²=id`, it permutes the `2^D` half-tilings). But:
- **It is only a GEOMETRIC symmetry — it does NOT respect tournament isomorphism** (verified n=5,6,7: two
  half-tilings of the same SC class can map under `V` to *different* classes). So `V` does **not** descend to the
  SC-class metagraph; it is **not** a fold and **not** a relabeling.
- **It is not the complement** (which is trivial on the SC world). `V` genuinely reshuffles the SC tournaments
  as tile-configurations.

So the answer to part (a): yes, there are meaningful operations beyond the complement-fold — `V` is the cleanest
— but `V` lives at the *tile/geometry* level, not the *iso-class* level. (Other non-fold operations, all
tile-level: span/range flips (flip all tiles of one span `d`), wiggly flips (single tile), vertex-star flips;
and one iso-level non-fold operation, the involutive **anti-automorphism** `τ` of each SC tournament.)

## R and V in span/anti-diagonal coordinates — R fixes span, V swaps the axes
Put `s = x+y` (anti-diagonal) and `d = x−y` (**span** — the arc's range). Then
```
    R : (s,d) → (2n+2−s,  d)        reflect the anti-diagonal, FIX the span
    V : (s,d) → (n+3−d,  n+3−s)      SWAP span ↔ anti-diagonal
```
So the two reflections are exactly the two natural mirrors of the `(s,d)` grid: `R` reflects along the
anti-diagonal direction keeping span fixed; `V` is the *transpose* mirror that exchanges the span and
anti-diagonal roles. In `(x,y)` they sit at **45° to each other** (R across slope −1, V vertical) — **not
perpendicular**.

## The grid (part b) — it's the coordinate lattice, not a family of mirror lines
The tiles occupy a lattice with **two perpendicular families of parallel diagonal lines**:
- constant **span** `d = x−y ∈ {2,…,n−1}` (slope +1),
- constant **anti-diagonal** `s = x+y` (slope −1),
on the **checkerboard sublattice** `s ≡ d (mod 2)` (since `s+d = 2x`). This *is* the grid the question senses —
perpendicular, parallel-within-family — and it is exactly the `(x,y)` integer lattice rotated 45° and scaled by
`√2`, restricted to the triangle `d ≥ 2` (and, for the half, `s ≤ n+1`).

**But these are coordinate lines, not lines of symmetry.** The actual reflection symmetry is only **D₁** — one
axis per region. There is **no family of parallel mirror lines**: no span-reflection (spans have unequal
multiplicities), and the two mirrors that do exist (`R`, `V`) are single and at 45°. So the honest answer is:
*the perpendicular parallel diagonal families are the span×anti-diagonal coordinate grid; the symmetry is a
single mirror per region — `R` running along an anti-diagonal grid line (span-fixing), and `V` the transpose
mirror swapping the two grid directions.*

## Status
- **Established (n=4..9):** full staircase symmetry `= D₁` (`R` = complement, across `x+y=n+1`, span-fixing);
  half-region symmetry `= D₁` (`V`, vertical `x=(n+3)/2`, exists because the half's row-centers are constant).
- **Characterized:** `V` swaps span↔anti-diagonal in `(s,d)`; it is a geometric involution of the half-tiling
  cube that does **not** respect tournament-iso (not a fold/relabeling) — verified n=5,6,7.
- **The grid:** the span×anti-diagonal checkerboard lattice (two perpendicular parallel diagonal families); the
  mirrors are single (D₁), `R` and `V` at 45°, not a parallel family.
- **Honest scope:** geometric/tile-level facts, verified small n; `V` being iso-blind means it is *meaningful as
  a shape symmetry*, not a tournament invariant.

Related: HYP-3810 (the half-tiling model), HYP-3811 (single terminal fold — this adds the half's *extra* mirror
`V`, geometric-only, so it doesn't contradict "no metagraph quarter-fold"), S18 (R = complement, span-fixing),
"everything is the triangle" (span = the arc range; the two legs vs hypotenuse). Scripts:
04-computation/mmg_halftiling_{V_reflection, symmetry_grid, V_characterization}_opus_20260701.py.
