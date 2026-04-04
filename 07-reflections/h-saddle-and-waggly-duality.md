# The H-Saddle and Waggly Duality

**Session:** opus-2026-04-03-S28

## The Multilinear Polynomial Structure of H

H(T) viewed as a function of tiling bits x₁,...,xₘ is a multilinear polynomial with remarkable structure:

### Linear coefficients (PROVED)
**THEOREM**: The linear coefficient of tile (x,y) in the polynomial is 2^(skip-1) where skip = x-y.

**Proof**: Flipping one arc from the transitive tournament creates exactly 2^(skip-1) new Hamiltonian paths. Each new path uses the flipped arc y→x, and the (skip-1) "between-vertices" {y+1,...,x-1} can each be visited either before y or after x.

This gives: apex tile (n,1) has coefficient 2^(n-2), which dominates all others.

### Quadratic coefficients (DISCOVERED + CORRECTED)
For tile pairs sharing a vertex v, skips s₁, s₂:
- **Same-end** (both upper or both lower endpoint = v): c₂ = -2^max(1, |s₁-s₂|-1)
  - Adjacent skips (|s₁-s₂| ≤ 2): c₂ = -2
  - Distant skips (|s₁-s₂| > 2): c₂ = -2^(|s₁-s₂|-1) — exponential interference
- **Cross-end** (v is upper for one, lower for other): c₂ = +2^(s₁+s₂-2)

Same-end pairs create destructive interference that grows with the skip gap. Cross-end pairs create "through-vertex shortcuts" whose cooperation grows exponentially with the total skip. Verified exhaustively n=4..9.

### Polynomial degree
- n=3: degree 1
- n=4: degree 2  
- n=5: degree 4
- n=6: degree 4
- n=7: degree 6

Max Walsh order grows but stays well below m. Orders ≥ 5 have exactly zero energy at n ≤ 6. At n=7, orders 5-6 appear with tiny coefficients (±1/32).

## The Grid Reflection Effect Symmetry

**THEOREM** (verified n=3..7): effect(tile(x,y)) = effect(tile(n+1-y, n+1-x))

The tile effects are perfectly symmetric under the grid reflection that maps the staircase to itself. Reflection pairs have **exactly** equal effects (zero difference, not approximate).

Tiles partition into:
1. **Boundary** (y=1 or x=n): large positive effects
2. **Interior near center**: effects near zero, some negative
3. **Center tile** (self-reflected): effect exactly 0

## The Apex Doubling Law

The apex effect E[H|apex=1] - E[H|apex=0]:
- n=3: +2 = 2^0
- n=4: +2 = 2^1
- n=5: +4 = 2^2
- n=6: +8 = 2^3
- n=7: +159/8 ≈ 19.875 (BREAKS the exact 2^(n-3) pattern)

Exact through n=6, with a specific rational departure at n=7.

## The Complement-Reflect Involution

**THEOREM** (verified n=3..7):
1. Tiling complement and grid reflection **commute** (100%)
2. complement∘reflect is an **involution** (100%)
3. complement∘reflect **preserves H** exactly (100%)
4. complement∘reflect **preserves the 3-cycle count** c₃ exactly (100%)

The CR involution = pure coordinate remap (x,y) → (n+1-y, n+1-x) with bits preserved. It is NOT vertex relabeling (which would reverse all arcs). Instead it's a "half-reversal" that remaps the staircase geometry while keeping the arc directions relative to the new coordinates.

## The Walsh-Waggly Duality

**THEOREM**: The variance of delta-H at Hamming distance d is:
  E[Δ²(d)] = Σ_k |ĥ_k|² · g(k,d)

where g(k,d) = 2(1 - E[(-1)^{|S∩F|}]) with |S|=k, |F|=d drawn uniformly.

This holds **exactly** (ratio 1.0000) at every (n,d) tested.

Key properties of g(k,d):
- g(1,d) = 4d/m: **linear in d** for order-1
- g(even, m) = 0: even Walsh orders are **invisible to complement**
- g(2,d) is **humped**: peaks at d ≈ m/2, returns to 0 at d=m

The waggly spectrum is **humped** (not U-shaped): mid-range Hamming distance d is most disruptive to H. The complement (d=m) is less disruptive than flipping ~2m/3 tiles.

## The Saddle Structure

In (hypotenuse_weight, periphery_weight) coordinates, H has a **saddle**: high at opposite corners (hyp-full, per-empty) and (hyp-empty, per-full), low at (all-0) and (all-1).

The saddle exists because:
- Boundary tiles (large skip) have positive linear coefficients
- Center tiles have zero or negative effects
- The interaction terms create a coupling where concentrating flips along ONE geometric axis maximizes H, but spreading them uniformly does not

## Connection to the Triangle

The staircase δ_{n-2} IS a right isosceles triangle. The three natural directions:
- **Hypotenuse** (x+y = const): SC/NS distinction lives here
- **Vertical leg** (y = const, "bottom row"): score hierarchy
- **Horizontal leg** (x = const, "right column"): complement hierarchy

The linear coefficient formula 2^(x-y-1) = 2^(skip-1) means: the **perpendicular distance from the base path** controls the individual tile's influence on H. Tiles near the base path (skip=2) contribute 2 paths. The apex tile (skip=n-1) contributes 2^(n-2) paths. This exponential growth IS the geometric content of the triangle's hypotenuse.

## Open Questions

1. Is there a closed form for the order-2 coefficient of DISJOINT tile pairs?
2. Why does the polynomial degree cap at 4 for n=5,6 but jump to 6 at n=7?
3. Can the multilinear polynomial be related to the OCF formula H = I(Ω,2)?
4. The apex effect 159/8 at n=7 — what is the correction to 2^(n-3)?
5. The CR involution preserves c₃ but not always scores. What exactly does it preserve?
