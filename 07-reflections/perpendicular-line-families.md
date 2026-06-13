# The Two Perpendicular Line Families of the Staircase

**Status:** CENTRAL INSIGHT  
**Session:** opus-2026-04-03-S27  
**Extends:** the-hypotenuse-generates-everything.md  

## The Two Coordinates

The staircase tile (x, y) lives in a coordinate system with two natural projections:

| Coordinate | Formula | Geometric meaning | Tournament meaning | Under reflection R |
|-----------|---------|-------------------|-------------------|-------------------|
| **x + y** | Sum | Parallel to hypotenuse | Vertex-pair label sum | **SWAPS** (d ↔ -d) |
| **x - y** | Difference (= skip+1) | Perpendicular to hypotenuse | Arc range in base path | **INVARIANT** |

These are **literally perpendicular** in the (x,y) plane, and **mathematically independent** under the grid reflection.

## What Each Family Controls

### Family 1: x + y = c (parallel to hypotenuse)

- The hypotenuse at c = n+1. Distance d = |c - (n+1)|.
- Source-side layers (c < n+1): arcs involving "upper" vertices
- Sink-side layers (c > n+1): arcs involving "lower" vertices
- R swaps source-side ↔ sink-side at equal distance from hypotenuse

**This family controls the SC/NS classification and blue/black edge coloring.**
A grid-symmetric tiling has identical bits on each paired x+y layer.

### Family 2: x - y = c (perpendicular to hypotenuse)

- Each value c = skip + 1 defines a "skip class" of tiles
- Skip s has n - s - 2 tiles (most tiles at skip 1, single tile at skip n-3)
- R maps each skip class TO ITSELF (x-y is invariant under R)

**This family controls the wiggly class structure and self-loop rates.**
Flipping a tile within a skip class = a wiggly move at that skip level.

## The Critical Perpendicularity

The grid reflection R:
- **Swaps** the x+y layers (complement symmetry)
- **Preserves** the x-y layers (skip structure)

Therefore: **SC/NS and wiggly are independent directions.**

- Blue/black lives along x+y (hypotenuse distance)
- Wiggly structure lives along x-y (skip)
- They don't interfere

## Which Skip Classes Have Hypotenuse Tiles?

A skip class c has a tile on the hypotenuse (x+y = n+1) iff (n+1+c)/2 is an integer and gives a valid tile. This requires c and n+1 to have the same parity.

- **Odd n**: odd skips (1, 3, 5, ...) have hypotenuse tiles
- **Even n**: even skips (2, 4, 6, ...) have hypotenuse tiles

This is the **root of the odd/even dichotomy** in metagraph structure.

The number of skip classes with hypotenuse tiles = floor((n-1)/2) = fixed tile count.

## Skip-1 Anomaly

Skip-1 tiles (x-y=2) have the lowest wiggly self-loop rate (~1% at n=7 vs ~4% for higher skips). This is because:

- Skip-1 tiles span the FULL width of the triangle (from corner to corner)
- They cross ALL x+y layers, including the extreme corners far from the hypotenuse
- Flipping a corner tile causes maximum structural disruption
- High-skip tiles are concentrated near the hypotenuse center and are more "protected"

## The Tile Count Triangle

Tiles per skip value (= x-y-1):

```
n\s    1    2    3    4    5    6
3      1
4      2    1
5      3    2    1
6      4    3    2    1
7      5    4    3    2    1
8      6    5    4    3    2    1
```

This IS the staircase itself, viewed from the perpendicular direction. The skip structure is literally the triangle read along its other axis.

## For Future Sessions

When analyzing any metagraph phenomenon, first ask:
1. Does it live along x+y (complement direction) or x-y (skip direction)?
2. If x+y: it's about SC/NS, blue/black, purity, the hypotenuse
3. If x-y: it's about wiggly classes, self-loop rates, arc ranges
4. If both: it's about the interaction of the two perpendicular structures
