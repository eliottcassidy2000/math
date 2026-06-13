# THM-280: Grid Reflection Induces Tournament Complement

**Status:** PROVED  
**Session:** opus-2026-04-03-S27  
**Verified computationally:** n=3 through n=7  

## Statement

Let T be a tournament on [n] represented by a tiling b of the staircase delta_{n-2} with fixed base path n -> (n-1) -> ... -> 1. Let b' be the grid-reflected tiling, defined by

    b'_{(n+1-y, n+1-x)} = b_{(x, y)}    for every tile (x, y).

Let T' be the tournament encoded by b'. Then **T' is isomorphic to T^op** (the tournament obtained by reversing all arcs of T).

Specifically, T' = sigma(T^op) where sigma is the vertex permutation v -> n+1-v.

## Consequence

The grid reflection on the staircase tiling space induces the complement involution T -> T^op on isomorphism classes. Therefore:

1. Transpose-paired classes = complement-paired classes (G_n/Z_2 quotient).
2. Self-complementary tournaments (T ≅ T^op) = transpose-self classes.
3. The two merge operations (grid reflection and T^op) are identical at the class level.

## Proof

### Setup

The staircase delta_{n-2} has m = C(n-1, 2) cells (tiles). Each tile is labeled by a pair (x, y) with 1 <= y < x <= n and x - y >= 2 (non-consecutive vertices in the base path). The base path arcs are n -> (n-1) -> ... -> 1, fixed for all tilings.

A tiling b assigns a bit b_{(x,y)} in {0, 1} to each tile:
- b_{(x,y)} = 0: arc from x to y (higher to lower)
- b_{(x,y)} = 1: arc from y to x (lower to higher)

The tournament T(b) is the directed graph on [n] with:
- Arc k -> (k-1) for k = n, ..., 2  (base path, always present)
- Arc x -> y if b_{(x,y)} = 0, or arc y -> x if b_{(x,y)} = 1  (for each tile)

### The grid reflection

The grid reflection R maps tile (x, y) to tile (n+1-y, n+1-x).

**Claim 1:** R is a well-defined involution on the set of tiles.

Proof: If (x, y) is a tile (meaning 1 <= y < x <= n, x - y >= 2), then the image (n+1-y, n+1-x) satisfies:
- 1 <= n+1-x < n+1-y <= n  (since x > y)
- (n+1-y) - (n+1-x) = x - y >= 2
So (n+1-y, n+1-x) is also a tile. Applying R twice: (x,y) -> (n+1-y, n+1-x) -> (n+1-(n+1-x), n+1-(n+1-y)) = (x, y). []

The reflected tiling b' is defined by: b'_{R(x,y)} = b_{(x,y)}, i.e., b'_{(n+1-y, n+1-x)} = b_{(x,y)}.

### The vertex permutation

Define sigma: [n] -> [n] by sigma(v) = n + 1 - v.

Note: sigma reverses the vertex ordering. It is an involution (sigma^2 = id). It maps the base path n -> (n-1) -> ... -> 1 to the reversed path 1 -> 2 -> ... -> n.

### Tile arcs: T' = sigma(T^op) for non-consecutive pairs

Consider a tile (x, y) with x - y >= 2. In T, the arc between x and y is determined by bit b_{(x,y)}.

**Case b_{(x,y)} = 0:** T has arc x -> y.

The reflected tiling has b'_{(n+1-y, n+1-x)} = 0, so T' has arc (n+1-y) -> (n+1-x), i.e., sigma(y) -> sigma(x).

In T^op: since T has x -> y, T^op has y -> x.
In sigma(T^op): sigma maps this to sigma(y) -> sigma(x).

So T' has sigma(y) -> sigma(x), and sigma(T^op) has sigma(y) -> sigma(x). They agree. []

**Case b_{(x,y)} = 1:** T has arc y -> x.

The reflected tiling has b'_{(n+1-y, n+1-x)} = 1, so T' has arc (n+1-x) -> (n+1-y), i.e., sigma(x) -> sigma(y).

In T^op: since T has y -> x, T^op has x -> y.
In sigma(T^op): sigma maps this to sigma(x) -> sigma(y).

Again T' and sigma(T^op) agree. []

### Base path arcs: T' = sigma(T^op) for consecutive pairs

For consecutive vertices k and k-1 (with 2 <= k <= n):

In T: arc k -> (k-1) (the base path, always present).
In T': the base path is unchanged: arc k -> (k-1).

In T^op: T has k -> (k-1), so T^op has (k-1) -> k.
In sigma(T^op): sigma maps this to sigma(k-1) -> sigma(k) = (n+2-k) -> (n+1-k).

We need T'(k, k-1) = sigma(T^op)(k, k-1) for all consecutive k.

T' has arc k -> (k-1) (base path).
sigma(T^op) at vertices k and k-1: we need the arc between sigma^{-1}(k) = n+1-k and sigma^{-1}(k-1) = n+2-k in T^op, then relabel.

Actually, sigma(T^op)(u, v) = T^op(sigma(u), sigma(v)) = T(sigma(v), sigma(u)).

So sigma(T^op)(k, k-1) = T(sigma(k-1), sigma(k)) = T(n+2-k, n+1-k).

Since n+2-k and n+1-k are consecutive vertices (they differ by 1), and n+2-k > n+1-k, the base path gives T(n+2-k, n+1-k) = 1 (arc present from higher to lower).

And T'(k, k-1) = 1 (base path arc present).

So sigma(T^op)(k, k-1) = T'(k, k-1) = 1 for all consecutive pairs. []

### Conclusion

For all vertices u, v in [n]:

    T'(u, v) = sigma(T^op)(u, v) = T^op(sigma(u), sigma(v))

Since sigma is a permutation, T' and T^op differ only by a vertex relabeling, hence **T' ≅ T^op**.

At the isomorphism class level, the grid reflection R on tilings induces the complement involution T -> T^op on tournament classes. QED.

## Remarks

1. The proof works for ALL n >= 3. There are no restrictions.

2. The key insight: the grid reflection (x,y) -> (n+1-y, n+1-x) on tile coordinates encodes the vertex permutation v -> n+1-v, which is the natural "reversal" involution. This reversal simultaneously:
   - Reverses the base path direction (handled by the base path being fixed)
   - Swaps the two legs of the staircase triangle (source column <-> sink row)
   - Produces T^op up to the relabeling sigma

3. The grid reflection is the UNIQUE involution of the staircase that corresponds to T^op. Other reflections of the triangle (e.g., reflecting along the hypotenuse) would correspond to different operations.

4. Self-complementary tournaments (T ≅ T^op) are exactly those whose isomorphism class is invariant under the grid reflection. A tiling in such a class is grid-symmetric iff it is a fixed point of R within its class fiber.
