# The Tiling Picture

## A Different Way to See a Tournament

One of the most useful insights in this project is that **every tournament can be represented as a pattern of black and white squares in a staircase-shaped grid**.

This "tiling picture" transforms the complex, arrow-laden tournament diagram into something you can literally look at — a visual pattern. And it turns out many deep properties of the tournament correspond to visual properties of the pattern.

---

## The Staircase

Start with a staircase-shaped grid. For a tournament on n players, the staircase has:
- n-2 steps
- A total of (n-1)(n-2)/2 squares

For 4 players, the staircase looks like this (3 squares):

```
[ ]
[ ][ ]
```

For 5 players (6 squares):

```
[ ]
[ ][ ]
[ ][ ][ ]
```

Each square in the staircase corresponds to one **non-consecutive pair** of players in a fixed ordering.

---

## Fixing a Reference Path

Here's the key idea. We fix one particular Hamiltonian path as a **reference point** — a baseline. Let's call it P₀.

P₀ is just the players listed in a fixed order: Player n → Player n-1 → Player n-2 → ... → Player 1. In our tournament, this means Player n beats Player n-1, who beats Player n-2, and so on down to Player 1. These "consecutive" arrows along P₀ are fixed — they never change.

All the other arrows in the tournament — the ones between NON-consecutive players in P₀ — are the **free variables**. They can point either way. Each one gets one square in the staircase.

---

## Coloring the Squares

For each free arc (non-consecutive pair in P₀):
- If the arc goes in the **P₀ direction** (upper-ranking player beats lower-ranking player), color the square **white (0)**.
- If the arc goes **against P₀** (lower-ranking player beats upper-ranking player), color the square **black (1)**.

The result is a unique black-and-white pattern in the staircase — a **tiling** of the staircase.

Every black-and-white pattern corresponds to exactly one tournament that contains P₀ as a Hamiltonian path.

---

## What the Tiling Tells You

The tiling picture is useful because:

### 1. Arc flips become single-square flips

If you want to ask "what happens if I reverse just one arc?", that's just flipping one square from white to black (or vice versa). In the tiling picture, this is a minimal local change.

This makes it easy to think about the **neighborhood** of a tournament — the set of all tournaments that differ from it by exactly one arrow.

### 2. Patterns correspond to properties

The staircase has three "sides" (like any right triangle), and each side of the staircase encodes a different aspect of the tournament:

- **The bottom row**: arcs involving the "sink" player (the one who loses most)
- **The left column**: arcs involving the "source" player (the one who wins most)
- **The diagonal**: arcs between players who are "close" in the P₀ ranking

Self-complementary tournaments (which tend to have high path counts) correspond to tilings with a specific **mirror symmetry**: the pattern looks the same when you reflect it along the diagonal.

### 3. Counting paths is a sum over tilings

The number of Hamiltonian paths H(T) can be expressed as a sum over features of the tiling. Each way of choosing a set of vertex-disjoint odd cycles in T corresponds to a specific subset of squares in the tiling pattern.

---

## The Hypercube

Here's a beautiful way to think about all possible tilings at once.

With m squares in the staircase, each square can be black or white — 2 choices. So there are 2^m possible tilings.

These 2^m tilings can be arranged as the vertices of a **hypercube**: an m-dimensional cube.

For 4 players: m = 3 squares, so 2³ = 8 tilings, arranged as the vertices of a 3-dimensional cube.
For 5 players: m = 6 squares, so 2⁶ = 64 tilings, arranged as the vertices of a 6-dimensional hypercube.

Two tilings are **neighbors** in the hypercube if they differ by exactly one square flip.

### The Metagraph from the Hypercube

Tournaments with the same abstract structure (same up to relabeling) get grouped into equivalence classes. These classes are the **vertices of the metagraph** — the map of all tournament types.

Flipping one square in the tiling corresponds to an edge in the hypercube. When you pass through that edge, you might stay in the same tournament class (the flip didn't change anything structurally) or move to a different class.

This is the "wiggly line" structure: the map of tournament types is obtained by collapsing the hypercube by relabeling symmetry, keeping track of which collapses go where.

---

## The Triangle is Everything

The staircase is a right isosceles triangle. And this triangle isn't just a visual convenience — it encodes deep structure:

- **The hypotenuse**: arcs between players who are "far apart" in the ranking. These arcs have the most freedom and carry the least information individually.
- **The vertical leg**: arcs from the top-ranked player. These encode the "score" hierarchy.
- **The horizontal leg**: arcs to the bottom-ranked player. These encode the complementary hierarchy.

Remarkably, four fundamental mathematical constants emerge from the geometry of this triangle:
- **√2**: from the ratio of hypotenuse to leg
- **π**: from the Wallis product formula counting tilings
- **e**: from the exponential growth rate of distinct tournament types
- **γ** (Euler-Mascheroni constant): from the asymptotic correction in counting formulas

These are the four most famous constants in mathematics, and they all appear from studying the shape of this staircase triangle.

---

## Key Words

- **Tiling**: a black-and-white pattern of squares in the staircase, encoding one tournament that contains the reference path P₀
- **Reference path P₀**: the fixed baseline Hamiltonian path used to set up the staircase
- **Square flip**: changing one square from white to black or back — corresponds to reversing one arc
- **Hypercube**: a high-dimensional cube whose vertices are all possible tilings (2^m of them)
- **Staircase Young diagram**: the triangular grid shape used to encode tournaments; a right isosceles triangle whose three sides encode three different aspects of tournament structure
