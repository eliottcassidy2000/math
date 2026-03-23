# The Principal Line of Symmetry

**Session:** kind-pasteur-2026-03-23-S20cp
**The idea:** Start from the transitive tournament (H=1). Follow its blue connections through the SC backbone. This traces a SPINE through the merged meta-graph. All other classes hang off this spine. Organize them symmetrically.

---

## The Structure

### The Transitive Class
At every n, the transitive tournament T_trans has:
- Score sequence (0, 1, 2, ..., n-1) — maximally spread
- H = 1 — unique minimum (Redei's theorem: exactly 1 Hamiltonian path)
- c3 = 0 — no directed 3-cycles
- |Aut| = 1 — no automorphisms (all scores distinct)
- Self-complementary: YES (complement = reverse ordering, isomorphic)
- Score variance = (n^2-1)/12 — maximum among SC classes

The transitive class is the MOST ORDERED tournament. Its complement is itself. It sits at one extreme of every ranking invariant.

### The SC Blue Spine
From the transitive, follow blue edges (SC-preserving) through self-complementary classes. This traces a TREE:

| n | SC classes | Spine size | Max depth | Blue edges in spine |
|---|------------|------------|-----------|---------------------|
| 3 | 2 | 2 | 1 | 1 |
| 4 | 2 | 2 | 1 | 1 |
| 5 | 8 | 8 | 3 | 7 (Hamiltonian!) |
| 6 | 12 | 12 | 7 | 11 |
| 7 | 88 | 88 | ? | ? |

At n=5, ALL 8 SC classes lie on a single Hamiltonian path through the blue subgraph! The blue SC subgraph is connected and contains a Hamiltonian path. At n=6, all 12 SC classes are reachable via blue from the transitive.

### The Principal Axis (Longest Path)
The longest blue path from the transitive through SC classes:

**n=3:** [transitive] --blue--> [3-cycle]
- H: 1 -> 3
- Score: (0,1,2) -> (1,1,1)
- Score variance: 0.67 -> 0.00

**n=4:** [transitive] --blue--> [4-cycle/regular]
- H: 1 -> 5
- Score: (0,1,2,3) -> (1,1,2,2)
- Score variance: 1.25 -> 0.25

**n=5:** The axis visits ALL 8 SC classes
- H along axis: 1, 9, 9, 11, 15, 13, 15, 3
- The path is NOT monotone in H — it oscillates!
- Score variance: 2.00, 0.80, 0.80, 0.40, 0.00, 0.40, 0.40, 1.60
- The regular tournament (H=15, variance=0.00) sits at position 4 of 8 (middle)

**n=6:** The axis visits 10 of 12 SC classes
- H along axis: 1, 5, 37, 33, 41, 45, 29, 17, 29, 45
- Score variance: 2.92, 2.25, 0.92, 0.92, 0.25, 0.25, 0.92, 1.58, 0.92, 0.25
- Two H-maxima (H=45, the regular tournament) at positions 5 and 9

### The Pattern: Variance Decrease Toward the Center

Along the principal axis, score variance GENERALLY DECREASES from the transitive (maximum variance) toward the regular tournament (zero variance). But the path is not a straight line — it meanders.

The regular tournament (maximum H, zero score variance) acts as a CENTRAL ATTRACTOR. The transitive tournament (minimum H, maximum variance) is the PERIPHERAL ORIGIN. The axis traces a path from periphery to center through the SC backbone.

---

## Bilateral Structure

### The Two Branches
At n >= 5, the transitive class has EXACTLY 2 blue neighbors, creating two branches:

**n=5:**
- Branch A: -> class 3 (H=3, score (0,2,2,2,4)) -> 5(H=15) -> 6(H=13)
- Branch B: -> class 4 (H=9, score (1,1,2,3,3)) -> 9(H=9), 7(H=11) -> 8(H=15)

Branch A leads to the "near-regular" H=15 class with score (1,2,2,2,3).
Branch B leads to the fully regular H=15 class with score (2,2,2,2,2).
The two branches MEET at the H-maximum level from opposite directions!

**n=6:**
- Branch A: -> class 6 (H=5, score (0,2,2,3,3,5)) — short, dead end
- Branch B: -> class 11 (H=17) -> 14(H=37) -> 19(H=33) -> ... (long chain)

One branch is stunted (single vertex), the other is the main trunk. This asymmetry GROWS with n.

### NS Class Attachment

NS (non-self-complementary) classes attach to the SC spine ONLY via black edges. They never appear on the spine itself.

At the transitive end of the axis, NS branches are ALL ABOVE (higher H than the axis vertex). At the maximum-H end, NS branches are ALL BELOW (lower H). There is a crossover point in the middle where the transition happens.

This creates a natural bilateral decomposition:
- **Left side** (from transitive toward the crossover): branches go UP
- **Right side** (from regular toward the crossover): branches go DOWN

### Score Staircase Pattern

EVERY SC class on the axis has a palindromic score sequence (self-complementary scores satisfy s_i + s_{n-1-i} = n-1). When visualized as a staircase:

```
n=3:  transitive: |.|#|   3-cycle: |||||
n=4:  transitive: |..##|  regular: |..##|
n=5:  ALL axis:   |..|##| (consistently 2 below, 1 at median, 2 above)
n=6:  ALL axis:   |...###| (consistently 3 below, 0 at median, 3 above)
```

At n=5 and n=6, the staircase visualization is IDENTICAL for every class on the axis! The binary representation (below/above median) doesn't change along the axis — only the specific scores within each half change.

---

## Cross-n Patterns

### 1. First Step Delta H
The H-jump from transitive to its first blue neighbor:
- n=3: Delta = 2 = 2^1
- n=4: Delta = 4 = 2^2
- n=5: Delta = 8 = 2^3 (to the H=9 branch) or 2 (to the H=3 branch)
- n=6: Delta = 4 (to H=5 branch) or 16 (to H=17 branch)

The LARGER branch jump follows 2^k: 2, 4, 8, 16 for n = 3, 4, 5, 6. This is exactly 2^(n-2)! This connects to the hypotenuse formula H = 1 + 2^d.

### 2. Axis Length / SC Count
- n=3: 2/2 = 1.00
- n=4: 2/2 = 1.00
- n=5: 8/8 = 1.00 (all SC on axis!)
- n=6: 10/12 = 0.83

The axis covers nearly all SC classes, but as n grows, some SC classes fall off the longest path.

### 3. Number of NS Classes Attached to Each Axis Vertex
At the transitive (H=1): many NS branches (2 at n=5, 6 at n=6).
At the regular (H=max): few or no NS branches.

The transitive is the most-connected-to-NS vertex on the axis. This makes sense: the transitive tournament is the most "borderline" SC class — one flip can easily create a non-SC tournament.

### 4. Mirror Fraction
- n=3,4,5: 0% (no off-axis vertex has a mirror with same H, same axis attachment, same distance)
- n=6: 16.7% (4 of 24 off-axis classes have mirrors)

Mirror symmetry is EMERGING as n grows. At small n, the structure is too sparse for mirrors. At larger n, the denser NS attachment creates mirror pairs.

---

## The Geometric Picture

Imagine the SC blue spine as a **vertical backbone** with the transitive at the bottom (H=1) and the regular at the top (H=max). The NS classes are **ribs** extending horizontally from the backbone.

At the bottom of the backbone, ribs point UP (NS classes with H > axis H).
At the top, ribs point DOWN (NS classes with H < axis H).
In the middle, ribs extend in both directions.

This creates a shape like a **double helix** or **fish skeleton**:

```
                  H=max (regular)
                     |
              NS --- SC --- NS
                     |
              NS --- SC --- NS
                     |
                    SC
                     |
              NS --- SC --- NS
                     |
                    SC
                     |
              NS --- SC
                     |
                  H=1 (transitive)
```

The "ribs" are longer near the middle of the backbone (more NS classes attach where the score structure is richest) and shorter near the extremes.

---

## Connection to the Triangle Foundation

The principal line IS the vertical leg of the triangle delta_{n-2}:
- **Bottom** = transitive = source column = score (0,1,...,n-1) = maximum hierarchy
- **Top** = regular = score (k,k,...,k) = no hierarchy
- The axis traces the path from maximum hierarchy to minimum hierarchy
- NS classes branching off represent INTERMEDIATE levels of hierarchy

The score variance along the axis decreases from bottom to top, exactly as the vertical leg of the staircase shortens as you move along the hypotenuse.

The two branches from the transitive correspond to the two "modes" of breaking the transitive order:
- **Mode A** (small Delta H): flip a single arc near the extremes of the score sequence
- **Mode B** (large Delta H = 2^(n-2)): flip an arc in the middle of the score sequence

Mode B creates a more dramatic restructuring (hence larger H jump). This is the hypotenuse-removal mode from the Triangle Foundation.
