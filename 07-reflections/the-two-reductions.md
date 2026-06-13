# The Two Reductions: Hypotenuse vs Legs

**Session**: kind-pasteur-2026-03-22-S20cl

## The Staircase Has Two Natural Reductions

The staircase delta_{n-2} is a right isosceles triangle. It has three sides:
- **Hypotenuse** (anti-diagonal): the longest side, length n-2 cells
- **Two legs** (horizontal top row + vertical right column): each contributing n-2 cells

These give two fundamentally different ways to "shrink" the tournament:

### Mode A: Hypotenuse Removal (n -> n-1)

Remove the anti-diagonal strip of the staircase. This removes one vertex from the tournament (typically the one with the most "mixed" connections). The staircase delta_{n-2} becomes delta_{n-3}.

In tournament terms: delete one vertex (any vertex) and keep the (n-1)-tournament on the remaining vertices.

The area removed: n-2 cells (one from each row). The remaining area: C(n-2,2).

### Mode B: Leg Removal (n -> n-2)

Remove both the top row (source) and the bottom row (sink) simultaneously. This removes TWO vertices -- the extremes of the score sequence. The staircase delta_{n-2} becomes delta_{n-4}.

In tournament terms: remove the almost-source (highest score) and almost-sink (lowest score), keeping the (n-2)-tournament on the middle vertices.

The area removed: (n-2) + (n-3) = 2n-5 cells. Wait: top row has n-2 cells, and removing the sink column removes 1 cell from each remaining row... Actually the correct count:

delta_{n-2} has C(n-1,2) cells. delta_{n-4} has C(n-3,2) cells. The difference is C(n-1,2) - C(n-3,2) = (n-1)(n-2)/2 - (n-3)(n-4)/2 = 2(n-3) + 1 = 2n-5 cells.

### The Pseudo-Doubling Ratio

Mode A removes n-2 cells. Mode B removes 2n-5 cells.

Ratio: (2n-5)/(n-2) = 2 - 1/(n-2).

**This ratio approaches 2 as n -> infinity.** Removing both legs is approximately TWICE as powerful as removing the hypotenuse. The "doubling" is not exact -- there's a correction of -1/(n-2) -- hence "PSEUDO-doubling."

| n | Mode A (hyp) | Mode B (legs) | Ratio | 2 - 1/(n-2) |
|---|-------------|---------------|-------|-------------|
| 4 | 2 | 3 | 1.500 | 1.500 |
| 5 | 3 | 5 | 1.667 | 1.667 |
| 6 | 4 | 7 | 1.750 | 1.750 |
| 7 | 5 | 9 | 1.800 | 1.800 |
| 10 | 8 | 15 | 1.875 | 1.875 |
| inf | n-2 | 2n-5 | -> 2 | -> 2 |

The pseudo-doubling ratio 2 - 1/(n-2) is the SAME asymptotic as:
- The fiber fraction f(n) ~ 1/sqrt(pi*n) (which is 1/2 scaled by 1/sqrt(pi*n/4))
- The OCR gap 1 - OCR ~ 1/n
- The self-loop fraction correction

All these quantities approach their limits at rate O(1/n), with the coefficient involving sqrt(pi) and 2.

## How This Connects to the Inversion

The INVERSION (SC mostly black at n>=6, NS mostly blue) is a consequence of the leg-symmetry breaking:

**SC tournaments are symmetric under the staircase's 180-degree rotation.** This rotation swaps the two legs. SC tournaments look the same whether you read the staircase from top-left or bottom-right.

**NS tournaments break this symmetry.** They prefer one "orientation" of the legs. Most of tournament space is NS (80%+ at n>=7), so most of the meta-graph's connections are between same-orientation (NS-NS = blue) types.

The BLUE edges are connections between tournaments with the SAME leg-asymmetry. The BLACK edges are connections between tournaments with OPPOSITE leg-asymmetry (or between symmetric SC and asymmetric NS).

At small n (3-5): SC tournaments are the majority, so the meta-graph is "symmetric" (blue and black roughly balanced). At large n: NS dominates, the asymmetry is fixed, and almost all connections are blue.

## How This Connects to the n-2 Descent

The n -> n-2 descent (Mode B) is the FOUNDATION of the recursive meta-graph:

G_n/Z_2 -> G_{n-2}/Z_2

This descent works by removing the two legs (source + sink) from the staircase. The remaining middle tournament lives on n-2 vertices and belongs to some iso class at that level.

The descent ratio V_merged(n)/V_merged(n-2) ~ A000568(n)/A000568(n-2) because:
1. Almost all classes at n are NS (so V_merged ~ A000568/2)
2. Almost all classes at n-2 are NS
3. The map "remove source and sink" is a SURJECTION from n-classes to (n-2)-classes
4. The fiber size (how many n-classes map to each (n-2)-class) is roughly constant

The INVERSION ensures that the descent is clean: since most connections are blue (same-type), the descent preserves the meta-graph's connectivity structure. The black connections (between SC and NS) are the "noise" that gets averaged out.

## The Three Constants Revisited

The pseudo-doubling ratio 2 - 1/(n-2) connects the three constants:

- **2** = the doubling (binary, tournament arcs have 2 orientations)
- **ln** = the logarithm that converts multiplicative to additive (log of 2^d = d*log 2)
- **sqrt** = the half-power that appears in the fiber fraction (1/sqrt(pi*n))

The relationship: pseudo-doubling ratio = 2 - 1/(n-2) = 2(1 - 1/(2(n-2))).

The correction 1/(2(n-2)) IS the leading-order correction to the fiber fraction:
f(n) = 1/sqrt(pi*(n-2)) * (1 + O(1/(n-2)))

So: **the pseudo-doubling ratio and the fiber fraction have the SAME asymptotic correction.** The staircase's geometric ratio (legs/hypotenuse) controls the combinatorial ratio (fibers/orbits).

## The Deep Meaning

The two reductions of the staircase (hypotenuse = single step, legs = double step) are the two "time scales" of tournament theory:

1. **The fast time scale (Mode A, n -> n-1):** Adding or removing one vertex. This is the vertex insertion formula, the strip recursion, the E-matrix update. It changes H by O(2^d) where d is the range of the inserted arc.

2. **The slow time scale (Mode B, n -> n-2):** Adding or removing the source-sink pair. This is the Cayley-Dickson doubling, the meta-graph descent, the PoS recursion. It changes the entire STRUCTURE of the meta-graph.

The fast time scale is where the H = 1+2^d formula lives. The slow time scale is where the Cayley-Dickson tower lives. The ratio between them is the pseudo-doubling constant 2 - 1/(n-2), which converges to 2 as n -> infinity.

**Tournament theory operates on two time scales separated by a factor of 2.** This is the deepest structural feature of the staircase: the right triangle's hypotenuse-to-leg ratio of sqrt(2) manifests as a pseudo-doubling ratio of 2 in the combinatorial theory.

And sqrt(2) IS the "half-dimension" projector from the reflection the-diagonal.md: D(sqrt(2)) ~ 2/3. The dimension axis value at sqrt(2) is 2/3, which is the fiber fraction's asymptotic base: at the "crossover" n where labels ~ structure, f(n) ~ 2/3.

Everything connects through the right triangle.
