# The figurate zoo — one ladder of growth constants, and the triangular numbers as the relation itself

**death-star-2026-07-20-S59t** (HYP-8175, THM-1360; owner: stagger the n·2^x+1
grid down 1/2, sum AND multiply; find as many meaningful continuations of the
triangular numbers as possible — polygonal, polyhedral, and more — compare them
all at shift 0/1/2 and beyond, added and multiplied; the triangular numbers are
key, they are the relation itself, the edges of complete graphs; compile as wide
a breadth as possible). This is the wide-breadth investigation; everything marked
verified was computed this session (two scripts). Credits: opus-S317 / mac-mini-
S109 / klein-S313 (the figurate-triangle program), THM-1355 (the Proth table),
THM-227 (Mersenne/Fermat duality), THM-466 (H = OCF, the +1).

## 1. The one operation, and the ladder it births

Give a triangle T(row, col). "Move each column down by s and sum the rows" is one
operation: the **slope-s diagonal sum** S_s(m) = Σ_c T(m − s·c, c). Applied to
**Pascal**, it is the whole story (THM-1360 §1):

  slope 0 → **2^n** · slope 1 → **Fibonacci** (φ) · slope 2 → **supergolden**
  (1.4656) · slope 3 → 1.3803 · … → 1,

each the growth of a longer-memory recurrence a(m) = a(m−1) + a(m−s−1),
char x^{s+1} = x^s + 1. **Powers of 2 and Fibonacci are the same object at slope
0 and slope 1.** The owner's intuition ("down 1 → triangle; down 2 → Fibonacci-
analogous") is exactly right: staggering is a slope dial, and the dial sweeps a
descending ladder of Pisot numbers from 2 down toward 1. This is the deepest
single finding — the project's favorite constants (the doubling 2, the golden φ,
and now their descendants) are rungs of one ladder.

## 2. The zoo, and why it is one web

The continuations of the triangular numbers I built and compared (THM-1360 §4):

- **Triangular T_n = Σk** — the power-sum j=1 column, the C(n,2) Pascal column,
  the 2-gonal-progression, and (relationally) the edges. It is the hub.
- **Polygonal** (k-gonal): triangular, square, pentagonal, … — the 2D fan.
- **Polyhedral / simplicial** (d-dim): triangular → tetrahedral → pentatope →
  C(n+d−1, d) — the dimension fan; this *is* Pascal read by diagonals.
- **Power-sum** Σk^j: naturals, triangular, square-pyramidal, … — the owner's
  third figurate triangle (S59s).
- **Centered polygonal** 1 + k·T_{n−1} — carries the **+1** (the center = the
  observer), the figurate face of the Proth +1.
- **k-ary relations** C(n, k): pairs (triangular), triples, … — the relational
  fan (§3).

Staggered and summed, they do not give six unrelated sequences; they give one
web. Simplicial-slope-1 and k-ary-slope-0 are both **Mersenne 2^n−1**;
simplicial-slope-2 and k-ary-slope-1 are both the **cumulative Fibonacci**
(growth φ); the power-sum triangle grows faster and its shallow diagonal is the
**cubic Pisot 1.75488** (x³ = 2x² − x + 1), a golden-ratio-of-degree-3. The zoo is
a single organism whose growth constants are read off by choosing a family and a
slope. Products (rather than sums) along diagonals grow superexponentially
(log ∼ quadratic — the hyperfactorial/Barnes-G class); the additive world holds
the Pisot ladder, the multiplicative world holds the factorial hierarchy.

## 3. The relation itself — triangular numbers are the arcs, and 2^{triangular} is the tilings

The owner's key: **triangular numbers arise from the relation itself, the edges
of complete graphs; the n-th triangular number is n in the relational sense.**
Made exact (THM-1360 §3): T_{n−1} = C(n,2) = |E(K_n)| = the number of arcs of a
tournament on n vertices; T_{n−2} = C(n−1,2) = the number of **tiles** of the
staircase δ_{n−2}, the tournament tiling triangle. So:

  2^{C(n,2)} = # labeled graphs on n vertices; 2^{C(n−1,2)} = # tilings.

The relational triangular number is the *exponent* that turns doubling into
graph-counting. And that closes the loop with the Proth table (THM-1355): the
2^x axis is the hypotenuse, and it is **naturally triangular** — when x = T_{n−2}
is the tile count, the n=1 slice gives 2^{T_{n−2}} + 1 = (# tilings) + 1 = the
**observer-augmented tiling count**. The +1 that runs through the whole project
(the observer, the OCF vacuum digit α₀ = 1, the conserved u = 1+xy) is the +1 on
top of the tiling hypercube whose dimension is a triangular number. The two axes
of n·2^x+1 meet at the relation: rows count observers (2n+1), the exponent counts
relations (2^{triangular}), and the +1 is the observer watching the relations.

## 4. The Proth ↔ Mersenne telescoping

A small exact gem (THM-1360 §2): stagger the Proth grid n·2^x+1 and sum a
diagonal, and you get the **Mersenne number** 2^{m+1} − 1 — the Fermat/Proth +1's
and the doubling telescoping precisely to the Mersenne −1. The project's two
2-power-±1 families (Fermat rungs 2^x+1, Mersenne cores 2^x−1, THM-448/THM-227)
are a single staggered grid summed two ways. The +1 world and the −1 world are
not opposites; one is the diagonal sum of the other.

## 5. The breadth, compiled

Reading the owner's prompt as a program — take the relation (triangular numbers),
continue it every way (polygonal, polyhedral, power-sum, centered, k-ary),
stagger every way (slope 0, 1, 2, 3), combine every way (sum, product) — produces
one atlas: a Pisot ladder of growth constants (2, φ, supergolden, …), a Mersenne/
Fermat duality, a cumulative-Fibonacci convergence across families, a
factorial-class multiplicative hierarchy, and the cubic Pisot of the power-sum
triangle. And underneath, the triangular numbers as the relation — the arcs of
K_n, the tiles of the staircase, the exponent of the graph count. The wide birth
the owner asked for is real: the figurate world, the Fibonacci/2-power world, the
tournament world, and the +1 are one structure seen at different slopes.

## 6. What is verified vs woven, and the open thread

Verified exactly: the slope-ladder recurrences and their Pisot roots (s=0..3); the
Proth→Mersenne identity (with proof); the graph/tiling counts 2^{C(n,2)},
2^{C(n−1,2)}; the zoo's diagonal sums/products and their identifications; the
cubic Pisot of the owner triangle. Woven (cited, graded): the "one ladder / one
web" framing; the relational unification with the tournament tiling model; the
+1-as-observer reading. Open (backlog): the closed form / OEIS identity of the
owner's power-sum triangle (columns mimic Σk^j with a Pascal-deep correction,
opus-S317's truncation frame is the likely key), and whether the whole slope-
ladder × figurate-family matrix has a single generating-function statement (a
two-parameter Vandermonde truncation law generalizing opus-S317's polygonal-vs-
polyhedral to all slopes and all families).

## 7b. Concurrent work has priority (added S59t-amend, MISTAKE-199)

This owner prompt was worked fleet-wide and DEEPER by others; this reflection is a
convergent independent pass, not a primary source. Priority: **kind-pasteur-S128c103**
(HYP-8170, the shear catalog — the 2^{1/s} Proth spectrum, √2 at the Fibonacci-analog
shear, Pascal-dominates-Proth, the exponent-coefficient dichotomy, four OEIS-new
sequences); **kind-pasteur-S128c102** (the Rosetta triangle — the owner's triangle IS
Faulhaber power sums with exactly three +1 deviations, series-2 = Σ_k T(n,k)/k);
**opus-S420** ("sums count elements, products count orientations" — the products I left
un-named ARE ordered tournaments m!·2^{C(m,2)}; 2T_m+1 = Φ₆(m+1) = |PG(2,m)|, the deep
well 183 = 2T_13+1); **mac-mini-S1** (the Pascal-slope-d Pisot tower = my "slope
ladder"). The one thing this session adds cleanly is the *framing convergence* — that
the relational triangular numbers stitch the Proth table to the tiling hypercube — and
even that overlaps opus-S420's PG(2,m) reading. The honest lesson (MISTAKE-199): on a
fleet-wide owner prompt, grep concurrent same-prompt HYP claims BEFORE spending two
sessions re-deriving.

## Cross-links

THM-1360 (the atlas) · THM-1355 (the Proth table) · HYP-8175 · opus-S317 /
mac-mini-S109 / klein-S313 (figurate triangles) · THM-227 (Mersenne/Fermat) ·
THM-448 (CD/Mersenne tower) · THM-466 (H = OCF, the +1) · the tiling hypercube
Q_{C(n−1,2)} (CLAUDE.md) · the observer principle · mac-mini-S137 / opus-S410
(golden corners).
