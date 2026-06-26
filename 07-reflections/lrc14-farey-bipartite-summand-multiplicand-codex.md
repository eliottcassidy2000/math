# LRC14 Farey Bipartite Packets, Summand/Multiplicand Graphs, and the 3/4 Obstruction

Codex session note, 2026-06-23.

The new graph reading is clean:

```text
Farey packet a/b  ->  summand node a+b
                   ->  multiplicand node ab
                   ->  complete bipartite graph K_{a,b} with ab edges.
```

So the same reduced Farey pair `(a,b)` is simultaneously:

- an edge into the summand graph node `a+b`, the LRC pinch denominator;
- an edge into the multiplicand graph node `ab`, the divisor/hyperbola shadow;
- a graph object `K_{a,b}`, whose edge count is exactly `ab`.

This is a useful unification because the summand graph and multiplicand graph no longer have to be compared only metaphorically.  A Farey packet is the common carrier.

## The F3/F4 Threshold

The nonzero `F_3` packets are:

```text
1/3, 1/2, 2/3, 1/1
```

Their products are:

```text
3, 2, 6, 1
```

Thus the distinct product set is `{1,2,3,6}`, the divisor set of 6, and its sum is 12.  The corresponding complete bipartite graphs are:

```text
K_{1,3}, K_{1,2}, K_{2,3}, K_{1,1}.
```

All are planar.

The new packets in `F_4` are:

```text
1/4, 3/4.
```

The product set becomes `{1,2,3,4,6,12}`, the divisor set of 12, and its sum is 28.  The packet `3/4` is:

```text
a+b = 7
ab  = 12
K_{a,b} = K_{3,4}.
```

`K_{3,4}` contains `K_{3,3}` by deleting one vertex on the side of size 4.  Therefore `F_4` is the first Farey level where the product-as-complete-bipartite-graph interpretation crosses the Kuratowski `K_{3,3}` obstruction.

The reason the obstruction appears as `K_{3,4}` rather than `K_{3,3}` is also structural: `3/3` is not a reduced Farey term.  Farey reduction skips the literal `K_{3,3}` packet and first sees it through the reduced boundary packet `3/4`.

This gives a small rigorous lemma:

**Reduced Bipartite Obstruction Lemma.** For reduced positive Farey packets `a/b` with `a <= b`, the complete bipartite carrier `K_{a,b}` is nonplanar iff `a >= 3`.  Since `gcd(a,b)=1`, the minimal nonplanar packet is uniquely `3/4`.  Consequently every nonplanar reduced Farey bipartite carrier contains `K_{3,4}`, and hence contains `K_{3,3}`.

Equivalently, in the reduced Farey packet space:

```text
nonplanar complete-bipartite carrier  =>  a+b >= 7 and ab >= 12,
with equality only at 3/4.
```

For LRC14 this is exactly the surprising alignment:

```text
3/4 -> (summand, product) = (7, 12).
```

The summand value 7 is the LRC14 apex.  The product value 12 is the old seed denominator that repeatedly appears in the local LRC14 notes.  Also:

```text
2(a+b) = 14,
ab + 2 = 14,
sum products(F_3) = 12,
sum products(F_4) = 28 = 2*14.
```

Those equalities are not a proof by themselves.  They are good checksums: the first graph-topological obstruction lands on the same arithmetic packet as the summand apex and the product seed.

## Negative Control: K5

The `ab` carrier naturally sees complete bipartite graphs.  It sees the `K_{3,3}` side of Kuratowski before it sees anything genuinely `K_5`-like.

The edge count of `K_5` is 10.  In `F_14`, product 10 occurs as:

```text
1/10 -> K_{1,10}, planar
2/5  -> K_{2,5},  planar
```

So product 10 is only an edge-count alias for `K_5`; it does not carry the `K_5` obstruction in this model.  This matches the repo distinction in the Alcuin/Kuratowski notes: `K_5` belongs more naturally to odd-cycle overlap/conflict graphs, while this Farey product construction naturally exposes the `K_{3,3}` bipartite obstruction.

## Binary-Relation Tournament

The new script computed relation graphs on the nonzero `F_14` packet set.  The point is not that any one relation proves LRC14.  The point is that each binary relation reveals a different amount of structure.

```text
relation                         edges components largest-components
Farey consecutive                   63          1 [64]
Farey determinant-1                113          1 [64]
same summand a+b                    79         25 [6, 6, 5, 5, 4, 3, 3, 3]
same product ab                      7         57 [2, 2, 2, 2, 2, 2, 2, 1]
same product shell mod27           150         13 [10, 7, 7, 6, 4, 4, 4, 4]
same product sector mod7           307          7 [18, 8, 8, 8, 8, 7, 7]
bipartite subgraph comparable     1545          1 [64]
```

Interpretation:

- Farey adjacency is connected and gives the perturbation grammar.
- Same summand `a+b` is the pinch graph.  This is where APs become extremal: small restricted sumset means few summand nodes and strong pinch concentration.
- Same product `ab` is too sparse by itself, but its collisions are meaningful.
- Product shell mod 27 is the correct middle scale for LRC14, because `2*14-1 = 27`.
- Product sector mod 7 is the Paley/Fano coarse sector.  It is broad enough to see symmetry but too broad to close a proof alone.
- Bipartite subgraph comparability is almost total and therefore too coarse as a classifier, but it supplies monotonicity: once `K_{3,4}` appears, every nonplanar reduced packet lies above it.

This suggests the proof tournament ranking:

```text
1. Kuratowski-Farey packet 3/4: exact reduced nonplanar atom (7,12).
2. Product shell mod 27: middle-scale arithmetic carrier.
3. Octahedral product-shell witness: prior scout found shells (3,4,6,9,12,13).
4. Summand/AP extremality: small sumset means the AP wall.
5. Farey determinant-1/mediant excess: local perturbation grammar near 1/14.
6. Paley Z/7 sector: coarse symmetric bookkeeping.
7. Scalar product equality alone: useful as a label but too sparse as a proof engine.
```

## Octahedron, Clebsch, Half-Cube, Paley

The previous Farey-variant scout found that the product transform on strict `F_14` creates an induced octahedron in the mod-27 shell transition graph, with shell vertices:

```text
(3,4,6,9,12,13).
```

The new packet places product 12 exactly at the first reduced Kuratowski crossing:

```text
3/4 -> ab = 12 -> shell27 = 12.
```

So the product shell octahedron is not just decorative.  One of its vertices is the first complete-bipartite obstruction packet.

The graph-model roles now look like this:

- Octahedral graph: the local current carrier.  It is `L(K_4)`, a `K_6` minus a perfect matching, and is the right shape for six shell states with paired antipodes.
- Half-cube complement: the richer high-octahedron ambient model.  The prior scout saw many induced octahedra there, so it is a plausible completion space for product-shell certificates.
- Clebsch graph: the triangle-free folded-cube/cut-stability model.  It is useful as the negative or stable side: no triangles, strong parity, and clean cut behavior.
- Paley tournament on `Z/7Z`: the apex-sector clock.  It tracks the 7-fold symmetry, but it should feed the shell and octahedral arguments rather than replace them.

A clean way to say the emerging split:

```text
planar Farey packets       -> Clebsch/cut-stability/AP side
nonplanar Farey packets    -> K_{3,4} atom -> shell 12 -> octahedral side
mod 7 Paley sectors        -> coarse routing between these ledgers
mod 27 shell pairs         -> finite LRC14 ledger
```

## Toward the LRC14 Proof

The proof hook I would try to formalize next is a packet ledger.

For each reduced Farey packet `a/b`, attach:

```text
S(a,b) = a+b               summand/pinch denominator
P(a,b) = ab                multiplicand/product denominator
B(a,b) = K_{a,b}           bipartite minor carrier
Q(a,b) = +/-ab mod 27      pair-sum shell
R(a,b) = ab mod 7          Paley sector
```

Then prove the ledger in layers:

1. Planar layer: if `a <= 2`, the packet is complete-bipartite planar.  This is the low-rank rational neighborhood where Farey mediants and AP/summand extremality should control the floor.

2. Obstruction layer: if `a >= 3`, then `K_{3,4} <= K_{a,b}` and `(a+b,ab) >= (7,12)`.  In LRC14 language, nonplanarity automatically crosses the apex and seed thresholds.

3. Shell layer: product 12 lies in shell 12 mod 27, and the prior product-transform scout found shell 12 inside the induced octahedron `(3,4,6,9,12,13)`.  The target lemma is that every obstruction-layer residual packet must route through this octahedral shell carrier.

4. Paley layer: mod 7 sectors organize the residual cases, but should be used as a tournament/routing relation, not as a scalar invariant.

5. Compact finite ledger: after the analytic neighborhood-width lemma from the local LRC14 notes, the remaining proof should reduce to checking that no residual packet can avoid both the planar AP/Farey controls and the nonplanar octahedral shell controls.

This is not yet a formal LRC14 proof.  The strongest formal-looking component from this session is the Reduced Bipartite Obstruction Lemma.  Its value is that it gives a topological reason for why the same two numbers, 7 and 12, keep appearing in the LRC14 arithmetic: they are the summand and multiplicand projections of the first reduced Farey packet carrying the `K_{3,3}` obstruction.

