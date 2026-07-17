# LRC(14) weighted-ratio strict-clique certificates (S53)

## Scope

This package makes the finite strict-clique part of the LRC(14) weighted-pair
route independently replayable.  It does **not** claim a complete LRC(14)
proof.  The remaining analytic/formal crux is the continuous-to-clean-grid
discrepancy lemma described below.

The executable is

```text
04-computation/lrc14_weighted_ratio_clique_certificates_codex_S53.py
```

and the seven compact proof DAGs are in

```text
05-knowledge/results/lrc14_weighted_ratio_clique_certificates_codex_S53/
```

Run either

```bash
python3 04-computation/lrc14_weighted_ratio_clique_certificates_codex_S53.py \
  check 05-knowledge/results/lrc14_weighted_ratio_clique_certificates_codex_S53
```

or regenerate a fresh directory with `generate` in place of `check`.

## Exact pair observable

For nonzero integer speeds `x,y`, put `g=gcd(|x|,|y|)`,
`a=|x|/g`, and `b=|y|/g`.  With the periodic Bernoulli polynomial

```text
B2(t) = t^2 - t + 1/6,
```

the continuous pair covariance is

```text
P(a,b) = [B2({(a-b)/14}) - B2({(a+b)/14})] / (a b).
```

Thus the observable depends on the primitive ratio, not the common scale.  Its
negative numerator has absolute value at most `12/49`, so

```text
max(0,-P(a,b)) <= 12/(49 a b).
```

Consequently every strict superlevel set is finite and can be reconstructed
exactly with integer arithmetic and `Fraction`; no floating-point decision is
used.

## Anchored ratio graph

For a positive threshold `t`, let `H_t` have one vertex for each oriented
primitive ratio `a/b` satisfying `-P(a,b) > t`.  Two vertices `r,s` are
adjacent exactly when the reduced quotient `r/s` is another allowed ratio.

If the strict weighted runner graph at level `t` has a `k`-clique, anchoring
one runner converts the other `k-1` runners into a `(k-1)`-clique of `H_t`.
Therefore a certificate that `H_t` is `K_(k-1)`-free proves that the runner
graph is `K_k`-free.

The strict thresholds and replayed graph fingerprints are:

| runner clique excluded | threshold | ratio vertices | ratio edges | proof nodes | branches |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | `2/441` | 14 | 0 | 1 | 0 |
| 4 | `5/2646` | 60 | 66 | 10 | 9 |
| 5 | `1/1764` | 272 | 861 | 52 | 51 |
| 6 | `1/3136` | 534 | 3048 | 137 | 137 |
| 7 | `5/37632` | 1454 | 15654 | 1173 | 1194 |
| 8 | `1/9996` | 1990 | 26217 | 1508 | 1513 |
| 9 | `1/14112` | 2960 | 51627 | 2326 | 2333 |

The top value is `tau_2=6/637`; no negative pair has larger weight.

## Certificate rule

Each JSON file stores only a compact integer-ID proof DAG.  The checker does
not trust or deserialize a precomputed graph.  It reconstructs `H_t` from the
Bernoulli formula, checks its vertex/edge fingerprint, and replays every proof
node.

At a state `(P,r)`, a deterministic greedy coloring partitions `P` into
independent color classes.  If its color bound is below `r`, the state contains
no `r`-clique.  Otherwise the checker reads the colored order backwards.  For
each required vertex `x`, it recursively verifies that `P intersect N(x)` has
no `(r-1)`-clique and then deletes `x` from `P`.  Replay also rejects:

- a node reused under a different `(P,r)` context;
- missing, extra, or reordered branches;
- cycles or unreachable nodes;
- a malformed size-one leaf.

This is a proof of nonexistence, not merely a report from a maximum-clique
search.

## Parallel-class circle layers

Primitive ratio colors are also the arithmetic labels of parallel pair-wall
families on the phase circle.  At the two highest nontrivial layers the exact
formula leaves very little structure:

- strictly above `5/588`, only ratio color `13` (and its reverse) remains;
- strictly above `4/539`, only ratio colors `12` and `13` (and reverses) remain.

For a fixed ratio `c>1`, orient an edge from `u` to `cu`.  Every vertex has at
most one predecessor and one successor, and magnitudes strictly increase along
an oriented edge.  Hence a fixed ratio color is a disjoint union of paths and
has at most `12` edges on thirteen runners.  The one-color and two-color caps
are therefore `12` and `24`.

This is the safe parallel-class-circle contribution.  Pairwise overlap of
circle wall families must not be promoted to a common phase by a Helly
assumption: three families can realize a boundary-of-a-triangle nerve, with
the three pairwise overlaps occurring at different phases and no triple
intersection.  The present bound uses only the fixed-color path property and
does not require such a common chamber.

## Layer-cake bound

Let

```text
W(v) = sum_{i<j} max(0,-P(v_i,v_j))
```

for the thirteen nonstationary speeds.  Write `T_k(13)` for the Turan edge
number of a `K_(k+1)`-free graph on thirteen vertices.  For `k=3,...,8`, these
numbers are `56,63,67,70,72,73`.  The strict-clique certificates, the two path
caps, and the trivial tail give

```text
W(v) <= 78 tau_9
      + sum_{k=3}^8 T_k(13) (tau_k - tau_(k+1))
      + 42 (4/539 - tau_3)
      + 24 (5/588 - 4/539)
      + 12 (tau_2 - 5/588)

     = 176738453/411675264
     = 0.429315211418678...
     < 13/30.
```

The exact margin is

```text
13/30 - W(v) >= 8270807/2058376320
                  = 0.004018121914656...
```

A stronger two-color grid-planarity argument would replace `24` by `22`, but
the package deliberately uses the simpler path-only cap.

## Zarankiewicz audit

Using the convention that `z(m,13;2,2)` is the maximum number of edges in a
`K_(2,2)`-free bipartite graph with part sizes `m` and `13`, the exact small
values audited in parallel are

```text
m                 1  2  3  4  5  6  7  8  9 10 11 12 13
z(m,13;2,2)      13 14 16 19 23 27 30 33 37 40 44 48 52
```

These values are useful diagnostics for incidence formulations, but a raw
uncolored Zarankiewicz bound is not sound here without an additional
`K_(2,2)` exclusion.  Multiplicative rectangles

```text
x, r x, s x, r s x
```

naturally create four-cycles.  Forgetting the ratio/parallel-class colors
therefore destroys precisely the information needed to distinguish a genuine
obstruction from a legal arithmetic rectangle.  A future refinement could use
a color-sensitive Zarankiewicz problem; the conservative S53 bound does not
assume one.

## Tournament and quotient audit

The pair observable is symmetric, so thresholding produces an undirected
graph.  Orienting edges by numerical order would be an arbitrary gauge: its
score sequence, directed cycles, SCCs, and tie Hamiltonian path would not
preserve the weighted-pair predicate.  The meaningful fingerprints here are
instead strict clique numbers, edge counts, threshold flips, and fixed-color
path counts.

The explored vertex sets were runners, primitive anchored ratios,
parallel-class circle walls, and wall-crossing events.  Primitive ratios were
chosen because the quotient preserves every pair-threshold and multiplicative
consistency condition needed by the anchored-clique implication.  It destroys
absolute scale, higher-support phase chronology, and the location of a common
circle chamber.  Those losses are harmless for the pair-mass bound but would
not be harmless for a direct many-runner intersection proof.

The challenged assumption is that runners or uncolored incidences must be the
vertices.  The ratio quotient is faithful for the present pair predicate;
uncolored Zarankiewicz and tournament quotients are not.

## Remaining bridge

The finite-grid formalization still needs the uniform discrepancy estimate

```text
abs(normalizedMass2(v,q) - C(v))
  <= (24 sum_i |v_i| + 78)/(q-1).
```

The intended proof counts at most `|v_i|+|v_j|` open circle components for a
pairwise bad intersection, charges at most `2/q` discrepancy per component,
and then charges at most `1/(q-1)` per pair when replacing the grid containing
the always-bad zero phase by the normalized nonzero grid.

At `q = cleanModulus v 534` and `sum_i |v_i| >= 13`, the target simplifies to

```text
5/1246 < 8270807/2058376320,
```

with residual room

```text
967423/183195492480.
```

Once this inequality is kernel-checked in the clean-grid development, the S53
pair-layer margin is large enough to absorb the discretization error.
