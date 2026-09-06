# Native Boolean flow: the one-tree obstruction and an orbit-uniform fractional map

Status: **PROVED ELEMENTARY REPRESENTATIONS / FINITE-EXACT CONTROLS**;
independent root and `nc2_seed` written audits PASS. An all-order polynomial lower bound for
the Boolean triangle Laplacian remains **OPEN**. No new theorem ID.

## 1. Inheritance and the exact target

The target is the combinatorial Laplacian of the simple graph whose vertices
are isomorphism classes of Eulerian graphs on exactly n vertices and whose
edges toggle a triangle. Its vertex measure is uniform on classes. The
[Boolean consumer note](overnight_hexagon_sep05_boolean_consumer.md) proves
the factorial comparison lower bound `12/[n(n-1)n!]`, explains both lost
multiplicities, and rules out the first weighted Fourier eigenfunction at
every n>=4. We do not repeat the weighted spectral census.

The closest constructive mechanism is
[THM-4069 — basis dependence and canonical envelope](../../01-canon/theorems/THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope.md):
an anchored-triangle basis supplies paths, but its edge image must not be
mistaken for the entire canonical orbit graph. The corrected near miss is
MISTAKE-495. The canonical hostile below is a bottleneck in a *chosen routing
tree*, not a bottleneck in the ambient graph. A recovered sidecar is the
[old Eulerian Burnside code](../../04-computation/even_graph_burnside_s260.py),
which computes fixed spaces by literal edge-orbit parity constraints. This
is the degree-even cycle space, not the different tournament/Royle meaning
of “even graph”; see the
[three-count audit](../../04-computation/even_graph_equinumerosity_three_counts_klein.py).

The live board is: uniform orbit measure; stabilizers; pairwise random paths;
congestion rather than diameter; fixed-point kernels; discarded multiplicity.
All maps below retain the actual endpoint measure. None identifies the
weighted and Boolean operators.

## 2. A single routing tree necessarily has enormous congestion

Let G be any connected simple graph with q>=2 vertices and maximum degree
Delta, and route every endpoint pair through one spanning tree T. For an
edge e whose deletion gives a side of size s, all s(q-s) unordered endpoint
pairs use e. Under the uniform-vertex identity

```text
sum_v (f(v)-mean(f))^2=(1/q) sum_{a<b}(f(a)-f(b))^2,
```

even the *unlengthened* load of e is s(q-s)/q. Choose a centroid of T.
Every branch has size at most q/2; one branch has size at least
(q-1)/deg_T(centroid)>=(q-1)/Delta. Its edge therefore has load at least

```text
rho_T >= (q-1)/(2 Delta).                            (1)
```

The usual path-length-weighted Cauchy--Schwarz congestion is no smaller.
For the Boolean triangle graph, `Delta<=binom(n,3)` and
`q>=2^binom(n-1,2)/n!`. Thus every one-tree all-pairs certificate has
superpolynomial congestion in n. This includes deterministic rootward
cycle-decomposition schemes when their paths form a single tree.

This is **not** a lower bound for optimal fractional congestion, and it is
**not** an upper bound or obstruction for the actual spectral gap. Random
pair-dependent paths or genuinely fractional flows are still available.

## 3. A precise pairwise fractional map with the right measure

Fix vertex zero and identify the cycle space with `F_2^m`, where
`m=binom(n-1,2)`, using the anchored triangles `(0,i,j)`. Let `[x]` be the
isomorphism class of x and `s(x)` its number of labelled representatives.
For each pair of distinct endpoint classes, independently choose uniform
labelled representatives x,y, and toggle the nonzero coordinates of
`D=x xor y` in a uniformly random order. Project each intermediate state
to its class. Every nonconstant step is a native Boolean triangle edge;
stationary steps contribute zero, and repeated quotient edges retain their
multiplicity. All stars are conjugate, and uniform representative choices
make the averaged class flow independent of the chosen anchor.

For a Boolean edge e, set

```text
rho_e = (1/q) sum_{x<y, [x]!=[y]} 1/[s(x)s(y)]
         * sum_{U subset D} sum_{i in D\U}
           1_{ { [x xor U], [x xor U xor {i}] } = e }
           / binom(|D|-1,|U|).                       (2)
```

Any fixed ordering of the labelled coordinate vectors may define `x<y`.
The indicator is zero for a stationary class step. Formula (2) follows
because a random order first visits U and then toggles i with probability
`1/[|D| binom(|D|-1,|U|)]`. Cauchy--Schwarz using the original path length
`|D|` cancels that first denominator. Summing the endpoint identity proves

```text
gap(L_B) >= 1/max_e rho_e.                           (3)
```

The representative weights in (2) are necessary: each endpoint *class*
has unit weight in the variance identity. Replacing them by uniform labelled
endpoints is precisely the earlier stationary-measure loss. Formula (2)
is an all-order representation, not a closed polynomial estimate.

The exact bounded probe compares subset-position aggregation with literal
enumeration of every ordering, for every labelled endpoint pair at n=3,4,5:

| n | classes q | native edges | max rho | certified gap lower bound |
|---|---:|---:|---:|---:|
| 3 | 2 | 1 | 1/2 | 2 |
| 4 | 3 | 2 | 2 | 1/2 |
| 5 | 7 | 8 | 50081/7875 | 7875/50081 |

This is stronger than the factorial comparison at these orders, but these
three rows are not evidence of an all-order bound by themselves. In fact
the maximum at n=4 and n=5 is on the C3--C4 edge, not the empty--C3 edge;
the tempting “only the isolated empty endpoint matters” shortcut already
fails in the cheapest cases.

## 4. The forgotten fixed-point sidecar has a closed rank

For a vertex permutation with cycle lengths `l_1,...,l_r`, the number of
orbits on edges is the inherited elementary formula

```text
b=sum_i floor(l_i/2)+sum_{i<j} gcd(l_i,l_j).
```

There is a useful closed form for the parity-constraint rank, eliminating
the old script's row reduction:

```text
rank A = r-1, if some l_i is odd;
rank A = r,   if every l_i is even.
dim Fix(cycle space,g)=b-r+1_{some odd l_i}.          (4)
```

Indeed invariant edge unions have constant degree parity on each vertex
cycle. Every even vertex cycle has an antipodal matching edge orbit,
providing its individual unit constraint column. Between two odd vertex
cycles, every cross-edge orbit gives the sum of their two unit columns,
spanning the even-sum subspace on the odd cycles. Odd/even cross-edge
orbits have zero odd-cycle coordinate and are already supported on the
even coordinates. Edges inside an odd cycle have even degree. No edge
orbit has odd total parity among the odd vertex cycles. These statements
give both spanning and the one possible missing rank, proving (4).

The checker independently compares (4) to the old literal edge-orbit
matrix for all 914 integer partitions of orders 1,...,16. The count agrees
with 2,3,7,16,54 classes at n=3,...,7 and with the inherited larger counts.
The matrix-rank mechanism and Eulerian Burnside count are inherited; the
closed rank formula is a transparent refinement, not a new equinumerosity
claim or a claim of literature novelty.

## 5. An exact uniform-orbit lift, not a mixing theorem

Consider all pairs `(g,F)` with `g in S_n` and F an Eulerian graph fixed by
g, and put uniform measure on these pairs. Any graph class O contributes

```text
sum_{F in O} |Aut(F)| = |O| * n!/|O| = n!
```

pairs. Hence projection to graph classes is exactly uniform; projection
to a labelled graph F has weight `1/[q s(F)]`, exactly the endpoint weight
needed in Section 3.

This also gives an exact finite sampling description: choose a permutation
cycle type with weight `class_size * 2^dimFix` using (4), then a uniform
permutation in that conjugacy class and a uniform vector in its fixed
cycle-space kernel. Output that labelled graph as a representative of a
uniform random class. This uses no rejection based on unknown graph
automorphism size. It is an identity for a fixed-pair mixture, not a claim
that a local Burnside Markov chain mixes rapidly, nor a polynomial-time
claim for enumerating all integer partitions.

The possible connection is now sharply typed: fixed-point kernels supply
the correct endpoint measure for a native fractional flow. What is still
missing is an all-order congestion bound for (2), or a better pairwise flow.
The rank formula alone does not supply it. A single rooted tree is ruled
out as a polynomial-certificate strategy by (1); the fractional map is not.

The [direct sampler continuation](eulerian-uniform-orbit-sampler-overnight-hexagon-sep05.md)
implements this exact mixture with `exp(O(sqrt(n)))` preprocessing and
expected polynomial per-sample bit complexity, without a graph-isomorphism
oracle. Its exhaustive distribution controls and n=30 witnesses do not
resolve the remaining local-flow congestion obligation.

## Reproduction

```bash
python3 -B 04-computation/overnight_hexagon_sep05_boolean_flow.py
python3 -B -O 04-computation/overnight_hexagon_sep05_boolean_flow.py
```

The [source](../../04-computation/overnight_hexagon_sep05_boolean_flow.py)
and [output](overnight_hexagon_sep05_boolean_flow.out) retain all 1,670
cross-orbit labelled endpoint pairs, every native edge and all random
orders at n<=5, plus the independent fixed-point rank bank. There are 928
gates; exceptions, not Python assertions, enforce them under optimization.
The script does not run another weighted spectrum census or allocate a
larger Boolean graph.

Normal and optimized runs agree, with semantic digest
`6ddd5dd6adaf696dc8d390d250d01d839135053dd16dca12dc3a904adddacc7e`.
Raw LF source SHA256 is
`b08a25da38a619cdde0c129bb50160e234882825e3e69ee2a6e89ac234241586`;
output SHA256 is
`706b9b89e393226d53ef9b58b622ae78199cd3a3a52c14449ede604fb370634b`.
