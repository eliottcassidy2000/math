# Square-sum graphs: two interpretations of the prime-group comparison

**Status:** PROVED elementary formulas and thresholds; FINITE-EXACT graph and
Hamiltonian censuses. No new moonshine correspondence is claimed.

## 1. Inputs, types, and the uncorrected overlap

The supplied groups are preserved verbatim:

- A = {2,3,5,7,11,13};
- B = {17,19,23,31,47};
- C = {29,41,59,71};
- D = {37,43,53,61,71};
- E = {73,79,83,89,97}.

There are 25 entries, 24 distinct primes, two occurrences of 71 and no 67.
Consequently these are not a partition of the first 25 primes. Replacing
**either** occurrence of 71 by 67 repairs the set-theoretic defect; the graph
data alone does not choose which occurrence was intended. In particular,
D:71->67 is a conditional proposed repair, not a silent correction.

Let p_i be the i-th prime, with p_1=2 and p_25=97. Define the **simple,
labelled square-sum graph** G_n on {1,...,n}, with distinct vertices x,y
adjacent exactly when x+y is a perfect square. There are two natural readings:

1. **Prime-order reading:** the label p_i selects G_{p_i}.
2. **Rank reading:** the label p_i selects G_i, the i-th graph in G_1,G_2,... .

Neither is identified with taking the first 25 *Hamiltonian-admissible*
orders or the first 25 terms of some other named square-sum sequence. Those
would require a third explicit indexing rule. A single square-sum permutation
is a Hamiltonian path in G_n, not the whole graph G_n.

The files `prime_order_first25.md` and `rank_order_first25.md` give both
25-row tables, including a MISSING row for the label 67 and C,D for label 71.
Under the prime-order reading, G_71 is duplicated and G_67 omitted; under
the rank reading, G_20 is duplicated and G_19 omitted.

## 2. Exact all-order arithmetic of the graphs

For every n>=1 and 1<=v<=n,

    deg_{G_n}(v) = floor(sqrt(n+v)) - floor(sqrt(v))
                  - 1_{2v is a square}.                         (1)

Indeed each root k with v<k^2<=n+v gives the unique prospective neighbor
k^2-v in [1,n]; remove the one forbidden loop if 2v=k^2.

On adjoining n to G_{n-1}, its earlier neighbors are k^2-n with
n<k^2<2n. Thus the **arrival-root set** and edge increment are

    R(n) = {floor(sqrt(n))+1,...,floor(sqrt(2n-1))},
    |E(G_n)|-|E(G_{n-1})| = |R(n)|.                              (2)

Here G_0 has no vertices or edges. A fixed root k appears in R(n) exactly
for floor(k^2/2)+1 <= n <= k^2-1. These square and half-square boundaries,
not primality, drive the natural arrival phases.

Edges of a fixed root k form a matching: each vertex has at most the one
partner k^2-v. Its size is

    e_k(n) = max(0, min(floor((k^2-1)/2), n-floor(k^2/2))).        (3)

To prove (3), orient each such edge as x<y. Its smaller endpoint ranges
from max(1,k^2-n) to floor((k^2-1)/2). Sum (3) over k>=2 for the total edge
count. Equations (1)--(3) are exact, not asymptotic. The natural maps
G_m -> G_n for m<=n are induced inclusions preserving labels and root-labelled
edges; they do not preserve Hamiltonicity or any supplied prime-group label.

## 3. Structural thresholds for every order

### 3.1 Connectedness

**Theorem.** G_n is connected iff n=1 or n>=14.

The graph G_14 is the tree with these 13 edges:

    1-3, 1-8, 2-7, 2-14, 3-6, 3-13, 4-5,
    4-12, 5-11, 6-10, 7-9, 11-14, 12-13.

The direct finite check of G_2,...,G_13 gives respectively
2,2,3,3,3,3,3,3,3,3,3,2 components. Every arrival n>=5 has an earlier
neighbor: for k>=3 the intervals

    [floor(k^2/2)+1, k^2-1]

cover every integer from 5 onwards. Consecutive intervals overlap or touch
because floor((k+1)^2/2)+1 <= k^2 for k>=3. Therefore every new vertex after
14 attaches to the connected old graph. This proves the infinite tail.

**Consequence.** Under the prime-order reading, A is **exactly** the set of
disconnected graphs among the first 25 prime orders. This is explained by the
composite threshold 14: A consists precisely of the primes below 14. Under
the rank reading this characterization fails: G_1 is connected and the
additional disconnected graphs G_7,...,G_13 lie outside A's six ranks.

### 3.2 Girth and odd girth

**Theorem.** G_n is a forest for n<=14; for 15<=n<=29 its girth and odd
girth are 5; for n>=30 both are 3.

The forest claim follows from the explicit G_14 tree and induced inclusion.
G_15 contains the 5-cycle [1,3,6,10,15]. The independent exact edge-deletion
shortest-path check of G_29 excludes cycles of lengths 3 and 4; inclusion
then excludes them at all smaller orders. G_30 contains the triangle
[6,19,30], with edge sums 25,49,36. Inclusion proves the infinite tail.
This theorem uses a fully checked finite base, not an unbounded search.

### 3.3 Two-vertex-connectedness

Use the convention that a 2-vertex-connected graph has at least three
vertices. **Theorem.** G_n is 2-vertex-connected iff n>=31.

The checker verifies a spanning open-ear decomposition of G_31, beginning
with triangle [6,19,30,6] and successively adjoining these paths:

    [6,3,1,8,17,19], [6,10,26,23,13,3],
    [23,2,14,11,5,4,12,13], [2,7,29,20,5],
    [7,9,16,20], [1,15,10], [7,18,31,5], [4,21,15],
    [3,22,14], [1,24,12], [11,25,24], [9,27,22], [8,28,21].

Every ear has distinct previously present endpoints, new interior vertices,
and square-sum edges. A cycle is 2-connected, and adjoining such an ear
preserves 2-connectedness: after deleting any one vertex, the old part
remains connected and every remaining portion of the ear attaches to it.
These ears cover all 31 vertices.

Every arrival n>=19 has at least two earlier neighbors. Check 19<=n<=24
in (2). For n>=25 set m=floor(sqrt(n))>=5; then

    2n-1 >= 2m^2-1 >= (m+2)^2,

because m^2-4m-5=(m-5)(m+1)>=0. Hence |R(n)|>=2. Adjoining a new vertex
with two distinct neighbors to a 2-connected graph preserves that property.
Induction proves the tail n>=31. All n<31 fail: the checker verifies that
they have fewer than three vertices, are disconnected, or have a degree-1
vertex. In particular G_30 still has the leaf 18, adjacent only to 7.

**Hostile boundary.** G_31 is 2-connected but has no Hamiltonian cycle;
G_32 has exactly one Hamiltonian cycle modulo rotation and reversal.
Both counts are independently exhausted. Do not promote 2-connectedness
to Hamiltonicity.

### 3.4 The inherited odd-root channels

Let H_n keep just edges whose square-sum root is odd. Session131 proved
that H_n has exactly the four residue components {0,1}, {2,7}, {3,6},
{4,5} modulo 8 for every n>=44. This packet independently verifies all
components through n=97, including the sharp transition at 44, but does
not claim that connectedness of these components makes them Hamiltonian.
The exact source and hash are retained in `source_pins.json`.

All prime-order graphs from p=47 onward consequently have the same four
odd-root components. Session131's even-root seam bound for an actual
Hamiltonian cycle at odd order is at least seven, also constant across
all the odd prime orders here where such cycles exist. These invariants
therefore cannot recover B,C,D,E.

## 4. The two exact 25-graph comparisons

### Prime-order reading p_i -> G_{p_i}

The 25 prime-order graphs fall into these verified regimes:

| Regime | Prime orders |
|---|---|
| Disconnected; no Hamiltonian path | 2,3,5,7,11,13 |
| Connected but no Hamiltonian path | 19 |
| Hamiltonian path but no Hamiltonian cycle | 17,23,29,31 |
| Hamiltonian cycle | 37,41,43,47,53,59,61,67,71,73,79,83,89,97 |

The obstruction at 19 is particularly transparent: deleting vertex 7
leaves three components, whereas deleting one vertex from a Hamiltonian
path leaves at most two. Their vertex sets are [9,16], [18], and the other
15 vertices. The positive cycle claims all have explicit checked cycles;
they do not depend on an unproved universal existence assertion.

B mixes all three connected regimes. C mixes path-only with cycles. D and E
are cycle-only over the stated range, but this does not distinguish them
from many vertices of B or C. The raw overlap 71 in C and D cannot be a
partition of any deterministic graph invariant's values.

### Rank reading p_i -> G_i

Among G_1,...,G_25, the nontrivial Hamiltonian-path orders are exactly

    i = 15,16,17,23,25,
    corresponding prime labels = 47,53,59,83,97,
    supplied groups = B,D,C,E,E,
    path counts modulo reversal = 1,1,1,3,10.

G_1 additionally has the trivial one-vertex path. No G_i with 1<=i<=25
has a Hamiltonian cycle. Thus the named groups do not classify path
existence in the first-25-order interpretation either.

A coarse Hamiltonian property cannot detect the duplicated label:
G_67 and G_71 both have Hamiltonian cycles, whereas G_19 and G_20 both
lack Hamiltonian paths. Exact vertex counts or edge lists do detect it.

## 5. What does and does not constitute a connection

The cleanest cross-group hostile probe is the arrival-root phase:

    R(41)=R(43)=R(47)={7,8,9},

although 41 belongs to C, 43 to D, and 47 to B. Their girth is 3, their
minimum degree is 2, and each has a Hamiltonian cycle. Likewise
R(53)=R(59)={8,9,10} crosses D and C. Therefore neither this phase nor the
listed coarse structural flags yields the supplied grouping. This refutes
those proposed classifiers, **not** every imaginable arithmetic relation.

Graphs at distinct orders cannot be isomorphic, since their vertex counts
differ. A comparison by group therefore needs a specified quotient or
feature map rather than literal graph isomorphism.

Connection contract:

- source: the prime labels with their supplied group memberships;
- target: either labelled initial-interval graphs G_p or rank graphs G_i;
- map: explicitly p->G_p or p_i->G_i;
- preserved: in the first map the label is the vertex count; in the second,
  prime order in the enumeration is the vertex count;
- lost: an arbitrary coarse graph signature generally discards the index,
  root labels, exact vertices, and any external arithmetic provenance;
- required sidecar: the chosen index convention, full membership list,
  and exact edge/root labels if a stronger transfer is proposed;
- decisive outcome: A matches the connectedness threshold only under the
  prime-order map; the remaining named groups fail the tested classifiers.

No abstract moonshine, genus, spectral, or sporadic-group relation has been
proved merely by placing these prime labels next to these graphs.
