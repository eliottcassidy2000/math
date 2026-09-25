# Halving, the doubles seam, and what a minor certificate must retain

**Date:** 2026-09-25. **Status:** PROVED elementary claims and an explicit
planarity threshold with exact certificates; CITED graph/tournament minor
theorems; OPEN global Collatz decoder. No convergence or equivalence of
prize problems is claimed.

## Inheritance, portfolio, and board

The closest proved mechanism is
[THM-371, first-doubling-unit-pair-seam](../../01-canon/theorems/THM-371-first-doubling-unit-pair-seam.md):
for odd m, phi(2m)=phi(m), whereas later dyadic lifts double phi; the
matching count creates its extra pair exactly at the first doubling.
The [four-core note](collatz_tournament_core_20260925.md) supplies the
substitution Q[H,H,H,1] and the even-order hostile with no pair module.
The [arithmetic operation note](arithmetic_seams_20260921_operations.md)
already proves that its full addition and multiplication shadows contain
arbitrarily large complete graphs. We therefore use a sparse operation
graph, rather than rediscover those saturated shadows.

Canonical hostile: orientation-forgetting sends every n-vertex tournament
to K_n. Corrected near miss: a two-dimensional incidence array is not a
planar embedding (2026-09-21 arithmetic-seams entry in MISTAKES).
Least-used relevant sidecar: the edge word and its arithmetic guards on
each contracted path; a branch set alone does not restore them.

Anchor: sound contraction certificates for a decoder. Niche: literal +2
versus row-scaled +2 on dyadic coordinates. Wildcard: the first nonplanar
prefix of the sparse +2/doubling graph.

| Concept | Preserved predicate / operation | Lost coordinate / decisive test |
|---|---|---|
| Dyadic rows | n=2^r q, q odd; doubling increments r | Literal +2 may change r; compute v2(n+2) |
| Pair-module halving | Uniform external orientations | An arbitrary even tournament need not admit pairs |
| Strong tournament minor | Directed routes through strongly connected branch sets | Pair modules are not strongly connected |
| Sparse operation graph | Explicit arithmetic adjacencies | Undirected paths forget direction and guards |
| Contracted certificate path | Endpoint route with stored word | Bare edge loses division and parity checks |
| Planarity obstruction | K3,3 subdivision or planar rotation | Numerical proximity to square-sum thresholds is not a map |

## 1. The exact arithmetic seam

Write n=2^r q with q odd. **PROVED:** literal addition by two has three
different behaviors:

```
r=0:  v2(n+2)=0;
r=1:  n+2=2(q+1),                 so v2(n+2)>=2;
r>=2: n+2=2(2^(r-1)q+1),         so v2(n+2)=1.          (S1)
```

Thus the odd chain is closed under +2; above it, +2 alternates between
the doubles seam and the sea of multiples of four. Doubling, by contrast,
is the uniform vertical operation r->r+1. This is a precise mechanism
behind the user's three-region picture, independent of any Collatz claim.

There is a second operation that must not be confused with literal +2.
In coordinates

```
n = 2^r(2k+1),                    (r,k) in N_0^2,
```

doubling sends (r,k) to (r+1,k), while **row-scaled** addition
`n -> n+2^(r+1)` sends (r,k) to (r,k+1). The undirected graph of these
two operations is exactly the quadrant square grid, hence planar.
Halving an even vertex is the partial inverse of the first operation.
This graph preserves its row coordinate; the literal +2 graph does not.

The elementary dilation audit makes the distinction visible another way:
on even labels, dividing an edge {n,n+2} by two gives {n/2,n/2+1},
not an edge of additive step two. Preserving the additive mode under a
halving quotient therefore requires storing its scale.

## 2. An exact K3,3 threshold for literal +2 and doubling

Define the finite **undirected simple** graph G_N on 1,...,N by edges
{x,x+2} and {x,2x}, whenever both endpoints are in range. The duplicated
edge {2,4} is counted only once. This is an explicitly selected model of
the two arithmetic operations, not a claim that the user had uniquely
specified this graph or that it is the square-sum graph.

**PROVED, with exact finite certificates:**

```
G_N is planar if and only if N<=15.                       (S2)
```

For G_15, the script stores a rotation at each vertex. Every dart occurs
once; its face permutation has six orbits and

```
V=15, E=19, F=6, V-E+F=2.
```

The graph is connected, so the rotation constructs a cellular embedding
in the orientable surface of genus zero. Every smaller G_N is its
subgraph and is planar.

For G_16 use the shores {5,8,12} and {6,10,14}. The following nine paths
give a subdivision of K_(3,3):

| Start | To 6 | To 10 | To 14 |
|---|---|---|---|
| 5 | 5,3,6 | 5,10 | 5,7,14 |
| 8 | 8,6 | 8,10 | 8,16,14 |
| 12 | 12,6 | 12,10 | 12,14 |

Their interiors are disjoint, using only 3,7,16. All path edges are
literal +2 or doubling edges. A planar drawing of this subdivision
would suppress to a planar K_(3,3), impossible because a simple
bipartite planar graph with six vertices has at most eight edges, while
K_(3,3) has nine. Every larger G_N contains G_16, proving (S2).

This is a genuine all-N result derived from two finite certificates.
NetworkX provides an independent check for every N=1,...,64. Its
nearness to the square-sum transition at 15 is currently **only a
numerical juxtaposition**: the graphs have different adjacency predicates,
and neither their Hamiltonian paths nor their obstructions have been
transported. The literal-operation graph is connected for every N:
odd n>1 descends by n-2, and even n descends by halving. Thus even
connectivity has no threshold resembling square-sum existence here.

## 3. K5 and K3,3 measure two different branching patterns

K5 needs five branch sets with every pair adjacent. K_(3,3) needs two
shores of three branch sets and all nine cross adjacencies. These are
the two planarity obstructions in Wagner's minor characterization;
Kuratowski's characterization uses subdivisions instead. The result in
Section 2 directly supplies a subdivision, so needs neither a search
for arbitrary minors nor the full characterization theorem.

For a tournament's underlying undirected graph these tests are automatic:
it is K_N, which contains K5 from N=5 and a K_(3,3) subgraph from N=6
(delete within-shore edges). Consequently these obstructions alone do
not see the tournament orientation, root, dyadic layers, or Collatz
route. In the sparse graph they detect a real arrangement obstruction;
in the complete graph they detect only its order.

Primary planarity source: Kuratowski,
[*Sur le probleme des courbes gauches en Topologie* (1930)](https://www.impan.pl/pl/wydawnictwa/czasopisma-i-serie-wydawnicze/fundamenta-mathematicae/all/15/0/92829/sur-le-probleme-des-courbes-gauches-en-topologie),
DOI 10.4064/fm-15-1-271-283. Wagner's original minor paper is
[*Uber eine Eigenschaft der ebenen Komplexe* (1937)](https://doi.org/10.1007/BF01594196);
the DOI was located, but its full text was not retrievable in this session.
No proof step here relies on an uninspected passage from that paper.

## 4. The relevant tournament analogue of graph minors

**CITED:** Robertson--Seymour's graph minor theorem says that finite
undirected graphs are well-quasi-ordered by minors. It implies a finite
excluded-minor basis for any minor-closed family. It does not assert that
a chosen arithmetic operation terminates or that its certificate
language is minor-closed. See the authors' institutional
[Graph Minors XX record](https://collaborate.princeton.edu/en/publications/graph-minors-xx-wagners-conjecture/),
DOI 10.1016/j.jctb.2004.08.001.

There is a closer positive theorem: **CITED**, Kim--Seymour,
[*Tournament minors* (2015)](https://doi.org/10.1016/j.jctb.2014.12.005),
[author manuscript](https://web.math.princeton.edu/~pds/papers/tourminors/paper.pdf),
proves tournaments well-quasi-ordered under taking a subdigraph and
contracting strongly connected subdigraphs. The branch-set certificate
uses disjoint strongly connected vertex sets and directed edges between
them representing the target arcs. A directed target path lifts because
inside each branch set one can travel from its entry point to its exit.

This is different from the pair-module quotient proposed for halving.
A two-vertex tournament has only one directed edge and is never strongly
connected. Therefore a prescribed partition into pairs is not a strong
minor contraction. A uniform pair module can still be represented by
deleting one vertex from each pair; that realizes the quotient as an
induced subtournament and hence a strong minor. But the minor model then
does not certify that both vertices formed one arithmetic unit. Retain
the pair partition and its uniformity certificate separately.

Likewise, Q[H,H,H,1] admits contraction of its three full H blocks by
this strong-minor rule if H is strongly connected. An arbitrary H need
not meet that hypothesis. Uniform substitution itself does preserve
quotient paths by choosing representatives, even without strong
connectivity; this is another, less information-preserving operation.

The 2011 Chudnovsky--Seymour paper
[*A well-quasi-order for tournaments*](https://doi.org/10.1016/j.jctb.2010.10.003)
instead treats **immersion**. Those two relations must not be merged.

## 5. A concrete sound arithmetic sidecar for contraction

In the inverse ordinary Collatz graph, retain letters

```
D(x)=2x;
O(x)=(x-1)/3, permitted exactly when x=4 mod6.
```

A contracted directed path is certified by its endpoints and its word,
with every intermediate value positive and every O guard checked. A
pure dyadic path can be compressed to its exponent k: its endpoints
satisfy y=2^k x and the omitted vertices are uniquely recoverable.
For example the compressed route from 4 to5 has word DDO and expands
uniquely as `4 -> 8 -> 16 -> 5`. A bare edge 4--5 contains none of that
information. A bare contraction to a four-vertex tournament contains
still less arithmetic information.

There is already a minimal obstruction to using ordinary word deletion
as the arithmetic minor relation: `DDO` is legal from 4, while its
subword `DO` attempts the forbidden O step at 8. Thus the language of
guarded root words is not closed under taking subwords. Retaining only
operation names while deleting intermediate steps does not preserve
the certificate predicate.

The connection contract is therefore:

```
source: a finite guarded arithmetic path or branching certificate;
target: its contracted path/branch-set skeleton;
map: collapse selected dyadic paths or certified blocks;
preserved predicate: existence of the represented directed route;
lost data: lengths, numerical endpoints, orientation if forgotten, O guards;
sidecar: endpoints + dyadic exponents / complete words + guard witnesses;
decisive test: expand and replay every stored word from the root4.
```

This supplies a sound verifier and a compressible class of certificates.
It leaves **OPEN** the coverage theorem that every positive input admits
such a rooted certificate and the construction of those certificates
from tournament structure without first assuming a Collatz orbit.
Well-quasi-ordering supplies a finite-basis method only after a relevant
closed family and a sound relation preserving the arithmetic sidecar
have been established; neither is supplied by cardinality alone.

## 6. Reproduction and scope

Run from the repository root:

```
python 04-computation/experiments/decoder_minors_20260925.py
```

See [script](../../04-computation/experiments/decoder_minors_20260925.py)
and [output](decoder_minors_20260925.out). Core verification is standard
library only: G15 rotation, all nine G16 paths, valuation and grid laws
on every integer1..256, and guarded replay of DDO. Optional NetworkX3.5
independently checks planarity for every N1..64. Positive controls are
the planar prefix/grid and legal dyadic contraction; hostiles are the
K_(3,3) subdivision, orientation-forgetting, and an inadmissible
two-vertex strong branch set. No hidden filters or sampled universes.

The new sparse obstruction changes the concept board as follows: row
scaling removes its crossing mechanism; literal +2 exposes the seam;
tournament completion hides that mechanism by making every pair an
edge; minor compression can retain it only with typed operation labels.
The next research question is a decoder whose chosen contraction rules
retain both the valuation scale and the actual inverse-Collatz guards.
