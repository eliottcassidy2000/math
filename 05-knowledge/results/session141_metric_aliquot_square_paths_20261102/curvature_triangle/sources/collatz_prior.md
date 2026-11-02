# Collatz on a hexagonal surface: exact representation and an unavoidable stretching obstruction

**Session140, 2026-11-02. Status:** PROVED elementary arithmetic and obstruction; CITED classical grid-minor input for the global representation; FINITE-EXACT reproduction controls. Collatz convergence remains OPEN. No literature-priority claim for the elementary growth bound or the graph-theoretic reduction.

## 0. Answer in one paragraph

The standard Collatz graph really can be drawn faithfully as a directed **subdivision** inside the honeycomb lattice, with one outgoing arrow and at most two incoming arrows at each used vertex. This is unconditional, but edges must acquire **unbounded routing lengths**. In fact each root prime to three has at least `2^r` distinct predecessors at exactly `5r` ordinary Collatz steps, whereas a planar lattice has only quadratic ball growth. Thus no injective bounded-step realization exists, not even in the triangular lattice of hexagon centers. The correct object is a port-labelled, non-uniformly stretched directed network, not a unique space-filling trajectory. Rigidity of the abstract graph does not give rigidity of its geometric embedding.

## 1. Inheritance and concept board

Read through `git show` because this worktree is sparse:

* THM-4523, `01-canon/theorems/THM-4523-collatz-functional-graph-is-rigid-backward-trees-identify-vertices.md`: inherited audited rigidity and backward-tree arithmetic; explicitly no connectivity conclusion.
* THM-4623, `01-canon/theorems/THM-4623-collatz-functional-graph-order-topology-endomorphisms-and-koopman-torsion.md`: corrected endomorphism residue restriction; reachability has a greatest element iff Collatz holds. Its standard and shortcut maps must not be interchanged.
* THM-4754, `01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md`: exact 3-adic affine/tournament dictionary, not a planar lattice realization.
* THM-4505, `01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md`: the powers-of-two-or-three **sum** graph is a different graph; the note explicitly does not address Collatz convergence.
* `05-knowledge/results/collatz_mod6_20260917_reverse_tree_pieces.md`, section 3: the pruned shortcut inverse tree, modulus-nine port types and classical growth context. Its averaged `(4/3)^k` is not used as a pointwise lower bound here.
* `05-knowledge/results/collatz_coherence_20260930_trichotomy_hexagon_motifs.md`, sections 2 and 5: the genuine six-cycle is multiplication by two on `(Z/9Z)^*`; local indegrees already recorded.
* `05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md`: inverse closure is exactly the root basin, not proof that the basin exhausts the positive integers. This is the hostile to a surjectivity claim based on a picture.

Anchor: make the proposed lattice model faithful. Niche: uniform inverse growth versus geometric capacity. Wildcard: the arithmetic six-cycle as a port selector.

| Live object | Preserved information | Missing coordinate / decisive test |
|---|---|---|
| Standard functional graph | Arithmetic successor and inverse ports | Does a proposed map preserve every original edge? |
| Honeycomb vertices / triangular centers | Local degree and planar coordinates | Odd-cycle and ball-growth tests |
| Subdivided routing | Incidence and direction | Original clock, branch markers, path lengths |
| 3-adic port types | Legality and survival of inverse branches | Division consumes precision; mod9 is not a closed quotient |
| Graph rigidity | Unique arithmetic labelling of an isomorphic graph | No metric, coordinates, or unique embedding follows |

## 2. Three maps, three different local pictures

All domains below are positive integers, except that `U` uses positive odd integers.

| Map | Forward rule | Exact inverse fibre over m | Known positive cycle |
|---|---|---|---|
| Ordinary `C` | even `n/2`; odd `3n+1` | `{2m}`, plus `(m-1)/3` iff `m=4 mod6` | `1 -> 4 -> 2 -> 1` |
| Shortcut `T` | even `n/2`; odd `(3n+1)/2` | `{2m}`, plus `(2m-1)/3` iff `m=2 mod3` | `1 -> 2 -> 1` |
| Accelerated odd `U` | `(3n+1)/2^a`, `a=v2(3n+1)` | `(2^a m-1)/3`, all `a>=1` with `2^a m=1 mod3` | `1 -> 1` |

For `U`, the fibre is empty if `3|m`, and infinite otherwise; the displayed admissible values are automatically positive odd integers. Thus the assertion "at most two incoming edges" applies to `C` and `T`, **not** `U`. An infinite-degree accelerated vertex cannot become a vertex of any finite-degree lattice by merely subdividing its incident edges. One must split the vertex into a port structure, for example recover the shared dyadic tower of even states. Re-inflating accelerated edges shares intermediate states; it is not merely independent subdivision of each accelerated edge.

For `C`, all incoming ports are distinct from each other and from the outgoing neighbour. There are no positive fixed points or two-cycles: the parity cases for `C^2(n)=n` give only `n=0` or negative values. The underlying simple graph therefore has degree two or three.

Honeycomb vertices have degree three and form a bipartite graph; the centers of the hexagons form the degree-six triangular lattice. The triangle `{1,2,4}` already excludes an edge-for-edge honeycomb embedding of `C`. It does not exclude a triangular-lattice embedding by itself, nor a honeycomb embedding after subdivisions. The growth theorem excludes bounded-step versions of both.

**A literal hexagon does encode the known cycle after subdivision.** In brick-wall coordinates take the six-cycle

```
(0,0) -> (1,0) -> (2,0) -> (2,1) -> (1,1) -> (0,1) -> (0,0).
```

Mark `(0,0)` as `1`, `(2,0)` as `4`, and `(1,1)` as `2`; the other three vertices are routing vertices. Each ordinary step around `1 -> 4 -> 2 -> 1` now takes two lattice edges. The unused port at `(2,0)` can receive `8` from `(3,0)`. This is an exact small survivor of the proposed picture, not a proof that all arithmetic vertices fit with this two-edge clock.

## 3. Sharp uniform five-step binary predecessor supply

Let `S={n>=1:3 does not divide n}`. For `n in S` choose `k=k(n)` by

```
n mod9:   1  2  4  5  7  8
k:        2  1  2  3  4  1.
```

Then `2^k n=4 or7 mod9`. The inverse moves are `D(x)=2x` and `O(x)=(x-1)/3`, with `O` legal for `x=4 mod6`.

### Theorem 1 (PROVED)

For every `n in S`, the two distinct integers

```
F0(n)=32n,
F1(n)=(16n-2^(4-k(n)))/3
```

belong to `S`, exceed `n`, and satisfy `C^5(Fi(n))=n`. Consequently

```
|C^(-5r)(n) intersect S| >= 2^r       for all n in S and integers r>=0.
```

**Proof.** The first predecessor uses five doublings. For the second, starting at `n`, take `k` inverse doublings, then `O`, then `4-k` inverse doublings. The `O` target is even, is `1 mod3`, and is `4 or7 mod9`; hence the odd predecessor is a positive integer not divisible by three. The resulting endpoint is exactly `F1(n)`. Both endpoints are positive and prime to three; `F1(n)<32n`. Also

```
F1(n) >= (16n-8)/3 > n.
```

At a fixed level, distinct parents have disjoint five-step fibres because `C^5` is a function. Inductively, selecting these two children at each parent produces `2^r` distinct endpoints after `r` blocks. This argument does not assume that any orbit converges. QED.

The selected macro-vertices also do not repeat across levels: any cross-level collision, pulled forward by the shorter number of blocks, would put the root at a positive depth of a strictly increasing selected chain. Nevertheless, the five ordinary steps of a macro-edge need not form a simple path near the known cycle. We claim a binary **macro-tree and fibre count**, not that these particular five-step walks are already an internally disjoint binary-tree subdivision.

The clock five is the smallest uniform fixed clock that guarantees two predecessors **remaining in S**: for root `7`, `C^(-4)(7)={112,18}`, with just one member in `S`; earlier fibres also have only one. This is a sharp statement about this uniform unit-surviving requirement, not an optimal asymptotic growth exponent.

### Shortcut companion (PROVED)

Put `j=k-1`, so `j` is in `{0,1,2,3}` and `2^j n=2 or8 mod9`. Using `j` doublings, the inverse odd move `E(x)=(2x-1)/3`, and `3-j` doublings gives the same `F1(n)`. Together with `16n` this proves

```
|T^(-4r)(n) intersect S| >= 2^r.
```

Four is sharp in the same sense: `T^(-3)(7)={56,18}` has only one unit modulo three.

### Independent concurrent check

The Hadwiger--Nelson lane independently uses odd `m in S`: precisely two exponents `a in {1,...,6}` satisfy `2^a m=4 or7 mod9`, yielding two odd `U`-predecessors in `S`, each within seven ordinary steps. The root coordinator's stricter version chooses `a in {3,5,7}` when `m=2 mod3`, or `{4,6,8}` when `m=1 mod3`; exactly two children survive and both strictly exceed `m`. This gives a genuinely increasing odd predecessor tree with at most nine ordinary steps per edge. These independent constructions validate the same mechanism; the equal-clock five-step bound above is stronger for the ball-growth application.

## 4. Every component has exponential balls; planar local coding is impossible

Use the ordinary undirected graph distance, retaining all standard edges.

### Theorem 2 (PROVED)

Every connected component of the standard Collatz graph has exponential ball growth. More precisely, for `n in S`,

```
|B_C(n,5r)| >= 2^r.
```

For arbitrary `x`, there is `z in S` at forward distance `t<=v2(x)+1` (or `t=0` if already `x in S`), and

```
|B_C(x,t+5r)| >= 2^r.
```

**Proof.** The first assertion follows from Theorem 1: each endpoint has an ordinary walk of length `5r` to its root. If `3|x`, divide out its powers of two to obtain an odd `m`, then take the odd step to `3m+1`, which is prime to three. This does not require eventual convergence. QED.

### Theorem 3 (PROVED; quantitative obstruction)

There is no injection of the vertex set of any entire Collatz component into either the honeycomb lattice or the triangular lattice for which all standard edges have image endpoints at lattice graph distance at most a common finite `L`. In particular, there is no edge-for-edge embedding, no bounded-length subdivision routing, and no bounded-multiplicity local coding of this form.

For an injection defined just on the finite witness formed by the selected predecessor walks of `r` blocks rooted at `n in S`, let `D_r` be the largest image distance across any ordinary edge of that witness. Then

```
D_r >= (2^(r/2)-1)/(10r),        r>=1.             (1)
```

For `T`, the denominator is `8r`.

**Proof.** In integer coordinates, every edge of either the brick-wall honeycomb or the triangular lattice changes each coordinate by at most one. Thus a lattice ball of radius `R` has at most `(2R+1)^2` vertices. The images of the `2^r` selected endpoints lie within distance `5r D_r` of the root image, so

```
2^r <= (10r D_r+1)^2.
```

This is (1), which diverges with `r`. Any uniform bound gives a contradiction. For an arbitrary component use Theorem 2's bounded root offset. For fibre multiplicity at most `M`, the same proof gives `2^r<=M(10r D_r+1)^2`; every fixed finite `M` still fails. QED.

**General capacity corollary.** The identical argument rules out a uniformly bounded-edge, bounded-multiplicity coding into *any* target with a uniform polynomial ball bound `|B(v,R)|<=A(1+R)^d`. Thus moving to a fixed higher-dimensional lattice does not repair local injectivity. More quantitatively, a uniformly delta-separated placement in Euclidean `R^d` obeys the packing bound `2^r<=(1+10r D_r/delta)^d`, hence `D_r>=delta(2^(r/d)-1)/(10r)`. The packing proof places disjoint radius-`delta/2` balls around the endpoints inside the root ball enlarged by `delta/2`. If one drops uniform separation and allows accumulation/arbitrarily small geometric scales, this obstruction no longer applies; it is not a prohibition on arbitrary topological plane drawings.

This is stronger than a planarity obstruction: the standard graph is planar. It is a **capacity-per-radius obstruction**. It survives removal of the known cycle: for example the full backward tree rooted at `5` contains no cycle (`5` reaches `1` and never returns to `5`), yet satisfies the same growth bound. It also applies inside hypothetical nonconvergent components.

## 5. Faithful positive repair: an unbounded subdivision lives in the honeycomb

### Theorem 4 (PROVED using the cited grid-minor theorem)

The standard Collatz digraph has a subdivision isomorphic to a directed subgraph of the infinite honeycomb lattice. Each used vertex has exactly one outgoing edge and at most two incoming edges. Original vertices are distinguished branch markers, and contracting the subdivided directed paths recovers exactly the original graph. In every such realization the edge-routing lengths are unbounded.

**Step A: generic functional-graph structure.** A finite undirected cycle in a functional graph must be a directed cycle: its `q` edges require `q` outgoing incidences from its `q` vertices, so each cycle vertex uses its unique outgoing edge on that cycle. Two different cycles cannot be in one weak component. Indeed adjacent vertices have forward orbits that meet, hence the same is true along any finite connecting path; different directed cycles never meet. Consequently each component is a tree or a tree attached to one finite cycle. For standard Collatz the underlying graph is countable, simple and subcubic. This uses no assertion about which components actually occur.

**Step B: an accumulation-free planar drawing.** Root a tree arbitrarily. Each graph-distance level is finite. Assign the children disjoint subintervals of their parent's interval inside `(0,1)`, place level-d vertices at height `d`, and join parent to child by straight segments in the corresponding strips. Edges in separate parent intervals are disjoint, so this draws the tree with no vertex accumulation point. For a unicyclic component, draw the cycle as a horizontal polygonal chain closed by a lower arc; attach its finitely many trees in disjoint upward strips. Put countably many components in strips `(2i,2i+1)`. Every bounded region then contains only finitely many vertices. No component has an infinite family of nested cycles.

**Classical input (CITED, primary text checked 2026-11-02).** Daniela Kühn, *Minor-universal planar graphs without accumulation points*, Theorem 2, characterizes minors of the infinite square grid. In particular, a graph with a vertex-accumulation-free plane drawing and no component containing infinitely many nested cycles is such a minor. This applies to Step B. The result is used as a theorem, not proved anew: [author-hosted primary PDF](https://web.mat.bham.ac.uk/D.Kuehn/unvacc3.pdf), Theorem 2, printed pages 4--5.

**Step C: explicit square-grid minor of the honeycomb.** Use the brick-wall honeycomb `H` with vertex set `Z^2`, all horizontal edges, and vertical edges `(x,y)--(x,y+1)` exactly when `x+y` is even. These vertical edges form a perfect matching. Contract each matching edge, and label it by its lower endpoint `(x,y)`, where `x+y` is even. Its four neighbours are `(x+-1,y+-1)`, independently. The coordinates

```
(a,b)=((x+y)/2,(x-y)/2)
```

turn these into `(a+-1,b)` and `(a,b+-1)`, exactly the square grid. Hence the Collatz graph is a minor of `H`.

**Step D: subcubic minor models are subdivision models.** In each connected branch set of the minor model choose the at most three attachment vertices used by incident model edges. Connect them by a finite minimal tree. With three attachments, its median gives a branch vertex with internally disjoint paths to the attachments; with one or two attachments choose the corresponding endpoint or point of their connecting path. Repeated attachments cause only zero-length arms and no problem. Branch sets are disjoint. These finite paths, joined across the model edges, give a subdivision of the original subcubic graph inside `H`. Orient each path in the original Collatz direction. Original vertices retain one outgoing edge; new subdivision vertices have one incoming and one outgoing edge. Theorem 3 forces unbounded path lengths. QED.

**Shortcut boundary.** The same representation can be made for `T` if its directed two-cycle is retained as two distinct oppositely oriented arcs, rather than collapsed to one undirected edge. First insert a separate subdivision vertex on each arc, making a simple four-cycle, then apply the same argument. Routing both arcs along the same subdivided path would incorrectly create internal vertices with two outgoing edges; that is not the construction here.

The theorem asserts a **subgraph**, not an induced or spanning subgraph, and not a computable canonical layout. Unused lattice edges are not Collatz transitions. Every original vertex is represented, but not every lattice vertex need be used. The proof preserves unknown components instead of silently deleting them.

## 6. What is unique, and what the hexagon really records

**Rigidity scope audit.** THM-4523's backward-tree separation is restricted to positive integers prime to three. It is false on all positive integers: every multiple of three has only the doubling inverse ray, so the rooted backward trees of `3`, `6`, `9`, and so on are isomorphic. The extension from restricted separation to whole-graph rigidity uses forward structure. An automorphism first fixes every 3-unit by the separation theorem. For an odd multiple `m` of three, `w=3m+1` and `2w` are 3-units and fixed; the two preimages of `w` are `2w` and `m`, so `m` is fixed. Every `2^j m` is then fixed along the unique inverse ray above `m`. This proves the extension, conditional only on the inherited audited 3-unit separation lemma; that long 3-adic lemma is not re-proved in this packet. None of Theorems 1--4 depends on it.

With this restriction and extension understood, THM-4523 gives uniqueness of the isomorphism between an abstract graph isomorphic to the standard Collatz graph and its arithmetic-labelled version. It does not prescribe coordinates, lattice positions, angular order, edge lengths or path choices. Even literal geometric uniqueness fails under translation of any realization. We make no new theorem about uniqueness modulo Euclidean symmetries.

No single forward Collatz orbit visits all positive integers: any orbit that ever visits `1` is thereafter trapped in `{1,4,2}`, so has only finitely many distinct earlier values; an orbit never visiting `1` plainly omits it. A spanning directed walk/space-filling orbit is therefore not the same object as the full branching functional graph. A geometric space-filling curve would require an additional traversal or interpolation rule, not supplied by the outdegree condition.

There is a genuine local arithmetic hexagon:

```
(Z/9Z)^* under doubling: 1 -> 2 -> 4 -> 8 -> 7 -> 5 -> 1.
```

The two safe ordinary odd-predecessor gates are residues `4,7`; this six-cycle is precisely what guarantees a surviving gate within four doublings. It is an arithmetic state cycle, not a planar honeycomb tile. Higher precision is essential: `1` and `10` have the same residue modulo nine, but `F1(1)=4` and `F1(10)=52=7 mod9`. Thus the chosen inverse move's next mod-nine state is not determined by the current mod-nine state. Dividing by three consumes a 3-adic digit; the full integer or a higher-precision sidecar must remain.

The faithful transport is:

```
source: C on positive integers
map: original vertices -> marked lattice vertices; edges -> directed finite paths
target: a directed honeycomb subgraph
preserved: incidence, direction, components, directed cycles, reachability
lost if markers/lengths are forgotten: one-step timing and original versus routing vertices
sidecar: branch markers, arithmetic labels or their rigid recovery, path/clock data
cost: routing lengths necessarily unbounded
not supplied: convergence, a spanning walk, a canonical or unique geometric embedding.
```

The transferable principle is **local degree compatibility is not global metric capacity**. One must test growth, clocks and quotient losses in addition to matching a three-port picture. This is a concrete instance of the repository's existing controlled-forgetting methodology, not a new claim equating the major conjectures.

## 7. Verification and stopping boundary

Run:

```
python3 verify.py > verification.json
python3 -O verify.py > verification-optimized.json
cmp verification.json verification-optimized.json
```

Both modes give identical PASS JSON. Checks: forward/inverse formulas through `n=10000`; independently enumerated inverse fibres through `n=500`; sharp-clock root `7`; selected macro-trees to depth fourteen for six roots in both maps (`16384` distinct leaves each); both independent odd-predecessor constructions below `m=2000`; exact square-grid matching contraction in a finite coordinate window. All assertions use explicit raising and remain active under `-O`.

The finite tests support implementation, not the all-height proof, and do not verify Kühn's theorem computationally. No assertion of Collatz convergence is made. The proof task stops at a faithful unconditional representation and a sharp fixed-clock local supply leading to an explicit nonlocality lower bound.
