# Complement, divisibility and the sharp exponent-three graph bridge

25 October 2026. **PROVED elementary:** the complement/monotonicity obstruction,
the explicit exponent-three bridge, and the resulting classification below.
**FINITE-EXACT:** independent controls in [graph_bridge_check.py](graph_bridge_check.py)
and [graph_bridge_receipt.json](graph_bridge_receipt.json). No literature novelty
claim or transfer of an open tournament problem is made.

## Inheritance: the eight-state repair was already solved

The original chart and graph/path coordinates are in
[`trijection_reset_20261020/arithmetic.md`](../trijection_reset_20261020/arithmetic.md)
and [`graph_audit.md`](../trijection_reset_20261020/graph_audit.md).
The strongest directly relevant predecessor is
[`black_detours_closed_blue_lines_klein_dictionary_20261020.md`, section 5](../black_detours_closed_blue_lines_klein_dictionary_20261020.md).
Its K1-K2 already prove the action obstruction and give the eight-state repair,
including 32 equivariant bijections, two monotone bijections, and two dictionaries
reading both signs. The no-map result and minimal twelve-state resolution are
also proved in
[`divisor_reflection_resolution_20261021/REFLECTION_RESOLUTION.md`](../divisor_reflection_resolution_20261021/REFLECTION_RESOLUTION.md),
while [THM-4733, divisor-plane versus tiling cube](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md)
gives the folded four-slot carrier. These are recovered results, not new
discoveries of this session.

The canonical hostile is `pq,pr`, fixed by complement followed by prime
exchange. The corrected near miss is to equate action-preserving encodings
with divisibility-preserving encodings. The sidecar retained below is every
ordered edge coordinate. The board is **fixed points / divisibility /
complement / threshold coordinates / marked-path closure**. The new question
is what extra obstruction appears when the divisor order is requested as well.

## 1. Objects and the recovered action obstruction

For `k>=1`, let `D_k` consist of the proper nontrivial divisors of `p^k q r`,
written `x=(a,b,c)`, with `0<=a<=k`, `b,c in {0,1}`, excluding bottom and top.
Divisibility is coordinatewise order. Let

```text
C(a,b,c)=(k-a,1-b,1-c),       Q(a,b,c)=(a,c,b),
f_k=2a-k(b+c),                g=b-c.
```

On graphs on the ordered vertices `0,...,m-1`, let `K` be full complement
inside `K_m`, and `R` be vertex reversal `i -> m-1-i`. We seek maps satisfying
`Phi C=K Phi`, `Phi Q=R Phi`. Here **monotone** means
`x divides y => Phi(x) subset Phi(y)`; an **order embedding** requires this
implication in both directions and injectivity.

For even `k`, `CQ` fixes `(k/2,1,0),(k/2,0,1)`. In the graph carrier `KR`
has no fixed point: reversal fixes the unordered pair `{0,m-1}`, whose bit
would have to equal its complement. Therefore no equivariant map exists,
even without injectivity or monotonicity. This is the inherited obstruction.
For odd `k`, action-only embeddings into sufficiently large ordered graph
charts exist by matching the two stabilizer types; the cited resolution
packet gives the exact capacity criterion.

For `k=2`, deleting `pq,pr` gives the inherited eight-state map. In the
backward-arc graph gauge, one of the two monotone maps is

```text
p -> {},          p^2 -> {02},       qr -> {01,12},
pqr -> {01,02,12}, q -> {01},         r -> {12},
p^2q -> {01,02},  p^2r -> {02,12}.
```

The other exchanges the two outer edge coordinates. This carries divisor
complement to full graph complement and prime exchange to order reversal.
It is not an order embedding: the empty image of `p` lies below the image
of `q`, although `p` does not divide `q`. The predecessor's tables use
*forward*-tile bits, whereas this packet uses graph edges for *backward*
arcs; closedness must be read in the stated gauge.

## 2. A general complement-order obstruction

**Proposition G1.** Let a poset carry an order-reversing involution `C`. If
`Phi` is monotone into subsets of a fixed universe and
`Phi(Cx)=E minus Phi(x)`, then every `x<=Cx` has `Phi(x)=empty`.

Indeed `Phi(x) subset Phi(Cx)=E minus Phi(x)`. A set contained in its own
complement is empty. This argument requires neither a graph nor an injection.

For `D_k`, the points satisfying `x<=Cx` are exactly

```text
p, p^2, ..., p^floor(k/2).
```

The binary coordinates must be zero, and `2a<=k`; bottom is excluded.
Consequently:

1. For every `k>=4`, **no injective monotone map** to any Boolean edge space
   can intertwine `C` with plain complementation. In particular all odd
   `k>=5` pass the action test but fail this stronger test: `p` and `p^2`
   must both map to the empty graph.
2. For every `k>=2`, **no order embedding** can intertwine those operations.
   The image of `p` is empty and hence is contained in the image of `q`,
   contradicting order reflection because `p` does not divide `q`.
3. These obstructions already use `C` alone. Changing the vertex reflection,
   increasing the graph size, or adding coordinates while keeping full edge
   complementation cannot repair them.

The first statement is about injective monotonicity; the second is strictly
stronger. The exponent-three construction below is a witness to that distinction.

## 3. The sharp positive case: fourteen divisors on four graph vertices

**Proposition G2.** For `k=3`, the following six edge predicates define an
injective monotone map with `C -> K` and `Q -> R`:

| Graph edge | Present exactly when |
|---|---|
| 01 | `a+2b >= 3` |
| 02 | `b=1` |
| 03 | `a>=2` |
| 12 | `a+b+c>=3` |
| 13 | `c=1` |
| 23 | `a+2c>=3` |

All coordinates are increasing in the exponents, proving monotonicity.
Complement sends `a` to `3-a`, sends each binary coordinate to its opposite,
and replaces each of `a+2b`, `a+2c`, `a+b+c` by five minus itself. Threshold
three is therefore complemented exactly. Prime exchange swaps edges `01,23`
and `02,13`, while fixing `03,12`; this is precisely reversal `i -> 3-i`.

For injectivity, read `b,c` from edges `02,13`. If `b!=c`, the two exterior
thresholds together with edge `03` record `a>=1,a>=2,a>=3`, recovering `a`.
If `b=c=0`, the allowed values are `a=1,2,3`, distinguished by `03` and
either exterior threshold. If `b=c=1`, the allowed values are `a=0,1,2`,
distinguished by `12` and `03`. Thus the map has an explicit inverse on its
image. Four graph vertices are minimal because three vertices carry only
eight graphs, fewer than the fourteen source divisors.

For `k=1`, the familiar squarefree map on three graph vertices is
`(edge01,edge02,edge12)=(b,a,c)`. It is an order embedding and respects
both actions. Combining this with G1 and the even-exponent fixed-point
obstruction proves the exact classification:

```text
For k>=1, D_k has an injective monotone ordered-graph encoding
with C -> full complement and Q -> order reversal iff k=1 or k=3.
It has an order embedding with those laws iff k=1.
```

The positive `k=3` map is not onto the full graph cube; its image contains
fourteen of the sixty-four graphs on four labelled vertices. It creates
extra comparisons, already `Phi(p)=empty subset Phi(q)`.
The proper-divisor domain is essential: evaluating the same formula at
the omitted unit gives the same empty image as `p`, and evaluating it at
the omitted top gives the same full image as `p^2qr`. It does not extend
to an injective map on the whole divisor lattice.

## 4. Predicates transported to marked five-vertex tournaments

Lift a graph on `0,1,2,3` to a tournament with directed marked path
`0 -> 1 -> 2 -> 3 -> 4`: graph edge `{i,j}` reverses the free arc
`i -> j+1`. The inherited graph/path trijection now yields exact arithmetic
readings for the exponent-three map:

```text
marked path closes to a cycle   iff a>=2;
tournament is strong           iff a>=2 or b=c=1;
marked path is blue            iff b=c;
edge02 minus edge13            =   b-c = g;
tournament is transitive       iff divisor is p.
```

Closure reads the apex edge `03`. For strength, apply the inherited
interval-cover criterion: each edge `{i,j}` covers path cuts `i,...,j`.
If `a>=2`, the apex covers all cuts. If `a<2`, cut zero requires `b=1`
and cut three requires `c=1`; when both hold, edges `02,13` cover every
cut. Blue means reversal-fixed; equivariance and injectivity give exactly
the `Q`-fixed divisors. The graph is empty exactly at `p`, which is the
transitive criterion in this gauge.

Thus this is more than a relabelling of a fourteen-element set: it preserves
both reflections and all divisibility implications, and it gives explicit
tests for three marked-tournament predicates. It still does not transport
divisibility as an equivalence or preserve arbitrary Hamiltonian paths.

There is a uniform limit on combining order with the old sign-reading
dictionary. For any `k>=2`, `p divides Cp`, while `f_k(p)=2>0` and
`f_k(Cp)=-2<0`. Any encoding that reads positive `f_k` as the presence of
the closing edge would require that edge to disappear along a divisibility
relation. Hence no monotone encoding in this backward-arc gauge can also
read `sign(f_k)` as closed/open. The other bit convention cannot evade the
conflict: `q divides p^kq` has `f_k(q)=-k` and `f_k(p^kq)=k`, so a bit
reading negative rather than positive `f_k` would disappear on that chain.
Thus neither sign convention is monotone. These two chains supply an
all-order mechanism for the predecessor's finite `k=2` sign-versus-monotonicity
conflict.

## 5. Operation typing and the positive incidence alternative

On marked tournaments, graph reflection is
`(T,P) -> (T^op, reverse(P))`. Graph complement followed by reflection is

```text
(T,P) -> (reverse just the arcs of P in T, reverse(P)),
```

followed by path normalization. Reversing the marked sequence alone while
keeping `T` fixed produces an illegal directed path. This is the operation
proved in THM-4733 Proposition A; the shorthand `(T,P)->(T,P^rev)` in the
earlier Klein note omitted the required arc reversal.

The companion ideal/star construction preserves divisor order in both
directions by representing the `p` exponent as an initial segment of a
marked chain of edge slots. There `C` becomes relative host complement
*followed by chain reflection*. This extra permutation is exactly what
G1 shows cannot simply be dropped while retaining the full order.

## 6. Reproduction and scope

```text
python -B 05-knowledge/results/reflection_bridge_20261025/graph_bridge_check.py --output 05-knowledge/results/reflection_bridge_20261025/graph_bridge_receipt.json
python -B -O 05-knowledge/results/reflection_bridge_20261025/graph_bridge_check.py
```

The checker uses the full `D_1,D_3` maps, every divisibility pair in those
universes, direct arc reachability for tournament strength, fixed and
forced-empty sets for `k=1..12`, and all `8!` possible bijections for the
inherited eight-state repair. Every check raises an explicit exception;
normal and optimized receipts agree. All-exponent conclusions rest on the
proofs, not the bounded checks. No external paper or historical novelty
claim is needed for these elementary propositions.
