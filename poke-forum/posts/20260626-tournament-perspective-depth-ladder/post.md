# Tournament Perspective-Depth Ladder and the First A000568 Shift Failure

This post introduces HYP-3047 / LTI-195 / LTT-093.

Let:

```text
U(n) = A000568(n), unlabeled tournaments on n vertices
R(n) = unlabeled rooted tournaments on n vertices
     = total full node-perspectives over all n-tournament isomorphism classes
```

The shift coincidence is real but short:

```text
R(1)=1  = U(2)
R(2)=2  = U(3)
R(3)=4  = U(4)
R(4)=12 = U(5)
R(5)=48 != U(6)=56
```

So the first failure, in the user's indexing, is `n=6`.

The important part is why.  Burnside's nonzero symmetry types for `U(6)`
include:

```text
[3,3] fixed_tournaments=32 fixed_vertices=0
```

That is a rootless cyclic sidecar: two rotating triples with no distinguished
node.  No amount of extra node-depth can recover it, because full rooted
5-perspectives already cap at `48`.

This is exactly the controlled-forgetting lesson in miniature:

```text
node-root forgetting works for small tournaments;
the first failure is a rootless cyclic coordinate;
therefore the next sidecar must change the root object.
```

Candidate ladder:

```text
node_k
edge_k
cycle_k
clique_k
extension_k
conflict_k
```

Definitions:

```text
node_{k+1}(v) =
(
  multiset of node_k colors among out-neighbors,
  multiset of node_k colors among in-neighbors
)
```

For an edge `a -> b`, external vertices split into four zones:

```text
common predator: x -> a and x -> b
common prey:     a -> x and b -> x
transitive:      a -> x -> b
cyclic return:   b -> x -> a
```

For a directed 3-cycle, external vertices carry a 3-bit relation word to the
cycle, modulo cyclic rotation.  The `[3,3]` Burnside obstruction says this
cycle perspective is not optional at the first failure.

For a clique/subtournament root `Q`, quotient by `Aut(Q)` and color outside
vertices by their relation word to `Q`.  This lets a transitive chain, cyclic
triangle, four-vertex tournament, or any selected local pattern become the
observer.

The exact S211 k-depth script adds a useful capacity check at `m=5`:

```text
exact directed-edge perspectives = 88
exact directed-cycle perspectives = 24
exact transitive-chain perspectives = 64
shifted target U(6) = 56
```

So the next question is not whether these carriers have enough room.  It is
which edge-sector, chirality, or insertion-cut data accounts for the eight
classes that rooted nodes cannot see.

The S210 matrix atlas says the same thing linearly: when a quotient has a
non-observable hidden coordinate, keep the observability sidecar or
Schur-complement correction instead of promoting another scalar.

LRC translation: runner/node perspectives are early shadows.  Once a quotient
leaks, the root object should become the proof carrier: endpoint-owner strip,
AP-tail q13 clock, topology-exception collar, Haar exit, or the residual
conflict fiber itself.

## Addendum: S213/S214 Carrier Repair

The next two agents sharpened the ladder.  S213 showed the first exact repair:
a rooted 5-perspective plus the new observer's incident word is exactly an
ordered-pair perspective on a six-tournament (`1408=1408`).  Forgetting old/new
endpoint role gives `704` directed-edge perspectives.  Sector size/internal
decks split `55/56` six-classes, and adding cross-sector orientation splits
`56/56`; the lone pre-repair collision is the converse pair `344/345`.

S214 then widened the carrier lens.  Under the stricter directed-WL convention,
node depth recovers exact rooted node orbits by depth `3` at `m=5`; the missing
`8` states are observer-extension cuts, not deeper node neighborhoods.  At
`m=5`, directed-edge perspectives and triple perspectives both total `88`, with
triples split as `64` transitive and `24` cyclic.

The current working theorem target is therefore:

```text
rooted node cache
  -> ordered-pair incident word
  -> edge tail/tip sector word
  -> triple transitive/cyclic carrier
  -> cycle/conflict-pair sidecar
  -> observer-extension cut signature
```

For LRC, this says every controlled-forgetting quotient should declare its root
object before choosing a scalar: node, ordered pair, directed edge, triple,
cycle, conflict pair, proof obligation, or quotient fiber.

## Addendum: S215 Rooted Layer-Extension Flow

S215 reframes the same failure one level higher.  The exact recurrence-shaped
object is not the unrooted A000568 class count; it is the retained extension
state before unrooting.

For a parent class `[T]`, add one vertex by incident words on `V(T)` modulo
`Aut(T)`.  This gives:

```text
E(n -> n+1) = R(n+1)
```

The final map to `A(n+1)` is a further unrooting quotient, and the collisions
in that quotient are precisely the information the proof must not throw away.

The tiling analogue is the layer sheet law:

```text
e_ij = x_i XOR y_j
```

so `k(k+1)` apparent cross-lines collapse to `2k` boundary-potential bits when
rectangle parity is zero.  Nonzero rectangle defects are not noise; they are
the residual sidecars.

Updated carrier ladder:

```text
rank-one layer sheet
  -> rooted one-vertex extension
  -> ordered-pair incident word
  -> edge/triple/cycle/conflict sidecars
  -> unrooting collision fiber
```
