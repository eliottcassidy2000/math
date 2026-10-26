---
id: THM-4069
title: "Even-graph tree-diameter classification and canonical all-cycle envelope"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The image on
  isomorphism classes of the fundamental-cycle cube depends on its spanning
  tree: already at n=4 the star gives P3 while the path gives K3. For every
  n>=3 the star image is connected bipartite with chromatic number 2, and a
  tree has only odd fundamental cycles exactly when it is a star. The
  complete basis-response law is E_(n,T)=union_(3<=k<=diam(T)+1)R_(n,k),
  where R_(n,k) toggles a k-cycle. Hence there are exactly n-2 distinct
  spanning-tree images, strictly nested by tree diameter. The canonical
  repair is the all-simple-cycle envelope. It is the union over all
  spanning-tree images and equals the path-tree image for every n. Thus the
  historical path-based E_n graph object survives
  canonically, and computations explicitly verified on that edge image keep
  their prior status; only the claim that an arbitrary spanning-tree
  fundamental-cycle image gives the same adjacency is retracted.
source: codex-frontier-synthesis-creative-20260825c / recovered even-graph dual lane
audit: >
  PASS. An independent quotient-level proof audit verifies the complete
  diameter classification and strictness witness; a fresh exhaustive scan of
  every labelled tree through n=6 confirms that the image depends only on
  diameter. The primary audit constructs cycle spaces from independent star and
  path bases, canonicalizes every Eulerian graph through n=6, compares both
  images with the all-cycle envelope, scans all 18,247 labelled trees through
  n=7 for the odd-basis/star equivalence, and verifies the generator-union
  identity through n=6. The independent audit instead scans all edge subsets,
  all (n-1)-edge tree subsets, and recognizes cycles by connected 2-regular
  support. Its n=3..6 tree counts are 3,16,125,1296 and its two constructions
  agree. Both normal/optimized output pairs byte-match; neither script uses
  Python assertions or floating-point literals.
depends_on: []
related:
  - THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay
  - THM-1430-graph-switching-is-exactly-E-n
  - THM-1945-the-invariant-monoid-orbit-dictionary
script: 04-computation/even_graph_cycle_envelope_thm4069.py
output: 05-knowledge/results/even_graph_cycle_envelope_thm4069.out
independent_audit_script: 04-computation/even_graph_cycle_envelope_thm4069_independent_audit.py
independent_audit_output: 05-knowledge/results/even_graph_cycle_envelope_thm4069_independent_audit.out
script_sha256: b99b6831212b9cbba3ace8afa8d80a62a4be978a62592ea65f319ef5e753b917
output_sha256: 4c211e8d5e438b6bb9a6e58cc342b881a540ea503a706b2358d59190957a5dfe
independent_audit_script_sha256: fa65cf189ba7d5e01a981f544e02c46965e6f88d13298ada098eb48a5d5fde41
independent_audit_output_sha256: 5951766166f681e4b2dac21080572a355100032eabbbeef6f264d69d2687b364
hash_basis: raw LF bytes
---

# THM-4069 -- the even-graph dual has a canonical envelope, not a canonical basis

**PROVED in the all-`n` scopes below.** The vertex set used by the historical
even-graph metagraph is canonical. Its path presentation also happens to give
a canonical adjacency, but for a subtler reason than previously stated:
quotienting a fundamental-cycle cube by graph isomorphism is not independent
of the spanning tree. The path works because it contains one generator from
every simple-cycle orbit.

## 1. The typed objects

For `n>=3`, let

```text
Z_n={F subset E(K_n): deg_F(v) is even for every v}.              (1)
```

Thus `Z_n` is the binary cycle space of `K_n`. Let `T` be a spanning tree.
For every chord `e=uv` of `T`, write

```text
C_T(e)=e union P_T(u,v),                                         (2)
B_T={C_T(e):e notin T}.                                          (3)
```

Here `P_T(u,v)` is the unique tree path. The fundamental cycles `B_T`
form a basis of `Z_n`, so the labelled Cayley graph

```text
Q_T=Cay(Z_n,B_T)                                                  (4)
```

is the cube `Q_m`, where `m=binom(n-1,2)`. This standard basis fact is also
immediate by chord coordinates: `C_T(e)` contains chord `e` and no other
chord, so the cycles are independent; subtracting the corresponding cycles
from an Eulerian graph leaves an Eulerian subgraph of a tree, hence zero.

Let `q:Z_n -> Z_n/S_n` send a labelled Eulerian graph to its graph-isomorphism
class. Define `E_(n,T)` to be the **simple edge image** of `(4)` under `q`:

```text
{[F],[F xor C]} is an edge
iff F in Z_n and C in B_T give distinct classes.                  (5)
```

Loops are discarded. Since `Q_T` is connected and `q` is surjective,
`E_(n,T)` is connected. Notice that the fixed set `B_T` need not be
`S_n`-invariant, so `(5)` is an edge image, not automatically an orbit graph
of an `S_n` action on `Q_T`. This is exactly the coordinate that the old
"canonical dual" wording hid.

## 2. Minimal basis dependence: `P_3` versus `K_3`

At `n=4`, `Z_4/S_4` has precisely three classes:

```text
0,       C_3 plus an isolated vertex,       C_4.                  (6)
```

Take first the star `T_star=K_(1,3)`. Its three fundamental cycles are all
triangles. From the empty graph a single move can therefore reach only the
triangle class. The symmetric difference of two distinct triangles in
`K_4` is a four-cycle, so the remaining non-loop transition is

```text
0 -- C_3 -- C_4.
```

Hence

```text
boxed: E_(4,T_star) ~= P_3.                                      (7)
```

For the path `T_path=0-1-2-3`, the three fundamental cycles have lengths
`3,3,4`. The four-cycle generator adds the missing edge from `0` to `C_4`,
while a triangle generator also connects the other two classes. Therefore

```text
boxed: E_(4,T_path) ~= K_3.                                      (8)
```

Equations `(7)` and `(8)` are the minimal witness: the vertex set is the
same, but the adjacency is spanning-tree dependent.

## 3. The exact star boundary

For every `n>=3`, all fundamental cycles of a star are triangles. Edge-count
parity

```text
[F] |-> |F| mod 2                                                (9)
```

is invariant under graph isomorphism and flips across every edge of
`E_(n,T_star)`. The graph is connected and has the nontrivial edge from the
empty graph to a triangle. Thus

```text
boxed: chi(E_(n,T_star))=2 for every n>=3.                       (10)
```

Stars are exactly the trees for which this all-odd-generator mechanism is
available. Indeed,

```text
|C_T(uv)|=dist_T(u,v)+1.                                        (11)
```

Let `A union B` be the tree bipartition. Every fundamental cycle is odd iff
every chord joins two vertices in the same bipartition part. Equivalently,
every cross-pair in `A x B` is already a tree edge, so `T=K_(|A|,|B|)`.
Because a tree has `n-1` edges,

```text
|A||B|=|A|+|B|-1
iff (|A|-1)(|B|-1)=0.                                           (12)
```

Thus one part is a singleton and `T` is a star. Conversely a star plainly
has only triangular fundamental cycles. Therefore

```text
boxed: every C in B_T has odd length iff T is a star.            (13)
```

This characterizes the standard parity proof, not the chromatic number of
every nonstar image; no such converse is claimed.

## 4. The canonical all-cycle envelope

Let `C_n` be the set of all labelled simple cycles of `K_n`, of every length
`3,...,n`, and define

```text
widehat(E)_n = simple orbit graph of Cay(Z_n,C_n) under S_n.      (14)
```

Unlike `B_T`, the generator set `C_n` is `S_n`-invariant, so `(14)` is a
canonical graph on the isomorphism classes of Eulerian graphs. Every
fundamental cycle is simple, hence

```text
E_(n,T) is a spanning subgraph of widehat(E)_n.                   (15)
```

The envelope is connected because it contains the connected image of any
fundamental-cycle cube.

### 4.1 Complete classification by tree diameter

For `3<=k<=n`, let `R_(n,k)` be the canonical edge relation obtained by
toggling one labelled `k`-cycle and then passing to isomorphism classes.
This is well-defined because `S_n` is transitive on the simple `k`-cycles.
For a tree `T`, put

```text
Lambda(T)={|C|:C in B_T}.                                        (15a)
```

Generators of the same length induce the same relation after quotienting,
so multiplicity and placement disappear and

```text
E_(n,T)=union_(k in Lambda(T)) R_(n,k).                           (15b)
```

If `D=diam(T)`, then

```text
Lambda(T)={3,4,...,D+1}.                                         (15c)
```

Indeed every chord cycle has length at most `D+1`. Conversely a diameter
path contains a pair at each distance `d=2,...,D`; that pair is a chord and
its fundamental cycle has length `d+1`. Thus

```text
boxed: E_(n,T)=E_n^(D):=union_(k=3)^(D+1)R_(n,k).                 (15d)
```

The inclusions

```text
E_n^(2) proper-subset E_n^(3) proper-subset ...
 proper-subset E_n^(n-1)                                        (15e)
```

are strict. The edge from the empty graph to the class of a single
`(D+2)`-cycle belongs to `E_n^(D+1)` and to no earlier layer: a move from
the empty graph is isomorphic to its generator, whose length is at most
`D+1` in `E_n^(D)`. Every diameter `D=2,...,n-1` occurs: take a path of
length `D` and attach any remaining leaves to its second vertex. Hence there
are exactly `n-2` spanning-tree images. In particular,
the full tree presentation is compressed losslessly to the one integer
`diam(T)` for this edge-image observable.

It also has two exact intrinsic descriptions.

### 4.2 Union over every spanning tree

Given a simple cycle `C`, delete one of its edges `e`. The remaining path is
a forest and therefore extends to a spanning tree `T` of `K_n` that does not
contain `e`. Then `C=C_T(e)`. Consequently

```text
C_n = union_(T spanning tree) B_T,
boxed: widehat(E)_n = union_T E_(n,T).                           (16)
```

The union in `(16)` is an edge union on the common canonical vertex set.

### 4.3 One path already realizes the envelope

Let `P_n` be the path `0-1-...-(n-1)`. For every `k in {3,...,n}`, its chord
`(0,k-1)` has a fundamental cycle of length `k`. The symmetric group is
transitive on labelled simple `k`-cycles. If an envelope edge is witnessed
by

```text
F -> F xor C,       |C|=k,                                      (17)
```

choose a permutation `pi` carrying `C` to the length-`k` path-basis cycle.
Then

```text
pi(F) -> pi(F) xor pi(C)                                        (18)
```

is a `P_n`-basis move with the same endpoint isomorphism classes. This gives
the reverse inclusion to `(15)` and proves

```text
boxed: E_(n,P_n)=widehat(E)_n for every n>=3.                    (19)
```

Thus the historical scripts based on the standard path did compute a
canonical object. What fails is the stronger, previously implicit assertion
that replacing that path by any spanning tree leaves the image unchanged.

For `n>=4`, the three classes in `(6)`, with `n-4` additional isolated
vertices, form a triangle in `widehat(E)_n`: `0` reaches `C_3` and `C_4`,
and a suitably overlapping triangle and four-cycle have triangular symmetric
difference. Hence the star/path separation persists for every `n>=4`:

```text
chi(E_(n,T_star))=2,       chi(widehat(E)_n)>=3.                 (20)
```

## 5. Exact audit and preserved consequences

The primary exact computation constructs `Z_n` from both star and path
bases. The independent computation instead scans all `2^binom(n,2)` edge
subsets and keeps the Eulerian ones. They agree on:

| `n` | labelled `|Z_n|` | iso classes | simple cycles | star edges | path/envelope edges |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 1 | 1 | 1 |
| 4 | 8 | 3 | 7 | 2 | 3 |
| 5 | 64 | 7 | 37 | 8 | 16 |
| 6 | 1024 | 16 | 197 | 31 | 90 |

The primary Pruefer scan tests all labelled trees through `n=7`; the numbers
of all-odd bases are `3,4,5,6,7`, exactly the labelled stars. The independent
subset scan reproduces all labelled-tree counts `3,16,125,1296` through
`n=6`, recognizes simple cycles without using cyclic orderings, and obtains
the same generator unions and path-envelope edge counts.

The correction boundary is therefore sharp:

- **survives:** the canonical vertex set of Eulerian graph isomorphism
  classes, its two-graph/switching interpretations where separately proved,
  and the graph object underlying the historical **path-basis** scripts;
- **retracted:** treating a general fundamental-cycle-basis image as
  independent of the spanning tree, or importing the path values to a star;
- **new canonical target:** `widehat(E)_n`, equivalently the historical path
  image by `(19)`, with basis-free definition `(14)` and union law `(16)`.

This theorem does not recertify historical numerical values beyond the exact
counts in the table. No tournament, LRC, asymptotic chromatic, perfection, or
two-graph transfer claim is made beyond this typed correction.
