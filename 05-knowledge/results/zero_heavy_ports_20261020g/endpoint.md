# Two exact ports in the closed blue endpoint layer (n-3,2)

**Status: PROVED** for the stated labelled three-black ledgers.
The three residual locations isolated by these ledgers are now **CLOSED**
by the incoming [complete endpoint-two theorem](../endpoint_two_repairs_20261020g/README.md),
integrated below in section9. The entire layer has three-black detours for
every n>=6, with a FINITE-EXACT head at orders6–8 and structural proof from9.
General black-core connectivity remains **OPEN**. Research: 2026-10-20.

## 1. Inheritance and the retained object

The closest mechanism is the three-K ledger and the intrinsic `(n-2,1)`
endpoint corollary in
[`kernels_games_detours_20261020f/detours.md`](../kernels_games_detours_20261020f/detours.md).
The incoming
[`strong_blue_games_20261020f/detours.md`](../strong_blue_games_20261020f/detours.md),
[`transport.md`](../strong_blue_games_20261020f/transport.md), and
[`AUDIT.md`](../strong_blue_games_20261020f/AUDIT.md) supply the two-block
repair, the exact relocation ledger, and the warning that a nonautomorphic
reordering need not be black. All nonstrong mixed blue edges already have
detours; the remaining global gate concerns strong mixed-to-mixed blue edges.
We do not reopen the proved nonstrong case.

The canonical hostile is a blue cyclic cut at even order: two reflections
compose to a shift by **two**, not necessarily a shift by one. The least-used
sidecar here is a Hamiltonian path with both endpoints prescribed in a
path-reversed induced tournament. Its existence is a real extra condition.

The object retained is a labelled tournament, its marked directed Hamiltonian
path, and its two exceptional endpoint arcs. The relation is the actual arc
orientation, with no ties or supplementary orientation choices. The map is
the marked complement K. We preserve the literal labelled target after three
black moves. Degree sequences, unweighted clone cores, and isomorphism classes
alone lose the locations of the exceptional arcs and the marked path.

Our concept board is: endpoint scores; two labelled ports; reflection
permutations; path-reversal endpoint availability; clone weights; and exact
edge-XOR ledgers. Every positive statement below keeps the required member
of this board explicitly.

## 2. Exact two-port normal form

Let `P=(p_0,...,p_(n-1))=(a,b,W,c,z)` be a blue directed Hamiltonian path,
closed by `z -> a`, with

`d_T(a)=n-3`, `d_T(z)=2`.

Write `V=(b,W,c)`, `r=n-2`, `W=(u,...,v)`, and let `rho` reverse P.
There is a unique middle vertex `x=p_i`, `2 <= i <= n-2`, such that

`N_T^-(a)={z,x}`, `N_T^+(z)={a,y}`, `y=rho(x)=p_(n-1-i)`.

Thus T is exactly the cyclic singleton extension
`a -> V -> z -> a`, with the two arcs `a-x` and `z-y` reversed.
The middle V is blue under its restricted path reversal. Conversely, every
such blue V and permitted x gives a source of this endpoint layer. The two
ports may coincide when x is the central vertex. They are still two distinct
arcs, and the formulas continue to apply.

**Proof.** The first score leaves exactly two inneighbours, one being z.
The first path arc excludes x=b. Blue reversal exchanges first-endpoint
inneighbours with last-endpoint outneighbours, determining y and all claimed
cross-arcs. The last path arc excludes y=c. The converse is immediate from
the paired flipped arcs. The closed path proves T is strong. This normal
form is an equivalence, not just a sufficient construction.

For a directed Hamiltonian path Q, use

`K_Q(T) = T^op XOR E(Q)`.

All edge sets here are sets of **unordered labelled edges**. Its score law is

`d_(K_Q(T)) = (n-1)1 - d_T + e_(first Q) - e_(last Q)`.

Blue means that positional reversal is an anti-automorphism. In particular,
blue Q has endpoint score sum n-1. The same path is blue in T iff it is blue
in K_Q(T). Three K moves end at `K_P(T)` exactly when their three path-edge
sets have XOR `E(P)`. Each ledger below also checks every path direction.

## 3. Two reusable blackness lemmas

### 3.1. The first full cyclic cut is black for n >= 7

The cut `Q=(V,z,a)` is directed. If it and P were both blue, their two
reflections would compose to the shift-by-two permutation of P, an
automorphism of T.

For odd n this permutation is transitive, so T would be regular. Its two
endpoint scores n-3 and 2 then force n=5. For even n=2h its parity orbit
containing a has h vertices. The induced tournament on that orbit is
vertex-transitive, hence regular; a has at most one inneighbour in that
orbit, because its only global inneighbours are z and x and z is in the
other parity orbit. Therefore `(h-1)/2 <= 1`, forcing n<=6. This proves the
claim for n>=7, with the small even boundary retained.

### 3.2. The common final word is black for n >= 8

Suppose V is closed (`c -> b`) and x!=c. Put `B=K_P(T)`. Then

`Q=(a,z,b,W,c)`

is directed in B. In B the outneighbours of a are exactly `{b,z,x}`.
The inneighbours of c include b and v: the original `c -> b` is
complemented, and the original path edge `v -> c` is retained.

If Q were blue, its reversal would take `N_B^+(a)` onto `N_B^-(c)`.
The images of b and z are respectively the second-last vertex w of W
and v. To include b, the image of x must be b. Consequently

`x=w=p_(n-4)`.

Compose Q's reversal with P's reversal, both anti-automorphisms of B.
In the p-labels the resulting permutation g is

`g(0)=n-3`, `g(1)=0`, `g(2)=n-1`,
`g(i)=i-2` for `3<=i<=n-2`, and `g(n-1)=n-2`.

For odd n it is one n-cycle. Regularity and `d_B(a)=3` would force n=7.
For even n=2h it has two h-cycles, with the orbit of a equal to

`{a} union {p_1,p_3,...,p_(n-3)}`.

A tournament automorphism cannot have an even cycle (its halfway power
would interchange a pair of vertices). If h is odd, this orbit induces a
regular tournament. Of a's three outneighbours b,z,x=w, only b is in that
orbit, so `(h-1)/2=1`, forcing n=6. Both alternatives contradict n>=8.

## 4. Open middle: one exceptional boundary port

Assume `b -> c` and n>=7. If `x!=v=p_(n-3)`, the three words

```
Q1 = (V,z,a)
Q2 = (z,V_without_b,b,a)
Q3 = (a,z,V_without_b,b)
```

give a three-black detour from T to `K_P(T)`.

**Directions and target.** Q1 is the full cyclic cut. In the first partner,
the join from z to the first vertex u of `V_without_b` exists exactly when
`y!=u`, equivalently x!=v; the other changed join uses `b -> c`.
After the second move the third word is directed. One can check that third
word in B itself because K retains its path edges: `a -> z`, `z -> u`,
all internal path edges, and `c -> b` hold there. The three path-edge sets
cancel to E(P), proving the literal labelled target.

**Blackness.** Lemma 3.1 handles Q1. After Q1 the scores of z,a are
n-3,1, so Q2 has endpoint sum n-2. For Q3 work in B. Its reversal exchanges
a with b and z with c. If x!=c, the three vertices a,c,z all point to b in
B, while `N_B^+(a)={b,z,x}`. Blue would force the image of x to be z,
hence x=c, a contradiction. If x=c, its image is z, but z does not point
to b in B: the port `z -> b` was complemented. Again Q3 is black.

### Endpoint-availability sidecar for the residual x=v

Let `C=K_V(T[V])`. If C has a Hamiltonian path H ending at b and starting
outside `{c,y}`, the same proof works with

```
Q1 = (V,z,a),   Q2 = (z,H,a),   Q3 = (a,z,H).
```

The new start condition is exactly the set of possible joins from z in both
the first partner and B. The internal H edges cancel between Q2 and Q3.
For Q3, write h for H's penultimate vertex. A hypothetical blue reversal
sends z to h. Since `{a,c,z}` all point to b in B, it would force h=c
and then x=c, impossible because x=v. No assumption that all such H exist
is being made.

This sidecar genuinely fails on an infinite family. Take C to be the
transitive r-tournament on `(b,u,...,c)` with only `b -> c` reversed.
It is blue and has the standard path; set `T[V]=K_V(C)` and choose x=v.
The **only** Hamiltonian path of C ending at b is `(u,...,c,b)`:
C-b is transitive and c is the only predecessor of b. Its start u=y is
forbidden. Thus the endpoint-availability sidecar is absent at every such
order. At n=7 the source rows are
`(46,100,8,82,103,76,5)` and the local C rows are `(14,28,24,16,1)`.
These are method obstructions, not counterexamples to black connectivity.

## 5. Closed middle: two standard ledgers

Now assume `c -> b`, n>=8, and set `epsilon=[b -> v]=[u -> c]`, the
equality following from blue reversal in V.

### 5.1. epsilon=0 and x!=c

Use

```
Q1 = (c,W,b,z,a)
Q2 = (z,c,W,b,a)
Q3 = (a,z,b,W,c).
```

The extra within-middle joins are `c -> u` and `v -> b`, supplied by
epsilon=0. The port-dependent join `b -> z` holds exactly when x!=c.
The remaining joins follow directly under K. Their edge XOR is E(P).

If Q1 were blue, `c -> b` and `c -> z` would force two different W
vertices (the second and first vertices of W) to point to a. But only x
among W can point to a. Thus Q1 is black. After Q1, z,a have scores
n-3,1, making Q2 black. Lemma 3.2 makes Q3 black.

### 5.2. epsilon=1 and x outside {c,u}

Use

```
Q1 = (c,b,z,a,W)
Q2 = (z,c,b,a,W)
Q3 = (a,z,b,W,c).
```

Here `b -> z` needs x!=c and `a -> u` needs x!=u. Together with
epsilon=1 these give all joins; the edge XOR is again E(P).

For completeness, the first word's blackness has a direct all-order
permutation proof, rather than an assumption about reordered paths. If Q1
were blue, `g=rho_Q1 rho_P` would be an automorphism of T. For n>=8,

```
g(0)=n-5, g(1)=n-3, g(2)=n-2, g(3)=1, g(4)=n-1, g(5)=0,
g(i)=i-4 (6<=i<=n-3), g(n-2)=n-4, g(n-1)=n-6.
```

Following the subtract-four chains and their displayed end jumps gives:

| Order | Cycle lengths of g |
|---|---|
| n odd | n |
| n=4h | 2h, 2h |
| n=4h+2 | h, h, h+1, h+1 |

For example, at n=4h one cycle is
`(0,4h-5,4h-9,...,3,1,4h-3,4h-7,...,5)`, of length 2h.
More explicitly, suppressing the ordinary subtract-four steps, the cycles
on the special vertices `{0,1,2,3,4,5,n-2,n-1}` are:

| n modulo 4 | Compressed cycles (insert all intervening subtract-four steps) |
|---|---|
| 0 | `(0,3,1,5)`; `(2,n-2,4,n-1)` |
| 1 | `(0,4,n-1,3,1,2,n-2,5)` |
| 2 | `(0,5)`; `(1,3)`; `(2,n-2)`; `(4,n-1)` |
| 3 | `(0,2,n-2,3,1,4,n-1,5)` |

This follows directly from the displayed end jumps and accounts for every
vertex in the four residue chains, proving the cycle-length table at all
orders. In the odd case transitivity contradicts the unequal
endpoint scores. In each even case there is an even cycle, impossible for
a tournament automorphism. Hence Q1 is black.

For Q2, let U1 be the first partner. Precisely
`N_U1^-(z)={b,y}`. The last vertex v points to b in U1. If x!=v it also
points to a. A blue Q2 would map `{b,y}` to `N_U1^+(v)`, but the image
of b under this reversal is a W vertex, different from a,b. This is
impossible. If x=v, then y=u. For n>=8 neither the image of b nor the
image of u under Q2's reversal is b, contradicting `v -> b`. Thus Q2
is black in both cases. Lemma 3.2 handles Q3.

## 6. Repairing the closed epsilon=1 port x=c

This boundary case has a different universal ledger, valid already for
n>=7. Here y=b, W=(u,Wtail), and

```
Q1 = (a,W,c,z,b)
Q2 = (z,Wtail,c,u,b,a)
Q3 = (b,z,Wtail,c,u,a).
```

One way to obtain it is the inherited two-block repair with singleton left
block (a), pivot b, and right block `R=(W,c,z)`. In the induced tournament
on R after reversing its displayed path edges, the word

`H_R=(u,c,reverse(Wtail),z)`

is an endpoint-preserving Hamiltonian path: `u -> c` comes from epsilon=1,
`c -> v` and the descending Wtail arcs are reversed path arcs, and the
last join to z is an untouched cross-arc. There is no assumption about
the remaining orientation inside W. Substituting H_R into the inherited
two-block identity gives the displayed words; their edge XOR is E(P).

All three are black. If Q1 were blue, `N_T^-(a)={z,c}` would imply
`N_T^+(b)={u,second(W)}`. But b also points to v, which is distinct
from those two vertices when |W|>=3. After Q1, z,a have scores n-3,3,
so Q2 has endpoint sum **n**, and is black. Finally evaluate Q3 in B.
Its reversal exchanges a,b and z,u. The arcs `a -> z` and `b -> u`
both hold in B, and the former would require `u -> b` under a blue
reversal. This contradiction makes Q3 black.

The same two-block construction can be retained as a sidecar when
epsilon=0: it needs an actual endpoint-preserving H_R and an actual black
first word. These requirements cannot simply be dropped. If V is the
full tiling F_r (only consecutive forward arcs, all other arcs backward),
then epsilon=0. In the path-reversed R tournament, u's sole outneighbour
is z. For |R|>=3 no Hamiltonian path can start at u and end at z.
This gives an infinite exact obstruction to that endpoint repair.

## 7. Exact clone mergers and the small boundary

In the **closed-middle case**, at x=c (so y=b), the following are equivalent:

1. a,b are a homogeneous directed pair a->b;
2. b dominates every vertex of W;
3. every vertex of W points to c;
4. c,z are a homogeneous directed pair c->z.

The pair criteria follow by comparing their relations to W,c,z or to
a,b,W. Conditions 2 and 3 are exchanged by the blue reversal. Equivalently,
`d_(T[V])(b)=r-2` and `d_(T[V])(c)=1`. In this case

`T = C3[TT2, T[W], TT2]`.

The two outer clone chains have exactly length two: a W vertex is
distinguished from either chain by the opposite cyclic block's vertices.
Internal clone chains of W are unrestricted. A canonical weighted-core
description must retain these outer weights 2,2 and the allocation of the
marked Hamiltonian path. Contracting them to an unweighted C3 does not
preserve the marked K operation.

For W=TT2, the six-vertex source has rows
`(14,12,56,48,35,3)`, alternating scores `(3,2,3,2,3,2)`, and its full
one-step cyclic cut is blue. This is an explicit boundary hostile to
extending Lemma 3.1 down to n=6. The source and its K partner are both
strong and mixed: explicit black paths are `(0,2,5,1,3,4)` and
`(0,1,4,5,2,3)` respectively.

There is also a closed-epsilon=0 six-vertex source with rows
`(14,4,40,50,39,3)`. Its source/partner black paths are
`(0,2,3,4,5,1)` and `(0,1,2,4,5,3)`.
For **each** six-vertex source, permuting the five fixed blocks
`(a),(b),(W=(2,3)),(c),(z)` gives no three-black ledger to the literal
labelled K target. This finite statement retains W contiguous and forward
and searches only these 120 words. It says nothing about longer detours,
split W words, or a final target reached only up to isomorphism.

## 8. Earlier template boundary and reproduction

For every n>=8, all sources in this endpoint layer are covered by the
displayed three-black ledgers **except potentially** these three positions:

| Middle condition | Remaining port | Retained extra obligation |
|---|---|---|
| open, b->c | x=v=p_(n-3) | an allowed-start HP of K_V(V) ending at b, or another ledger |
| closed, epsilon=0 | x=c=p_(n-2) | an endpoint-preserving R repair plus blackness, or another ledger |
| closed, epsilon=1 | x=u=p_2 | a different first join or a different ledger |

The first two displayed endpoint-repair routes have the infinite hostile
families proved above. These were open obligations for this particular bank
of ledgers; section9 records their subsequent closure by different
whole-block words. The infinite method obstructions remain valid and explain
why a different repair was needed. No bridge is asserted by their failure.

`endpoint_check.py` is self-contained, standard-library-only, and uses
explicit RuntimeError checks. It enumerates every blue fixed-path middle
of orders 4 through 8 and every permitted x: **32,140 labelled sources**,
without an isomorphism filter. It checks the normal form, scores, port
locations, every applicable three-black ledger and its literal target,
clone equivalences, and the available open sidecars at n=7,8. Additional
controls test both reflection cycle tables through n=100, both infinite
endpoint-repair obstruction families through r=24, and the two 120-word six-vertex hostile
universes. The finite checks are controls for the proofs, not replacements
for them.

Reproduce from the repository root:

```
python 05-knowledge/results/zero_heavy_ports_20261020g/endpoint_check.py
python -O 05-knowledge/results/zero_heavy_ports_20261020g/endpoint_check.py
```

The checked normal/optimized output is retained in `endpoint_check.out`.

## 9. Incoming completion of all three residuals

Commits 8a97a7e82e and f3e658d7cc arrived during final integration. They
prove the [complete endpoint-two theorem](../endpoint_two_repairs_20261020g/README.md).
The present investigators independently checked the three new mechanisms
below; they retain physical ports even on our infinite hostile families.
Write W=(u,D), preserving this note's convention that x is the exceptional
incoming port. The incoming notes use different letters for those roles.

* **Open middle, x=v:** the arbitrary-block repair in
  [probe.md, Theorem B](../endpoint_two_repairs_20261020g/probe.md) is
  `(a,b,c,z,u,D)`, `(b,c,a,z,u,D)`, `(b,u,D,c,a,z)`.
  It works at every n>=6 and needs no prescribed-endpoint path in K_V(V).
* **Closed middle, epsilon=0, x=c:** the block repair in
  [lift.md, Theorem A](../endpoint_two_repairs_20261020g/lift.md) is
  `(a,u,D,c,z,b)`, `(D,a,u,c,z,b)`, `(D,a,b,u,c,z)`.
  It works at n=7 and all n>=9; the remaining small orders have the
  separately replayed finite certificates. It applies to the full-tiling
  hostile, without requiring the absent endpoint-preserving R path.
* **Closed middle, epsilon=1, x=u:** write W=(u,U,v), with U nonempty.
  [theory.md, section5](../endpoint_two_repairs_20261020g/theory.md) gives
  `(a,b,u,U,z,v,c)`, `(a,u,U,z,v,b,c)`, `(a,u,U,v,b,c,z)`.
  Three explicit arc contradictions prove blackness for every n>=7.

Together with the preceding lemmas this covers all n>=9. The incoming
finite head exhausts all12,64,320 marked sources at n=6,7,8 and retains
the five exceptional literal certificates. Independent arc-set and bit-row
replays also cover all31,744 sources at n=9,10 without strongness filters.
The large-order sample checks supplement, rather than replace, the proofs.

Hence the earlier three-location frontier is superseded. The next direct
endpoint target is the closed score pair (n-4,3), which has **two** extra
incoming middle ports and their reflected outgoing partners. These marked
port defects are not the heavy-coordinate deficiency or black-orbit deficit.
The completed layer does not imply that every strong mixed blue edge has
such a presentation, so all-order black-core connectivity is still open.
