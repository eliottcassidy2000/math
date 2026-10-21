# Order-ten closure and an adaptive fourth-layer cut

20 October 2026. **FINITE-EXACT:** all 86,016 literal order-ten closed blue
sources with endpoint degrees `(n-4,3)` admit three-black detours.
**PROVED:** an open-middle branch of the next endpoint layer `(n-5,4)`
admits the same kind of detour for every `n>=7`; a second guarded
end-pair branch, including closed middles, does so for every `n>=14`.
These are labelled, source-bound statements. Full fourth-layer coverage and
black-core connectivity remain **OPEN**.

## 1. Inheritance, carrier and scope

The immediate dependencies are the
[endpoint-three packet](../endpoint_three_repairs_20261020h/README.md),
its [movable-cut proof](../endpoint_three_repairs_20261020h/theory.md),
[open-middle and cyclic-cut proofs](../endpoint_three_repairs_20261020h/root.md),
and [end-pair proof](../endpoint_three_repairs_20261020h/lift.md).
That packet proves the third endpoint layer structurally for `n>=11`
and checks the retained finite cases through order nine. It leaves order ten
outside the structural cutoff. Its sharp cyclic-cut hostile at `n=4k-2`
is a failure of that proposed path, not a disconnection theorem.

The reset for the representation is the
[marked tournament trijection](../trijection_reset_20261020/README.md), with
the [session-17 scope synthesis](../metagraph_trijections_reset_session17_20261020.md).
A graph without its vertex order, or a tournament without its marked path,
does not retain the source of a repair.

Write the marked Hamiltonian path as

```
P=(0,1,...,n-1)=(a,b,x,3,...,y,c,z),
a=0, b=1, x=2, y=n-3, c=n-2, z=n-1.
```

For a Hamiltonian word `Q`, let `E(Q)` be its undirected consecutive-pair
set and let `K_Q(T)=T^op xor E(Q)`: reverse exactly the arcs outside Q.
The word is **blue** when its positional reversal is an anti-automorphism
of its actual source, and **black** otherwise. Every word below is checked
in the tournament reached by the preceding words.

Assume P is blue and closed, so `z->a`. For endpoint layer k,

```
d_T(a)=n-1-k,  d_T(z)=k,
N^-_T(a)={z} union I,  N^+_T(z)={a} union r_P(I).
```

Here I is a literal subset of positions `2,...,n-2` of size `k-1`.
The words “open middle” and “closed middle” mean `b->c` and `c->b`,
respectively; they do not name a strong-component decomposition.

The live concepts are the source path, exceptional endpoint ports, actual
complement masks, anti-reversal witnesses, and ordered graph charts.
The new operation adapts the cut to the last incoming port. It changes
neither the source nor the desired labelled endpoint.

## 2. Complete closure at order ten

**FINITE-EXACT theorem.** Every closed standard-blue tournament of order ten
with `d(a)=6,d(z)=3` has a three-black path from T to the literal
`K_P(T)`.

The enumeration is complete without an isomorphism, strong-target, or mixed
filter. Delete a,z. The remaining eight-position Hamiltonian tournament is
blue; its 21 free pair bits form 12 reflection orbits, giving 4,096 literal
middle presentations. Choose the two positions of I from `2,...,8`.
Thus the universe has exactly

```
4096 * binom(7,2) = 86,016
```

sources. Conversely every such middle presentation and I restores exactly
one source satisfying all endpoint and blue guards.

The inherited deterministic selected bank succeeds on 86,013 of them.
There are exactly three failed selected rows. The following supplied
certificates resolve all three; no search is performed by the production
consumer. A row integer encodes its outgoing neighbours by set bits.

| Source bit rows | I | First word | Second word | Third word |
|---|---|---|---|---|
| 438,724,808,339,676,842,149,812,594,73 | 3,6 | 0123458679 | 0123586794 | 0123567894 |
| 414,628,616,466,420,841,657,806,582,25 | 5,6 | 0123458679 | 0123987465 | 0123974685 |
| 126,44,120,752,994,704,898,775,559,7 | 7,8 | 0123467859 | 0123859467 | 0123894567 |

All digits here are individual vertex labels. During discovery the inherited
complete first/last-path search found these after 10,21,2 candidate pairs,
respectively. These small search counts are not completeness evidence;
the final 86,016-source census and literal replay are.

[receipts.json](receipts.json) retains each exceptional source, target,
three words, explicit anti-failure witnesses and all three ordered graph
charts. The whole deterministic certificate stream has its own SHA-256.
An independent directed-arc-set replay accompanies the bit-row computation
for every one of the 86,016 repairs.

Combined with the inherited all-order theorem, this closes every `n>=10`.
The earlier [finite producer](../endpoint_three_repairs_20261020h/probe_complete.py)
explicitly asserts that both full and retained universe counts are
96,640,7680 at orders 7,8,9; its
[receipt](../endpoint_three_repairs_20261020h/probe_complete.json) and
[scope statement](../endpoint_three_repairs_20261020h/probe.md) remove no
source at those orders. Consequently the combined closed blue third layer
is repaired at every `n>=7`. This uses the equality of those counts,
not an unqualified removal of the earlier filter. The separate order-six pure
`C3 -> C3` target hostile is unchanged. This is a layer theorem, not a
proof of unrestricted black-core connectivity.

## 3. The fourth-layer adaptive end-pair theorem

**Theorem (PROVED).** Suppose `n>=7`, P is a closed blue Hamiltonian path,
the middle is open `b->c`, and

```
I={h,y,c},  h in {2,...,n-4}.
```

In particular the marked endpoint degrees are `(n-5,4)`.
Put `e=n-4` and choose

```
k=e if h=e, otherwise k=y;
L=(3,...,k-1),  R=(k,...,y).
```

Empty L is allowed. R has one or two vertices. Then

```
Q1=(a,b,c,z,x,L,R),
Q2=(b,c,z,x,L,a,R),
Q3=(b,x,L,a,R,c,z)
```

are actual black Hamiltonian paths in succession and

```
K_Q3 K_Q2 K_Q1(T) = K_P(T).
```

### Directedness and exact endpoint

Q1 uses the supplied middle arc b->c, the original c->z, and z->x,
the last following from y in I by blue reversal. Every remaining join is
an original consecutive path edge.

Let p be the last vertex of `(x,L)` and r the first vertex of R. The cut
was chosen so that `p notin I` and `r in I`. Thus T has `a->p`
and `r->a`. Neither pair is an edge of Q1, so in `S=K_Q1(T)` the
new joins are `p->a->r`. All other edges of Q2 were retained by Q1.

At `U=K_Q2(S)`, the new joins `b->x` and `y->c` have been reversed
twice, whereas `c->z` was retained twice. The remaining Q3 edges were
retained from Q2. Hence Q3 is directed.

Deleting the common interiors of L and R from the path-edge ledger gives

```
E(Q1) xor E(Q2) xor E(Q3) = E(P).
```

Three complements therefore give the exact desired labelled target.
The proof also covers the empty-L boundary.

### The second and third words are black

At S, `a->b` remains true. If R has one vertex, Q2 reversal pairs this
arc with `y->c`, which is false at S. If R has two vertices, the paired
arc is `y->z`. In this case h=e and x is not in I, so y->z is true in T;
Q1 omits that pair, making it false at S. Both cases contradict blue Q2.

At U, `b->c` is true and `x->z` false: Q1 and Q2 both retained both
pairs. Q3 reversal pairs those two arcs. Therefore Q3 is black.

### The first word is black, including order nine

Set `t=n-5`. Q1 reversal pairs `a->t` with `c->y`. The latter is
false, while the former is true unless h=t. It remains to treat h=t.

If Q1 were blue, `g=r_Q1 r_P` would be a tournament automorphism.
Its orbit through a is

- all vertices for odd n;
- the even positions for `n=2 mod4`;
- `{0,2,6,10,...,n-6}` for `n=0 mod4`.

These follow by composing the two displayed reversals. A tournament
induced on a transitive automorphism orbit is regular.

For odd n, a's full outdegree n-5 would equal (n-1)/2, forcing n=9.
For `n=2 mod4`, the only in-neighbour of a in its even-position orbit
is c; regularity forces `n/2=3`, outside n>=7.
For `n=0 mod4`, a has no in-neighbour in its orbit, forcing that orbit
to have size one, again outside the range.

At the remaining n=9,h=4, the exact g-cycle is

```
(0,3,1,4,7,5,8,6,2).
```

Both 4 and 8 beat a, but they are `g^3(a)` and `g^(-3)(a)`.
Applying `g^3` to `g^(-3)(a)->a` yields `a->g^3(a)`,
contradicting the other incoming arc. Hence Q1 is black here as well.
This proves the all-order theorem without extrapolating a finite census.

## 4. A second structural branch, including closed middles

**Theorem (PROVED).** Suppose `n>=14`, P is closed and blue,
`I={h,y,c}`, `h notin {x,e}` with `e=n-4`, and `c->x->y`.
There is no condition on b versus c. Put `W=(3,...,e)`. The words

```
Q1=(a,x,W,y,c,z,b),
Q2=(z,W,y,x,b,a,c),
Q3=(b,z,W,y,x,a,c)
```

give three black moves to the same literal K_P endpoint. This extends
the inherited E2 xy ledger to a third incoming port; its additional
blackness argument below handles that port explicitly.

**Legality and XOR.** The exclusions h!=x,e give a->x and 3->z.
Q1 is therefore directed. In S=K_Q1(T), the new joins of Q2 are
z->3,y->x,x->b,b->a,a->c, all complements of the appropriate original
nonpath arcs. W and its join to y are retained. After Q2, the new Q3
joins b->z and x->a reverse the previously retained Q1 arcs, while its
other joins are retained. As in E2, the complete path-edge XOR is E(P).

**First word.** Its reversal pairs 3->z with x->c. The former is true
and the latter false in T, so Q1 is black.

**Second word.** At S its reversal sends a to 3 and vertex5 to x.
If h!=5, then S(a,5)=0 while S(x,3)=1, contradicting blue Q2.
If h=5, suppose Q2 blue. Its anti-equation on (a,4) gives T(1,3)=1:
both relevant pairs were complemented by Q1, and T(a,4)=1.
P reversal then gives T(e,c)=1. But Q2 reversal also sends vertex7
to e and z to c. Since n>=14, vertex7 is not an out-neighbour of z
in T: the latter are exactly a,b,x,n-6. Thus T(7,z)=1. Q1 omits
both pairs {7,z} and {e,c}, so the anti-equation on (7,z) instead
gives T(e,c)=0. This contradiction proves Q2 black.

**Third word.** Check it at the opposite endpoint B=K_P(T), where
blackness is equivalent to blackness before its own K move. P is still
blue there. If Q3 were blue, g=r_Q3 r_P would be an automorphism of B.
It fixes a,b,c,z, and advances by two cyclically on the ordered middle
set {2,...,y}. Its cycles have length n-4 for odd n and (n-4)/2 for
even n. In the stated range all have at least three vertices. However,
g fixes a and must preserve its middle out-neighbour set {h,y}, of size
two. Such a set cannot be a union of those cycles. Hence Q3 is black.

This proof is structural. The checker separately replays deterministic
guarded samples at orders 14 through 30,40,60. The stated range avoids
the special order13 incidence in the second-word argument; no failure of
reachability at that order is inferred. The omitted ports x and e are
actual legality boundaries of this ledger, not merely numerical cases.

## 5. A sharp guard hostile and the remaining fourth-layer branches

The extra port cannot be ignored. At n=7 take

```
rows=(6,36,8,19,103,77,15), I={3,4,5}.
```

The old unadapted end-pair word inserts a immediately before y:

```
Q1=(0,1,5,6,2,3,4),
old Q2=(1,5,6,2,3,0,4).
```

Its join `3->0` is false after Q1. The adaptive cut instead gives

```
Q2=(1,5,6,2,0,3,4),
Q3=(1,2,0,3,4,5,6),
```

which both pass. The first failed implication was “a beats the vertex
immediately before y”; the new incoming port can be precisely that vertex.
The sidecar is the literal last exceptional port, not the port count alone.

The present coverage board for k=4 is:

| Guard | Status |
|---|---|
| Open middle, y not in I, n>=15 | Inherited all-k cyclic-cut repair |
| Open middle, y,c both in I, n>=7 | New theorem above |
| Closed middle, c in I and y not in I | Inherited degree-free movable-cut repair |
| Closed middle, I={h,y,c}, c->x->y, h notin {x,n-4}, n>=14 | New second structural branch |
| Other fourth-layer guards | Not closed by this packet |

The inherited degree-free branch requires c->b,z->b,a->y,d(z)>=3,
which the displayed k=4 guard supplies. The cyclic-cut threshold is sharp:
at n=14 there are blue cyclic cuts, but this says nothing by itself about
the existence of other black detours. A productive next obligation is the
open branch `y in I,c notin I`, where two additional incoming ports
can obstruct the single-port repair at different joins.

## 6. Exact transport to ordered graphs

For a marked source `(T,Q=(q_0,...,q_(n-1)))`, put an edge
`{i,j-1}` in an ordered graph on `0,...,n-2` iff
`j>=i+2` and `q_j->q_i`. This is the inherited trijection,
retaining Q as its map from positions to literal source labels.

Under this same Q-chart, K_Q is graph complement. Blue is invariance under
the specified order reflection; black is failure of that invariance.
Thus each repaired move has an exact ordered-graph interpretation.
To concatenate moves one must decode the current labelled tournament and
re-encode in the next Q-chart. Complementing the old graph three times
without those chart transitions is a different operation.

Source: the labelled closed-blue marked tournament.
Target: its labelled K_P endpoint.
Map: the three authenticated complements, or their ordered graph charts.
Preserved predicate: actual directed Hamiltonian joins, blackness and the
literal endpoint.
Lost data under an unmarked quotient: which vertices carry the endpoint
ports, which graph cuts represent the next path, and which labelled target
is required.
Needed sidecar: each ordered chart and its literal Q.
Cheapest hostile: the inherited three isomorphic one-edge graphs on an
ordered three-vertex carrier lift to different tournament types.

This interface turns the new endpoint repairs into valid moves in the
ordered-graph representation. It does not infer global connectivity from
the existence of a local repair.

### Incoming compression of the same order-ten boundary

Incoming commit4ec382bd66, read during the final clean sync, supplies an
independent [168-cell affine certificate](../endpoint_three_order10_20261021/README.md)
for the same order-ten theorem. It splits the21 port sets by three actual
arc bits, proves all4,536 directed joins symbolically, and authenticates
501 inconsistent blue-step systems plus three singleton exceptions using
explicit GF(2) equation witnesses. Its three literal exceptional sources
and replacement words agree with section2 here.

This strengthens the representation of that boundary: our complete source
replay and its affine parameter-space proof reach the same unfiltered
conclusion by different verification routes. The incoming proof imports
none of the inherited Python implementations. Its affine contradiction
certificates are a useful next method for the remaining fourth-layer
guards, where another large census should not replace a structural
partition. This packet does not claim to have already produced that full
fourth-layer partition. Incoming results and their attribution are retained
separately rather than counted as newly derived here.

## 7. Reproduction and audit boundary

Run from the repository root:

```text
python -B -X utf8 05-knowledge/results/endpoint_four_20261020h/endpoint_probe.py
python -B -O -X utf8 05-knowledge/results/endpoint_four_20261020h/endpoint_probe.py
```

Both emit the same deterministic [receipts.json](receipts.json).
The [script](endpoint_probe.py) consumes the inherited selected bank and
independent directed-set replay; the three exceptional rows are supplied
literal certificates. Production performs no path search.

The complete fourth-branch finite controls are 16,96,1024,10240 sources
at orders 7,8,9,10, respectively. Orbit and XOR formulas are additionally
checked at every order 7 through 128. These support the written proof;
they do not replace it.

The second structural branch has 687 deterministic guarded replays at
orders 14 through 30,40,60. Those samples include the exceptional port h=5;
the written proof supplies the unbounded quantifier.

The consumers reject noninteger/bool source aliases, wrong source
identity, wrong endpoint layer, wrong literal target, malformed paths and
off-scope fourth-layer ports before accepting a receipt. All checks remain
active under optimization.

Final normal and optimized output agree on **98,340 active checks** and
the exact saved JSON, LF SHA-256
`9a0b806e9397eacdf61626bd1535ff36772f130edde46ede2ea8f64ab3688ebb`.
The full order-ten stream hash is
`d04f612042ba52edd7907516b84f45a09204e6bf50083f9374365bf4f60b019e`.

Independent peer audit: the family-algebra reviewer passed the complete
order-ten census/API and the all-order open-branch proof, independently
reproducing the original 97,652-check receipt. The parent independently
audited the added section4 closed-branch proof, including all actual
joins, the h=5 contradiction and the fixed-a middle-orbit argument.
The final 98,340-check normal/optimized pair additionally reproduces all
687 new guarded samples. No unresolved correction remains in these
stated scopes.
