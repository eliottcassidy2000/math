# Session 125 — an actual rooted-forest / even-Fibonacci numeration bridge

**Status:** PROVED elementary bijections and determinant identities;
FINITE-EXACT implementation checks in the accompanying receipts. Classical
matrix-forest and Fibonacci results are credited below. No novelty claim.

## 1. Objects, conventions and the exact claim

Let `P_m` be the path with ordered vertices `1,...,m`, with `m>=1`. A rooted
spanning forest is a subset of its edges with exactly one chosen root in
each component. It is uniquely a sequence of blocks `(k,r)`, where `k>=1`
is the component's number of consecutive vertices, `1<=r<=k` is its root's
position inside that block, and the block lengths sum to `m`.

Write `M` for a monomer covering one cell, and `D` for a domino covering two
cells. A tiling is an ordered sequence of these tiles, without gaps.

Set `F_0=0,F_1=1`, and `A_i=F_(2i+2)` for `i>=0`. These weights begin
`1,3,8,21,55,...`. A length-`L` Chung–Graham word is an MSB-first string
`d_(L-1)...d_0` over `{0,1,2}` having no factor `2 1* 2`; leading zeros
are included. Its value is `sum_i d_i A_i`. The empty word has length zero.
The classical numeration theorem, proved independently in
`../session125_origin32_numeration_20261031/numeration/PROOFS.md`, says these words represent precisely
`0,...,A_L-1`, once each and in lexicographic order. The classical attribution
is also independently confirmed by [Burns, Theorem 2](https://arxiv.org/html/2502.18870v1).

**Bijection theorem.** For every `m>=1`, putting `L=m-1` gives explicit
mutually inverse maps

```
rooted spanning forests of P_m
             ↕
monomer/domino tilings of 2m−1 cells
             ↕
legal Chung–Graham words of fixed length m−1.
```

They preserve the refinement

```
number of forest edges = number of dominoes = number of NONZERO digits.
```

The arithmetic digit sum is NOT this statistic: a digit `2` contributes one
forest edge, not two. The graph here is a path with ordered vertices, not
an arbitrary square-sum graph. The last line is a code for its forests,
not a relabeling of the vertices by Fibonacci weights.

## 2. Rooted forest ↔ odd tiling, with inverse

Send a block `(k,r)` to

```
D^(r−1) M D^(k−r),
```

and put one extra `M` between consecutive blocks. A forest with `c`
components produces `m−c` dominoes and `2c−1` monomers, so its tiling has
length `2(m−c)+(2c−1)=2m−1`.

Conversely an odd-length tiling has an odd number of monomers, since
dominoes cover an even number of cells. Read its monomers left to right:
the first, third, fifth, ... are root markers; the second, fourth, sixth,
... are separators. Split at the separators. Each remaining segment has
one root-marker monomer, say `a` dominoes before it and `b` after it.
Recover `(k,r)=(a+b+1,a+1)`. The segment lengths give total order `m`.
Both procedures visibly undo each other, including `m=1` (`M`, one
singleton root). Every domino represents a forest edge.

## 3. Legal word ↔ tiling: fully oriented two-state transducer

Take a length-`L` word. Number the `2L+1` cells `0,...,2L` from left to
right. Read digits in **LSD-first** order `d_0,d_1,...,d_(L-1)`; at step
`i`, process cells

```
a = 2(L−i)−1,    b = a+1.
```

Thus the two-cell blocks are processed right to left. The final cell `0`
is not a block. There are two states:

* **Free:** neither cell `a,b` is already occupied.
* **Locked:** cell `b` is occupied by a domino started in the block to its
  right; cell `a` remains free.

Start Free. Here `D(x,y)` occupies adjacent cells `x,y`, and `M(x)`
occupies cell `x`.

| State before | Digit | Place now | State after |
|---|---:|---|---|
| Free | 0 | `M(a), M(b)` | Free |
| Free | 1 | `D(a,b)` | Free |
| Free | 2 | `M(b), D(a−1,a)` | Locked |
| Locked | 0 | `M(a)` | Free |
| Locked | 1 | `D(a−1,a)` | Locked |
| Locked | 2 | illegal | — |

At the final single cell, place `M(0)` if Free; if Locked, cell `0` is
already covered and no tile is placed. Cross-boundary dominoes are placed
in full when first created; the Locked state records their occupied cell.
Consequently there are no overlaps, gaps, or protrusions.

The legality condition is exactly the one forbidding `2 1* 2`. It is
unchanged on reversing a word, so LSD processing accepts the same words
as the usual MSB legality automaton. Every nonzero digit creates one
and only one domino; zeros and finalization create none.

**Reading direction remains part of the state.** Reversal preserves the
accepted language, not the terminal state. For example, `201` finishes Free
when read MSB-first but Locked when read LSD-first; the latter is the forest
transducer's state. Its tiling is `DMMMD`, whereas `102` gives `MDMDM`.
Do not identify this pending-domino state with an MSB append-digit or Beatty
eligibility state without reversing the reading convention.

**Inverse.** Scan any tiling using the same blocks from right to left.
If Free and `b` belongs to `D(a,b)`, output `1`. Otherwise `b` is a monomer:
if `a` is also a monomer, output `0`; if `a` is joined to `a−1`, output `2`
and enter Locked. If Locked, inspect the only unoccupied cell `a`: a
monomer outputs `0` and enters Free; a domino to `a−1` outputs `1` and
stays Locked. These cases exhaust all possibilities. Output digits were
read LSD-first, so reverse them to display the word MSB-first. The
single-cell finish is forced. This inverse both proves surjectivity and
rules out two words encoding the same tiling.

### All eight forests of P_3, not merely the number eight

Blocks are written `(length, root-position)`.

| Value | Word | Tiling | Blocks |
|---:|---|---|---|
| 0 | 00 | MMMMM | `(1,1),(1,1),(1,1)` |
| 1 | 01 | MMMD | `(1,1),(2,1)` |
| 2 | 02 | MMDM | `(1,1),(2,2)` |
| 3 | 10 | MDMM | `(2,1),(1,1)` |
| 4 | 11 | MDD | `(3,1)` |
| 5 | 12 | DDM | `(3,3)` |
| 6 | 20 | DMMM | `(2,2),(1,1)` |
| 7 | 21 | DMD | `(3,2)` |

`22` is forbidden. The canonical numeral `100` for integer eight belongs
to the next length and codes a forest of `P_4`, not a ninth forest of `P_3`.

## 4. Weighted refinement, recurrence and the origin boundary

Give every domino weight `t`, equivalently every forest edge weight `t`,
and every nonzero digit weight `t`. The transition matrix, with rows and
columns ordered `(Free,Locked)`, is

```
T(t) = [[1+t, t],
        [  1, t]].
```

Both final states finish with weight one. Therefore

```
f_m(t) = [1,0] T(t)^(m−1) [1,1]^T.
```

The classical matrix-forest theorem identifies the same polynomial as
`det(I+t L(P_m))` (weighted graph version of
[Chebotarev–Shamis, Theorem 5](https://arxiv.org/html/math/0602575v2)).
Since the trace of `T(t)` is `1+2t` and its determinant is `t²`, its
characteristic equation gives

```
f_1(t)=1,  f_2(t)=1+2t,
f_m(t)=(1+2t)f_(m−1)(t)−t²f_(m−2)(t),  m>=3.
```

At `t=1`, `T(1)=[[2,1],[1,1]]`, the square of the ordinary Fibonacci
transfer matrix. Thus `f_m(1)=F_(2m)`, also the classical path formula
[Chebotarev, Theorem 1](https://arxiv.org/html/math/0612792v4).

There is a direct coefficient formula: a tiling of `2m−1` cells with `k`
dominoes has `2m−1−2k` monomers and `2m−1−k` tiles in total. Choosing
which `k` tiles are dominoes gives

```
[t^k] f_m(t) = binomial(2m−1−k,k),  0<=k<=m−1.
```

Hence this same binomial coefficient counts legal length-`m−1` CG words
with exactly `k` nonzero digits. For `P_3`, the polynomial is
`1+4t+3t²`; for `P_4`, it is `1+6t+10t²+4t³`.

**Do not invent an empty object at formal index zero.** Extending the
recurrence to `m=2` requires formal `f_0(t)=0`. The empty graph really
has exactly one rooted forest and determinant one. Our bijection begins
at `m=1`, with the empty word mapping to the one-cell tiling and `P_1`.
It does not assign the empty word to the empty graph.

### What operations are, and are not, transported

Prepending one leading zero to a word prefixes two monomers to its tiling,
regardless of the old final state. In the forest, this prefixes one rooted
singleton component and shifts all old labels by one. This is a genuine
commuting operation of the three models. But leading zeros preserve the
integer value while changing the forest's order: fixed word length is
load-bearing and must not be quotiented away without a size sidecar.

At fixed length, numerical successor orders forests of the **same** path;
it is not path-vertex addition. The `00 -> 01` successor changes roots/edges
of a forest of `P_3` without adding a vertex. The structural bijection
therefore neither establishes nor claims that ordinary arithmetic successor
is the square-sum graph growth operation.

### A genuine operation-preserving bridge: right-singleton extension

Appending `0` at the **right/LSD end** of a displayed word is different
from prepending a leading zero. The new initial LSD step is Free-to-Free
and places two new rightmost monomers; all remaining cells are exactly the
old tiling. Consequently

```
word w -> w0    <->    tiling T -> TMM    <->    forest F -> F followed by (1,1).
```

If `n=V(w)` and `P(n)=sum_i d_i A_(i+1)`, the companion
[Beatty theorem](../session125_origin32_numeration_20261031/numeration/BEATTY_BRIDGE.md)
gives `P(n)=ceil(beta*n)`, where `beta=(3+sqrt(5))/2` and `P(0)=0`.
Thus right-singleton extension transports this exact numeric shift, whereas
left-singleton extension transports leading-zero padding and preserves value.
These are two different typed operations, not two descriptions of successor.

For a **nonempty canonical positive word**, no leading zero is equivalent
to the first forest component not being a singleton: the proved left-padding
bijection gives both directions. Such a forest with exactly one edge must
have first component `(2,1)` or `(2,2)`, followed only by singleton components.
Its words are exactly `10^k` and `20^k`, respectively, for `k>=0`. Their
ranks are `A_k` and `2A_k`. With the companion invariant
`Q(n)=P(n)^2-3nP(n)+n^2`, they have Q=1 and Q=4 because
`A_(k+1)^2-3A_k A_(k+1)+A_k^2=1`; this identity follows from the recurrence
and the initial pair `(1,3)`. This classifies the **one-edge forest subclass**,
not all integer solutions of a norm equation.

In particular, `1 -> 10 -> 100` has ranks `1 -> 3 -> 8` and appends rooted
singletons to the same rooted first edge. The path orders here are 2,3,4,
whose **total** forest counts are 3,8,21. The rank of an individual encoded
forest is not the total count of forests of that path. No operation here
is square-sum Hamiltonian gluing or square-sum vertex arrival.

## 5. Scalar counts forget the marked port: a minimal clean hostile

For any finite simple graph `G`, let `R(G)=det(I+L(G))`, and put
`A=I+L(G)`. This matrix is positive definite because
`x^T A x = sum_v x_v² + sum_(uv in E)(x_u−x_v)²`.

Adjoin one new vertex `z`, with neighbor set `S` consisting of distinct
old vertices. Let `s=|S|`, `b=1_S`, and `M=A+diag(b)`. Then

```
I+L(G+z) = [[ M,  −b],
            [−b^T, 1+s]].
```

The Schur determinant identity gives the exact response

```
R(G+z) = det(M) (1+s−b^T M^−1 b)
       = (1+s)det(M) − b^T adj(M)b.                 (1)
```

The second expression is an integer formula requiring no divisions.
It includes `S=empty`: adding an isolated vertex leaves `R` unchanged.
The formula counts all rooted forests; it is not a Hamiltonicity test.

For a single neighbor `v`, the determinant lemma and rank-one inverse
formula reduce (1) to

```
R(G+leaf_at_v) = R(G) (2+(A^−1)_(vv))
              = 2R(G) + adj(A)_(vv).              (2)
```

On `P_3`, in its natural vertex order,

```
A = [[2,−1, 0],             A^−1 = (1/8) [[5,2,1],
     [−1,3,−1],                           [2,4,2],
     [ 0,−1,2]],                          [1,2,5]].
```

Its rooted-forest count is eight. Attaching a leaf to an endpoint gives
`2*8+5=21`, the count of `P_4`. Attaching it to the middle gives
`2*8+4=20`, the count of `K_(1,3)`. The same old graph, old scalar count,
and one-edge birth have different responses because the marked port is
different. No function of the old scalar count and arrival degree alone
can govern all graph growth. Retain a port-dependent boundary response.

## 6. The exact boundary sidecar for one-vertex queries

Let `U` have columns the coordinate vectors of the vertices in `S`, and
let `C=U^T A^−1 U`, the Green matrix restricted to these named ports.
Since `M=A+UU^T`, the determinant lemma and Woodbury identity give

```
det(M)=R(G)det(I+C),
U^T M^−1 U = C(I+C)^−1.
```

Substituting into (1) and using
`I−C(I+C)^−1=(I+C)^−1` yields

```
R(G+z) = R(G)det(I+C) [1 + 1^T(I+C)^−1 1]
       = R(G) [det(I+C)+1^T adj(I+C)1].            (3)
```

Thus the scalar count together with the **named-port Green matrix** is
sufficient for this query; the scalar alone is insufficient. To answer
future queries at arbitrary old vertices, retain the full inverse, or a
boundary that already includes every future attachment port. Retaining
only the currently queried submatrix does not license later queries at
forgotten vertices. General inverse/Schur updates can propagate a fixed
retained boundary as vertices are added.

This is a different sidecar from the full connectivity/endpoint-pairing
state needed for Hamiltonian continuation. A determinant and Green matrix
solve the rooted-forest counting query, not the Hamiltonian query.

## 7. Actual square-sum graphs and what the new origin hides

For `G_n`, vertices remain the ordinary integers `1,...,n`, with an edge
between distinct `u,v` exactly when `u+v` is a square. Exact Bareiss
Laplacian determinants, cross-checked by the Schur arrival calculation:

| Graph | Rooted spanning forests |
|---|---:|
| G_12 | 9,240 |
| G_30 | 4,337,712,683,055,095 |
| G_31 | 21,517,231,556,064,823 |
| G_32 | 98,839,029,392,469,636 |
| G_33 | 795,810,222,452,889,290 |

The arrival ports are `{5,18}` for vertex 31, `{4,17}` for vertex 32,
and `{3,16,31}` for vertex 33. The `G_12` count is also independently
`R(P_3)R(P_4)R(P_5)=8*21*55` from its three path components.

These large counts do not say whether a spanning cycle exists. In
particular, the independently inherited `G_31/G_32` Hamiltonian threshold
must use its named-port and global-connectivity proof, not these counts.
Declaring `G_32` coordinate zero (or one) transports graph order and
successor; it does not erase the prehistory of the boundary response.
The numeral system is a coordinate for that growth index. The forest
bijection above is a separate, real combinatorial meaning of the same
classical weights, with an explicit limit on what it transports.

## 8. Evidence, independent controls and provenance

Reproduce with Python standard library only:

```sh
python test_forest_bridge.py verification.normal.json
python -O test_forest_bridge.py verification.optimized.json
cmp verification.normal.json verification.optimized.json
```

The test uses explicit checked exceptions, not optimization-removable
`assert`. Its universes are all ternary words of lengths 0 through 11;
all odd tilings of lengths 1 through 23; all rooted forests of paths
`P_1` through `P_12`; every labelled simple graph of orders 0 through 5;
and every neighbor set for every old graph of orders 0 through 4.
No graph-class, connectedness, or Hamiltonicity filter is used.

Independent paths: words are checked against a regular-expression legality
test; tilings are recursively generated by first-tile branching; forests
are generated by cut masks and Cartesian products of root choices.
Every graph's forest count is also computed by edge-subset enumeration,
union-find cycle rejection, and the product of its component sizes.
Determinants are cross-checked by the full signed permutation expansion.
The one-apex Schur formula is checked against those independent forest
counts, and the port-restricted Green formula is checked separately.
Weighted path determinants are evaluated at `t=0,1,2,3`, in addition to
exact coefficient histograms. The leading-zero operation, rank values,
invalid words/tilings, `m=1`, and formal-versus-empty boundary are explicit
controls. Receipts record the exact totals and source hashes.

Repository inheritance (read through `git show HEAD:<path>` from sparse
source HEAD `59cd2ec44b9bf5f2598b1271a927eae383fdb5df`):
`05-knowledge/results/session32_square_monster_20261023/CAT_MAP_PATH_FORESTS.md`.
The inherited recurrence and warning were retained, and upgraded here to
explicit object-level maps with a weighted statistic.

Primary sources directly inspected on 2026-10-31:

1. Chebotarev–Shamis, *Matrix-Forest Theorems*,
   [arXiv:math/0602575v2](https://arxiv.org/html/math/0602575v2),
   Theorems 5–7 (weighted rooted forests and inverse response).
2. Chebotarev, *Spanning Forests and the Golden Ratio*,
   [arXiv:math/0612792v4](https://arxiv.org/html/math/0612792v4),
   Theorem 1 and its direct combinatorial proof (path Fibonacci counts).
3. Burns, *Chung-Graham and Zeckendorf representations*,
   [arXiv:2502.18870v1](https://arxiv.org/html/2502.18870v1),
   Theorem 2 and the stated LSD convention (classical numeration attribution).

The odd-tiling Fibonacci recurrence is elementary and classical; its
specific implementation here is fully proved rather than claimed novel.
No source-worktree edits or remote publication were performed by this lane.
