# Short-cycle determination of collinear-triple moments

**Status: PROVED, by the covering and partition-inversion argument below;
numeric variance coefficients are separately FINITE-EXACT.**

**Independent audit:** root and wildcard_portfolio checked the covering lift,
same-shore partition inversion, non-induced copy normalization, geometric
multiplicities, and the weighted cycle bound. Both reviews passed. The
finite compiler and its optimized run have identical output; its finite
checks corroborate the formulas but are not the proof of the all-size result.

Let `G` be a finite simple bipartite 2-regular graph with `n` vertices in
each specified shore. Let `c_L(G)` count its connected components that are
cycles of length `L`. Thus `L` is even, `L>=4`, and
`sum_L (L/2)c_L=n`. Independently assign uniform random bijections from
the two shores to `{0,...,n-1}`. An edge becomes the grid point with those
two labels. Let `X_3` be the number of unordered three-element edge sets
whose grid points are collinear on a nonaxis line.

For every `n>=6`, there are rational numbers `A_n,B_n,D_n`, depending only
on `n`, such that every such skeleton satisfies

```text
Var(X_3)=A_n+B_n*c_4(G)+D_n*c_6(G).                 (1)
```

The mean depends only on `n`. No sign assertion about `B_n,D_n`, geometric
density estimate, or existence of a labeling with `X_3=0` follows merely
from (1). All cycles of length at least eight are invisible to the first
two moments once `n,c_4,c_6` are fixed.

The inherited exact joint-incidence compiler is
[overnight_20260906_no3line_pairprofiles.py](../../04-computation/overnight_20260906_no3line_pairprofiles.py).
It enumerates the grid geometry and every skeleton cycle type for
`n=4,...,8`. The argument here supplies its universal structural
justification; it does not extrapolate those finite checks. No external
literature theorem is invoked.

## 1. Homomorphism and injection normalization

All graphs and maps preserve the two specified shores. For a finite
simple bipartite graph `H`, define:

- `hom(H,G)`: the number of vertex maps preserving shores and edges;
- `inj(H,G)`: the number of such maps injective on each shore;
- `C_H(G)=inj(H,G)/aut(H)`, where `aut(H)` counts shore-preserving
  automorphisms.

For graphs with no isolated vertices, `C_H(G)` counts edge subsets of `G`
whose graph, on the incident vertices, is isomorphic to `H`. These are
**non-induced** copies: other edges of `G` between the selected vertices
are irrelevant. If `H` has `r` left and `s` right vertices, then

```text
C_H(K_n,n)=(n)_r*(n)_s/aut(H).                     (2)
```

Here `(n)_r` is the falling factorial. A graph with a vertex of degree
greater than two has `inj(H,G)=0` for the present target.

For a partition `pi` of a finite set, write

```text
mu(pi)=product_(blocks B of pi) (-1)^(|B|-1)(|B|-1)!.
```

Partition inversion on the two shores separately gives the exact formula

```text
inj(H,G)=sum_(pi_left,pi_right)
    mu(pi_left)*mu(pi_right)*hom(H/(pi_left,pi_right),G).       (3)
```

The quotient identifies only vertices in the same shore. It has no loops;
parallel edges can be collapsed without changing its homomorphism count.
Its number of distinct edges is therefore at most `|E(H)|`. Formula (3)
follows by sorting every homomorphism according to its two fiber
partitions and applying Mobius inversion to the product partition lattice.

## 2. Rooted cycle-cover lemma

Let `F` be connected and bipartite, with `m>=1` edges, and fix a vertex
in its left shore. For even `L>=4`, let `h_L(F)` count homomorphisms from
`F` to the cycle `C_L` sending this root to a specified left vertex.
Let `h_infinity(F)` count homomorphisms from `F` to the integer path
`...,-1,0,1,...`, with parity as the bipartition, sending the root to zero.
The latter count is finite because `F` is finite and connected.

**Lemma.** If `L>m`, then `h_L(F)=h_infinity(F)`.

**Proof.** Projection modulo `L` sends every integer-path homomorphism
to a rooted cycle homomorphism. It is injective: the lift of a neighboring
vertex is uniquely determined by its residue and the known lift, since
the two residues obtained by adding `1` and `-1` differ for `L>=4`.

Conversely, choose a spanning tree of `F` and lift a cycle homomorphism
uniquely along that tree from root zero. For every edge outside the tree,
the associated fundamental cycle is simple and has length at most `m`.
The sum of its oriented increments, each `+1` or `-1`, is a multiple of
`L` and has absolute value at most `m<L`. It is therefore zero. Thus
every remaining edge also lifts, giving the required homomorphism to
the integer path. This proves the bijection. The same proof shows that
trees have `h_L=h_infinity` for every `L>=4`. QED.

Cycle vertex transitivity now gives

```text
hom(F,G)
 =sum_(even L>=4) (L/2)c_L(G)*h_L(F)
 =n*h_infinity(F)
   +sum_(4<=L<=m, L even) (L/2)c_L(G)
       *(h_L(F)-h_infinity(F)).                    (4)
```

If a connected component is a single isolated vertex, its homomorphism
count is simply `n`; it adds no cycle variable.

## 3. The general short-edge subgraph theorem

**Theorem.** For each finite bipartite graph `H` with `m` edges,
`inj(H,G)` is an integer polynomial in `n` and the variables
`c_4,c_6,...,c_(2 floor(m/2))`. Every monomial satisfies the cycle-weight
bound

```text
sum_(even L>=4) L * exponent_of(c_L) <= m.          (5)
```

The exponent of `n` is not counted in this cycle weight. For graphs with
no isolated vertices, the analogous polynomial for `C_H(G)` has rational
coefficients and the same bound.

**Proof.** Homomorphism counts multiply over connected components. In
(4), any selected factor containing `c_L` comes from a component with at
least `L` edges. Hence in the product for an arbitrary quotient in (3),
the sum of the lengths of the selected cycle factors cannot exceed the
sum of the component edge counts. That sum is at most `m`. Formula (3)
is an integer linear combination of these products. Dividing by the fixed
integer `aut(H)` gives the copy-count assertion. QED.

In particular, for every `H` with at most six edges,

```text
C_H(G)=a_H(n)+b_H(n)c_4(G)+d_H(n)c_6(G).            (6)
```

There can be no quadratic short-cycle term, because its minimum possible
weight would be `4+4=8`. For at most three edges there is no cycle term
at all. Neither statement assumes `H` is connected or a matching.

## 4. Exact passage to the geometric event

A nonaxis collinear triple of distinct grid points has three distinct row
labels and three distinct column labels, so its graph is a three-edge
matching. For two such triples, their union is a simple bipartite graph
with at most six edges and degree at most two. The triples can share zero,
one, or two edges when they are distinct; diagonal pairs can also be
included to compute `E[X_3^2]` directly.

The independent shore permutations act transitively on all grid copies
of a fixed shore-preserving graph type `H`. Consequently, for each fixed
grid edge set `U` of type `H`,

```text
P(U is contained in the randomly labeled skeleton)
    =C_H(G)/C_H(K_n,n).                            (7)
```

One can also obtain (7) by averaging the indicator over all copies of
`H`: every labeling contains exactly `C_H(G)` of them. Extra skeleton
edges on the same row/column vertices are allowed, which is why induced
copy counts would be the wrong observable.

Let `W_n(H)` count ordered pairs of **distinct** nonaxis collinear grid
triples whose union has type `H`. Multiplicity is essential: different
pairs can produce the same union. Then

```text
E[(X_3)_2] =sum_H W_n(H)*C_H(G)/C_H(K_n,n),         (8)
```

where `(X_3)_2=X_3(X_3-1)`. There are finitely many relevant `H`, and for
`n>=6` each prototype can be interpreted in the same complete graph.
Substituting (6) into (8) proves affine dependence on `c_4,c_6`.
The analogous one-triple expression involves just a three-edge matching,
whose copy count depends only on `n`. Thus `E[X_3]` depends only on `n`,
and `Var(X_3)=E[(X_3)_2]+E[X_3]-E[X_3]^2` proves (1).

This is the exact source-to-target map: retain the union's bipartite
incidence type and the multiplicity of geometric event pairs. Replacing
that type by only its edge count or overlap cardinality destroys
information needed by (7). Replacing it by the three statistics
`n,c_4,c_6` is justified only after the short-edge theorem and its weighted
geometry sum. No independence between the two line events is assumed.

For the degree-at-most-two types used by the compiler, the sorted list
of component tuples `(edge_count,left_vertices,right_vertices,is_cycle)`
determines the shore-preserving type. Components are paths or even cycles.
Their shore-preserving automorphism counts are `L` for an `L`-cycle, `2`
for an even-edge path, and `1` for an odd-edge path; repeated isomorphic
components add the usual multiplicity factorial. These are precisely
the denominator factors in (2).

## 5. Boundary controls and higher moments

Two cheap hostile controls locate the edge-budget boundary. The graphs
`C_16` and `C_8` disjoint-union `C_8` have the same `n=8,c_4=c_6=0`,
but have respectively zero and two copies of the eight-edge graph `C_8`.
Also, the shore-preserving injection count of two disjoint four-cycles is
`16*c_4*(c_4-1)`: choose two distinct target components and one of four
shore-preserving identifications for each. Thus both a new cycle length
and a quadratic cycle term are actually needed at edge budget eight.
For the short-cycle normalization itself, rooted closed-walk counting gives
`hom(C_4,G)=6n+4c_4` and `hom(C_6,G)=20n+24c_4+6c_6`, while their
injection counts are respectively `4c_4` and `6c_6`. These identities
also check the factors `L/2` in (4) and the distinction between homomorphisms
and injections.

For `n>=6`, all three affine coefficients can be recovered without an
all-cycle-type census. Evaluate the statistic on:

```text
G_0=C_(2n),
G_4=C_4 disjoint-union C_(2n-4),
G_6=C_6 disjoint-union C_(2n-6).
```

The `G_4-G_0` contrast is the `c_4` coefficient. The `G_6-G_0` contrast
is the `c_6` coefficient for `n>=7`; divide it by two at `n=6`, where
`G_6` consists of two six-cycles. This applies to each copy profile and
to their variance-weighted sum separately. The compiler's `n=4,...,8`
checks and numerical coefficients remain finite calculations; the
all-`n` conclusion follows from Sections 1--4.

More generally, fix `k>=1` and, for a uniform full prototype list, assume
`n>=3k`. An ordered tuple of `k` distinct geometric triples has a union
with at most `3k` edges. Define its multiplicity weight as in (8), now
for the whole ordered tuple. Types of degree greater than two have zero
copies in `G`; all other types are covered by the same theorem. Hence

```text
E[(X_3)_k] is a polynomial in c_L for even 4<=L<=3k,
with every monomial satisfying sum L*exponent(c_L)<=3k.      (9)
```

Its coefficients may depend on `n` through the exact grid geometry;
no polynomial or asymptotic formula for that dependence is asserted.
Ordinary moments have the same cycle-length and weight bounds by their
finite expansion in factorial moments. For smaller `n`, the same proof
works after omitting prototypes with too many vertices in a shore.

The next factorial moment can therefore see `c_8` and `c_4^2`, as well
as `c_4,c_6`; the theorem allows these dependencies, without asserting
that their geometry-weighted coefficients are nonzero. A universal sign
or a zero-defect labeling requires further geometric information.

## 6. Exact small-grid coefficients

The [compiler transcript](overnight_20260906_no3line_pairprofiles.out)
gives the following coefficients in (1), checking every possible cycle
profile in each displayed size:

| n | A_n | B_n | D_n |
|---|---:|---:|---:|
| 6 | 247199/22500 | 16091/8100 | -8/225 |
| 7 | 4823011/264600 | 26219/22050 | 1117/264600 |
| 8 | 3635551/144060 | 33133/35280 | 23/31360 |

The sign change of `D_n` rules out a size-independent sign inference from
six-cycle incidence alone. The 38 union types observed for `n=6,7,8`
are a finite observation, not a claim that every possible geometric union
type has already appeared. These distinctions matter for the next test:
derive the geometry-weighted coefficients, rather than infer a tail law
from the existence of the short-cycle expansion.
