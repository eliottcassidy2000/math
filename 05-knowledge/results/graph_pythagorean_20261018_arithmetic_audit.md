# Graph/Pythagorean arithmetic audit, 2026-10-18

**Status: PROVED elementary statements below; FINITE-EXACT implementation controls.**
These are local deductions and a synthesis of inherited canon, with no novelty
claim. No proof of Kuratowski's theorem or a Tutte theorem is claimed.
No new external result is imported as a CITED input; literature provenance
remains with the existing core-paper references attached to inherited canon.

## Inheritance and concept board

- Closest proved mechanism: [THM-3333, Gaussian-square Farey/Pythagorean light
  cone](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
  especially sections 7--8. It already realizes triangular polarization by
  complete graphs and missing complete bipartite cuts. This note does not claim
  that construction as new.
- [THM-3334, Berggren parabolic spine](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
  sections 1--2, proves that the hypotenuse/even-leg unit-gap spine converges to
  the degenerate null ray `(0,1,1)`.
- [THM-3335, square-triangular Pell selector](../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
  sections 1--2, separates the consecutive-parameter condition from the extra
  square-even-leg condition. Those two Pell constructions must not be confused
  with the consecutive-*legs* branch below.
- [THM-3357, Berggren three-branch Walsh collapse](../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md)
  already retains the middle Pell branch and warns that Gaussian squaring does
  not commute with summing spinors. We use a single spinor at each step.
- Canonical hostile: `(5,12,13)` has hypotenuse/even-leg difference one but is
  further from isosceles than `(3,4,5)`. Corrected near miss: the preceding
  theorems' explicit separation of raw spinors, primitive triangles, and
  aggregate moments. Least-used sidecar here: two-clique relative graph
  homology in THM-3333 section 8; it preserves cut rank, not planarity.

The working board is (i) Euclid's rank-one cone; (ii) unit-gap parabolic versus
consecutive-leg Pell motion; (iii) clique/cut polarization; (iv) edge-disjointness
graphs; (v) Seidel switching. The anchor is a faithful graph/Pythagorean map,
the niche is an exact bridge among the three named graphs, and the wildcard is
whether Seidel spectrum itself recovers `3,4,5`.

## 1. A faithful marked-clique correspondence

For a simple graph write `A(G)` for its adjacency matrix. Take a complete graph
`K_c` with a distinguished partition into nonempty shores of sizes `x >= y`.
Compress `A(K_c)+I=J` to the orthonormal shore-indicator vectors. Its matrix is

```
B = [[x, sqrt(x*y)], [sqrt(x*y), y]],
det(B)=0,                 tr(B)=x+y=c.
```

Therefore

```
(a,b,c)=(x-y, 2*sqrt(x*y), x+y),         a^2+b^2=c^2.       (1)
```

The cut graph is `K_(x,y)` and its nonzero adjacency eigenvalues are
`+-sqrt(x*y)`: this follows by the same two-dimensional compression, with zero
diagonal. Thus the even coordinate of (1) is exactly twice its spectral radius.
All other adjacency eigenvalues of the cut are zero.

The output is integral iff `x*y` is a square. To justify the parity subtlety:
if `2*sqrt(x*y)` is an integer, its square is divisible by four, so that integer
is even and `sqrt(x*y)` is an integer. Equivalently, write

```
x=d*m^2,             y=d*n^2
```

with the same positive squarefree `d`. The output is then
`d*(m^2-n^2,2mn,m^2+n^2)`.

The output is a positive primitive triple **with this coordinate `a` odd** iff

```
x=m^2 > y=n^2,       gcd(m,n)=1,       m-n odd.             (2)
```

Indeed primitivity first forces `d=1`, then `gcd(m,n)=1`; odd `a` forces
opposite parity. Conversely these conditions give the usual primitive triple:
an odd prime common to two coordinates divides both `m,n`, while parity rules
out a common factor two. Thus primitive triples with the odd leg distinguished
are in bijection with isomorphism classes of marked clique partitions (2).
This is the ordinary Euclid correspondence packaged as a graph construction.

At the root `(m,n)=(2,1)`, the complete graph is `K_5`, the partition is `4+1`,
and the cut is `K_(4,1)`. A partition yielding `K_(3,3)` instead has `x=y=3`,
and gives the degenerate triple `(0,6,6)`. It cannot represent the root under
this map. A general graph is not determined by these two shore sizes.

There is also an edge-count construction consistent with inherited THM-3333.
Put `F(G)=|V(G)|+2|E(G)|`. For nonempty simple graphs of orders `m,n`, set

```
a=F(G)-F(H),       b=2mn,       c=F(G)+F(H).
c^2-a^2-b^2=4*(F(G)F(H)-m^2*n^2) <= 0.                 (3)
```

Equality holds iff both graphs are complete, because `0<F(G)<=m^2` and
`0<F(H)<=n^2`. When they are complete, `b=2|E(K_(m,n))|` is the doubled cut
term and `a=m^2-n^2` is the difference of their loop-completed relation counts.
This is a real algebraic role for complete and bipartite graphs **inside one
construction**, not an assignment of whole graph classes to separate legs.
Deleting one edge of `G=K_m` gives exact defect `-8n^2` in (3).

## 2. Two meanings of difference one, with different limits

For the Euclid lift `a=m^2-n^2, b=2mn, c=m^2+n^2`,

```
c-b=(m-n)^2,             c-a=2n^2.                         (4)
```

Hence `c-b=1` iff `m=n+1`. This branch is

```
(3,4,5), (5,12,13), (7,24,25), (9,40,41), ...
```

and satisfies `a/c -> 0, b/c -> 1`. It approaches a degenerate triangle.
This is precisely the inherited parabolic branch, not a new finding.

Isosceles means `a/b -> 1`, which is equivalent for `m>n>0` to

```
m/n -> 1+sqrt(2).                                         (5)
```

The central operation `(m,n)->(2m+n,m)` preserves coprimality and opposite
parity and flips the sign of

```
Q(m,n)=m^2-2mn-n^2=a-b.
Q(2m+n,m)=-Q(m,n).                                        (6)
```

Starting at `(2,1)` therefore gives

```
(3,4,5), (21,20,29), (119,120,169), (697,696,985), ...
```

with `|a-b|=1` throughout. Because `b` grows unboundedly, `a/b -> 1` and
`a/c,b/c -> 1/sqrt(2)`. The underlying operation even has a finite-set model:
`(U,V)->(U disjoint_union U disjoint_union V,U)`. None of these operations
has yet supplied a cubic graph or preserved a snark predicate.

## 3. An exact structural bridge among K5, Petersen and K3,3

Let `D(G)` be the graph whose vertices are edges of `G`, adjacent exactly when
the two edges are disjoint; equivalently `D(G)` is the complement of `L(G)`.
Then `D(K_5)` is the Petersen graph `P`. Explicitly, on labels modulo five,

```
outer_i={i,i+2},          inner_i={i+3,i+4}
```

give the standard outer 5-cycle, inner star, and five spokes. This explicit
map, including every edge and nonedge, is checked in the companion script.

There is a second exact step which reaches the bipartite obstruction. In the
two-subset model of `P`, take `v={0,1}` and `T={2,3,4}`. Its neighborhood `N(v)`
is the three two-subsets of `T`. Perform **Seidel switching** on `N(v)`: toggle
every edge crossing this set and its complement, leaving all other pairs
unchanged. Then `v` becomes isolated. Label the other nine vertices by
`(a,i) in {0,1,2} x {0,1,2}`, with `t_i` enumerating `T`, as follows:

```
(0,i) -> {0,t_i},     (1,i) -> {1,t_i},     (2,i) -> T-{t_i}.
```

After switching, two such vertices are adjacent iff `a!=b` and `i!=j`.
For the first two rows this is ordinary disjointness; between the third row
and either other row, switching reverses the original predicate `i=j`.
Within a row there are no edges. Consequently

```
switch_(N(v))(P) = K_1 disjoint_union D(K_(3,3)).           (7)
```

Moreover `D(K_(3,3))` is isomorphic to `L(K_(3,3))`. On edge labels in `F_3^2`,
the bijection `(a,i)->(a+i,a-i)` takes `a!=b and i!=j` to sharing exactly one
coordinate: two nonzero elements of `F_3` are equal or negatives. Therefore

```
K5 --D--> P --switch N(v), remove isolated v--> L(K3,3).    (8)
```

This is a genuine chain of graph operations connecting the three named
graphs. It is not an isomorphism between them. Switching destroys cubicity:
the image has one isolated vertex and nine degree-four vertices. Nor does
these operations' existence preserve planar embeddability or edge coloring.

The obvious extension of the first operation leaves snarks immediately:
`D(K_c)` has `binom(c,2)` vertices and degree `binom(c-2,2)`. It is cubic
exactly when `c=5` (among `c>=5`). The next unit-gap clique `K_13` gives degree
55; the next central-branch clique `K_29` gives degree 351.

## 4. Seidel spectral radius gives the literal 3,4,5 root

For a simple graph on `N` vertices, define `S=J-I-2A` and let `rho(S)` be its
spectral radius, the largest absolute value of an eigenvalue. Then

```
rho(S(P))=3,       rho(S(K5))=4,       rho(S(K3,3))=5.      (9)
```

Here is a proof without numerical diagonalization.

- `S(K_N)=I-J` has eigenvalues `-(N-1)` once and `1` with multiplicity `N-1`.
- For `K_(p,q)`, let `s` be `+1` on one shore and `-1` on the other. Then
  `S=ss^T-I`, with eigenvalues `N-1` once and `-1` with multiplicity `N-1`.
- In the Petersen two-subset model every adjacent pair has zero common
  neighbors and every distinct nonadjacent pair has one. Thus
  `A^2=2I-A+J`. Since `AJ=JA=3J`, expansion gives `S(P)^2=9I`.

Equation (9) is an exact common-invariant realization of `3^2+4^2=5^2`.
It is stronger than choosing unrelated graph statistics for the three
numbers. It remains a statement about these three graphs, rather than a
classification or a planar-obstruction proof.

The operation (7) explains why a conference Seidel matrix of order ten and the
3-by-3 rook graph occur together. Switching is conjugation `S -> D S D` by
a diagonal sign matrix, so it preserves the entire Seidel spectrum. At a
fixed order, complete and complete bipartite graphs have opposite Seidel
spectra and identical Seidel spectral radii. Switching carries the complete
bipartite graph to the empty graph, whose Seidel matrix is the negative of
the complete graph's. Thus absolute spectrum cannot distinguish their roles.

The normalization also has a sharp hostile family. For every simple cubic
graph on `N>=14` vertices,

```
rho(S)=N-7.                                                (10)
```

Indeed the all-ones vector has eigenvalue `N-7`. On its orthogonal complement
`S=-I-2A`; all adjacency eigenvalues lie in `[-3,3]`, as follows from the
maximum-coordinate bound for a 3-regular matrix. Hence all restricted Seidel
eigenvalues are in `[-7,5]`, while `N-7>=7`. This proves (10), regardless of
whether the graph is bipartite, edge-colorable, a snark, or disconnected.

More sharply, the Petersen graph is the **unique simple cubic graph** whose
Seidel matrix satisfies `S^2=(N-1)I`. To prove this, its constant-vector
eigenvalue forces

```
(N-7)^2=N-1,             so N=5 or N=10.
```

Cubic order must be even, leaving `N=10`. Expansion then forces
`A^2=2I-A+J`, so the graph has no triangles or 4-cycles. From any vertex its
three neighbors and six distinct grandchildren exhaust the graph. Each
grandchild has two neighbors among those six; they form a 6-cycle, since
triangles are forbidden and a simple 2-regular graph on six vertices has no
other triangle-free decomposition. Sibling grandchildren must be opposite
on that cycle: adjacency makes a triangle, distance two makes a 4-cycle.
This uniquely reconstructs Petersen.

For every simple graph, `tr(S^2)=N(N-1)`, so
`rho(S)>=sqrt(N-1)`, with equality iff `S^2=(N-1)I`. Petersen achieves this
spectral balance exactly. By (10), every growing cubic family instead has
`rho(S)/sqrt(N-1)=(N-7)/sqrt(N-1) -> infinity`. Consequently this particular
balance is **not** an asymptotic signature of large snarks.

## 5. An infinite extension through conference Seidel matrices

The scalar identity (9) extends naturally, although its Petersen coordinate
continues as a **regular graph with a conference Seidel matrix**. The following
construction is proved here; no claim that the construction is new is intended.
Here a conference Seidel matrix means `S^2=(N-1)I`; this does not call its graph
a conference graph in the standard strongly regular graph terminology.

Let `t` be an odd prime power, `q=t^2`, and let `chi` be the quadratic character
of `F_q`, with `chi(0)=0`. Define the Paley graph on `F_q` by adjacency
`chi(x-y)=1`. This is undirected because `chi(-1)=1`. Put

```
Q_xy=chi(x-y),           Q*1=0,           Q^2=qI-J.         (11)
```

For completeness, the off-diagonal character sum is `-1`. Translation and
scaling reduce it to `sum_u chi(u(u-1))`. The number of solutions of
`v^2=u(u-1)` equals both `q+sum chi(u(u-1))` and `q-1`: after multiplying by
four, factor `(2u-1-2v)(2u-1+2v)=1` and choose any nonzero first factor.
The diagonal sum is `q-1`, proving (11).

Adjoin an isolated vertex. The resulting Seidel matrix is

```
S = [[0,1^T],[1,-Q]],             S^2=qI.                  (12)
```

Thus its Seidel radius is `t`. One can also choose a regular graph in this
switching class, as follows. Regard `F_q` as an affine plane over `F_t`.
Choose a nonzero direction `alpha` with `chi(alpha)=-1`. All nonzero elements
of `F_t` are squares in `F_q`, so every line parallel to `alpha F_t` is an
independent set. Every vertex on a different parallel line has exactly
`(t+1)/2` neighbors on it.

Here is the last count's proof. The identity
`chi_q(z)=chi_t(Norm(z))` reduces the line sum to the character sum of the
quadratic `Norm(beta+u alpha)` over `F_t`, where `beta` is not in `alpha F_t`.
Its leading coefficient is `Norm(alpha)` and its discriminant is
`(beta*alpha^t-beta^t*alpha)^2`, nonzero (and nonsquare in `F_t`). Any quadratic
with nonzero discriminant has character sum minus the character of its
leading coefficient: complete the square and count
`(x-y)(x+y)=D` for the resulting nonzero `D`, just as above. Here that sum is
`-chi_q(alpha)=1`, yielding `(t+1)/2` neighbors.

Switch the union of any `k=(t-1)/2` such parallel lines in the graph (12), and
call the resulting graph `H_t`. It has `t^2+1` vertices and degree

```
d=t(t-1)/2.                                                (13)
```

Indeed the isolated vertex acquires `kt=d` neighbors. An unselected field
vertex changes degree from `(t^2-1)/2` to that number minus `k`, also `d`.
A selected vertex has `(k-1)(t+1)/2` neighbors within the selected set;
toggling its other incidences changes its total degree to
`(t^2+1)/2+k-t=d`. Switching preserves (12), so `rho(S(H_t))=t`.

Define `h=(t^2+3)/4`, an integer because `t` is odd. The three graph radii are
then exactly

```
(rho(S(H_t)), rho(S(K_((t^2+1)/2))), rho(S(K_(h,h))))
  = (t, (t^2-1)/2, (t^2+1)/2).                            (14)
```

They form a primitive Pythagorean triple, since the last two entries are
consecutive and the difference of their squares is `t^2`. Its Euclid
parameters are `((t+1)/2,(t-1)/2)`. The complete graph order in (14) is
precisely the hypotenuse, agreeing with the marked-clique construction (1).
The balanced bipartite graph is a different graph from that construction's
square-shore cut, and that distinction is essential.

| t | H_t order | H_t degree | Complete graph | Bipartite graph | Radii |
|---:|---:|---:|---|---|---|
| 3 | 10 | 3 | K5 | K3,3 | (3,4,5) |
| 5 | 26 | 10 | K13 | K7,7 | (5,12,13) |
| 7 | 50 | 21 | K25 | K13,13 | (7,24,25) |
| 11 | 122 | 55 | K61 | K31,31 | (11,60,61) |

At `t=3`, the uniqueness proved in section 4 identifies `H_3` as Petersen.
For `t>3`, degree (13) exceeds three, so this family cannot consist of
snarks. In fact its normalized radii converge to `(0,1,1)`, precisely the
parabolic limit, rather than to an isosceles triangle. This is an exact
infinite-family repair of the proposed three-graph correspondence with a
precise boundary to the snark interpretation. This construction covers odd
prime-power `t`; no assertion is made here about every other odd `t`.

At the root, two operations happen to coincide: `H_3=D(K5)`. The family
`D(K_c)` from section 3 and the family `H_t` from (13) are different away
from that root. Their vertex counts already disagree: `binom(c,2)` versus
`2c`, which are equal for positive `c` only at `c=5`.

## 6. Transfer and stopping ledger

| Source | Target/map | Preserved | Lost / necessary sidecar |
|---|---|---|---|
| Marked clique with square shores | (1) | Exact Euclid triple, shore cut radius | Arbitrary graph topology; retain marked partition |
| Two complete graphs | Relation count and doubled bipartite cut | Identity (3) | Planarity and coloring; those require incidence data |
| Euclid spinor | Parabolic or central iteration | Unit hypotenuse gap or unit leg gap, respectively | Confusing the two destroys the limiting shape |
| K5 | Edge-disjointness graph P | A concrete incidence construction | Original graph not preserved as a minor certificate |
| P with distinguished vertex | Switched isolated vertex plus D(K3,3) | Seidel spectrum and explicit pair map | Degrees and cubic edge-coloring predicate |
| The named triple of graphs | Seidel radii (3,4,5) | One common spectral invariant | Complete versus bipartite role at fixed order |
| Square-order Paley switching class | Graph triple (14) | Entire primitive unit-gap family for odd prime-power t | Cubicity beyond t=3; limit is parabolic |
| Cubic graph | Large-order Seidel radius N-7 | Order | Snark and edge-colorability information |

The successful correspondence is therefore threefold: a faithful arithmetic
map on **marked** complete graphs, an exact graph-operation chain among the
three named examples, and a common spectral `3,4,5` identity extending to
the infinite family with conference Seidel matrices (14). The
obstruction to the proposed theorem proof is now explicit: neither graph
operation nor the scalar radius preserves the required forbidden-subdivision
or cubic-coloring predicate. A theorem transfer would require such a sidecar.

## 7. Exact controls and reproduction

```
python 04-computation/graph_pythagorean_20261018_arithmetic_probe.py
```

The companion output is
[`graph_pythagorean_20261018_arithmetic_probe.out`](graph_pythagorean_20261018_arithmetic_probe.out).
The script uses standard-library integer arithmetic and explicit `check`
calls, so optimized Python does not remove the tests. It checks:

- all positive shore pairs `x>=y` with `x+y<=200`: 279 integral triples,
  including 32 positive primitive triples;
- all 2040 primitive Euclid pairs with `1<=n<m<=100`, including 99 unit-gap
  pairs, and eight central-branch states;
- every edge/nonedge of the Petersen model and both maps in (7)--(8);
- `S(P)^2=9I` and complete/bipartite minimal-polynomial identities at every
  order 2--20 and every nonempty shore split;
- explicit finite-field construction of the regular graphs `H_t` for
  `t=3,5,7,11`, checking every degree and every entry of `S(H_t)^2=t^2 I`;
- direct disjointness degrees at clique orders 5, 6, 13, 29;
- four direct clique/cut edge-count identities and missing-edge hostiles.

These finite tests audit the implementation. The displayed algebra and
explicit maps prove the unrestricted claims in their stated domains.
