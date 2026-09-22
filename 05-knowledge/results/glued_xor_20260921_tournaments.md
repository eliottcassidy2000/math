# Four vertices: two diagonal bits, three matching terms, and three distinct quotients

**2026-09-21. Status: PROVED elementary statements; FINITE-EXACT complete
order-four census and two named order-five controls.** No novelty claim.
The interpretation is conditional: a fixed directed four-cycle is meant
by the surrounding path. A fixed Hamiltonian path instead leaves three
bits, and its answers differ.

The two-diagonal idea is exact with that four-cycle fixed. Its four
completions are all isomorphic, each with five Hamiltonian paths. XOR
describes changes of the two labelled diagonal orientations, not four
different unlabelled tournament structures. Vertex switching is another
operation: it changes whole cuts, and can turn a transitive tournament
into a strongly connected one. Three signed perfect-matching terms give
an exact invariant linking the six arcs to a three-part description.

## Inheritance and the intrinsic relation

The closest proved mechanism is the path-and-chord model in
[the arithmetic-seams tournament note](arithmetic_seams_20260921_tournaments.md),
building on [THM-781, Hamiltonian-path inversion](../../01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md).
Its hostile is the transitive tournament: a full undirected cycle space
does not entail any directed cycle. The corrected near miss is identifying
a path-dependent chord operation with a base-independent quotient.
The least-used coordinate is which labelled vertices carry the degrees,
triangles, and selected Hamiltonian path.

The switching loss is already inherited from
[THM-1415, canonical star quotient](../../01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md),
whose tournament quotient differs from graph switching. The matching
action and its parity boundary are inherited from
[THM-2753, six-edge parity and three-matchings](../../01-canon/theorems/THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration.md).
The acting-group caution is shared with
[THM-3173, six-state actions and pointed frame cube](../../01-canon/theorems/THM-3173-six-state-free-factor-actions-and-pointed-frame-cube.md).
The present contribution is a complete small model specifying which
binary operation is being used and exhibiting its preserved and lost data.

Vertices are four labelled objects, with one actual chosen orientation
on every distinct pair. The orientation itself is the pairwise observable;
there are no ties. This is a study of tournaments already supplied as
oriented complete graphs. No tournament is manufactured from equal
arithmetic values, prime frequencies, or an unspecified dominance rule.

| Live object / operation | Preserved predicate | Required extra coordinate |
|---|---|---|
| Fixed directed cycle and diagonal XOR | The four cycle arcs and unlabelled type | Labels and individual diagonal bits |
| Vertex switching over F2 | Cycle parity, signed-matrix determinant | Original vertex signs to restore directions |
| Vertex relabelling | All unlabelled directed structure | Which vertex was distinguished |
| Three perfect matchings | Pfaffian formula with signed products | Orientation and vertex-order signs |
| Fixed-path presentation | A chosen directed Hamiltonian order | Presentation multiplicity H(T) |

Anchor: make the two-edge claim exact. Niche: determine the switching
quotient. Wildcard: explain six arcs versus three matchings without
discarding sign information. The research move is to specify the action
before identifying two cubes or two Klein four-groups.

## 1. A fixed directed square really does leave two independent bits

Fix

```text
0 -> 1 -> 2 -> 3 -> 0.
```

Write `a=1` if `0->2`, and `b=1` if `1->3`. Otherwise each diagonal points
the other way. The four states `(a,b)` form `F2^2`; independent diagonal
reversals are the commuting translations

```text
A(a,b)=(a XOR 1,b),  B(a,b)=(a,b XOR 1).
```

**PROVED.** All four completions are isomorphic. Indeed the vertex
rotation `i -> i+1 mod4` acts by

```text
R(a,b)=(1-b,a):  00 -> 10 -> 11 -> 01 -> 00.          (1)
```

This is an affine operation on the bit square, not merely an exchange of
its coordinates. In particular `R^2=A B`, the simultaneous reversal of
both diagonals, and `R A R^-1=B`. The flips and R together act as the
eight symmetries of the bit square; the two flips alone form `C2 x C2`.

There is no single vertex relabelling that implements A on all four
states. Such a relabelling must preserve their four common directed arcs,
hence be a rotation of the directed square. None of the four powers of R
is A. A single flip is isomorphic to a suitable rotation at each state,
but the suitable rotation depends on that state. This is the difference
between a uniform action and a pointwise existence of an isomorphism.

The labelled outdegrees are

```text
(1+a, 1+b, 2-a, 2-b),
```

and the cyclic-triangle indicators on `012,123,013,023` are

```text
(1-a, 1-b, b, a).                                   (2)
```

Thus all four have degree multiset `(1,1,2,2)` and exactly two directed
triangles, but the vertices and triples carrying those properties change.
Each has five Hamiltonian paths, one Hamiltonian cycle up to cyclic
rotation, and trivial automorphism group. Triviality of the automorphism
group follows because the two vertices of each given outdegree cannot
be exchanged: their mutual directed arc would reverse. The remaining
counts are proved by the eight-state table below and independently
checked by exhaustive path enumeration and dynamic programming.

In particular the square rotation in (1) is an isomorphism between
different labelled completions; it is not a nontrivial automorphism of
any one completion.

## 2. A fixed Hamiltonian path instead leaves three bits

Prescribe only `0->1->2->3`, and keep a,b as above. The third free bit is
`c=1` for `0->3`; `c=0` closes the directed square. The complete table is:

| a | b | c | Type | Number H of Hamiltonian paths |
|---:|---:|---:|---|---:|
|0|0|0|strongly connected|5|
|1|0|0|strongly connected|5|
|0|1|0|strongly connected|5|
|1|1|0|strongly connected|5|
|0|0|1|strongly connected|5|
|1|0|1|source plus directed triangle|3|
|0|1|1|directed triangle plus sink|3|
|1|1|1|transitive|1|

The triangle indicators directly give

```text
c3=2-c(a+b),  H=1+2c3=5-2c(a+b).                    (3)
```

The H entries can be checked by the 24 possible vertex orders. This
small instance also agrees with the inherited odd-cycle-collection
formula in the arithmetic-seams note; no longer directed odd cycle can
occur on four vertices.

Every tournament has a Hamiltonian path: inductively insert the new
vertex before the first path vertex it beats, or at the end if it beats
none. When insertion is internal, the preceding vertex beats it by
the choice of the first position. Thus after relabelling, the table
contains every order-four tournament. Its five strong entries all
contain a directed four-cycle, so Section1 identifies their type. The
other three are distinguished by degrees. This proves the four-type
classification without presuming strong connectivity implies a cycle
in arbitrary directed graphs.

The complete labelled census is:

| Type | Sorted outdegrees | c3 | H | Automorphism order | Labelled count |
|---|---|---:|---:|---:|---:|
|transitive|0,1,2,3|0|1|1|24|
|source plus directed triangle|1,1,1,3|1|3|3|8|
|directed triangle plus sink|0,2,2,2|1|3|3|8|
|strongly connected|1,1,2,2|2|5|1|24|

The counts sum to `2^6=64`. Source and sink types are exchanged by
reversing every arc; they are distinct directed isomorphism classes.

## 3. Vertex switching is a three-dimensional gauge quotient

Now encode every edge by `x_ij=1` for `i->j`, for `i<j`. Choose vertex
bits `t_i`. Reversing the cut between the chosen and unchosen vertices
acts by

```text
x_ij -> x_ij XOR t_i XOR t_j.                        (4)
```

Adding one to all four t_i changes nothing. Consequently the cut space
has dimension three, and the quotient of the six-dimensional edge space
has dimension three. There are eight labelled switching classes, each
of size eight. With the labels and reference orientation fixed, one
complete invariant is

```text
(x01 XOR x02 XOR x12,
 x01 XOR x03 XOR x13,
 x02 XOR x03 XOR x23).                               (5)
```

Every cut flips zero or two edges on each triangle, so (5) is invariant.
Conversely, set `t0=0` and determine `t1,t2,t3` from the edges incident
to zero; matching the three triangle parities then matches the remaining
three edges. This proves completeness. A parity coordinate is not the
indicator that its triangle is directed cyclically. Its reference edge
ordering also matters under relabelling.

**PROVED.** Prescribing any orientation on a fixed spanning tree gives
exactly one representative of each switching class. Starting at a tree
root, (4) recursively determines all vertex bits, uniquely up to the
irrelevant common bit. Hence the eight fixed-path states in Section2
form a transversal of the switching classes, as well as a conditional
family of genuinely directed Hamiltonian paths.

For the four fixed-cycle states, (5) becomes

```text
(a,1 XOR b,1 XOR a).                                (6)
```

They occupy four different switching classes, an affine plane in the
three-bit quotient. This can also be seen without coordinates: a
nontrivial cut in K4 has three or four edges, whereas the difference of
two such completions has only one or two diagonals. Therefore the two
diagonal flips are not vertex switches. Isomorphism and switching are
different equivalence relations, even on these four examples.

**Hostile.** Start with the transitive order `0->1->2->3` and all its
forward chords. Switch vertex1. The result contains the directed cycle
`0->2->3->1->0`; it is strong, and H changes from1 to5. Thus switching
does not preserve transitivity, strong connectivity, directed triangles,
or Hamiltonian-path count. Fixing a path by switching is not a harmless
replacement of the original tournament for these questions.

The complete switching-orbit census is six classes containing four
transitive and four strong tournaments each, and two classes containing
four source-triangle and four sink-triangle tournaments each. Up to
relabelling as well as switching there are two classes, recovering the
order-four part of THM-1415. Each labelled switching class has mean H=3,
consistent with the general inherited mean `n!/2^(n-1)` in
[THM-3729, rooted Pfaffian response](../../01-canon/theorems/THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average.md).

## 4. Six oriented edges enter three signed matching terms

Let S be the skew sign matrix: `s_ij=+1` for `i->j`, `-1` for `j->i`,
and zero on the diagonal. Its order-four Pfaffian is

```text
pf(S)=s01*s23 - s02*s13 + s03*s12.                    (7)
```

The terms are the three perfect matchings `01|23,02|13,03|12`. The minus
sign and the directed signs are essential data. Direct determinant
expansion gives `det(S)=pf(S)^2`, which here is either1 or9.

Vertex switching is `S -> D S D`, with diagonal signs D. Every matching
term is multiplied by the same product of all four diagonal signs;
thus pf changes by that common sign and its square is unchanged.
Relabelling vertices by p gives the separate law

```text
pf(p.S)=sign(p) pf(S).                               (8)
```

The ambient permutation of the six **unoriented** edges always has even
sign, as established in THM-2753. Hence (8) cannot be recovered from that
six-edge permutation sign: the directions and the matching ordering
restore information the unoriented statistic discarded.

For the fixed directed square, the first and third terms of (7) cancel:

```text
pf(S)=-s02*s13=-(-1)^(a XOR b),  det(S)=1.             (9)
```

This is a literal XOR of the diagonal bits expressed as a signed
matching product. Reversing one diagonal changes the Pfaffian sign;
reversing both does not. The Pfaffian sign is not an unlabelled invariant
and is not switching-invariant. Its square loses this binary bit.

Across all64 tournaments,

```text
det(S)=1+8*(c3 mod2).                                (10)
```

It separates the two switching-plus-isomorphism classes exactly:
det1 for transitive/strong, det9 for source-triangle/sink-triangle.
It does not separate all four isomorphism classes, and in particular
cannot recover H. Formula (10) follows from the four-type table or
direct expansion, and is independently checked for every orientation.

There are several different groups of order four here. The two diagonal
flips depend on the chosen square. The normal Klein group in
`S4 -> S3` acting on the three matchings consists instead of vertex double
transpositions. A single diagonal flip has no uniform vertex-relabelling
realization, so these actions cannot be identified merely because both
groups are abstractly `C2 x C2`.

## 5. The chosen path biases what one sees

The 64 labelled tournaments have total H equal to192 and total H squared
equal to768. Therefore uniform sampling of labelled tournaments has
mean H=3. Uniform sampling of path-marked presentations instead weights
a tournament by its number H of Hamiltonian paths, giving mean

```text
sum H^2 / sum H=768/192=4.
```

Equivalently the eight states with one fixed Hamiltonian path have mean
H=4, and five of the eight are strong. With a directed four-cycle fixed,
all four states are strong and mean H=5. These are three different
sampling spaces, not contradictory counts. This is the concrete
presentation-multiplicity sidecar inherited from THM-781.

The one-type collapse after fixing a cycle is special to order four.
Two named order-five controls both contain `0->1->2->3->4->0`:

* the transitive order with only `0->4` reversed has sorted degrees
  `(1,1,2,3,3)` and H=9;
* the cyclic regular tournament `i->i+1,i+2 mod5` has degrees
  `(2,2,2,2,2)` and H=15.

They are not isomorphic. The script checks these two examples only; it
does not claim an exhaustive order-five classification. In general a
fixed Hamiltonian path leaves `(n-1)(n-2)/2` free arcs, whereas a fixed
directed Hamiltonian cycle leaves `n(n-3)/2`. Their exceptional values
three and two at n=4 explain the two different small cubes.

## Reproduction and exact scope

```powershell
python 04-computation/experiments/glued_xor_20260921_tournaments.py
python -O 04-computation/experiments/glued_xor_20260921_tournaments.py
```

The standard-library script exhausts all64 orientations, all24 vertex
relabellings, all8 cuts, and all4096 pairs for completeness of the
switching invariant. It compares direct and degree triangle counts,
permutation and dynamic-programming path counts, and determinant and
Pfaffian evaluations. The JSON contains all64 rows, the eight fixed-path
states, all switching classes, explicit hostiles, source LF hash, and
scope. Normal and optimized decoded JSON must agree. There is no forced
numeric tournament, no general cycle-classification claim, and no
implication for Collatz convergence.
