# Native triangle parity, a failed Peck source, and forced Boolean zero modes

Status: **PROVED ELEMENTARY / FINITE-EXACT / INDEPENDENTLY AUDITED**.
Stanley's theorem below is explicitly **CITED**.
No exact Boolean nullity or Laplacian-gap formula is asserted. No new ID.

## 1. The source hypothesis is stronger than rank symmetry

The native graph is the Eulerian isomorphism-class triangle-toggle graph,
not an anchored-coordinate Boolean lattice. The closest mechanisms are
[THM-4069 — basis dependence](../../01-canon/theorems/THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope.md)
and the audited
[native Boolean flow note](overnight_hexagon_sep05_boolean_flow.md).
The discarded coordinate to recover is an invariant rank with the *actual*
cover relation, not the Hamming weight in a non-invariant triangle basis.
The live board is: cover rank; unit incidence; invariant subspaces;
edge parity; twisted fixed spaces; forced zero modes.

[Stanley, *Quotients of Peck Posets*, Theorem 1, pp. 29–32](https://math.mit.edu/~rstan/pubs/pubfiles/60.pdf)
proves that a finite-group quotient of a **unitary Peck** poset is Peck.
The source must have invertible complementary-rank products of its
unit-coefficient upward incidence maps. His two following cautions matter:
ordinary Peck alone is insufficient, and a quotient need not be unitary.
The paper does not supply a spectral-gap theorem for a Booleanized operator.

Here the direct source model fails already at n=4. The eight labelled
Eulerian graphs are empty, four triangles, and three four-cycles. Toggling
any of the four triangles gives the complete bipartite graph K_(4,4), with
even-edge side empty plus the three four-cycles. The three S_4-orbit sizes
are `(1,4,3)` in the order empty,triangle,C4.

Suppose an S_4-invariant grading makes **every native triangle edge** a
cover. Its orbit ranks a,b,c satisfy `|a-b|=|c-b|=1`. After normalizing
the smallest rank to zero, exactly four cases remain:

```text
(a,b,c)=(0,1,0), (1,0,1), (0,1,2), (2,1,0).
```

The three-rank cases have profiles `(1,4,3)` or `(3,4,1)` and fail rank
symmetry. The two-rank cases have profile `(4,4)`, but their upward map
is the all-ones 4-by-4 matrix, of rank one rather than four. Thus **no
S_4-invariant unitary-Peck grading has this full labelled cover graph**.
This does not rule out a different source object, a properly retained
additional grading, or a construction at other individual orders.

The Boolean quotient is nevertheless the three-vertex path, which itself
is unitary Peck with ranks `(1,1,1)`. Positivity of that quotient does not
repair a missing hypothesis on its labelled source. Even replacing the
unit entries in the two-rank source by S_4-equivariant weights cannot help:
the empty--triangle and C4--triangle edge orbits each have a constant
weight, so the four incidence rows are still identical.

The intrinsic distance from empty is another cheap hostile. At n=5 its
rank profile is `(1,1,2,2,1)`, not rank-symmetric. This does not prove that
no quotient regrading works. In fact putting empty,C4 at rank zero; all
odd-edge classes at rank one; and bowtie,K5 at rank two gives `(2,3,2)`
and an upward-square determinant of absolute value one. It is a positive
small quotient control, not an all-order construction or the missing
unitary source in Stanley's theorem.

## 2. The invariant parity does survive, and it has an exact index

For every n>=3, a triangle toggle reverses the parity of the number of
edges. Relabelling preserves that parity, so the native graph is bipartite.
Let `q_+` and `q_-` count its even- and odd-edge classes. In that order,

```text
B = [ 0  C ; C^T  0 ].
```

Consequently its adjacency spectrum is symmetric about zero and

```text
nullity(B)>=|q_+-q_-|.                               (1)
```

This concerns Boolean **adjacency**, not its combinatorial Laplacian.
The difference can be computed exactly without building the native graph.
For a vertex permutation g, let V_g be its fixed Eulerian cycle space and
let `epsilon(F)=|E(F)| mod2`. A character sum on this binary vector space is

```text
sum_(F in V_g) (-1)^epsilon(F)
 = 2^dim(V_g), if epsilon is identically zero on V_g;
 = 0, otherwise.                                    (2)
```

Weighted orbit averaging therefore gives the twisted Burnside identity

```text
Delta_n:=q_+-q_- = sum_(lambda admissible) 2^f(lambda)/z(lambda) >=0, (3)
```

where `z(lambda)=prod_l l^m_l m_l!` and the already-proved fixed dimension is

```text
f(lambda)=sum_i floor(l_i/2)+sum_(i<j)gcd(l_i,l_j)
          -r+1_{some l_i odd}.
```

The **complete** admissibility criterion is:

```text
(i) every vertex-cycle length is even; or
(ii) there are exactly one or two fixed points, and every other
     cycle length is divisible by four.                         (4)
```

Here is a full parity-constraint proof. An odd vertex cycle of length at
least three supplies a fixed odd Eulerian cycle; three fixed vertices
supply a triangle. Either makes the sum zero. We may thus retain only
even cycles and at most two fixed points. On an even vertex cycle of
length l_i, its antipodal matching orbit has degree-parity column e_i and
edge-count parity `l_i/2 mod2`. Any covector expressing total edge parity
from degree parity must therefore have coefficient `c_i=l_i/2 mod2`.

All other internal edge orbits have even cardinality and even degree.
Between two even cycles of lengths l_i,l_j, put g=gcd(l_i,l_j). The degree
column has entries `l_j/g mod2` and `l_i/g mod2`, while its edge-orbit
length is even. If the two lengths have equal 2-adic depth, both entries
are odd and c_i=c_j; if their depths differ, only the deeper cycle has
an odd entry and its c is zero. Thus all even/even constraints hold.

Between a fixed vertex and an even cycle, the column is e_i with even
edge-orbit length, forcing c_i=0. This is possible precisely when every
even cycle has length divisible by four. With two fixed vertices, choose
their two coefficients to sum to one; this handles their single edge.
With one fixed vertex its coefficient is free. These choices establish
sufficiency, and the forced coefficients establish necessity. In
particular `(2,1)` is excluded: its invariant K3 is the minimal hostile
to forgetting the interaction with fixed points.

There are no admissible partitions when `n=3 mod4`. At even n the single
n-cycle is admissible; when `n=1 mod4`, `(1,n-1)` is admissible. Therefore
for every n>=3,

```text
q_+=q_- iff n=3 mod4;
q_+>q_- otherwise.                                  (5)
```

In the balanced odd case complementing by K_n independently gives the
parity-reversing involution on classes, since K_n is Eulerian and has an
odd number of edges. The stronger expression (3) retains the exact index
in the other cases. It forces native Boolean zero modes at every order
outside `n=3 mod4`, regardless of any multiplicity-weighted spectral law.

## 3. A different object counts the index at n=4k+1

For n=4k+1 the admissible cycle types are exactly `(1,4 mu)`, with mu a
partition of k having r parts. If `S=sum_(i<j)gcd(mu_i,mu_j)`, their
contribution simplifies to

```text
2^f(1,4mu)/z(1,4mu)=2^(2k-2r+4S)/z(mu).             (6)
```

Under a permutation of cycle type mu, the number of orbits on ordered
distinct vertex pairs is `k-r+2S`: within each vertex cycle there are
`mu_i-1` directed differences, and between two cycles there are two
directions, each with gcd-many orbits. An **ordered pair of loopless
digraphs** therefore has `2^(2k-2r+4S)` fixed assignments. Burnside applied
to this different object proves the exact identity

```text
Delta_(4k+1) = number of simultaneous-isomorphism classes of
               ordered pairs of loopless digraphs on k vertices. (7)
```

The two digraphs are distinguished; equivalently every ordered vertex
pair carries one of four colors, including absence in both. Loops are
forbidden, opposite arcs are independent, and the same vertex permutation
acts on both digraphs. This is not the tournament/Royle “even graph” count.
It follows in particular that

```text
nullity(B_(4k+1)) >= Delta_(4k+1)
                 >= 2^(2k(k-1))/k!.                  (8)
```

The connection is a cycle-type map `(1,4mu) -> mu` preserving a fixed-point
weight. It does not construct individual graph pairs from signed Eulerian
classes, or supply actual Boolean kernel vectors. Targeted repository
search found no earlier typed statement; no external novelty is claimed.

## 4. Exact controls and the remaining wall

The standalone
[checker](../../04-computation/overnight_hexagon_sep05_boolean_peck.py)
tests (4) against every free vector of an independently built fixed-space
basis for all 914 partitions through n=16. It then evaluates (3) through
n=30 without enumerating graph states. The first values are

```text
n:       3 4 5 6 7  8  9  10 11    12  13      14 15         16
Delta_n: 0 1 1 4 0 19 10 360  0 25112 720 6975976  0 7030565576.
```

Literal Eulerian orbit enumeration at n=3,...,6 independently verifies
both parity counts, the bipartition and (1); nullity equals the index at
those four orders only. The two Peck hostiles and the alternative n=5
quotient grading are checked by exact matrices. An independent literal
four-color ordered-arc enumeration gives 1,10,720 digraph-pair classes
at k=1,2,3, agreeing with (7).

Exact Boolean nullity at arbitrary order is still **OPEN**; a zero index
does not force adjacency invertibility. No Laplacian gap bound follows
from an adjacency zero mode. The representation-theoretic survivor is
the parity index and its twisted fixed-space sidecar, not an unjustified
Peck or weighted-spectrum transfer.

```bash
python3 -B 04-computation/overnight_hexagon_sep05_boolean_peck.py
python3 -B -O 04-computation/overnight_hexagon_sep05_boolean_peck.py
```

The companion [output](overnight_hexagon_sep05_boolean_peck.out) records
all controls and the exact semantic digest; no assertion is disabled by -O.

Root independently derived criterion(4), checked the complete written
proof and digraph-pair weight identity, and replayed all1,028 gates.
Normal and optimized outputs agree. Raw LF source SHA256:
`07f1f991ec6a8a491134a0812c996580fbba3b115747d07afae91f880125451b`;
output SHA256:
`efdb354c11bd996f5157b21d6c0fd43f8ed44dae94735eedc38585822fea1523`.
