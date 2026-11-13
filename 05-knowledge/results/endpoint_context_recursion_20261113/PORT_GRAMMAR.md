# Changing the output port closes three local obstructions

**PROVED, scoped, by the arguments and finite symbolic certificates below;
FINITE-EXACT for the census. Chan–Pak remains OPEN.** No literature-priority
claim is made. Local identifiers PG1–PG4 refer only to this file.

The inherited starting point is
[TRANSPORT](../surplus_selector_20261113/TRANSPORT.md), with its actual paired
graphs in `transport_words.json.gz`. Its fixed-root context rules left slopes
1/3, 4/9 and 6/11 without a next consumer. The present change is to specify
an entry vertex and a possibly different output vertex of the context.
The use of boundary-conditioned path tuples is inherited methodology:
[THM-4879, two-port transitive seam gadgets](../../../01-canon/theorems/THM-4879-two-port-transitive-seam-gadgets-closed-transfer-matrices-with-a-three-power-entry-and-the-seam-parity-r-equals-s-mod-2.md)
already gives transfers for a different, disjoint single-seam operation.
Here the operands share a vertex, and deletion count is retained explicitly.

## PG1. A four-coordinate compiler with two specified ports

Write z(B,r)=(h,s,q,F): total Hamilton paths; those starting at r; H(B−r);
and ordered two-path covers whose first path is nonempty and ends at r,
with the second possibly empty. H(empty)=1. Let A have entry r and exit u.
Identify the root of B with r, and orient every other cross edge
A−r → B−r. Denote the resulting tournament A⋆B.

For each fixed (A,r,u), with |A|>=2, there is an integer matrix M(A,r,u), independent of B with |B|>=2,
such that

    z(A⋆B,u) = M(A,r,u) z(B,r).

The singleton context is identity. For all other contexts the number of
vertices increases by |A|−1. These are actual tournament constructions.
Any word of these operations therefore has its exact matrix product and
its literal graph word; equality of a numerical profile alone is not an
inverse construction certificate.

Proof and coefficient algorithm. Put Y=B−r. In a full path, or a pair of
paths whose first path ends at u∈A, the vertices of Y form at most two
nonempty consecutive runs. Indeed a Y-run can return to A only through r,
which is used once. In the two-path case a run in the first path must
return before its final vertex u, so at most one further run can lie in
the second path. Enumerate permutations of A and one or two Y-markers,
with no consecutive markers within a path. For F, cut immediately after
u. A marker is a nonempty path on some subset of Y. Its internal size
is not fixed. The ordered run weights are as follows, using temporarily
e=end-root paths, G=the converse F, and P=P_Y(2):

| Runs | Adjacencies to r inside a path | Weight |
|---|---|---|
| one | neither | q |
| one | r→Y | s |
| one | Y→r | e |
| two | Y→r→Y | h−s−e |
| two | Y→r only | F−q−e |
| two | r→Y only | G−q−s |
| two | neither | P−2q |

Every full path or qualifying cover has a unique such run decomposition.
For the three observers H, start-u, and first-end-u, every valid two-run
word has a Y→r adjacency: the earlier run cannot end a path at u and
cannot return to A elsewhere. Thus the last two rows never contribute.
The e terms cancel bijectively. In a one-run word with Y→r, insert a
second Y-marker at the end of the path (H/start-u), or at the end of the
second path (first-end-u, including an initially empty second path).
This produces exactly one of the two-run words with coefficient −e;
deleting that final marker is the inverse. The start-u restriction is
preserved. Therefore only h,s,q,F remain.

For deletion of u≠r, apply the H observer to A−u with entry r. For u=r,
deletion leaves (A−r)→(B−r), so q'=H(A−r)q. This proves the fourth row and
the entire integer matrix law. The proof handles nonempty B−r; B singleton
can either be checked directly or included separately as the identity
attachment. The research language below uses |B|≥2 throughout.

`port_compiler.py` performs this finite symbolic enumeration. It is not a
fit to a finite family of B: the seven run weights are universal. Its
independent physical controls cover all labelled A of orders2–3, every
entry/exit, and every labelled B of orders2–4 (5,920 controls).

## PG2. Exact escape rules for all three stranded slopes

A repair pair has difference d(a,b)=(0,a,b,a), b≠0. It preserves both
h and the next rooted-birth target E=h+F−s. Its slope is a/b. A port matrix
is a consumer when its image is another such nonzero difference.

Rows below are adjacency bit masks; bit j in row i means i→j. Vertex
labels start at zero. Multipliers act on the displayed primitive direction.

| Input | Output | Context rows | Entry | Exit | Image of primitive difference |
|---|---|---|---:|---:|---|
| 1/3 | 1/2 | (2,20,9,35,45,7) | 2 | 5 | (0,−3,−6,−3) |
| 4/9 | 1 | (2,12,33,21,39,11) | 5 | 3 | (0,−2,−2,−2) |
| 6/11 | 1/6 | (66,84,73,67,109,79,0) | 0 | 3 | (0,1,6,1) |
| 1/2 | 1/2 | (2,4,9,3) | 3 | 2 | (0,−1,−2,−1) |

The exact matrices and all 22 observed transition types for the four input
slopes are in `port_search.json.gz`. The first three rows close the local
fixed-root composition failures; they do not solve an exact-target inverse
problem. In particular the output h generally differs from the input h.
The 1/6 output has an inherited fixed-root consumer taking it to zero;
1 has an inherited self-consumer. Their usefulness for target selection
requires more than mere nonzero composition.

Census scope: every ordinary tournament through8, every entry satisfying
current-H neutrality for one of {1/3,4/9,6/11,1/2}, every output vertex in
the context. There are 1,310 eligible entries and 10,139 entry/exit pairs.
Four independent realizable input profiles determine each matrix; every
reported first representative is then certified by the symbolic compiler.
Minimum context sizes for those three escape directions are6,6,7 in this
complete head. These are not all-order classifications of all consumers.

## PG3. A renewable four-vertex context

For A=(2,4,9,3), entry3, exit2, the matrix is

    M = [ 3 −2 −2  6 ]
        [ 1 −1 −1  2 ]
        [ 1 −1 −2  3 ]
        [ 3 −2 −3  7 ].

It sends d=(0,1,2,1) to −d. Consequently EVERY actual repair pair with
difference t d, t≠0, can be iterated indefinitely, with difference
(−1)^k t d. No additional arithmetic entry test is needed. The port is
part of each literal context, so recursive entry is physical, not inferred
from a divisibility calculation. Each step adds three vertices.

All outputs satisfy F+s=h+q. Their output root has indegree one, from
context vertex1, regardless of B. This is weaker than unit readiness:
its defect q−s=F−h equals the INPUT F−q and need not be1.

On F+s=h+q, put E=h+F−s. The observed target pair obeys

    (h',E') = (h+4E, 2h+7E).

This matrix has determinant −1 and characteristic polynomial x²−8x−1.
Hence both scalar sequences satisfy X_(k+2)=8X_(k+1)+X_k. The full matrix
has characteristic polynomial x(x+1)(x²−8x−1): the −1 mode transports
the exact repair difference while the targets grow.

The smallest seed is already H5: states(5,2,3,6) and(5,1,1,5), with
common E9. It gives H=5,41,333,2705,... under the return. See
[LANGUAGE](LANGUAGE.md) for a second strong context, its physical
factorization, and the exact finite-language inverse.

The inherited H605 pair can also be re-rooted at vertex3, with states

    (605,208,319,716),  (605,202,307,710).

They have E=1113 and differ by −6d. The resulting target pairs begin

    (605,1113), (5057,9001), (41061,73121), ... .

This proves an infinite family of actual renewable pairs. It is a sparse
family, and no assertion that these targets exhaust the odd integers is
made. The larger eight-vertex return discovered first is retained in
`port_compiler.json.gz`; it has target recurrence X_(k+2)=144X_(k+1)+5X_k
on this initial hyperplane and a determinant−5 target map. Its targets
E after the first step are divisible by5. The smaller context removes
that particular determinant obstruction, not the global coverage problem.

## PG4. What target-decreasing inversion still requires

For the four-vertex return the arithmetic inverse is

    (h,E) = (4E'−7h', 2h'−E').

If 7h'/4<E'<2h', both are positive and h<h'. But arithmetic inversion
does not prove that this state is realized by a tournament, nor that a
given tournament admits this context decomposition. A complete proof
needs a cover of target/endpoint states by physically invertible rules,
with every selected predecessor realized and a well-founded decrease.
The present rule proves renewable forward entry and gives a precise
inverse test; it does not supply that cover.

Every context word retains its input as an induced subtournament. The
incoming all-order13→27 containment obstruction therefore remains in
force even after changing ports; see[OBSTRUCTIONS](OBSTRUCTIONS.md).

The attained H=53,084,231 still has its inherited all-order exclusion of
unit-ready representatives. Nothing here repairs it at the same target.
The exclusion concerns indegree one AND defect one, not all indegree-one
roots. Indeed two unit-ready states with equal h and E have identical
(h,s,q,F): F=h+1, s=2h+1−E, q=s+1. A nonzero repair difference cannot
remain inside that narrow class on both sides. This explains why varying
the readiness type is essential.
