# A minimal finite incidence menu, a one-bit repair, and its all-order limit

October 23, 2026. **PROVED:** the rooted strong-four identity and the
doubly-regular obstruction to every depth of the stated refinements.
**FINITE-EXACT:** all 192 specified feature menus on every ordinary
tournament class of orders 4 through 8, and a guarded one-bit repair at
order 8. **OPEN:** a small universally complete metagraph invariant.
Minimality below is within the explicit menu, not among all invariants.

## Inheritance and the target

The closest mechanism is the association of vertex and neighbor incidence
in [the bucket packet](../metagraph_bucket_obstruction_20261023/README.md).
Its full first-round tags separate 6878 of the 6880 order-eight classes;
one more round separates all. The canonical hostile is its two pairs of
colliding actual tournaments, not an abstract array that may be unrealizable.
The corrected near miss is that equal target signatures need not leave a
transport ambiguity, since capacities may force their placement. Conversely,
class separation alone is not a replay of the source/path/target ports.
The least-used coordinate here is **which vertex owns two outgoing edges**.

Our anchor is the smallest useful finite signature; the niche is a cheap
rooted motif formula; the wildcard is the old skew-Hadamard doubling tower.
The concept board is: motif incidence, ownership versus marginal counts,
adaptive measurements, port transport, uniform equitable cells, and
individualized neighborhoods. Each construction below retains an intrinsic
orientation, with no artificial ordering of tied vertex tags.

The all-order lossless fallback remains the root-curl and root-score
representation in [the port-rigidity packet](../metagraph_port_rigidity_20261023/PROOF.md).
This packet asks how much less information suffices in small orders, and
where that compression must fail. No universal classification or new
A000568 enumeration complexity bound follows.

## 1. Rooted strong-four incidence in cubic time

Let c_e count cyclic triangles through an unordered edge e, c_v the number
through vertex v, and b_v the number of induced strong four-vertex
tournaments containing v. Then, for every tournament,

    b_v = sum_{e containing v} binom(c_e,2)
          + sum_{cyclic triangles {v,a,b}} (c_{ab}-1).                 (1)

Proof. Two cyclic triangles with a common edge induce the unique strong
four-vertex tournament, and each such four-set has exactly those two
cyclic triangles. Mark v in this bijection. If v is on the common edge,
it is counted by the first sum. If v is one of the two outer vertices,
its triangle is fixed and the other triangle can be chosen in c_ab-1
ways. Each rooted four-set appears once, in exactly one case. This also
gives sum_v b_v = 4 M2, where M2 = sum_e binom(c_e,2).

Two passes over cyclic triples compute all b_v in O(n^3) arithmetic
operations, instead of scanning all four-subsets. The consumer independently
compares (1) with literal four-subset enumeration on every inherited
representative and both larger hostile examples. This is a formula-level
cost reduction; no end-to-end enumeration speedup has been benchmarked.

## 2. The exact finite menu

Choose any subset of the three seed features

    bit 1: s_v = outdegree(v),
    bit 2: c_v = cyclic triangles through v,
    bit 4: b_v = strong four-sets through v.

Let rho_0(v) be the selected coordinate tuple, possibly empty. Refine by

    rho_{r+1}(v) = (rho_r(v), sorted multiset of rho_r(w) over v -> w).

The graph record is the sorted multiset of its vertex records. Independently
append any subset of the scalar moments

    bit 1: M2 = sum_e binom(c_e,2),
    bit 2: M3 = sum_e binom(c_e,3),
    bit 4: M4 = sum_e binom(c_e,4).

The search uses all 8 seed subsets, all r in {0,1,2}, and all 8 moment
subsets: exactly 192 menus per order. A menu (f,r,m) is weaker than
(f',r',m') when f is a subset of f', r <= r', and m is a subset of m'.
Refinement retains its input, so this is a genuine information order.
An entry is minimal when it separates every class and no strictly weaker
entry in this menu does. Each failed entry has an explicit first collision
and a digest of its full collision partition in the receipt.

| Order | Ordinary classes | All minimal menus (feature bits, rounds, moment bits) |
|---:|---:|---|
| 4 | 4 | (0,1,0), (1,0,0) |
| 5 | 12 | (0,1,3), (0,2,0), (1,0,3), (1,1,0), (2,1,0), (3,0,0), (4,1,0), (5,0,0) |
| 6 | 56 | (0,2,1), (1,1,1), (2,1,1), (3,0,1), (4,1,0), (5,0,0) |
| 7 | 456 | (2,2,1), (4,1,0) |
| 8 | 6880 | (2,2,1), (4,2,0) |

In particular, **only b_v and one neighbor round** suffice through order
seven. At eight, **only b_v and two rounds** suffice. The alternative
minimal menu at eight is c_v, two rounds, and M2. Degree and cyclic-triangle
counts can therefore be removed from the first route's initial tags.

For b_v alone the zero-round/one-round class counts are respectively
24/56 at order six, 188/456 at seven, and 2793/6878 at eight. This is
class separation, not a theorem that these invariants characterize arbitrary
integer arrays or that all larger tournaments remain distinguishable.

## 3. An adaptive one-bit replacement for the second round

Write the first-round b-only record in its simple form

    R(v) = (b_v, sorted multiset of b_w over v -> w).

Only two ordinary order-eight collision pairs remain. Their canonical
upper-triangle integer codes in the inherited n8 TSV are

    (27467784, 27467793),     (6427781, 6427778).

Converse exchanges the two pairs; the members of either pair are not
converses of one another. The inherited Hamiltonian counts are 523 and
527. Each pair has not only the same multiset of R(v), but the same
**complete matrix of arc counts between R-record cells**. Thus remembering
additional marginal edge totals between these cells still loses the
decisive information.

For record types R,B,C put

    d_B(v) = number of w with v -> w and R(w)=B,
    W(R;B,C) = sum_{v with R(v)=R} d_B(v) d_C(v).                    (2)

The first pair is separated by

    R = (20,(16,18,20)),
    B = (16,(16,16,20)),
    C = (20,(14,16,20)).

The second pair is separated by

    R = (20,(14,16,16,20)),
    B = (16,(14,16,18,20)),
    C = R.

In the displayed pair order the respective measurements are 0 and 1.
The repeated R-cell has two vertices. The first moments say how many
edges leave it toward B and C; (2) says whether the **same member** owns
both incidences. This ownership correlation survives neither cell totals
nor separate one-neighbor marginals. It is an invariant under relabelling.

The public `adaptive_wedge_signature(A)` returns the first-round record
multiset and one bit. It measures the sum of these two selectors only if
n=8 and the record multiset equals one of the two exact guarded multisets
in `GUARD_KEYS`; elsewhere it returns bit zero. The receipt verifies
6880 distinct resulting signatures and checks both identical arc-count
matrices. The extra measurement is used on only four classes. One bit
is optimal **inside each fixed two-candidate fibre**, since zero extra
information cannot distinguish its two candidates. There is no unrestricted
bit-compression optimality claim.

The API's domain is an actual tournament 0/1 adjacency matrix, with zero
diagonal. It is an invariant evaluator, not an acceptance verifier for
untrusted arbitrary matrices. Its guard is representation-independent;
it does not select a hidden vertex labelling or compare canonical codes.
The codes above name reproducible witnesses only.

### From class signatures to the metagraph

For every literal source/path port, its path-complement target is an actual
tournament. A complete ordinary-class signature gives its unique ordinary
target class. One can then merge converse classes, retaining the port's
closed/open color, endpoints, and multiplicity. This logical implication
does not replace the actual enumeration of ports.

The companion [order-eight transport packet](../metagraph_order8_20261023/README.md)
performs that enumeration and independently tests the public API. Its
weaker one-round signature genuinely leaves alternative transportation
tables; adding our bit or the full second round fixes the target choices.
The present receipt proves the class-separation ingredient only. A
converse-merged class name can be formed from the unordered pair of ordinary
signatures of T and its converse; merging that early must not erase the
orientation gauge needed by the source/path port.

## 4. Why indefinitely more rounds and moments cannot finish the job

Let T be a doubly regular tournament on n=4t+3 vertices: every vertex has
outdegree 2t+1, and any two vertices have t common outneighbors. Then

    c_e = t+1 for every edge,
    c_v = (2t+1)(t+1) for every vertex,
    b_v = (n-1)t(t+1) for every vertex.                              (3)

For an edge u->v, the two vertices have t common outneighbors and t
common inneighbors. The remaining intermediate types have t vertices
on u->w->v and t+1 on v->w->u. The latter count is c_uv. Counting at v
gives 2c_v=(n-1)(t+1). Formula (1) then gives

    b_v = (n-1) binom(t+1,2) + c_v t = (n-1)t(t+1).

Thus all seed tags are identical and depend only on n. A refinement round
gives each vertex the same old tag and 2t+1 copies of that tag. Induction
proves that **every finite number of rounds** agrees between any two DRTs
of the same order. Including both in- and outneighbor multisets does not
help. Moreover every edge-incidence moment of every order is fixed:

    M_k = binom(n,2) binom(t+1,k).

This rules out the whole finite menu even after allowing arbitrarily many
rounds and arbitrarily high edge-incidence moments.

There is a further boundary. Start ordered pairs with the three colors
`equal`, `arc`, `reverse arc`, and refine a pair (u,v) by the multiset of
color pairs `(color(u,w),color(w,v))` over w. The four intermediate counts
just displayed, together with the two diagonal cases, depend only on the
original relation and n. The diagonal profile also depends only on n.
Consequently no relation splits at the first round, and induction makes
the rank-three partition stable forever. This is a statement about this
precisely defined, unindividualized pair refinement; it does not say that
all possible edge invariants or rooted constructions fail.

### An explicit realized pair at order 31

Take Paley P31, with i->j when j-i is a nonzero square modulo 31, and
D(D(P7)), where D is the skew-Hadamard doubling from THM-4557,
[`doubling-a-doubly-regular-tournament-keeps-its-automorphism-group`](../../../01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md).
This is an inherited construction, not a new DRT classification theorem.
The consumer reconstructs both adjacency matrices and verifies:

| Quantity | P31 | D(D(P7)) |
|---|---:|---:|
| Vertex tag (s_v,c_v,b_v) | (15,120,1680) | (15,120,1680) |
| Every edge's c_e | 8 | 8 |
| (M2,M3,M4) | (13020,26040,32550) | (13020,26040,32550) |
| Vertices whose outneighborhood induces a DRT | 0 | 9 |

The last line is an isomorphism invariant, proving these are different
ordinary classes. Negation modulo 31 is an explicitly checked converse
isomorphism of P31, so the two examples also belong to different
converse-merged nodes. Both have the same stable rank-three pair profiles.
The finite check establishes a realized witness; (3) and the induction
establish its obstruction at every refinement depth. No order-minimality
claim is made, and no 31-vertex Hamiltonian-path enumeration is required.

The simultaneous [order-thirteen packet](../metagraph_allsizes_20261023/README.md)
already gives all-depth vertex-refinement failures, repaired on its pair
by M3=234 versus 221. Our DRT pair explains why that repair cannot be
universal: all M_k agree. Conversely, its source-prefix extension proves
vertex-refinement failure at every n>=13; this packet does not improve
that smallest-order or every-order statement.

### A single five-vertex motif repairs this realized obstruction

The outneighborhood test can be reduced to a natural small count. Suppose
T is doubly regular of order n=8h+7, and put U=N^+(v). Then U has 4h+3
vertices and is regular of degree 2h+1: this is the common-outneighbor
count of v and each u in U. For a pair u,w in U let k_uw be their number
of common outneighbors inside U. Define

    B_v = sum_{pairs u,w in U} binom(k_uw,2),
    J_v = sum_{pairs u,w in U} (k_uw-h)^2.

Count the sum of k_uw by their common target. Every target has indegree
2h+1 in U, so the mean of k_uw over the binom(4h+3,2) pairs is exactly h.
Expanding the square gives the all-order identity in this family

    J_v = 2 B_v - h(h-1) binom(4h+3,2).                             (4)

Crucially, **B_v is the number of transitive five-subsets having v as
source**. A term in B_v chooses two vertices in U and two common targets
in U, with all four top-to-bottom arcs present. The arcs inside the two
pairs order them, making this a transitive four-set in U. Conversely each
transitive four-set has exactly one top pair and bottom pair. Adjoining v
gives a transitive five-set with its unique source marked. Thus the
unrooted transitive-five count is sum_v B_v; no automorphism division is
needed.

The variance J_v vanishes exactly when all k_uw equal h, which is exactly
double regularity of U. It is nonnegative for every v. Hence (4) also gives
the general bound

    #transitive five-sets >= n h(h-1) binom(4h+3,2) / 2,

with equality exactly when every outneighborhood is doubly regular. This
is an elementary consequence of the retained local variance, with no
priority claim.

At order 31, h=3 and J_v=2(B_v-315). The independent exact results are:

| Tournament | J_v histogram | B_v histogram | Total transitive five-sets |
|---|---|---|---:|
| P31 | 60 repeated 31 times | 345 repeated 31 times | 10695 |
| D(D(P7)) | 0^9, 36^7, 48^7, 84^8 | 315^9, 333^7, 339^7, 357^8 | 10395 |

Here x^m denotes multiplicity m, not an exponentiation. A **single global
transitive-five count** therefore distinguishes these two examples even
though every triangle-edge moment and every stated unindividualized
refinement round agree. This repairs this pair, not all DRT pairs.

This is also the first possible size of an **unrooted induced motif**
separating this pair. In any DRT of order 4t+3, all motif counts through
size four depend only on t. Vertex and edge counts are immediate; cyclic
triple counts follow from regularity. Among the four four-vertex classes,
the strong count is binom(n,2) binom(t+1,2). A source plus a cyclic triple
is counted by summing the cyclic-triple count in every outneighborhood,
which is regular of order 2t+1 and degree t. Its total is

    n [binom(2t+1,3) - (2t+1) binom(t,2)].

The converse class has the same formula using inneighborhoods. The
remaining class is transitive and is obtained by subtraction from
binom(n,4). At order 31 the four respective counts are

    strong: 13020; source+C3: 4340; C3+sink: 4340; transitive: 9765.

The sidecar checks all four counts by literal four-subset score sequences.
Thus no unrooted induced motif of size at most four can separate these
two DRTs, while the transitive-five motif does. This restricted minimality
does not assert minimality among rooted motif records or arbitrary
invariants.

`neighborhood_defect.py` independently reconstructs the tournaments and
checks four descriptions: direct squared deviations, common-target
overlap counts, literal transitive four-subsets in every outneighborhood,
and literal five-subsets of the whole tournament. Paley P7 is the positive
zero-defect control. Its default and optimized modes reproduce
`neighborhood_defect.json`; it neither imports nor changes the frozen
adaptive-signature API.

## 5. What this changes in the proof strategy

The useful map is from an oriented tournament to rooted motif records,
then to neighbor multisets. It preserves class invariants, but destroys
the joint ownership of several incidences at one vertex. At order eight
one adaptive mixed moment repairs exactly that loss. In doubly regular
families even unindividualized pair refinement remains homogeneous; the
successful sidecar here retains a transitive five-set whose four nonroot
vertices all share one source. The two repairs act at different scales,
but both retain a common owner that marginalization had discarded.

The next precise question is whether a small family of **rooted five-vertex
incidence records**, refined only in still-ambiguous fibres, separates the
relevant port-target classes before an entire canonical form is required.
Test it against both the n13 pair and the n31 DRT pair, retaining
source/path/target information. The uniform
cell obstruction must be an input to this design, not a reason simply
to add more rounds of the same operation. No claim toward universal
Collatz coverage or an arbitrary-size tournament count follows from this
finite signature calculation.

## Reproduction and universe

```text
python -B 05-knowledge/results/metagraph_minimal_signatures_20261023/check.py
python -B -O 05-knowledge/results/metagraph_minimal_signatures_20261023/check.py
```

Default execution compares the frozen receipt; `--write` intentionally
regenerates it. All checks use explicit exceptions and remain active under
`-O`. Standard library only. The finite universe is the complete inherited
ordinary class lists with sizes 4,12,56,456,6880. The consumer checks their
hashes against the prior audited bucket receipt, and evaluates all 192
menus independently. Completeness of those lists is inherited from that
packet's exhaustive labelled/orbit controls; this consumer does not
repeat its canonical enumeration. Each motif value has an independent
four-subset check, each failed menu retains a concrete collision, the
adaptive repair retains both literal adjacency matrices and identical
cell-arc matrices, and the larger hostile matrices are reconstructed
without trusting a stored class name. No random sampling is used.

The inherited fixture hashes deliberately use the exact original bytes,
as does the earlier receipt they authenticate; a line-ending-only rewrite
of an inherited fixture is rejected and would require an explicit receipt
normalization change. This packet does not silently weaken those inherited
hash checks. Its own Python source hash is LF-normalized.

Normal and optimized replays both pass 28,622 checks, including the frozen
receipt comparison. Source LF SHA256:
`623fd524e2ec5d4b2d6707fefc88f15c7ce50f1dfc2dd2fa36c8041c6217479d`.
Receipt LF SHA256:
`f2586741801dd07b5cfdc5d6143827ac96f67d8f67d858c87d31e6090ebce1fd`.
The verifier compares canonical JSON encodings, since tuple-valued exact
records become arrays when serialized. Its first default replay exposed
that representation mismatch; the comparison was repaired without changing
any mathematical calculation or the public adaptive API.
