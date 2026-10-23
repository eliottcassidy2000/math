# An all-depth obstruction to vertex-incidence refinement

October 23, 2026. **PROVED:** the stated refinement collision at every
order n >= 13 by an explicit context extension. **FINITE-EXACT:** two
strong circulants at order 13 and the bounded checks below. **OPEN:** the
smallest order of failure among all tournaments and a generally complete
small incidence signature. This is not a same-fibre metagraph counterfeit.

## Inheritance and the precise target

The prior [incidence packet](../metagraph_bucket_obstruction_20261023/README.md)
separates every ordinary class through order 8 by at most two rounds of
vertex-neighbor refinement seeded with

    tag(v) = (outdegree(v), cyclic triangles through v, strong four-sets through v).

The finite success does not imply that more rounds always suffice. The
closest hostile is its pair of first-round collisions at order 8; the
corrected near miss is the [Session37 functional-versus-affine recovery
boundary](../session37_inheritance_audit_20261023/REPORT.md). A kernel in an
abstract linear coordinate space is not enough: this packet supplies two
actual tournaments with different Hamiltonian counts. The least-used
coordinate is the third moment of edge-triangle incidence.

Board: uniform vertex tags; actual realized pairs; refinement depth;
ordinal contexts; edge-overlap moments; literal path counts. The anchor is
the all-order boundary of the proposed signature, the niche is its cheapest
repair on this pair, and the wildcard is the same distinction between
diagonal and bilinear information in the dyadic electrical packet.

## Two explicit strong sources

On Z/13, put u -> v when v-u belongs to S. Take

    S_A = {2,7,8,9,10,12},
    S_B = {1,2,4,7,8,10}.

Each set selects exactly one of d and -d. Translation acts transitively,
and the step-2 cycle visits all vertices, so both tournaments are strong.
Negation is an isomorphism to converse, so both are self-converse.

| Quantity | A | B |
|---|---:|---:|
| Uniform vertex tag | (6,21,112) | (6,21,112) |
| Cyclic triples | 91 | 91 |
| Induced strong four-sets | 364 | 364 |
| Hamiltonian paths H | 3683797 | 3704857 |
| Automorphism order | 13 | 13 |
| Ordinary tilings H/Aut | 283369 | 284989 |

The H inequality proves nonisomorphism. The automorphism computation
exhausts every permutation fixing vertex 0 by adjacency-constrained
backtracking; both rooted stabilizers are trivial. Translation transitivity
then gives order 13 for the whole group.

**All-depth failure theorem.** Starting from these tags, every finite round
of vertex refinement that retains the old tag and the multiset of tags of
outneighbors produces identical signatures for A and B. This remains true
if the multiset of inneighbor tags is included as well.

Proof: initially all 13 vertices in both objects have one common tag. Each
vertex has exactly six outneighbors and six inneighbors. Thus its next tag
is the same function of one old tag and six copies of that tag, in both
objects. Induction proves the claim at every depth. No limit in the number
of rounds repairs this collision. This proof applies to any deterministic
injective serialization of the stated refinement records.

This does **not** show that every bounded-size invariant fails, or that
these classes have the same already-known node capacities. Their H and
tiling counts differ. It is a failure of the proposed transported invariant
to identify the target, not a counterfeit preserving all source margins.

## A collision at every larger size

Let E_k(T) be the ordinal sum of a transitive k-vertex tournament followed
by T: every new vertex beats every old one. The old core tags remain
(6,21,112), because no cyclic triple or strong four-set can cross a
uniformly directed component cut. The successive new vertices have tags

    (13+k-1-i, 0, 0),    i=0,...,k-1.

Use the matching partition consisting of those k singleton positions and
one cell containing the 13 core vertices. Its tag records match between
E_k(A) and E_k(B), and every vertex has the same neighbor-counts into every
matched cell. The partition therefore stays matched at every refinement
round. Within the core no split ever occurs.

Every Hamiltonian path must traverse the transitive prefix in its unique
order before entering the core. Consequently H(E_k(T))=H(T), so the two
extensions remain nonisomorphic. This proves a realized failure at every
n=13+k. The base pair is strong; the k>0 extensions are not claimed strong.
More generally the same reasoning preserves indistinguishability under
common ordinal contexts whenever their cells and incidence tags are
retained. No claim about arbitrary substitution operations is required.

## One more edge-overlap coordinate repairs this pair

For an oriented edge u -> v, record

    (common outneighbors, common inneighbors, cyclic triangles through uv).

The complete multisets are

    A: (1,1,5)^13, (2,2,4)^13, (3,3,3)^52;
    B: (2,2,4)^52, (3,3,3)^13, (4,4,2)^13.

Every vertex has the corresponding outgoing-edge multiset with all
multiplicities divided by 13. Thus even a local edge-incidence sidecar
distinguishes the pair without enumerating Hamiltonian paths.

Let c_e count cyclic triangles through an unordered edge. Their first two
binomial moments agree:

    sum_e c_e = 273,
    sum_e binom(c_e,2) = 364.

The second is exactly the strong-four count already inherited by the seed.
But the next moment differs:

    sum_e binom(c_e,3) = 234 for A, 221 for B.

This counts three distinct cyclic triangles sharing their marked edge. It
retains a five-vertex overlap that repeated averaging of uniform vertex
tags cannot recover. Adding the third moment repairs this particular pair;
it is not asserted to complete the invariant on all tournaments. These
moments also survive unchanged under source-prefix extension, since all
new edges belong to zero cyclic triples.

## Exact verification and scope

The consumer checks every vertex tag by literal triple/four-set enumeration,
cross-checks strong-four counts by the independent edge-overlap identity,
and computes H in two directions: forward all-start/end subset DP and a
backward remaining-set recursion from vertex 0, multiplied by the proved
translation orbit size. It also includes direct-permutation controls and
an exhaustive root-stabilizer search. Both the old colors and all in/out
neighbor colors are retained while checking stable partitions on the
extensions k=0,...,8. The theorem for all k is the argument above, not an
extrapolation of those finite controls.

`scout.py` is separate bounded discovery. It tested the 16 circulants at
order 9 and 32 at order 11 without finding the requested tag/H collision,
then found this order-13 pair. Those earlier searches are not exclusions
of arbitrary tournaments, nor a proof of minimality.

```text
python -B 05-knowledge/results/metagraph_allsizes_20261023/check.py
python -B -O 05-knowledge/results/metagraph_allsizes_20261023/check.py
```

Both modes reproduce the same receipt with 297 checks including its
comparison. Receipt LF SHA256:
`5a3413d4159279dbeeee833f3e76bc57fa12905aa48d153322976db5c418a1e1`.
