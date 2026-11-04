# One-probe invisibility of balanced modules, and the circulation left by a collision

**Status: PROVED** by the finite combinatorial arguments below; **FINITE-EXACT**
controls are specified separately. No historical novelty, all-order two-port
reconstruction, Hamilton-spectrum completeness, or complexity lower bound is claimed.

## Inheritance and the new step

The closest mechanisms are Session176's positive directed-cut insertion formula
([circuits/PROOF.md](../../session176_tournament_information_20261104/circuits/PROOF.md))
and its exact uniform-context path-cover interface
([information/PROOF.md](../../session176_tournament_information_20261104/information/PROOF.md)).
The canonical hostile is the pair of opposite labelled triangles: identical
complete one-insertion responses, separated by the same second insertion.

Recovery matters here. Session27 already proves full labelled reconstruction
using n-1 external vertices, even from the second response bit
([MODULAR_FAITHFULNESS.md](../../session27_repunits_lattices_20261022/session27-response/MODULAR_FAITHFULNESS.md)).
Its [response proof, sections 6-8](../../session27_repunits_lattices_20261022/session27-response/proof.md)
already gives global vertex-transitive converse collisions and two-port
tomography of endpoint and ordered two-path-cover data. Its universal two-port
injectivity question remains open. Session176's ordered sums of reversed
triangles are therefore an extension of a recovered information-loss mechanism,
not the first such mechanism.

The new step below allows a balanced module to be reversed **inside any ambient
tournament quotient**, while preserving every one-vertex insertion count.
The useful sidecar is endpoint balance separately at every path-cover run count.
It is stronger than regular degree. A second result describes every possible
one-insertion collision as a positive circulation on its flipped arcs.

Live board: uniform modules; ordered path covers; endpoint balance;
positive cut capacities; successive nonuniform probes. The first three make
the invisibility proof work; the fourth constrains all collisions; the fifth
explains why this quotient is still not a sufficient state for future work.

## 1. Definitions

All vertices are named, and every path is a directed vertex word. For a
nonempty tournament M on r vertices, let q_M(k) count ordered lists of k
nonempty vertex-disjoint paths partitioning V(M), and put q_M(k)=0 for k>r.
For a fixed slot i among those k paths define

    B_k(u) = number of these covers whose slot i starts at u,
    E_k(u) = number of these covers whose slot i ends at u.

The counts do not depend on i, since permuting the ordered paths is a
bijection. Call M **endpoint-balanced at every run count** if

    B_k(u)=E_k(u) for every u and 1<=k<=r.                (B)

For two distinct slots i,j define D_k(u,v) to count covers whose slot i
ends at u and whose slot j starts at v. These slots are distinguished;
D_k(u,u)=0 because distinct paths are vertex-disjoint. For k>=2,

    sum_v D_k(u,v)=E_k(u),
    sum_u D_k(u,v)=B_k(v).                              (1)

Each cover counted by one side contributes exactly once to the other.
Slot-permutation symmetry means D_k is independent of the particular ordered
pair i,j. Reversing every path, then swapping these two distinguished slots,
shows that M^op has matrix D_k^t. It also exchanges B_k and E_k and preserves q.

## 2. The contextual module-reversal theorem

**Theorem C1.** Let Q be any finite tournament with a distinguished vertex z,
and replace z by a nonempty tournament M satisfying (B). All edges between
each vertex outside M and M are uniform, inherited from Q. Let T=Q[M] and
T'=Q[M^op], keeping all names. For every bit vector b on their common vertex
set, add one vertex x with x->u iff b_u=1. Then

    H(T+x_b)=H(T'+x_b).                                (2)

No restriction is imposed on Q, the orientations among its other vertices,
or any of the probe bits. H(T)=H(T') as well.

**Proof.** Fix a permutation of all vertices outside M, including x, and
specify which gaps of that word are occupied by nonempty M-runs. This is
an exterior skeleton. Every Hamilton word has exactly one such skeleton
and one ordered k-path cover of M assigned to its occupied gaps. A gap
whose endpoint is an old outside vertex has a legality condition independent
of the individual vertex of M, by the module hypothesis. Unoccupied gaps
also have legality conditions independent of the orientation inside M.
Discard skeletons failing any of these conditions.

There are only three possible kinds of restrictions involving x:

* No adjacent M-run: the contribution is q_M(k).
* One adjacent M-run: the contribution is either
  sum_u (1-b_u)E_k(u) or sum_u b_u B_k(u).
* M-runs immediately before and after x: their slots are distinct, and the
  contribution is sum_{u,v}(1-b_u)b_v D_k(u,v).

The first case is unchanged by converse. The second is unchanged by (B).
For the third let J=D_k-D_k^t. Equation (1) and (B) say every row and
column sum of J is zero. Since J is skew-symmetric,

    sum_{u,v}(1-b_u)b_v J(u,v)
      = sum_v b_v sum_u J(u,v) - sum_{u,v} b_u b_v J(u,v)
      = 0 - 0.                                        (3)

Thus each skeleton contributes equally. Summing proves (2). With no x,
only unrestricted path-cover counts occur, proving H(T)=H(T'). Endpoints,
empty exterior prefixes/suffixes, and paths that become orphans after deleting
x are all included in the skeleton bijection. No insertion-only assumption
is being made. QED.

**Corollary C2.** Every vertex-transitive M satisfies (B), so C1 applies.
Indeed automorphisms act transitively on u and preserve all fixed-slot
counts. Summing B_k or E_k over u gives q_M(k), hence each equals q_M(k)/r.

**Corollary C3.** Suppose a tournament has s pairwise disjoint modules,
each endpoint-balanced at every run count. Independently reversing any
selection of them leaves the entire one-insertion response unchanged.
Reverse them one at a time and apply C1: each remaining module is still a
module with the same internal orientation until its turn. If each chosen
module has at least two vertices, its converse changes its named arcs, so
the 2^s choices are distinct labelled bases. This does not claim 2^s
nonisomorphic bases; directed triangles, for example, are self-converse up
to relabelling.

The theorem works for arbitrary fixed ambient quotients, including strong
ones, and for different balanced modules at different holes. It subsumes
ordered sums of triangles but does not depend on an ordered-sum factorization.

## 3. Why ordinary degree balance is insufficient

Condition (B) implies regularity: for k=r-1, the unique nonsingleton path
has two vertices, and

    B_(r-1)(u)-E_(r-1)(u)
      = (r-2)! (d_out(u)-d_in(u)).                    (4)

When the distinguished slot is a singleton u, its contributions agree.
When it is a two-vertex path, the starting/ending possibilities are exactly
the out-/in-neighbors of u, with (r-2)! orders for the remaining singletons.
For r=1 the balance and regularity statements are immediate.

The converse fails. The regular seven-vertex tournament with outgoing masks

    M=(22,88,50,69,104,11,37)

has all outdegrees three, H(M)=171, but B_1(0)=24 and E_1(0)=25.
Add x with x->0 and every other old vertex pointing to x. The counts for
M and M^op are respectively **403 and 407**. Thus degree balance cannot
replace the run-count endpoint condition, even with no ambient vertices.
This is an exact witness, not a minimal-order classification.

The module hypothesis is also load-bearing. Start with opposite directed
triangles and adjoin the same nonuniform vertex with outbits (0,0,1).
Both resulting four-vertex tournaments have H=5, but a further vertex
with outbits (0,1,0,0) gives 13 versus 11 (the inherited Session176 witness).
The original triangle is not a module after the first attachment.

## 4. A strong ambient collision and its two-step separation

Take Q to be a directed triangle, with one vertex inflated to C3. The
two choices of orientation inside that module give strong five-vertex
tournaments with outgoing masks

    A=(10,12,9,16,7),  A'=(12,9,10,16,7).

Both have H=15. They have identical full one-insertion responses by C1.
Add x with outmask 1, then y with outmask 2 on vertices (0,1,2,3,4,x).
The intermediate counts are both 31; the final counts are **81 and 87**.
All first-stage and second-stage attachments use the same names and bits
in both constructions. This separates a contextual one-probe collision
after the second probe, even when the original ambient tournament is strong.
The checker verifies all 32 first-stage assignments and this displayed
two-step witness; it does not claim a smallest strong witness.

## 5. Every one-insertion collision leaves a positive circulation

Let F_T(b)=H(T+x_b). Session176 proves

    F_T(b)=H(T)+2 sum_(u->v in T) w_uv b_u(1-b_v),

where each symmetric w_uv is a positive integer, determined by F through
mixed second differences. The same polynomial determines each weighted
outdegree d_w(u)=sum_(u->v) w_uv.

**Theorem C4.** If F_T=F_T' on the same named vertex set, orient the set
D of differing edges as in T and give each weight w_uv. At every vertex,

    sum_(u->v in D) w_uv = sum_(v->u in D) w_uv.        (5)

Consequently every differing edge belongs to a directed cycle in D.
In particular no differing edge joins distinct strongly connected components
of T. The same holds with T' in place of T, and their ordered SCC
decompositions agree.

**Proof.** Equality of response polynomials gives equality of H, all
weights, and all weighted outdegrees. Subtract the two weighted outdegree
equations at u; unchanged edges cancel, giving (5). Replace each edge
by w_uv parallel copies. The resulting finite directed multigraph is
balanced. Its edges decompose into directed cycles by repeatedly following
outgoing unused edges until a vertex repeats and removing the resulting
cycle. Hence every differing arc lies in a directed cycle of T, so stays
within a strong component of T. Repeat with T'. Since no cross-component
edges change, each strong component of either tournament is a union of
components of the other; the two partitions and their orders therefore
agree. Equivalently, Session176's zero-response suffix sets recover this
common ordered partition directly. QED.

At least three named arcs must change in a nontrivial collision, and opposite
triangles attain three. A three-arc collision reverses a directed triangle
whose three common weights are equal. These are necessary conditions;
the weights themselves are endogenous to the base tournament.

For example, start with A=(24,1,3,6,14) and reverse the directed triangle
on labels {0,1,4}, obtaining B=(10,16,3,6,13). All three old triangle
weights are four, so the proposed flip is a balanced circulation in the
old weighted graph. Nevertheless H(A)=11 and H(B)=15, and the weight
matrix changes. Circulation using the old weights alone is therefore
not sufficient to preserve the response.

## 6. Exact maps, losses and stopping boundary

The source of C1 is a module with its full ordered path-cover endpoint
profiles; the operation is internal converse inside a fixed quotient;
the target predicate is equality of every one-insertion Hamilton count.
The map forgets orientation circulation that a single probe sees only
through a directed cut. The retained sidecar is B_k=E_k at every k; it
makes the skew cut contribution vanish. Vertex-transitivity supplies that
sidecar but is not declared necessary.

C4 goes in the opposite direction: from an observed collision to a
positive balanced flow on the changed named arcs. It loses the demand that
the new tournament have the same endogenous weights, so it is a necessary
structural obstruction, not a synthesis equivalence.

The cheapest hostile tests are respectively the regular order-seven witness,
the nonmodule 13/11 witness, and the balanced old-weight triangle flip 11/15.
The strong 81/87 witness marks the stopping boundary: contextual one-probe
invisibility does not yield a continuation-closed state. General two-port
reconstruction remains OPEN. The known n-1-port reconstruction is retained,
not overwritten by this narrower theorem.

## 7. Reproduction

Run `python3 check.py` and `python3 -O check.py` from this directory. They
print identical JSON receipts. No packages, approximate arithmetic, random
filters, or unchecked assertions are used. `scout.py` supplies the exact
subset-DP primitives and preserves the initial bounded discovery experiment;
`check.py` states and checks the full committed universe.

Finite evidence corroborates the proofs and is not their all-order basis.
The receipt records all universes, totals, named hostile outputs, and the
SHA256 hashes of both executable source files. Independent checks use literal
Hamiltonian-word enumeration on the smaller controls.
