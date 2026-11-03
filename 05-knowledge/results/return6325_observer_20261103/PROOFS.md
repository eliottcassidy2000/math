# Observer completion, response filtration and structural halving

**PROVED scoped elementary mechanisms / FINITE-EXACT coefficients and atlases.**
Same-agent proof review; independent algorithms are identified in AUDIT.
Local labels O1--O6. No priority claim for standard observability or set cover.

## O1 Completion is relative to a target and a legal query menu

Let X be a declared source universe, b its baseline observation, and t its
desired target distinction. A context is a specified legal program followed
by an observation; include a failure symbol if legality depends on the source.
For a menu C, define kappa_C(b;t) as the least number of contexts whose
joint responses, together with b, determine t on X. Infinity is allowed.
This is nonadaptive, unless explicitly stated otherwise. Contexts are queried
on the same original source: the newborn read heads below are deleted after
their responses are recorded, so each query has one birth and one deletion.

For finite X, form pairs (x,y) with b(x)=b(y) and t(x)!=t(y). Each context
c covers the pairs with different c-responses. The completion number is
exactly the minimum set cover of these pairs. This proves the finite search
used here, including its infinity outcome. Distinguishing one fixed pair
needs either zero, one, or infinitely many contexts; the richer invariant
concerns recovery of a partition or a family, not one already chosen pair.

For a linear source span V, scalar queries are columns in V*. If W is the
baseline query span and W_infinity the span of every legal future scalar
query, then exactly dim(W_infinity)-dim(W) added scalar contexts suffice to
recover all future responses **linearly on V**. A basis extension by actual
query columns attains the lower bound. This does not count words with a
multi-output readout as single scalar queries, nor imply a nonlinear lower
bound on a smaller physical or integral subset. A particular target subspace
can require a different sparse selection; dimension alone need not solve it.

## O2 Response-holonomy filtration needs both kernels and source partitions

For total linear legal operators A_l and initial query columns Q, set

    W_d = span{A_w q : |w|<=d, q in Q},
    K_d = {v : v W_d=0}.

Then W_(d+1)=W_d+sum_l A_l W_d. Once equality occurs, closure holds for
every future word. Before stabilization each step increases dimension, so
depth is bounded by dim(V)-dim(W_0). The annihilators K_d descend.
In the unital operation algebra A define

    E_d = {D in A : D W_d=0}.

These are descending left ideals: (aD)W_d=0. If D is the difference of
two legal word operators, membership in E_d says their discrepancy is
invisible on every source through depth d. The algebraic units satisfying
(g-I)W_d=0 form a descending subgroup filtration. Arbitrary units of A are
not thereby legal attachment words; the inherited singular-word obstruction
continues to apply.

The actual persistent-spare language is the seven macros 3p,5p,N,J0,...,J3
from [JE4](../joint_port_geometry_20261103/JOINT_EFFECTS.md), on the rational
27-dimensional physical source span and exact 95-dimensional algebra
[A1--A4](../return701_algebra_20261103/ALGEBRA95.md). Exact rational closure gives:

|Initial scalar readouts|Depth|Response rank|Invisible source dimension|dim E_d|
|---|---:|---:|---:|---:|
|H|0|1|26|85|
|H|1|4|23|62|
|H|2 and all later|10|17|51|
|H, s, E1, S2|0|4|23|58|
|H, s, E1, S2|1|14|13|12|
|H, s, E1, S2|2 and all later|27|0|0|

Here s counts starts at the active root; E1 counts ends there; S2 counts
one/two-component covers with a path starting there, with the inherited
cap conventions. All matrices, query words and exact independence/closure
certificates are in `observer_filtration.json`. H-only has completion nine
for its own future responses and can never recover all 27 linear coordinates.
The four-readout menu needs exactly 23 additional scalar contexts, all with
at most two macros. This refines the old closed observer; it does not assert
closure under arbitrary native complement or recovery of the whole tournament.
An explicit decoder also recovers the capped observer of the actual saved
6325 source from the 27 basis-query responses. Independent elimination modulo
1009 confirms the rational rank lower certificates; rational closure supplies
the upper bounds. This reconstructs the declared cap, not the labelled source.

There is a separate finite, nonlinear definition. For a group G acting on X,
let H_C={g : f(gx)=f(x) for every x and every observed response f}. It is a
subgroup, by composition and inversion, and shrinks when contexts are added.
It need not be normal for an arbitrary menu. Crucially H_C is not a substitute
for the partition of X into equal-response sources. O3 gives an actual example
where the subgroup stops shrinking before the partition becomes discrete.

## O3 A complete one-birth observer for the hostile eight-source fibre

Use the literal 81/27 pair and the 54-element affine triple-symmetry group
of [G2](../digit_cycle_ternary_20261103/SYMMETRY_AND_PORTS.md). Its labelled
tournament orbits have sizes two and six. Their union X has eight actual
sources, each with H=462915 and the same 67 labelled cyclic triples.
The source and group permutations are retained in `completion_hostile.json`.
The two full automorphism orders are freshly recounted, not inferred from 54.

For a newborn w, let x_v=1 mean v->w. Its Hamilton count is the quadratic
multilinear response

    H(T+w)=H(T)+sum_v L_v x_v-sum_(v<u)(K_vu+K_uv)x_v x_u,
    K_vu=sum_(empty!=S!=V) end_v(T[S]) start_u(T[V-S]),
    L_v=end_v(T)-start_v(T)+sum_u K_vu.

Proof: a path either starts at w, ends at w, or splits into an old prefix
ending at v, then w, then a suffix starting at u. The latter contributes
K_vu x_v(1-x_u). The endpoint cases give the displayed constant and linear
terms. This retains the shared source partition and joint endpoint incidence.
The contemporaneous [S172 undirected one-vertex response theorem](https://github.com/eliottcassidy2000/math/blob/2791e01b10/05-knowledge/results/session172_memory_dynamics_20261103/squares/NOTE.md)
uses the same endpoint/two-path deletion mechanism. Here tournament incidences
are x_v and 1-x_u, giving the directed constant/negative-quadratic form.
Its undirected sharp-rank theorem is not transferred to this atlas.

All 4096 incoming cuts are evaluated exactly. For full reconstruction of X,
the optimal completion numbers by allowed incoming-set size are:

|Maximum incoming-set size|Completion number|Witness cuts|
|---:|---:|---|
|1|infinity|Four unresolved pairs remain|
|2|2|{0}, {5,8}|
|3|2|same|
|4 or unrestricted|1|{0,1,5,8}|

The last cut, mask 291, has responses, in the saved state order,

    2450493,2421909,2665485,2343867,2753631,2633253,2329119,2772027.

All are distinct. No menu confined to any fixed two old labels separates
all eight. The two-query witness uses three old labels across its two cuts.
Every selected newborn produces a strong tournament. The optimality statements
are finite exhaustive set-cover results for the declared atlas and menu.

The baseline / all-singleton / all-one-birth response partitions have sizes
1 / 4 / 8, while their invisible-group orders are **54 / 9 / 9**. Singleton
responses already give a faithful action of the effective quotient G/N, of
order six, while retaining four pairs of physically different sources.
Thus a faithful symmetry action need not be a faithful source observer.
The common order-nine kernel N is the affine reduction kernel at modulus three;
the 81/27 distinction and the orientation gauge are never discarded.

These group actions describe actual relation transformations, not free legal
commands. A specified transformation can be compiled by isolated arc repairs
with temporary helpers; none of that native cost is included in a read-head
query. Ports and old vertex names remain fixed during the birth queries.

## O4 Adjacent twins give a strictly monotone legal readout

Let u->v be a two-vertex module: their relations to every other vertex agree.
Write C for its contraction, with resulting vertex r. Let F_x count covers
by a distinguished nonempty path ending at x and a second possibly empty
path, and s_x count spanning paths starting at x. Then

    F_v(T)-F_u(T)=F_r(C),
    s_u(T)-s_v(T)=s_r(C).

For the first identity, swapping u and v maps the covers ending at u
bijectively to those ending at v without the arc u->v. The remaining covers
contract that arc bijectively to covers ending at r. For the second, swap
paths starting at v into paths starting at u; the omitted ones start u,v,
and contract to paths starting at r. The module condition preserves every
other arc. These are bijections, not count fits.

The legal newborn losing only to x has count R_x=H+F_x-s_x: paths starting
at the newborn give H-s_x; every other path has an old prefix ending at x,
giving F_x. Therefore

    R_v(T)-R_u(T)=F_r(C)+s_r(C)>0.                 (O4)

Strictness follows from the cover consisting of the singleton r and any
Hamilton path of C-r. In a transitive fibre, consecutive vertices are twins.
Their one-old-port readouts therefore strictly recover their positions.

## O5 Exact completion budgets after structural halving

Fix a labelled quotient Q, B named fibres, each of known size m=2^a, a>=1,
and consider all tournaments Q[TT_m,...,TT_m] obtained by varying the internal
orders of the old names. The query menu is **single-old-vertex incoming cuts**
as in O4; each query's integer response is read on the original source.
All unmarked invariants are identical across these internal reorderings.

If only quotient and unordered fibre membership survive, the exact completion
number for the full labelled tournament is **B(m-1)**. O4 makes each named
vertex's rank recoverable; m-1 ranks determine the last one. Conversely,
if two names in a fibre remain unqueried, interchanging them preserves every
queried rank but changes the tournament. The lower bound holds adaptively too.

If the entire unordered binary contraction hierarchy also survives, the
number is **Bm/2**. Query one leaf in every bottom pair. Its rank identifies
which side of every ancestor it occupied; collectively these paths visit
every internal node and recover all orientations. If one bottom pair has no
queried leaf, exchanging its leaves remains invisible. If all internal
orientations are retained, no additional queries are needed.

For the previous two-level normalization, m=4 and B=3^(k-1): the budgets are
**3B, 2B, 0**, depending on the retained data. Its matching records retain
the unordered hierarchy, while ordered original fibres restore everything.
These are operational query counts, not bit counts. No phase or marked
response is inferred from the bare odd core alone.

Literal controls at k=2 give

    H(original)=462915, H(normalized)=421425,
    H(first contraction)=45, H(odd core)=3.

Both hostile sources give these same scalars. Structural halving is not
division of H; the ordered-cover substitution kernel can reconstruct H
from the decorated quotient. Exact finite query atlases of sizes 2,8,24,128
and literal twelve-vertex read heads test O4/O5 independently.

## O6 The new return has an exact dyadic carry and a readiness obstruction

Let A be the saved 6325 source, B its repaired 6325 target, with the same
old root name 31 and sole in-neighbour 21. Clone that root by TT_m. The root
crossing matching has size one on both sources. The old spare name 6 has
matching size two on both, with different saved covers.

The inherited matching/run theorem [BC1](../boundary_carry_returns_20261103/PROOFS.md)
restricts the scalar clone response to c1+c2(2^m-2). Exact coefficients give

    H(A_m)=6325+4065(2^m-2),
    H(B_m)=6325+5407(2^m-2).

Thus this particular neutral return has carry debt **1342(2^m-2)**. One
clone-birth context distinguishes its endpoints: at m=2 the counts are
14455 and 17139. Equal matching capacity is not equality of its coefficients.
The same-root-star helper compiler can transport the relation-level repair
over all m, but does not cancel this scalar difference.

Readiness also fails uniformly along this clone extension. For m>=2, retain
the first clone as the root and put x=2^(m-1). The returned profile is

    (H,s,q,F)=(10814x-4489, 5407x-4516,
               5407x-4489, 10814x-4462).

Consequently F-H=27, rather than 1 or 2. O4's first identity implies that F
strictly increases along the clone order, while H is constant. **No clone
position is 1-ready or 2-ready for any m>=2.** At m=1 the actual profile is
(6325,891,893,6327), which is 2-ready; the empty-clone boundary is separate.

Proof of the displayed all-m profile: keep the first root exterior and
vary its remaining transitive twin block of size t=m-1>=1. Its crossing
matching is one. [BC7](../boundary_carry_returns_20261103/PROOFS.md) bounds
each of these one/two-cover queries by the three modes 1,2^t,3^t; deletion q
is encoded by a singleton-root component. The three exact samples t=1,2,3
determine all coefficients by an invertible Vandermonde matrix. The 3^t
coefficients vanish in this first-root profile; t=4 is an extra check.
`matching_carry.json` retains both endpoint-mark profiles and all finite
clone-position defects. This excludes the stated clone family, not exterior
re-marking, other parameterizations, or further native corrections.
