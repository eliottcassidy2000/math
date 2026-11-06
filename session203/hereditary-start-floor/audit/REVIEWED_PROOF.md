# The entire four-choice interface has a hereditary positive-start floor

Status: **PROOF CANDIDATE, pending independent audit.** This is an all-order
obstruction to using one specified graph interface in the endpoint budget of F(8,137).
It is not a counterexample to F(8,137), and it proves no spectrum exclusion.

## 1. Rooted path insertion through a strong induced kernel

For a tournament T and vertex r, write s_r(T) for the number of directed
Hamilton paths beginning at r. Let S be a nonempty induced subtournament.
Suppose r reaches every vertex of T. Choose a directed path P from r to S
whose only vertex in S is its last vertex v. If r is already in S, use P=(r)
and v=r. Then

    s_r(T) >= s_v(S) >= min_(x in S) s_x(S).                 (1)

Here the second inequality is useful when S is strong, so its minimum is
positive. No assumption that r belongs to S is needed.

To prove the first inequality, concatenate P with any Hamilton path of S
starting at v, listing v only once. This is a directed path starting at r
on the fixed vertex set K=V(P) union V(S). While K is not all of V(T), choose
an edge x->w with x in K and w outside K. Such an edge exists because r is
in K and reaches every vertex of T: a directed r-to-outside path has a
first edge leaving K. Fix a deterministic choice of w, depending only on K.

A vertex w having an incoming arc from K can be inserted into EVERY directed
Hamilton path p_1,...,p_t of T[K] without placing it before p_1. If p_t->w,
append w. Otherwise the indicators [p_i->w] contain a1 and end in0, so
some consecutive pair is1,0; insert w in that gap. Resolve any choice by a
fixed deterministic rule. This preserves the first vertex r. Repeat until
all vertices are included.

Deleting every vertex outside S from the final word recovers the original
Hamilton word on S. Therefore distinct initial words give distinct outputs,
proving (1). In particular, if s_r(T)>0, then r reaches every vertex and
(1) applies. This is a rooted version of the elementary path-insertion
injection; no literature-priority claim is made.

## 2. The four explicit strong kernels

Use the actual four-choice construction from
`session203-four-choice-release/PROOF.md`. The old core has outgoing rows

    (102,0,10,83,135,30,182,47).

Replace old vertex5 by a transitive spine of length5+2k. Old vertex3 beats
spine positions precisely according to

    1^k c 0^k,

where c is respectively01010,01110,10001,10101 in the increasing spine
order. Every other outside incidence is inherited uniformly from the old
core. The old vertex1 remains a global sink.

At k=0, delete that global sink to obtain four11-vertex tournaments S_c.
Direct integer Hamilton-path counting and directed reachability give:

    middle c   H(S_c)   min_x s_x(S_c)   S_c strong
    01010       7213          119           yes
    01110       7221          127           yes
    10001       7227          133           yes
    10101       7235          141           yes

`KERNELS.json` gives their complete outgoing rows and every start count.
These four literal finite objects are the only counted inputs needed for
the all-order conclusion.

For every k>=0, retaining the seven old outside vertices and the five middle
spine vertices induces exactly the k=0 member with the same middle word.
Removing its old sink therefore exhibits S_c as an induced subtournament of
the k-th member. This is a literal vertex selection, not a profile identity.

## 3. Consequence for the bounded endpoint menu

Let T be ANY tournament containing ANY one of the four-choice members at ANY
scale k as an induced subtournament. For every root r of T, either

    s_r(T)=0,

or, by (1) applied to the retained strong kernel,

    s_r(T)>=119.                                          (2)

For a specified middle word the sharper bound is119,127,133, or141 from
the table. This allows arbitrary additional vertices, arbitrary cross arcs,
and arbitrary rerooting. It is stronger than retaining the original root
under O/Z padding or a common rooted amalgamation.

The F(8,137) criterion permits only

    s=r+3j, r in {1,2}, 0<=j<=8,

hence1<=s<=26. Equation (2) rules out every such endpoint for every T in the
stated hereditary class. Thus an attempt to use the large four-choice fiber
for this particular bounded-start criterion must destroy its induced strong
kernel, for example by deleting vertices or changing internal arcs. Merely
rerooting it or attaching further vertices cannot do so.

This does not rule out other physical fibers, non-induced rewrites, larger
endpoint budgets, the original staggered criterion, or unrestricted odd
Hamilton-path coverage. The four-choice theorem itself remains valid.

## 4. Stronger result for arbitrary colored bodies

The previous conclusion does not depend on the transitive spine or its four
particular color words. Retain the same old eight-vertex core and replace
vertex5 by ANY tournament body Q of order at least3. Give every old vertex
other than3 its inherited uniform incidence to Q. Allow the incidences
between old vertex3 and individual body vertices to be completely arbitrary.
Delete the old global sink1. Call the resulting graph K(Q,c), where c records
those arbitrary incidences.

There are exactly eight labelled tournaments Q on three vertices and eight
colorings c, so exactly64 labelled three-vertex-body kernels to check. Every
one is a strong nine-vertex tournament and has minimum root-start count at
least27. The bound is sharp. The twelve colored isomorphism representatives
have the following exact minima (bit i of c says old3 beats body vertex i):

    body rows    colors c        corresponding minimum start counts
    (6,4,0)      0,1,2,3,4,5,6,7    27,35,31,39,29,37,33,41
    (2,4,1)      0,1,3,7            45,51,57,63

The first row is transitive; the second is a directed triangle, whose cyclic
symmetry identifies the remaining colors. `COLORED_KERNELS.json` supplies
all64 full adjacency lists and every actual root-start count, so the proof
does not depend on the representative classification.

For a larger arbitrary body Q, retain ANY three of its vertices together
with the six remaining old nonsink vertices. This induces exactly one of
those64 kernels. Therefore the insertion lemma proves the following
all-order statement:

**Interface floor theorem.** If a tournament T contains any member of this
arbitrary-colored-body interface with body order at least3 as an induced
subtournament, then every root r of T satisfies

    s_r(T)=0 or s_r(T)>=27.                               (3)

Again, arbitrary additional vertices, cross arcs and rerooting are allowed.
Since the largest endpoint allowed by F(8,137) is26, this entire interface
cannot provide any of its required positive-start witnesses. The four-choice
floor119 is a sharper statement for its special subfamily, but the floor27
already excludes every arbitrary-colored or nontransitive body in the same
interface. Neither assertion excludes unrelated graph constructions.

The finite kernels are independently replayed with a separately labelled
constructor, root-by-root forward DP, strong-connectivity search, and a
permutation count on a sharp kernel.

## 5. Exact controls

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The standalone verifier builds all four kernels, checks every tournament
edge and strong connectivity, counts every root by subset dynamic programming,
and verifies the retained induced-kernel map for k=0 through20. It also checks
the root-preserving insertion lemma for every labelled tournament through
order4 and every relevant root/kernel choice. Exceptions rather than
removable assertions guard every check. Finite controls support the literal
objects and implementation; the arbitrary-scale and arbitrary-extension
quantifiers follow from Sections1-3.


The stronger finite-kernel certificate is independently reproduced by:

    python3 verify_colored_floor.py > COLORED_FLOOR_NORMAL.json
    python3 -O verify_colored_floor.py > COLORED_FLOOR_OPTIMIZED.json
    cmp COLORED_FLOOR_NORMAL.json COLORED_FLOOR_OPTIMIZED.json

This checks all64 labeled colored bodies, all their endpoints, and strong
connectivity. A factorial permutation check independently confirms all
root counts of one sharp kernel. The arbitrary-body and arbitrary-extension
claims follow from induced restriction and the insertion proof, not an
extrapolation from body order3.
