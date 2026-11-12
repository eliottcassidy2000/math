# Ternary hosts: exact births, cyclic sections, and hierarchical protection

**Status: PROVED IN NOTE; independently reviewed.** The claims
below concern one explicit infinite hereditary family and its marked
extensions, not all tournament isomorphism classes. No external priority is
claimed. Chan–Pak remains open.

## 0. Inheritance and the cross-topic map

Let `T_0=K1`, `T_(e+1)=C3[T_e,T_e,T_e]`, with top arcs `0→1→2→0`.
The balanced ternary tree addresses its vertices, first differing digit
most significant first. This is the same tower as
[THM-3134, section4](../../../../01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md).
The ordered-cover transform, its exact values, and wreath divisibility are
inherited, not discoveries of this note. The all-order ordinal/C3 species
and response inventory are also inherited from
[Session268](../../session268_composable_tournaments_20261111/modular/PROOF.md).

The distant inspiration is [THM-3539](../../../../01-canon/theorems/THM-3539-fixed-keller-newest-prime-decomposition-centralizer-and-lca-packet-floor.md):
a large action can be understood through its occupied ancestor tree.
[THM-3511](../../../../01-canon/theorems/THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile.md)
provides the complementary lesson that sections, not shallow scalar
portraits, retain continuation. We do NOT transfer any Keller or Rule30
conclusion to tournaments. We transfer the operation of retaining sections.
Here the local group is `C3`, not `S3`: reflections reverse arcs and cannot
be silently added as automorphisms.

Closest mechanism: the tower's cyclic wreath subgroup. Canonical hostile:
identical size/multiset of component costs can have different embedding
births (§2). Corrected near miss: neither arbitrary permutations of children
nor a Kraft-like unordered capacity determines an ordered tournament.
Least-used sidecar: the actual ancestor occupancy of named endpoints.

Live board: cyclic sections; canonical strong components; alphabetic tree
height; setwise versus pointwise stabilizers; named path endpoints.

## 1. A complete hereditary grammar inside the tower

Let `D` be the smallest class containing `K1` and closed under nonempty
ordinal sums and `C3[A,B,C]` with nonempty children. Then

    D = union_(e>=0) {nonempty induced subtournaments of T_e}.       (1)

The inclusion from right to left follows by inspecting the occupied top
blocks. Three occupied blocks give a cyclic substitution; two give an
ordinal sum in the direction of their arc; one gives its own subtournament.
Empty blocks are deleted, not treated as fake vertices.

For the converse, embed the children in sufficiently large equal-depth
towers. A cyclic triple uses three top blocks. An ordinal sum of two objects
uses two top blocks in the required orientation; an arbitrary finite sum
uses a binary parenthesization. Padding to a larger depth is legitimate:
`T_e` is induced in one top block of `T_(e+1)`. This proves (1).

There is an unambiguous structural form:

* a singleton;
* an ordered list, of length at least two, of strong objects;
* a cyclic triple of arbitrary nonempty objects, modulo rotation, for a
  nontrivial strong object.

Adjacent ordinal-sum nodes are flattened. Only rotations identify cyclic
triples. To justify uniqueness, for `C3[A,B,C]` every proper module is
contained in one top block: a module meeting exactly two blocks is
separated by every vertex of the third; a module meeting all three and
omitting a vertex is separated by that vertex, which sees the other two
blocks in opposite directions. Thus the three children are precisely the
maximal proper modules. Ordered strong components are intrinsic. Recursing
establishes both uniqueness and a complete isomorphism test in `D`.

Consequently this is a genuine full grammar on an infinite family. It does
not exhaust all tournaments: the prime five-vertex tournaments already lie
outside. The inherited ordinary class generating functions are

    A=S/(1-S),       S=x+(A(x)^3+2A(x^3))/3.                       (2)

Equation (2) is explicitly credited to the earlier modular inventory.
Our new addition is the exact host-birth refinement, not the existence of
this species equation.

## 2. Exact first appearance is an alphabetic minimax code

Define `b(U)=min{e: U induced in T_e}` for `U in D`.

    b(K1)=0,
    b(C3[A,B,C])=1+max(b(A),b(B),b(C)).                         (3)

If `U=C_1→...→C_r` is its strong-component list, then

    b(U)=min_B max_i( depth_B(i)+b(C_i) ),                     (4)

where `B` ranges over full ordered binary trees on the component list.
Equivalently, with `d[i,i+1]=b(C_i)`,

    d[i,j]=1+min_(i<k<j) max(d[i,k],d[k,j]),      i+1<j.        (5)

This is an exact finite dynamic program, not an estimate extrapolated from
small levels. It computes a minimum depth and recursively constructs an
embedding achieving that depth.

**Proof of (3).** An embedding can be descended through single occupied top
blocks until it first occupies more than one. Since the source is strong,
it cannot occupy just two top blocks. All three are occupied, and uniqueness
of the source's maximal modules identifies its three children with the three
intersections. This gives the lower bound; embedding the children supplies
the upper bound.

**Proof of (4).** At the first multi-block branching of a nonstrong image,
exactly two top blocks are occupied. The two intersections are consecutive
nonempty lists of whole strong components: a strong component cannot meet
two blocks with only forward cross arcs. Repeat for the two lists. This
produces an ordered binary tree and the lower bound in (4). Conversely place
each binary split in two cyclic top blocks with the correct direction and
embed the individual components at their required depths. Hence equality.
Extra unused initial levels cannot improve a minimum. Equations (3)–(5)
therefore classify every birth in (1), at all sizes.

Examples:

    b(TT_r)=ceil(log_2 r),
    b(C→...→C) = b(C)+ceil(log_2 r)       (r identical strong C).

In particular the largest transitive subtournament of `T_e` has exactly
`2^e` vertices, and an ordinal list of `r` triangles first appears in the
balanced host of order `3^(1+ceil(log_2 r))`.

Order matters, even with the entire multiset of child heights retained:

    b(C3→K1→K1)=2,      b(K1→C3→K1)=3.                        (6)

The common component-cost multiset is `{1,0,0}`. In the first case the
expensive component takes one side of the top split. In the second it sits
between two sides and must receive an additional level. This is the
alphabetic, not ordinary Huffman, constraint. Both examples have five
vertices; their smallest balanced hosts have respectively 9 and 27.

### A one-register dyadic compiler (no cubic DP required)

For the ordered costs `c_i=b(C_i)`, initialize P=0 and process them in order:

    w_i=2^c_i,       P := w_i*(ceil(P/w_i)+1).                (6a)

Then the same exact answer is

    b(U)=ceil(log_2 P).                                     (6b)

This uses one integer register and one alignment operation per component.
It avoids constructing the entire ordered binary tree or trying every split.
The bit complexity of growing integers is retained; this is not a constant-
time claim.

**Proof.** An ordered prefix code in a binary tree of height h identifies
each codeword with a disjoint aligned dyadic interval in `[0,2^h)`, in leaf
order. A leaf for cost c_i must have interval length at least `2^c_i`.
If larger, choose an aligned subinterval of exactly that length within it;
this preserves order and disjointness. Conversely intervals of these exact
lengths are prefix-free codewords at depths `h-c_i`. Complete their binary
trie, deleting unused branches and suppressing unary nodes. This gives an
ordered full binary tree whose weighted leaf depths do not exceed h.
Therefore the minimax feasibility test is exactly the existence of ordered
disjoint aligned intervals of lengths w_i in `[0,2^h)`.

The earliest possible start of interval i is the first multiple of w_i
not less than the previous endpoint. Its endpoint is precisely (6a).
Greedy placement is coordinatewise no later than any legal placement,
by induction and monotonicity of the rounding operation. Hence feasibility
is exactly `P<=2^h`, proving (6b). The example (6) is the smallest visible
alignment-cost distinction, not a shortage of total unaligned capacity.
The verifier compares (6a) against independent split-DP for every cost list
of lengths1–7 with entries0–3 (21,844 lists).

The verifier generates all induced *types*, not all 2^27 subsets, at depth3
by the child grammar and checks constructive minimal-depth embeddings.
The proof, rather than that finite run, supplies all future levels.

## 3. Full marked-extension symmetry grammar

The maximal-module argument also proves

    Aut(T_e)=W_e=W_(e-1)^3 semidirect C3,
    |W_e|=3^I_e,              I_e=(3^e-1)/2.                 (7)

For a subset `S` of vertices, attach one distinguished new vertex with its
out-neighborhood (or in-neighborhood, consistently) equal to `S`. Two such
marked extensions are isomorphic iff their subsets are in one `W_e` orbit.
Their root-fixing automorphism group is the setwise stabilizer of `S`.
No claim is made about forgetting the new root: unmarked extensions can
collide through a different possible choice of deleted vertex.

Let

    O_e(x,z)=sum_(subset orbits [S]) x^|S| z^v3(|Stab(S)|).

Then the exact symmetry-refined recurrence is

    O_0=1+x,
    O_e=(O_(e-1)(x,z)^3
           +(3z-1) O_(e-1)(x^3,z^3))/3.                   (8)

Indeed an ordered triple of child orbit types is taken modulo cyclic
rotation. Nonconstant triples have orbit length3 and stabilizer exponent
the sum of the three child exponents. A constant triple contributes one
orbit and an extra factor3 to its stabilizer. The diagonal triple
polynomial is `O_(e-1)(x^3,z^3)`, so subtract it before dividing the
nonconstant triples by3, then add `z` times it. This proves (8), including
its integer coefficients.

Thus one compact rule gives *every* rooted one-vertex extension family of
every tower, with all symmetry orders, not just its total number. At x=z=1
it yields counts 2,4,24,4624,32955830624 for e=0,...,4. These count extensions
of these specific hosts, not A000568.

Every exponent `0,1,...,I_e` occurs. For e=0 this is immediate. Every child
orbit has a distinct complementary orbit with the same exponent: the host
has odd order, so complement changes subset cardinality. By induction any
sum of three integers in `[0,I_(e-1)]` can be realized with a nonconstant
triple of child types, even if the three exponents coincide (choose two
complementary types). This gives `[0,3I_(e-1)]`. A constant triple of maximum
exponent gives `3I_(e-1)+1=I_e`. Consequently the first balanced host depth
at which the root-fixing group can have order `3^a` is the least e with
`3^e>=2a+1`. This is a family-relative birthday, not a minimum among all
tournaments of every possible order.

The same recurrence works with q colors and a multivariate base
`O_0=x_1+...+x_q`, replacing each variable by its cube in the diagonal
term. Taking q=2^k exactly classifies old-to-new attachments for k named new
vertices with their mutual arcs fixed. The structural representative is a
nested cyclic necklace, not merely the pair (selected size, stabilizer).

The older P12 in
[modular-centers TOURNAMENTS.md](../../modular_centers_collatz_pell_20261111/TOURNAMENTS.md)
already classified the two free binary-subset orbits. Equation (8) retains
every nonfree orbit as well; the free case is an inherited control.

## 4. Endpoint occupation gives a hierarchy of protected factors

Setwise and pointwise stabilizers are different. For a named set `S` of
leaves, let `u(S)` be the number of internal tree vertices with at least one
marked descendant. To fix every marked leaf, the rotation at each occupied
internal vertex must be trivial; rotations in all other subtrees are free.
Therefore

    |Fix_(W_e)(S)|=3^(I_e-u(S)).                            (9)

Among all s-element subsets, the largest occupancy is

    max u(S)=sum_(d=0)^(e-1) min(3^d,s),
    E_e(s):=min v3|Fix(S)|=sum_(d=0)^(e-1) max(3^d-s,0).    (10)

The upper bound at depth d is immediate. It is simultaneous: distribute the
s marks as evenly as possible among the three children and recurse. At
levels with at least s cells this makes the marks occupy distinct cells;
at levels with at most s cells every cell is occupied. Induction or the
balanced recursive allocation proves attainability. For s=0 all I_e nodes
are unoccupied. For s beyond the number of leaves take `E_e(s)=0`.
Thus (10) is the **exact endpoint-stabilizer floor**, not merely a bound.

For every tournament U obtained by adding k vertices, with completely
arbitrary old-new and new-new orientations,

    3^E_e(2k) divides H(U).                                (11)

**Proof.** Group candidate Hamilton paths by the ordered new vertices,
positions of the new vertices, and the exact old labels adjacent to new
vertices. There are at most2k distinct such old labels. For a given group,
`Fix(S)` acts on all old labels, preserving old arcs and every used new-old
arc since its endpoint is fixed. It preserves the grouping, and acts freely
on ordered Hamilton paths: a fixed sequence fixes every old vertex. Its
order is divisible by `3^E_e(2k)` by (9)–(10). Every group size is divisible
by that common power. Summing proves (11). No assumption that U itself has
any nontrivial automorphism is used.

The endpoint-action principle is inherited from P11 in the same
modular-centers note. The new content is its all-level exact pointwise
stabilizer floor and its application, rather than only a single untouched
bottom block. The sharper triangle-chain content formula in
[TRIANGLE_ALL_K_PROOF](../../berggren_target_return_20261111/arithmetic/TRIANGLE_ALL_K_PROOF.md)
is a different source family and is not replaced here.

For one added vertex, e>=1 gives

    E_e(2)=(3^e-1)/2-(2e-1).

Thus every one-vertex extension of the 27-tower has count divisible by
`3^8`; at 81 vertices the divisor is `3^33`. These are universal statements
about extensions of a specified source, not minimum values or exact gcd
claims. The floor is sharp as a **pointwise-group** bound; equality for
Hamilton-count gcd at every e,k has not been proved. At e=2,k=1 the literal
512-extension gcd is3, agreeing with the floor.

This explains a structural obstruction for Chan–Pak: strong self-similarity
can preserve many layers of arithmetic factors even after arbitrary
attachments destroy visible global symmetry. An unrestricted source
selector must leave or sufficiently puncture this family. Neither exact
births nor all rooted extension types automatically supplies a prescribed
absolute Hamilton count.

### Prime-layer extension (same mechanism, no full-Aut claim)

For any odd prime p, orient the regular cyclic tournament R_p by
`i→j iff j-i mod p is in {1,...,(p-1)/2}`. Its translations give C_p.
In the iterated tower R_(p,e), retain only the explicit cyclic wreath
subgroup G_(p,e), whether or not it is the full automorphism group. Then

    |G_(p,e)|=p^((p^e-1)/(p-1)),
    E_(p,e)(s)=sum_(d<e) max(p^d-s,0),
    p^E_(p,e)(2k) divides H(U) for every k-vertex extension U.       (12)

The occupation proof and free path action are unchanged. This gives,
for example, universal one-vertex divisors 5^3 for the 25-vertex tower
and 7^5 for the 49-vertex tower. These are not claims about all tournaments
at orders25 or49, nor about prime-field Paley tournaments. The generic
subset-orbit recurrence under this explicit subgroup is

    O_e=(O_(e-1)^p+(pz-1)O_(e-1)(x^p,z^p))/p.                    (13)

Primality matters: every nonconstant p-tuple has cyclic orbit length p.
For composite local arity, proper periodic necklaces require additional
divisor terms. For p>3 (13) is a subgroup-orbit classification, potentially
finer than actual marked isomorphism classes. Only the ternary maximal-
module argument above promotes it to the full rooted-isomorphism grammar.

## 5. Reproduction and limits

Run `python3 verify.py` and `python3 -O verify.py` from this directory.
Both transcripts agree. Controls include exhaustive greedy/split-DP comparisons, full adjacency-based enumeration of all9! permutations at order9, literal binary color orbits on
1,3,9 vertices; all 512 subsets for the pointwise floor at order9; all512
one-vertex Hamilton counts over the order9 tower using subset DP; exact
orbit-size mass through order81; full induced grammar at depth3; constructive
minimal embeddings; and the ordered-cost hostile (6).

This is a complete recursion for `D`, its host births, and tower marked
extensions. It intentionally does not claim an exhaustive catalogue of
prime quotients of unbounded arity, or that a scalar symmetry profile alone
is a complete extension state. That distinction is precisely why this
family can be understood fully without declaring all tournaments understood.
