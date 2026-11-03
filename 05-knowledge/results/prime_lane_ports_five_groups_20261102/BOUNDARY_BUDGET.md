# A sharp symmetry budget for arbitrary multi-vertex attachments

**PROVED:** the boundary-divisibility theorem below, with a sharp family.
**FINITE-EXACT:** complete small attachment controls. No prime-lane
completeness or asymptotic prime production is claimed.

## Inheritance and the new operation

The closest mechanism is the inherited [weighted one-apex valuation gate,
BRANCHES §2](../prime_inverse_ports_20261102/BRANCHES.md): in an ordered block
word, at most two matrix transitions can avoid the p-divisibility of its
diagonal H/E factors. Its unweighted r-2 version is credited there to
Session124, and it already gives a sharp strong one-apex family.

Here there may be **any r new vertices**, with arbitrary mutual arcs and
arbitrary attachments to the old vertices. The old source need not be an
ordered block word. The retained structure is a product of independently
supported source automorphism groups. This is a different sufficient
divisibility certificate; it does not replace the older, sometimes stronger
path-cover-content valuation. The hostile to a stronger uniform claim is
the seven-vertex H41 seed, which removes two triangle factors with one birth.

## B1. General boundary-orbit divisor

Let T be a finite directed graph on a named vertex set V. Suppose pairwise
disjoint subsets M_1,...,M_k have groups G_i of automorphisms of T supported
on M_i: each element fixes every vertex outside M_i. Assume these actions
are faithful and write g_i=|G_i|. Their supports are disjoint, so the product
of any subcollection acts faithfully on T. Modules with internal
automorphisms are a standard source of such actions.

Let U be any extension of T by r named new vertices X with **U[V]=T**:
no old edge is removed, added, or reoriented. Define

    D_r = gcd_{S subset {1,...,k}, |S|<=2r} product_{i notin S} g_i.

Then D_r divides the number H(U) of directed Hamiltonian vertex orders.
The theorem also applies to undirected graphs, counting ordered directed
traversals of their Hamiltonian paths. Empty path sets have count zero.

Equivalently, for every prime ell, let w_i=v_ell(g_i), W=sum_i w_i and
w_(1)>=...>=w_(k). Then

    v_ell(H(U)) >= W - sum_{j=1}^{min(2r,k)} w_(j),       (B1)

with the conventional infinite valuation when H(U)=0.

**Proof.** For each Hamiltonian word P of U, let S(P) be the set of old
blocks containing a vertex adjacent along P to a vertex of X. Each new
vertex has at most two neighbors in the word, so |S(P)|<=2r. Partition the
path set by the exact set S(P).

On the stratum S, act by product_{i notin S}G_i on the old vertex names and
fix all other names. Every old-old step remains an edge because these are
source automorphisms. Every step incident to X remains unchanged, since
all its old endpoints lie outside the acted-on blocks. Thus the new word
is still a Hamiltonian word of U. Its touched-block set remains S.

This action on full vertex orders is free: an element fixing the entire
word fixes every old vertex, and faithfulness then makes it the identity.
Consequently the size of stratum S is divisible by product_{i notin S}g_i.
The common divisor D_r divides their sum H(U). Prime-wise minimizing the
surviving weight removes at most the2r largest w_i, yielding(B1). QED.

**What is not assumed:** the groups G_i need not extend to automorphisms
of U. Attachments can break all source symmetries globally. The action is
on each actual path stratum, using the absence of a new-vertex adjacency
at that block. This is why merely inspecting Aut(U) misses the divisor.

## B2. Triangle sources and sharpness for every budget

Suppose T has k disjoint independently rotating directed-triangle modules.
Choose G_i=C3. Every extension by r arbitrary vertices satisfies

    3^max(k-2r,0) divides H(U).                           (B2)

In particular a construction targeting an odd prime p other than3, or7p
with p an odd prime other than3, must have

    r >= ceil(k/2).

This is a necessary cost for this particular retained source structure.
It does not say a source with sufficient r realizes every prime, nor that
all tournament constructions start with such a source. For a general
target h it gives r>=ceil((k-v_3(h))/2) when the right side is positive.

The exponent is sharp for every k>=2r. Let W be the following fixed source
extension: old vertices0,1,2 form0->1->2->0; old vertices3,4,5 form
3->4->5->3; every first triangle vertex beats every second triangle vertex.
Adjoin vertex6, which beats precisely0 and3 among the old vertices.
The inherited exact two-triangle table gives H(W)=41; the independent
checker also enumerates all literal vertex orders. W is strong.

Take the ordered join of r copies of W and k-2r further directed triangles.
Deleting the r designated new vertices leaves exactly the ordered join
of k directed triangles, with every old arc retained. Its path count is

    H = 41^r 3^(k-2r),

so v_3(H)=k-2r. This proves the all-budget lower bound cannot be improved
without additional hypotheses. These joins need not be strong when there
is more than one component; the sharpness claim is for all tournaments.
The older one-apex source separately supplies strong sharp examples.

## B3. The square-sum connection retains actual incidences

The proof uses only path adjacency, source automorphisms, and the number
of new vertices. It therefore transports verbatim to a square-sum graph
when its **actual** source automorphisms with disjoint supports are supplied.
Deleting r vertices from a path exposes at most2r old-new incidences in
either setting. This is an exact shared boundary mechanism, not an equality
of Hamiltonian counts or a conversion of a symmetric relation to a tournament.

The owner's five groups of primes below100 are not automatically modules or
symmetry orbits of their square-sum graph. Group membership alone supplies
none of the G_i required by B1. The [prime-group audit](prime_groups.md)
retains the named prime ports and tests the group quotient directly. An
arbitrary partition into five color classes cannot be substituted into B1.

This also refines the order-five lesson: a constructive extension must keep
which old vertices meet a new vertex **along the chosen paths**, not only
block sizes, aggregate H, or abstract symmetry names. A count of touched
blocks gives the divisor; realizing a desired value additionally requires
the exact surviving stratum counts and compatibility of their shared cuts.

## B4. New-edge congruence and sharpness with a strong target

The boundary budget becomes stricter for a path required to use an edge
between two new vertices. Such an edge uses two of the at most2r incidences
at new vertices, leaving at most2r-2 old-new adjacencies. This gives a
second all-size consequence of the same supported-group action.

**New-edge congruence (PROVED).** Assume r>=2, and let U and U' differ only
in the orientation of one edge between two of the new vertices. In the
notation of B1,

```
H(U')-H(U) = 0 (mod D_(r-1)).                             (B4)
```

For k independently rotating triangle modules this specializes to

```
H(U')-H(U) = 0 (mod 3^max(k-2r+2,0)).                    (B5)
```

**Proof.** Paths not using the changed pair survive literally and cancel
in the difference. Partition each of the two remaining path sets by its
exact touched old-block set. Every affected path contains the prescribed
new-new edge, so at most2r-2 old blocks can be touched. The product of groups
on untouched blocks preserves this edge condition because it fixes every
new vertex, and acts freely exactly as in B1. Each side's affected-path
count is divisible by D_(r-1), and so is their difference. QED.

Therefore (B4) also holds after any sequence of changes to edges wholly
among the new vertices. The old source and all old-new ports remain fixed.
This is stronger than merely knowing that H is divisible by D_r. In the
triangle case, its first nonzero3-adic digit remains fixed under new-edge
rewiring whenever the original count attains the B2 exponent.

**Strong sharpness (PROVED).** For every r>=1 and k>=2r there is a **strong**
extension of the ordered join of k triangles by r vertices, with all old
edges retained, such that

```
v_3(H)=k-2r.                                             (B6)
```

Put t=k-2r. Let W_t be the one-apex extension of t+2 ordered triangles in
which the apex beats exactly one vertex of the first triangle, exactly
one of the second, and no vertex of the remaining t triangles. It is
strong: the apex reaches the first old component and every last component
returns to the apex. The inherited one-apex identity in
[BRANCHES §2](../prime_inverse_ports_20261102/BRANCHES.md) gives

```
H(W_0)=41,
H(W_t)=3^t(45*4^t-1), t>=1.                              (B7)
```

The t=0 exception is essential: substituting zero into the second formula
would incorrectly give44. Thus v_3(H(W_t))=t in every case.

Start with the ordered join of r-1 copies of W_0 and one copy of W_t.
Its old source is the ordered join of k triangles; its r named apices are
the only new vertices. Its count is41^(r-1) H(W_t), of exact3-adic
valuation t. If r=1 it is already strong. If r>=2, reverse the one edge
from the first apex to the last apex, so that it points from the last
apex back to the first.

All component interiors remain strong. Forward connections between each
pair of components remain available through their old vertices; the new
back edge closes a directed route from the last component to the first.
Thus every vertex reaches every other vertex and the target is strong.
By(B5) the changed count differs by a multiple of3^(t+2), so its valuation
remains exactly t. No old edge or old-new port changed. QED.

No all-strong statement is made at r=0: an ordered join of several triangles
is then the fixed nonstrong source. Also strong sharpness does not imply a
prime target. For example, the r=2,k=4 construction changes1681 to
7009=43*163, retaining the nonzero residue modulo9 but producing a composite.
The relevant positive statement is sharpness of the necessary symmetry
budget inside strong tournaments, not universal realization at that budget.

Independent full-adjacency controls use(r,t)=(1,0),(1,1),(1,2),(2,0),(2,1):

|r|t|vertices|before|after reversing the new edge|
|---:|---:|---:|---:|---:|
|1|0|7|41|41 (no reversal needed)|
|1|1|10|537|537 (no reversal needed)|
|1|2|13|6471|6471 (no reversal needed)|
|2|0|14|1681|7009|
|2|1|17|22017|107985|

The last differences are5328, divisible by9, and85968, divisible by27.
The independent checker verifies strongness, exact valuation, and every
retained old edge. Its no-producer-import implementation and controls are
in [boundary_independent.md](boundary_independent.md).

## Controls and next interface

Run `python3 -I -B boundary_budget.py` and repeat with `-O`. It generates
every one-vertex attachment to k ordered triangles for k=1,2,3 (8,64,512
cuts). At k=2 all Hamiltonian orders and their exact touched-block strata
are enumerated literally and checked against an integer subset DP. At
(k,r)=(1,0),(2,1),(3,1),(4,2),(5,2), full adjacency DP verifies the sharp
join formula and literal deletion recovers the original named source.
No strongness, favorable cut, or prime-output filter is used in the
complete attachment tests. The sharp family is a separate construction.

The useful remaining question is now quantitative: for r at or above this
necessary budget, which genuine boundary-stratum inventories admit every
required prime or7p? The [controlled lane exchange](realization.md) gives
explicit strong families crossing between some prime and7p values; neither
that family nor B1 makes sufficient budget into universal realization.
