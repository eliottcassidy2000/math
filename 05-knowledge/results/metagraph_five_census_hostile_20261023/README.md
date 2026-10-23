# Full rooted-five records still have an all-round hostile pair

2026-10-23. **FINITE-EXACT:** an explicit strongly connected pair at order
23 with equal full induced and rooted censuses through five vertices,
separated by transitive six-sets. **PROVED:** its blindness under every
round of vertex neighborhood refinement, its extension to every order
n>=23, and the structural impossibility of any fixed singly rooted motif
cutoff followed by arbitrary ordinary vertex-refinement rounds. **OPEN:** completeness of
actual-path flags, more elaborate joint incidence observers, or an all-size
metagraph-wiring classifier.

This is the next hostile after the [order-31 strong-four failure](../metagraph_minimal_signatures_20261023/README.md).
Rooted five-vertex records are useful on actual ports, but they cannot be
promoted to a complete all-size tournament signature. The missing coordinate
in this pair first appears in a six-vertex induced count.

## Explicit pair and finite result

On Z/23Z orient i->j when j-i belongs to the following difference set:

    S1 = {1,2,3,4,6,10,12,14,15,16,18},
    S2 = {1,2,3,4,7,9,10,12,15,17,18}.

They are encoded by bitmasks 382046 and 431774. Each contains exactly one
member of every pair {d,-d}, so each gives a tournament. Translations act
transitively on vertices. Both contain the step 1 and hence the directed
Hamiltonian cycle 0,1,...,22,0; they are strongly connected and 11-regular.

An independent Python adjacency/canonicalization census gives:

| Observer | First tournament | Second tournament |
|---|---:|---:|
| Full unrooted five-type vector (12 entries) | equal | equal |
| Full rooted-five vector at each vertex (48 entries) | equal | equal |
| Every induced and rooted census of size<=4 | equal | equal |
| Transitive six-sets | 1173 | 1265 |
| Transitive six-sets containing a specified vertex | 306 | 330 |

The 12-entry vector, in increasing canonical five-tournament word order,
is

    (2921,851,2691,1127,5520,2691,3795,851,5934,1403,4071,1794).

All 48 rooted-five types occur. The total masses are C(23,5)=33649 and
C(22,4)=7315 for the rooted vector. The different transitive six counts
alone prove nonisomorphism, independently of the multiplier argument used
to classify the search universe. Six is the smallest induced motif size
that can separate this pair, in either the unrooted census or complete
vertex-rooted incidence sense. This is pair-specific minimality.

## Why every ordinary vertex-refinement round remains blind

Define the initial color of v to be its full tuple of counts of all rooted
induced tournament types of sizes at most five. Additional degree or
cyclic-triangle coordinates are already contained in this data. Refine a
color by retaining its old value and the multisets of colors in its
out-neighborhood and in-neighborhood, using the same rule in both graphs.

Vertex transitivity makes the initial color constant within each graph;
the exact rooted census makes those two constants equal. At each round
every vertex sees eleven copies of the common color in each direction.
Induction proves equality forever. The assertion remains true if only the
out-neighbor multiset is used, or if equal global <=5 census data are
adjoined to every color.

This is a hostile to rank-one vertex refinement. It does not defeat a
joint multi-vertex observer, a path-marked flag, or all possible graph
algorithms. In particular, the actual path forest in the new flag
transport packet is information not present in these 48 counts.

## Extension to every n>=23

For k>=0 form

    T_i(k) = TT_k -> T_i,

where the same transitive k-vertex prefix dominates every vertex of the
23-vertex core. These tournaments are not asserted strongly connected
when k>0.

The complete rooted <=5 census agrees at corresponding prefix vertices
and at any chosen core vertex. To prove this, partition each selected
motif into its prefix and core parts. Its orientation between the parts
is forced. For a root in the core, the core portion is counted by the
corresponding rooted core census; for a root in the prefix, it is counted
by the unrooted core census. The prefix subobject and root position are
identical on both sides. Summing over the common prefix choices therefore
gives equal initial colors.

All core vertices remain monochromatic: each sees eleven core neighbors
in either direction and the same list of prefix colors as incoming
neighbors. Each prefix vertex sees the corresponding earlier/later prefix
colors and 23 copies of the core color as outgoing neighbors. The same
induction proves equality at every refinement round.

Finally the cores are the unique nontrivial strongly connected components
and are nonisomorphic. Any isomorphism of the padded tournaments would
identify their cores, impossible. Thus for every n>=23 there are
nonisomorphic tournaments invisible to this entire rooted-five vertex
refinement hierarchy. The construction preserves the obstruction by a
common prefix; it is not an infinite strongly connected family.

The independent audit adds a uniform six-motif separator. With TT_r(T)
the transitive r-set count and TT_0=1,

    TT_6(TT_k -> T) = sum_(r=0)^6 C(k,6-r) TT_r(T).

The core counts agree through r=5, so the transitive-six difference stays
exactly 92 at every padded order. This independently proves nonisomorphism
throughout the family and shows precisely what the five-level observer lost.

One could also recover all smaller unrooted motif counts by double
counting extensions of a k-set into a five-set: each occurs in
C(23-k,5-k) five-sets. Vertex transitivity then distributes the marked
counts uniformly. The executable audit computes all smaller rooted and
unrooted vectors directly rather than assuming this consequence.

## No fixed singly rooted motif cutoff can suffice at all orders

The explicit pair is not an isolated failure of the number five. An
elementary counting argument proves a boundary for **every fixed k**.
This is an application of the pigeonhole principle and the prime-cyclic
isomorphism lemma below, not a claim of novelty for those mechanisms.

**Theorem.** Fix k>=1. For every sufficiently large odd prime p there are
nonisomorphic strongly connected p-vertex tournaments which have identical
complete induced and singly vertex-rooted motif-count records of every
size <=k. Every round of ordinary vertex neighborhood refinement started
from these records also agrees. Consequently this observer fails at
infinitely many strong prime orders. By a common dominating transitive
prefix, it fails at every integer order beyond some threshold depending
on k, with no assertion of strong connectivity for the padded objects.

**Proof.** Let t_k be the finite number of unlabelled k-vertex tournaments.
For an odd prime p>=k there are exactly

    2^((p-1)/2)

skew connection sets in Z/pZ. By the proved multiplier-isomorphism lemma,
each tournament isomorphism class contains at most p-1 such sets. Hence
there are at least 2^((p-1)/2)/(p-1) isomorphism classes. Each full
unrooted k-census vector has t_k integer coordinates between zero and
C(p,k), so there are at most

    (C(p,k)+1)^t_k

such vectors. Thus two distinct isomorphism classes have equal vectors
whenever

    2^((p-1)/2) > (p-1) * (C(p,k)+1)^t_k.             (1)

For fixed k the left side is exponential in p while the right side is
polynomial, so (1) holds for every sufficiently large p. This is an
effective finite-search criterion, although the resulting search can
be expensive. For k=5, t_5=12, the prime p=1063 satisfies (1); the actual
23-vertex pair supplies a much smaller explicit witness.

Equality of the k-census implies equality of every lower unrooted census.
For an l-vertex tournament H and l<=k, double-count pairs consisting of
an induced H and an enclosing k-set:

    sum_F N_F(T) N_H(F) = N_H(T) * C(p-l,k-l),         (2)

where F ranges over all unrooted k-types and N denotes induced-subset
counts. The positive extension coefficient in (2) can be cancelled.

Now fix a rooted l-type (H,h), and let O be the orbit of h under Aut(H).
Across all ambient roots, each induced H contributes exactly |O| copies
of that rooted type. Translations make all p roots equivalent, so its
count at any particular root is exactly

    |O| * N_H(T) / p.                                (3)

Thus the two complete singly rooted records also agree, without losing
the root-orbit multiplicity. Their initial vertex colors are constant
and equal. Both tournaments are (p-1)/2-regular, so the same constant
color and neighborhood multisets persist at every refinement round.

Each prime-cyclic tournament is strongly connected: any a in its nonempty
connection set generates Z/pZ, and

    0, a, 2a, ..., (p-1)a, 0

is a directed Hamiltonian cycle. Distinct multiplier classes in the
collision are therefore nonisomorphic strong tournaments. Finally, choose
one such prime p for this k and add the same dominating TT_(n-p) prefix.
The preceding motif-partition and refinement argument works with k in
place of five. The original cores remain the unique nontrivial strongly
connected components, proving nonisomorphism for every n>=p. QED.

The theorem excludes a uniform solution based solely on a **fixed** list
of singly rooted induced motif counts, even if it includes every type up
to that cutoff and permits infinitely many ordinary refinement rounds.
It says nothing against pair-rooted records, joint contractions, actual
Hamiltonian-path markings, adaptive motif depth growing with the order,
or other operations retaining more incidence information. Those are the
precise surviving directions for the all-size metagraph program.
The theorem is an obstruction to complete tournament signatures of the
specified kind. It does **not**, by itself, produce a same-margin
counterfeit of the full metagraph or prove that global port-capacity
constraints must remain ambiguous: those are additional information.

### Quantitative boundary when the motif cutoff grows

The parent lane sharpened the counting argument to a growing cutoff.
For every fixed epsilon>0 and all sufficiently large odd primes p, a
cutoff

    k(p) <= (sqrt(2)-epsilon) sqrt(log_2 p)

still forces a collision of the stated observer. Indeed
t_k<=2^(k(k-1)/2), and

    log_2[(C(p,k)+1)^t_k] <= 2^(k(k-1)/2) (k log_2 p+1).

For k<=c sqrt(log_2 p) with c<sqrt(2), the right side is at most
p^(c^2/2) O((log p)^(3/2))=o(p). In contrast, the base-two logarithm of
the lower bound on isomorphism classes is (p-1)/2-log_2(p-1), asymptotic
to p/2. Thus the same pigeonhole collision follows. If an observer in
this exact family is complete at every sufficiently large prime order,
its cutoff must satisfy

    liminf_(p prime -> infinity) k(p)/sqrt(log_2 p) >= sqrt(2).

This is a necessary lower bound for complete singly rooted motif counts
followed by ordinary vertex refinement on prime cyclic tournaments. It is
not a bound for actual-path flags, joint pair queries, or global port
capacity reconstruction.

## Constructive repair: keep one actual vertex identity

The failure of a root *inside each counted motif* does not imply failure
after distinguishing one actual vertex of the ambient tournament. These
are different operations. `marked_probe.py` gives the selected actual
vertex a unique color and all others the common rooted-five color. It
then refines using both directed neighbor-color multisets, retaining the
old color at each step.

All 23 choices of one marked vertex become discrete after three strict
refinement rounds in both order-23 hosts. Their cell-size profiles already
differ at round two:

    S1: (1,1,1,2,2,2,2,3,3,3,3),
    S2: (1,1,1,1,1,2,2,3,3,4,4).

Both have (1,11,11) after round one and 23 singleton cells after round
three. The exact incidence matrix read in the final color order gives
different canonical words for the two hosts. Thus retaining the identity
of even one endpoint repairs this concrete five-record collision.

The exhaustive ordered-two-mark experiment gives:

| Host | Ordered distinct pairs | Discrete at round 2 | Discrete at round 3 |
|---|---:|---:|---:|
| S1 | 506 | 230 | 276 |
| S2 | 506 | 322 | 184 |

All final partitions are 23 singleton cells. The final marked canonical
words of the two hosts are disjoint. The enumeration includes every
possible pair of actual Hamiltonian-path endpoints, but does not assert
that every ordered pair can be such an endpoint pair.

The hostile control prevents another all-size promotion: in the Paley
tournament P_31, each of its 31 choices of one marked vertex stabilizes
after one round at (1,15,15), not a discrete partition. Therefore the
single-mark repair is a finite property of this pair. The experiment
supports retaining actual endpoints and adding joint marked records where
needed; it proves no universal single-mark theorem.

    python -B marked_probe.py
    python -B -O marked_probe.py

The complete marked universe has 1,058 p23 cases, plus all 31 one-mark
Paley31 controls. The receipt has 7,494 generating checks and 7,495 on
replay. It checks monotone partition refinement, every final singleton
partition, all marked universes, and the independent hostile control.

## Complete bounded prime-cyclic search

The C++ discovery program enumerates every skew difference set at each
listed prime p, then keeps its least numerical representative under all
nonzero multipliers modulo p. It computes the complete five-tournament
census using an independently generated 1024-word canonicalization table.

| p | Cyclic tournament classes | Distinct five-census vectors | Collisions |
|---:|---:|---:|---:|
| 13 | 6 | 6 | 0 |
| 17 | 16 | 16 | 0 |
| 19 | 30 | 30 | 0 |
| 23 | 94 | 93 | 1 |
| 31 | 1096 | 1094 | 2 |

The two order-31 collision pairs, retained as finite discovery outputs,
are masks (24985982,126262814) and (239508878,167952302). The independently
audited theorem in this packet uses only the explicit order-23 pair.
The smaller rows give no claim about arbitrary tournaments of those orders.

**Why multiplier orbits are full isomorphism classes at prime order.**
Let P be the translation subgroup of order p. Any isomorphism phi between
two such tournaments carries P to another order-p subgroup of the target
automorphism group. The target group embeds in S_p, whose order contains
exactly one factor of p, so both subgroups are Sylow p-subgroups. Some
target automorphism beta conjugates phi P phi^(-1) to P. Then beta phi
normalizes P, hence is an affine permutation x->a x+b: the normalizer of
a regular p-cycle has exactly this form, as follows by conjugating the
translation x->x+1 to x->x+a. Such an isomorphism sends S to aS.
Conversely every such multiplier gives an isomorphism. No computational
assumption about the absence of extra non-affine automorphisms is used.

## Reproduction and scope of computation

    g++ -O3 -std=c++17 search.cpp -o search
    ./search
    python -B check.py
    python -B -O check.py

Run from this directory (choose a temporary executable path on Windows).
`search.cpp` is the complete bounded discovery scan. `check.py` imports
none of it: it builds actual directed adjacency matrices, independently
canonicalizes every small tournament over all permutations, counts every
induced <=5 set and every <=5 set containing root 0, checks the translation
automorphism and strong directed cycle, and enumerates all six-sets via
their actual score vectors. Scores 0,1,...,5 characterize a transitive
six-tournament. Vertex transitivity transports root-0 counts to every
vertex, with no sampling.

The frozen Python receipt has 26 aggregate acceptance checks (27 on replay)
over these complete finite enumerations. Its LF SHA256 is
`e07d815325e112874961e6cd6943baf122052accc06231569092ac3e8527cbac`.
Equal five-census data are a positive equality control; differing
transitive-six data and disjoint multiplier orbits are independent hostile
controls. A separate peer audit accompanies this packet.
