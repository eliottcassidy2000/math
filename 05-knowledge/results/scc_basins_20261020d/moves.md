# Releasing endpoints at four- and five-vertex SCC boundaries

Status: PROVED for the bounded black detour theorem below; FINITE-EXACT
controls are labelled separately. All-order black-core connectivity is OPEN.
Scope: the converse-merged complement-line metagraph, not single-arc flips.

## 1. Inheritance and target

The [all-order descent theorem](../black_descent_20261020c/cycles.md) puts
every black component in contact with a nonstrong tournament. The
[terminal-basin reduction](../black_descent_20261020c/scc.md) therefore
leaves the problem of joining those SCC boundary states. The concurrent
[prefix theorem](prefix.md) removes the old two-singleton guard from its
four-step interchange whenever both endpoints are black-active.

The closest complementary mechanism is the exact common target in
[the nonstrong rotation lemma](../defect_transitions_20261020/blue.md),
Section 2. Its fixed-boundary implementation fails on

    R4 -> Z -> R4^op,

with Z pure blue: the standard first and last vertices force every path
between them to be blue. That is the canonical hostile, proved in Section
5 of that note. The least-used coordinate is the choice of endpoint inside
the terminal R4. The present construction releases it explicitly.

Concept board: triangle basins; ordered SCCs; a retained middle path;
endpoint degrees; path-edge XOR; and black detours around blue edges.
The new result joins basins across a previously obstructed family of
blue lines. It does not infer connectivity from degree inventories.

Write K_P(T) for complementing every arc outside a directed Hamiltonian
path P. Thus

    K_P(T) = T^op triangle E(P),
    d_(K_P(T)) = (n-1)1-d_T+e_first(P)-e_last(P).       (1)

For an odd number q of successive path presentations P_1,...,P_q,

    K_(P_q)...K_(P_1)(T)
       = T^op triangle (E(P_1) xor ... xor E(P_q)).    (2)

Every line is unchanged by replacing this operation with full path
reversal and then taking converse. A presentation is blue precisely when
its positional reversal is an anti-automorphism.

## 2. Uniform detours for every four-vertex outer SCC

**Theorem (PROVED).** Let U be a nonstrong self-converse tournament whose
first and last strong components have order four. Every blue line
incident to U has a black-only detour of length at most five between its
endpoints. In particular every such blue line is nonbridging, and its
endpoints belong to the same black component.

There is a unique strong four-vertex tournament R4. In path coordinates
one representative has outgoing-neighbor bitmasks

    R4=(2,4,9,3).

Thus U has the form R4^op -> Z -> R4, where Z may be empty. We do not
require Z to be pure blue. Fix the blue path defining the line under
consideration. Its restrictions to the two outer SCCs are paired by
positional reversal, and its middle restriction is blue.

The five Hamiltonian paths in this representative of R4 are

    0123, 1230, 2301, 3012, 3120.

Reading R4 along these paths gives the following five right-hand path
presentations B. The left R4 presentation is its reversed converse.
Put A=(the left R4 path, the middle Z path). In every word below A is
kept contiguous and in the same forward order. Write 0,1,2,3 for the
four right-hand vertices; these labels are disjoint from those in A.
The blue path whose complement is the target is

    P_0=(A,0,1,2,3).

Apply K successively along the words in the indicated row:

| B outgoing bitmasks | successive black path words |
|---|---|
| (2,4,9,3) | (A,1,2,3,0); (2,A,1,3,0); (2,A,0,1,3) |
| (2,12,9,1) | (A,2,0,1,3); (1,3,A,2,0); (1,2,3,A,0) |
| (6,12,8,1) | (A,1,2,3,0); (3,A,1,2,0); (0,3,A,1,2); (0,2,3,A,1); (2,3,A,0,1) |
| (6,4,8,3) | (A,1,2,3,0); (2,A,1,3,0); (2,A,0,1,3) |
| (10,4,9,2) | (A,1,2,0,3); (2,0,3,A,1); (2,3,A,0,1) |

### Path validity and exact target

These are finite boundary identities that extend to an arbitrary A path.
Its consecutive internal path arcs are retained by every K operation,
so the forward A word is always available. Its other internal arcs are
complemented at every step; the number of steps is odd, exactly as in
K_(P_0). Every arc from an interior A vertex to B that is never a path
boundary likewise has the correct odd complementation parity.

The only remaining information is the first and last A vertices and
the four vertices of B. Replace A temporarily by a two-vertex path
(L,R), initially dominating B. For each of the five rows, direct
six-vertex substitution checks every displayed path arc and

    E(P_1) xor ... xor E(P_q) = E(P_0).                (3)

The checker below performs these five finite identities with no
isomorphism search. Since every visit to A has the same direction and
the only outside attachments are at L or R, the same boundary identities
hold for A of any length at least two. Here |A|>=4. Equations (2) and
(3) give the literal labelled target K_(P_0)(U).

### Blackness

First consider P_1. It changes only the path through the final R4,
leaving the initial R4 and Z words fixed. If both P_0 and P_1 were blue,
the composition of their positional reversals would induce a nonidentity
automorphism on the final R4. The restriction is nonidentity because
the two terminal words differ. Every automorphism preserves each SCC,
and R4 has trivial automorphism group. For completeness, its score
classes have sizes two and two, so an automorphism has order dividing two;
a tournament has no nontrivial involutory automorphism, since it would
reverse the arc between a swapped pair. This contradiction proves that
P_1 is black in all five rows.

For the subsequent steps it is enough to test endpoint degree sums.
A blue presentation would require that sum to be n-1. Formula (1)
gives the following differences from n-1, in chronological order after
the first step:

| B | endpoint-sum differences for P_2,...,P_q |
|---|---|
| (2,4,9,3) | n-5, 6-n |
| (2,12,9,1) | n-4, 4-n |
| (6,12,8,1) | n-5, 4-n, n-4, 4-n |
| (6,4,8,3) | n-5, 5-n |
| (10,4,9,2) | n-4, 4-n |

All are nonzero because n>=8. These computations use only the four
terminal scores and the initial score of the first A vertex. The latter
is n-1-d_B(3), by the original blue reversal; no hidden hypothesis on Z
or its internal degrees is used.

Thus every displayed move is black, and the endpoint is the original
blue neighbor. Removing repeated vertices from the walk gives a black
path of length at most three or five. This proves the theorem.

## 3. Regular five-vertex outer SCCs need only three black moves

**Theorem (PROVED).** Let the first and last strong components of a
nonstrong self-converse U both be the regular five-vertex tournament R5.
Every blue line incident to U has a black-only detour of length at most
three. The middle Z is arbitrary; it need not be pure blue.

The mechanism is the same, with A=(the first R5 path, the Z path).
R5 on residues modulo five with positive increments {1,2} has 15
Hamiltonian paths in three automorphism orbits. Representatives are
01234, 02341, and 02413. Reading along these words gives the three
right-hand presentations in the table. The target blue word is
(A,0,1,2,3,4).

| B outgoing bitmasks | successive black path words |
|---|---|
| (6,12,24,17,3) | (A,3,0,1,2,4); (1,2,4,A,3,0); (1,2,3,4,A,0) |
| (18,12,9,17,6) | (A,0,1,3,4,2); (4,A,0,3,2,1); (2,4,A,0,3,1) |
| (10,20,9,18,5) | (A,0,1,4,2,3); (4,3,A,0,2,1); (3,A,0,2,4,1) |

The completeness of these three path presentations is a finite
five-vertex fact: direct enumeration gives precisely these three
induced matrices, each five times. The accompanying checker verifies
this over all 120 permutations. Equality of two induced matrices means
exactly that the map between their ordered words is an automorphism.

The seven-vertex boundary ledger, using a two-vertex A proxy, proves
path validity and the odd edge-XOR identity exactly as in Section 2.
The first move changes the right-hand Hamiltonian path to a different
automorphism orbit: in each row, reading B along that new terminal word
gives a matrix different from B. If this move and the original path
were both blue, their reversals would induce an automorphism between
those two paths in the last SCC, which is impossible. Thus the first
move is black even though its endpoint degree sum is balanced.

Every right-hand vertex initially has score 2, and the first A vertex
has score n-3. Formula (1) gives endpoint-sum differences from n-1
equal to n-5 and 5-n at the second and third moves. Both are nonzero
since n>=10. Hence those moves are black too, and their labelled target
is precisely the original blue partner.

## 4. What this adds to SCC normalization

### Extension to every strong five-vertex outer component

**Theorem (PROVED by finite boundary certificates).** The conclusion of
Section 3 holds with R5 replaced by any strong five-vertex tournament S.
Thus every blue line at S^op -> Z -> S has a three-black detour, for
arbitrary middle Z and every blue path presentation.

The exact finite certificate set is
[moves_five.json](moves_five.json), produced and replayed by
[moves_five.py](moves_five.py). There are 2^6=64 tournaments on five
vertices with the fixed path 01234; exactly 50 are strong. The checker
visits all 64, using literal reachability for strong connectivity, and
supplies one three-word certificate for each of the 50 strong matrices.
Each word retains A=(the initial S path, the middle Z path) contiguously
and forward, and no word ends in A. The first word changes the final S
path to a different automorphism orbit, certified by unequal induced
matrices. This is precisely the first-step blackness argument in
Section 3; it does not assume S has trivial automorphism group.

For every certificate the seven-vertex two-port proxy verifies all path
arcs and the exact odd edge-XOR identity. The next two steps have symbolic
endpoint-sum differences among

    n-j, j-n (4<=j<=8), or one of -2,-1,1,2.

They are nonzero for all n>=10. The score of the first A vertex is
n-1-d_S(4), supplied by the original blue reversal; no degree of an
interior A vertex is used. The same lifting proof as Section 2 therefore
proves the unbounded theorem. This is a complete finite boundary proof,
not an extrapolation from sampled middle tournaments. As an additional
control, all 25 blue middle presentations of orders 0,...,5 yield 1,250
literal detours and 3,750 independently checked black steps.

The theorem repairs the entire infinite hostile family
R4 -> Z -> R4^op, including every blue path presentation, not merely
the standard one. Its moves change global endpoints while preserving
the chosen middle path as a contiguous block. A fixed-endpoint search
cannot discover this route.

The middle may be mixed. Consequently the theorem also handles examples
outside the incoming mixed-center detour theorem's padding argument:
here blackness can come solely from choosing incompatible outer R4 paths,
even when the middle has no black path at all.

The theorems identify the two endpoints of every such blue line in the
black-component quotient. They do not prove that every SCC basin can be
transformed to one with an outer SCC of order four or five. Larger outer
components, and strong mixed blue sources, still need additional moves.

## 5. Reproduction and exact scope

Run from any directory:

    python 05-knowledge/results/scc_basins_20261020d/moves.py
    python -O 05-knowledge/results/scc_basins_20261020d/moves.py

The self-contained checker verifies five six-vertex R4 boundary
identities and three seven-vertex R5 identities, the complete path
presentation lists, R4's trivial automorphism group, and the affine
endpoint-degree expressions. It then tests every fixed-path
blue middle presentation through order seven, with all five outer path
presentations for R4 and all three for R5, checking each actual path,
its blackness, and literal
equality to the original blue partner. The resulting counts and receipts
are in [moves.json](moves.json). These finite tests audit the symbolic
argument; the unbounded conclusion follows from the boundary identities
and the degree formulas above.

There are 601 blue middle presentations in the stated range. They give
3,005 complete R4 detours and 1,803 complete R5 detours, comprising 15,626
individually checked black steps. No isomorphism search is used to verify
their final targets.

## 6. General endpoint reductions

The [inherited nonstrong rotation lemma](../defect_transitions_20261020/blue.md),
Section 2, gives a common black neighbor for two complement partners whose
source paths have the same first, second and last vertices. Its internal
shared-edge version is in Section 5 of that source. The following uses
retain those exact coordinates.

**Outer endpoint criterion (PROVED).** Fix a blue path at
U=S^op -> Z -> S. If its terminal S path has a Hamiltonian replacement
ending at the same vertex but in a different Aut(S)-orbit, the blue
line has a three-black detour.

Indeed change only that last component path. The global first, second
and last vertices are unchanged (a nontrivial strong S has order at
least three). The new whole path is black: otherwise the composition of
the two positional anti-reversals would identify the two path orbits
inside the terminal SCC. Its black complement partner and the original
blue complement partner have the inherited common black hub. Therefore
the remaining outer-component problem concentrates on **endpoint-rigid
presentations**, in which every terminal path ending at the specified
vertex lies in one automorphism orbit. This is a sufficient criterion,
not a claimed classification of those presentations.

**Nonstrong middle criterion (PROVED).** A nonstrong self-converse mixed
tournament C has a black Hamiltonian path sharing the first vertex of
any prescribed blue path. Consequently every blue line at 1 -> C -> 1
has a three-black detour when C is nonstrong and mixed.

Write C=S -> D -> S^op according to its first and last SCCs. If the
path-orbit count h(S)=H(S)/|Aut(S)| is greater than one, change only the
last component path to a different orbit; the first vertex stays fixed,
and the composition argument above proves blackness. If h(S)=1, the
SCC purity law makes D mixed: otherwise C would be pure blue. Substitute
any black path of D, retaining both outside paths; if the whole path
were blue, its positional reversal would restrict to an anti-automorphism
of that black D presentation, a contradiction. The first vertex again
stays fixed. In 1 -> C -> 1 this preserves global first, second and last,
so the common-hub lemma applies.

This also explains the paired C3 case: C3 -> D -> C3 is mixed precisely
when D is mixed, since h(C3)=1. Its paired outer blocks have length three,
so the inherited mixed-center hub lemma already handles every one of
its blue lines. There is no additional C3-only gap.

## 7. Singleton padding and the endpoint-hostile family

The remaining singleton-padded problem is U=a -> C -> z with C strong
and mixed. Fix the blue central path P. If C has a black path starting
at the same vertex as P, Section 6 again gives a three-black detour.
The following two families repair every failure of that criterion found
in the complete fixed-blue-path universe through central order nine.

### A cyclic blow-up with a blue middle

**Theorem (PROVED).** Let C have three blocks v,B,w with arcs
v -> B -> w -> v. Let B have a specified blue Hamiltonian path of
length m>=2. For U=a -> C -> z the blue word

    P0=(a,v,B,w,z)

has a three-black detour, applying K along

    (a,w,v,B,z), (B,z,a,w,v), (B,w,z,a,v).             (4)

Here B is retained contiguously forward. The first path is valid because
w->v->B. Its second vertex w and penultimate vertex R=last(B) have
outdegree sum d_B(R)+4 <= m+2, whereas n-1=m+3. The inequality uses
d_B(R)<=m-2, because R has the preceding path vertex as an inneighbor.
Thus this first path is black.

Both source paths share the internal arc v->L, L=first(B). Their
complement partners contain the closing endpoint arc z->a. Rotate their
Hamiltonian cycles at v->L to obtain the last two words in (4). They
lead to the same labelled hub U triangle {az,vL}; this is also checked
directly by path-edge XOR. The endpoint-score differences of these two
moves are -d_B(L) and +d_B(L), since

    d_U(v)+d_U(L)=(m+1)+(d_B(L)+2)=n-1+d_B(L).

Because m>=2 and L begins a Hamiltonian path, d_B(L)>=1. Both moves are
black. These identities prove (4) for arbitrary B, not only for a pure
blue B. If B is pure blue, all central Hamiltonian paths starting at v
are blue: the last vertex must be w, whose only outneighbor is the already
used v, and the remaining path is exactly a path through B. This gives
an infinite genuine hostile to the fixed-start criterion, now bypassed.

### The reverse-transitive-path presentation

Let C_h have path vertices 0,...,h-1, all consecutive arcs i->i+1,
and all other arcs directed from the larger label to the smaller.
Its standard path is blue. For h>=4, the corresponding blue line at
U=a -> C_h -> z also has a three-black detour.

Replace the central path by (1,2,...,h-1,0). Its central endpoint
degrees are both 1, so in U the second and penultimate vertices have
score sum 4, unequal to n-1=h+1. This source path is therefore black.
It shares the arc (h-2)->(h-1) with the original word. Both vertices
have U-score h-1, so their sum differs from n-1 by h-3>0. The inherited
shared-edge rotations are therefore black and give the same hub.
This is another uniform three-move proof. For h=3 the source U is pure
blue and the blackness inequality vanishes, exactly as required.

### Complete finite probe and remaining scope

[moves_endpoint_probe.py](moves_endpoint_probe.py) enumerates every
fixed-path blue matrix of central orders 3 through 9. It enumerates
all Hamiltonian paths from the prescribed first vertex, and when those
are all blue, checks all paths to separate mixed from pure-blue centers.
Every mixed hostile is in exactly one of the two families above:

| central order | all blue matrices | strong matrices | mixed fixed-start hostiles | cyclic-blow-up presentations | reverse-transitive-path presentations |
|---|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | 0 | 0 | 0 |
| 4 | 4 | 3 | 2 | 1 | 1 |
| 5 | 16 | 14 | 3 | 2 | 1 |
| 6 | 64 | 59 | 2 | 1 | 1 |
| 7 | 512 | 495 | 6 | 5 | 1 |
| 8 | 4096 | 4026 | 3 | 2 | 1 |
| 9 | 65536 | 65017 | 7 | 6 | 1 |

Thus every blue line at 1 -> C -> 1 with C strong and mixed of order
at most nine has a three-black detour (FINITE-EXACT, using the proved
uniform constructions). This is a structured family through total order
eleven, not a full metagraph census at those orders. The assertion that
these are all fixed-start hostile families at arbitrary order is OPEN.

### Coverage of every nonstrong mixed source through order eleven

**Corollary (FINITE-EXACT with the stated uniform reductions).** Every
blue line incident to a nonstrong mixed source of total order at most
eleven has a black-only detour.

To check exhaustiveness, let S and S^op be the first and last SCCs.
Their common order is at most five. If it is four or five, use the
uniform boundary theorems above. If it is three, S=C3 and the central
remainder must be mixed, so the mixed-center hub applies. If it is one,
the central remainder is mixed and has order at most nine. A nonstrong
remainder is handled by Section 6; a strong one by the complete probe
and the two uniform repairs in this section. There is no strong
two-vertex SCC. This covers every possible outside SCC size.

This is larger than the complete order-nine black-connectivity census,
but concerns only blue lines with a nonstrong mixed endpoint. It leaves
blue lines with two strong endpoints and all-order connectivity open.

## 8. No extra blue bridge at a mixed paired-black degree-one vertex

The completed [deficit-two center theorem](../deficit_two_centers_20261020/README.md)
proves D(z)=z^4. Combined with the
[paired-black leaf classification](../paired_black_leaf_reduction_20261020/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md),
every mixed leaf of the paired-black quotient is uniquely

    M_W = W -> R4 -> W^op,

where W is a rigid word of singleton/C3 SCCs. This classifies leaves of
the paired-black quotient, not leaves of the whole blue/black graph.

**Corollary (PROVED).** No mixed-to-mixed blue bridge is incident to any
mixed degree-one vertex of the paired-black multigraph. For nonempty W,
every incident blue line
has a black-only detour of length at most three.

For |W|>=2 use the inherited mixed-center hub lemma with center R4.
The only rigid word of length one is a singleton. For W=1 there are
three blue R4 path presentations; with outer singleton labels 0 and 5
and R4 labels 1,...,4 their explicit black K paths are:

| R4 matrix in its four path positions | three successive black words |
|---|---|
| (2,4,9,3) | (0,2,3,4,1,5); (4,1,5,0,2,3); (4,5,0,1,2,3) |
| (6,12,8,1) | (0,2,3,4,1,5); (4,1,5,0,2,3); (4,5,0,1,2,3) |
| (10,4,9,2) | (0,1,4,2,3,5); (1,4,2,3,5,0); (1,2,3,4,5,0) |

These are the shared-edge hub construction, at edge 3->4 in the first
two rows and 0->1 in the third. All six-vertex identities and colors are
checked in [moves_leaf.py](moves_leaf.py). For empty W, direct enumeration
of the order-four metagraph gives one blue edge from R4 to the pure-blue
transitive node, and one blue loop at R4. Neither is a mixed-to-mixed
bridge. Thus the empty-word boundary is also covered.

Here degree one is measured after dividing each literal black-line
multiplicity by two, with loops counted twice. It is equivalent, for a
mixed node, to black presentation deficit beta=2. A vertex of degree one
in this multigraph is also a leaf in its simple support, but the converse
need not hold: a simple-support leaf may carry several parallel pairs.
No classification of all simple-support leaves is used or asserted.

This removes the entire all-order mixed beta=2 family from the search
for a genuine mixed-to-mixed blue bridge. It does not imply black-core
connectivity or exclude a bridge between other regions.

## 9. Score-changing cycles and highly connected components

### A five-black detour from one noncomplementary cycle edge

**Theorem (PROVED).** Let C be a strong mixed tournament of order h.
If C has a Hamiltonian cycle containing an arc xy with

    d_C(x)+d_C(y) != h-1,

then every blue line incident to 1 -> C -> 1 has a black-only detour
of length at most five.

Fix its original blue central path P starting at v. Cut the supplied
Hamiltonian cycle to start at the same v, obtaining Q. If Q is black,
the fixed-first/second/last hub immediately gives the three-step detour.
Suppose Q is blue. Its positional anti-reversal sends an internal arc
xy to a different internal arc with the opposite nonzero endpoint-score
imbalance. The closing edge of Q is balanced because its endpoints are
the endpoints of a blue path; hence the initially supplied imbalanced
cycle edge is internal. A fixed middle arc under reflection would be
balanced, so the two imbalanced internal edges are indeed distinct.

Cut the Hamiltonian cycle at one of those edges to obtain a path Q'.
Its endpoint-score sum is unbalanced, making Q' black. It still shares
the other imbalanced edge with Q. In the singleton-padded tournament
both shared-edge degrees increase by one while the complementary sum
increases by two, so the imbalance is unchanged. The inherited
shared-edge hub gives a three-black detour from the source to the blue
partner of Q. The partners of P and Q themselves have a two-black
connection through the fixed-first/second/last hub. Concatenation gives
the asserted length-five detour.

If every edge of a Hamiltonian cycle is balanced, the degree sequence
around that cycle alternates between two complementary values. At odd
order those two values must coincide. Consequently the theorem covers
every nonregular strong center of odd order, and every even-order strong
center with at least three distinct degrees. The remaining score case
at this stage is an even center with two complementary degree values
and only alternating Hamiltonian cycles. Section 11 below excludes it.

### Hamiltonian connectivity eliminates fixed-start and endpoint rigidity

Call a tournament Hamiltonian-connected when there is a directed
Hamiltonian path for every ordered pair of distinct endpoints.

**Lemma (PROVED).** In a mixed Hamiltonian-connected tournament C, every
vertex v begins some black Hamiltonian path.

Otherwise choose, for every w!=v, a path from v to w; all are blue by
the contrary assumption, with positional anti-reversals r_w. Fix one
w0. Then r_w r_w0 is an automorphism sending w0 to w, so one
Aut(C)-orbit O contains every vertex other than possibly v. An
anti-automorphism normalizes Aut(C), hence maps its vertex orbits to
orbits of equal size. If v lay outside O, its orbit would be a singleton,
but r_w0 maps v to w0, contradicting |O|>=h-1>1. Therefore C is vertex
transitive. Move the initial vertex of any black path to v by an
automorphism; its image is black, a contradiction.

Thus every blue line at 1 -> C -> 1 has a three-black detour when C is
Hamiltonian-connected and mixed. There is also an outer-component
consequence: every blue line at S^op -> Z -> S has a three-black detour
if S is Hamiltonian-connected. For the specified terminal vertex w,
choose a new initial vertex with the opposite adjacency to w from the
old initial vertex. Strongness supplies such a vertex. A Hamiltonian
path between the chosen endpoints exists, and no automorphism fixing w
can identify it with the original ordered path. Apply Section 6.

**CITED inputs.** Thomassen proves that every 4-strong tournament is
Hamiltonian-connected, and that every arc of a 3-strong tournament lies
on a Hamiltonian cycle: C. Thomassen, *Hamiltonian-connected tournaments*,
Journal of Combinatorial Theory B 28 (1980), 142–163,
[DOI 10.1016/0095-8956(80)90061-1](https://doi.org/10.1016/0095-8956(80)90061-1).
These exact imported statements were checked against the publisher's
primary abstract on 2026-10-20. Here k-strong means that deletion of fewer
than k vertices leaves a strong tournament. **Source-access limit:** the
publisher article/PDF returned HTTP 403, and a bounded author-copy/API
search did not recover the original full text. These remain explicitly
CITED inputs; the original proof was not independently audited in this
session. The derived theorems below are proved relative to those exact
published statements, not to a stronger unverified version.

Therefore both preceding detour conclusions hold for 4-strong
components. Any outer-component obstruction not covered by the earlier
small-block theorems must in particular have vertex connectivity at most
three, or fail the stronger Hamiltonian-connectivity condition itself.

### All odd-order singleton-padded strong centers are covered

**Theorem (PROVED using the cited input and the explicit small-order
controls).** If C is strong, mixed, and of odd order, every blue line
at 1 -> C -> 1 has a black detour of length at most five.

The nonregular case was just proved. For regular C of order h, its
vertex connectivity is at least ceil(h/3), by the following elementary
cut calculation. If removing k vertices leaves a nontrivial cut A->B,
write a=|A| and b=|B|, so a+b=h-k. Summing outdegrees over A and
indegrees over B, with every original degree equal to (h-1)/2, gives

    h >= a+2b,     h >= 2a+b.

Adding yields 2h>=3(h-k), or k>=h/3. A regular tournament is strong:
the same calculation with k=0 excludes a separating cut. Thus every
regular tournament of order h>=11 is 4-strong and the cited theorem plus
the preceding lemma gives a three-black detour.

The remaining odd orders are 3,5,7,9. The complete fixed-blue-path probe
in Section 7 checks every regular presentation at these orders. It finds
no mixed regular fixed-start hostile; the regular pure-blue controls
are C3 and C5. These finite checks close the finite boundary, while the
cut calculation and Thomassen's theorem handle every larger order.
This proves the all-odd conclusion without assuming any unproved
classification of the endpoint-hostile families.

This theorem concerns the singleton-padded source 1 -> C -> 1. It does
not identify all strong centers with each other, does not cover every
large outer SCC of connectivity at most three, and does not establish
all-order black-core connectivity.

## 10. The empty-context dissolution boundary has its own four-black ledger

**Theorem (PROVED).** Let W be any nonempty tournament whose last SCC
is a singleton. Then

    R4 -> W   and   C3 -> 1 -> W

are connected by four black moves, restoring every arc of W and every
arc from the displayed four-vertex block to W. Equivalently, after
converse, every tournament A beginning with a singleton admits

    A -> R4   <->   A -> 1 -> C3

by four black moves.

Use B=(6,4,8,3), the same R4 presentation as in the concurrent dissolution
proof, and D=(10,0,3,6), whose first three strong vertices form C3 and
whose singleton sink is vertex 1. Retain a chosen Hamiltonian path W
as a contiguous forward block and apply K along

    (0,1,2,3,W),
    (3,W,2,0,1),
    (3,1,W,2,0),
    (2,0,3,1,W).                                     (5)

For |W|>=2 the six-vertex proxy W=(L,R) checks every path arc and the
literal identity B->W to D->W. Every internal W path edge is retained
four times; every other internal W edge is complemented four times.
Every exterior attachment uses only L or R, so the proxy lifts to an
arbitrary W. For |W|=1 the five-vertex identity is checked separately;
no two-distinct-ports argument is applied to a collapsed port.

The endpoint R=last(W) is its singleton sink, initially of outdegree
zero in the whole source. The initial degrees at vertices 0,1,2,3 are
n-2,n-3,n-3,n-2. Repeated use of equation (1) gives the four endpoint-sum
differences from n-1:

    -1, 4-n, n-3, 1.

They are nonzero for every n=|W|+4>=5, so all four steps are black.
These formulas do not involve any other degree or arc of W.

This specifically closes the awkward boundary of the earlier right-end
dissolution. If A->C3->1 is pure blue, its SCC palindrome forces
A=1 or A=1->C3->V, where V is a pure-blue rigid palindrome. In either
case A begins with a singleton, and the theorem instead reaches the
active target A->1->C3, whose first and last SCCs have different types.
In particular the smallest case A=1 has a literal five-vertex certificate.

The formal proxy, the singleton-port boundary, all endpoint formulas,
and every fixed-path context W with a sink and order at most seven are
checked by [moves_end.py](moves_end.py), with receipts in
[moves_end.json](moves_end.json). This ledger is different from simply
erasing an empty A in a two-context formula: that erasure can destroy
a required path arc. The present words have been verified independently.

## 11. The even score obstruction closes: all singleton-padded mixed centers

**Theorem (PROVED relative to the CITED Thomassen arc theorem, with
complete finite boundary certificates).** Every nonregular strong
self-converse tournament has a Hamiltonian cycle containing an arc xy
whose endpoint degrees satisfy

    d(x)+d(y) != n-1.                                  (6)

Consequently every blue line incident to a mixed source 1 -> C -> 1
has a black-only detour of length at most five, at every order. For a
nonstrong mixed C use Section 6; for strong regular C use Section 9;
for strong nonregular C combine (6) with its five-move construction.

### First reduction: the only possible connectivity is two

Suppose a nonregular strong self-converse T has every Hamiltonian cycle
balanced at every edge. One Hamiltonian cycle then makes the score
sequence alternate between two distinct complementary values delta and
Delta=n-1-delta. Thus n=2m is even, and the low- and high-score classes
L and H each have m vertices. Both have at least two vertices.

An internal edge of L has endpoint sum 2delta!=n-1. If T were 3-strong,
the CITED Thomassen arc theorem would put this edge on a Hamiltonian
cycle, contrary to assumption. Hence kappa(T)<=2.

Nor can T have a cut vertex. Every vertex is extremal in score. If a
minimum-score vertex v has T-v nonstrong, take a Hamiltonian cycle cut
starting at v. Its predecessor u lies in the final strong component B
of T-v. Put t=|B|. If t>1, strongness of B gives d(u)<=t-1, while any
vertex of an earlier component has degree at least t, so d(u)<Delta.
If t=1, u has degree one, still less than Delta for n>=4. Therefore
d(v)+d(u)<n-1, contradicting the balanced-cycle assumption. The
maximum-score case follows by converse. This reproduces the elementary
extremal-deletion mechanism from
[cycles.md, Section 7](../black_descent_20261020c/cycles.md).

All vertices are minimum- or maximum-score, so a cut vertex is
impossible. The supposed residual therefore has kappa(T)=2 exactly.

### A two-vertex cut is incompatible with alternating cycles for n>=10

Choose a two-vertex cut S and a nontrivial partition of T-S into A->B,
obtained by splitting its ordered strong components. Write a=|A| and
b=|B|, so a+b=n-2. Removing S from an alternating Hamiltonian cycle
leaves at most two directed paths. Each such path visits A before B,
because no arc goes from B to A. Hence its A segment is alternating,
and summing at most two such segments gives

    | |H intersect A| - |L intersect A| | <= 2.         (7)

The same inequality holds for B. It is a necessary consequence of the
assumed alternating cycles, not an arbitrary restriction on the cut.

First suppose Delta-delta=1. Then Delta=m. Summing outdegrees over A
and indegrees over B, each bounded above by m, yields

    2m >= a-1+2b,       2m >= b-1+2a.

Thus 4m>=3(a+b)-2=6m-8, forcing m<=4. This excludes n>=10.

Now suppose Delta-delta>=3; the gap is odd, so these are the only two
possibilities. If A contained a low-score vertex and B a high-score
vertex, then their degrees would give

    delta >= b,        Delta <= b+1,

contradicting the gap. Therefore A is entirely high or B is entirely
low. Passing to the converse, which exchanges the two degree classes
and the two sides of the cut, lets us assume B is entirely low.

If A were entirely high too, at least m-2 of the m high vertices would
lie in A. But (7) would give a<=2, contradicting m>=5. Thus A contains
a low-score vertex, and that vertex beats B, so b<=delta. Averaging the
outdegrees of the low vertices of B gives

    b*delta <= b(b-1)/2+2b,
    b >= 2delta-3.

Together these force delta<=3. Strongness supplies at least one arc
from L to H; otherwise H->L is a separating cut. Hence

    m*delta > m(m-1)/2.

For m>=5 this forces delta=3 and m in {5,6}. The two bounds on b then
give b=3. Every low vertex of A spends all three outarcs on B, so it
loses to every other vertex of A and S. There can consequently be at
most one low vertex in A.

Let t=|L intersect S|. Counting low vertices gives t>=m-4>=1. Since B
contains exactly three lows and no highs,

    |H intersect A|=m-2+t,
    |L intersect A|=m-3-t,
    |H intersect A|-|L intersect A|=1+2t>=3.

This contradicts (7). Thus no supposed residual exists at n>=10.

### Complete small even boundary

The only remaining even orders are 4,6,8. The exact checker
[moves_residual.py](moves_residual.py) enumerates every fixed-path blue
presentation at those orders, retains the strong two-score ones, and
finds an explicit Hamiltonian cycle and an edge satisfying (6) in each.
Its receipts are [moves_residual.json](moves_residual.json):

| order | complete blue presentation universe | retained strong two-score presentations | certificates |
|---|---:|---:|---:|
| 4 | 4 | 3 | 3 |
| 6 | 64 | 21 | 21 |
| 8 | 4096 | 537 | 537 |

The fixed-blue universe covers every self-converse isomorphism class,
not merely classes selected by a heuristic. Indeed H(T) is odd by
Redei's path parity theorem, and |Aut(T)| is odd because a tournament
has no nonidentity involutory automorphism. Therefore the set of
Hamiltonian-path orbits has odd size. Converse followed by a chosen
anti-automorphism induces an involution on that orbit set; its square
is an automorphism. It has a fixed orbit. Normalizing a path in such
an orbit gives a positional anti-reversal, hence a blue presentation.
This is also the inherited odd-fibre/blue-fixed-point law.

For every retained matrix the checker explicitly validates the full
cycle, its closing arc, and the claimed degree imbalance. The search
starts cycles at vertex zero without loss, since every Hamiltonian
cycle contains that vertex. Both normal and optimized Python runs pass.
These finite certificates close only the finite boundary; the
two-vertex-cut argument proves the unbounded exclusion.

Odd nonregular orders were already covered in Section 9, so this
finishes (6). Together with the regular-center argument, it closes
every mixed singleton-padded source. The remaining all-order
nonstrong-source frontier concerns larger outer SCCs not covered by
Sections 2–6 or the Hamiltonian-connectivity criterion, rather than
any residual even singleton-padded center. All-order black-core
connectivity itself is still open.
