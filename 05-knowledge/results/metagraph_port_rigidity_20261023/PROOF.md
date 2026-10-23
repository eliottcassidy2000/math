# The missing cut in a tournament's signed curl

2026-10-23. **PROVED, elementary:** labelled reconstruction and port transport,
with the proof below. **FINITE-EXACT:** all labelled tournaments of orders
2 through 6, all fixed-path ports of orders 3 through 7. **CITED:** the broader
orientation/connected-subgraph connection is established Tutte-polynomial
theory. This is not an all-order classification from scalar bucket counts,
a new A000568 algorithm, or a claim of literature priority.

## 0. Inheritance and the distinction that matters

The closest mechanism is [THM-4739, Kirchhoff and Hodge
laws](../../../01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md):
partner scores plus the per-vertex cyclic-triangle profile, attached to
actual source port types and endpoint-count constraints, determine the
metagraph at order 5 but leave alternatives at order 6. Node profiles alone
are not asserted sufficient. The canonical hostile is the
48,000 L3 transportation solutions at order 6; it is not an enumeration of
48,000 distinct merged graphs. The corrected near miss is confusing a
signed alternating triangle curl with the unsigned set of cyclic triples.
The least-used sidecar here is one distinguished vertex's score, or the
actual directed Hamiltonian path with its position labels.

Live board: signed curl / the SCC cut / one scalar sidecar / authenticated
ports / cyclic versus linear component order / operation-relative memory.
Anchor: exact port acceptance. Niche: the lost SCC cut. Wildcard: connected
labelled graphs. The active research move is controlled forgetting with an
explicit restoration sidecar.

The fixed observer and path convention are the lossless ordered-graph
coordinates of [the marked trijection
packet](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md).
No unsigned graph invariant is silently treated as a tournament invariant.

## 1. Exact curl fibres

Let V be a fixed labelled set, |V|=n, with n>=2. Write a_ij=1 when i beats j,
a_ji=-a_ij, and a_ii=0. Define the alternating signed curl

    c(i,j,k)=a_ij+a_jk+a_ki in {-3,-1,1,3}.

Fix a root r. It suffices to retain the root triangles
kappa_ij=c(r,i,j): every other triangle satisfies

    c(i,j,k)=kappa_ij+kappa_jk-kappa_ik.                 (1)

**Theorem 1 (exact lost coordinate).** Two labelled tournaments have the same
root curls iff one is obtained from the other by cyclically rotating the
linearly ordered list of its strong components, retaining each labelled
component and its internal arcs. A tournament with k strong components
therefore has exactly k labelled tournaments in its curl fibre.

**Proof.** Equal root curls imply equal full curls by (1). Let b be the second
arc array and d=b-a. Its triangle sums vanish. Put u_i=d_ir, u_r=0; then

    d_ij=u_i-u_j.

Every difference is 0 or +/-2. Consequently u has at most two values: three
different values would have pairwise difference at least 2 and extreme
difference at least 4. If all values agree, a=b. Otherwise the two level sets
form a nonempty proper cut. Across this cut d is +/-2, so a and b have
opposite arcs and all a-arcs point uniformly from one part to the other.
Within each part the arrays agree.

A uniformly directed cut of a tournament lies between two consecutive
strong components. Reversing it moves the initial block of components to
the end. Conversely the difference caused by that reversal is a gradient,
so every triangle curl is unchanged. There are exactly k places to cut the
cyclic list. The k resulting labelled arrays are distinct, because each
nonempty proper component prefix reverses some arcs. QED.

In particular, full signed curl determines a **strong labelled tournament**
without any score sidecar. This is much stronger than an energy or profile,
and much more information is retained.

## 2. A single score restores the cut

Let s_T(r) be the outdegree of the distinguished root.

**Theorem 2 (one-score reconstruction).** The signature

    S_r(T)=( (c(r,i,j))_{i,j != r}, s_T(r) )             (2)

is injective on all labelled tournaments on V. The curls are a labelled
alternating tensor, not a sorted list. Relabelling moves the root, the
tensor indices, and the score together.

**Proof.** If a nontrivial uniformly directed cut S -> V\S is reversed,
the root loses n-|S| wins when r is in S, and gains |S| wins otherwise.
Neither change is zero. By Theorem 1 every distinct member of the curl
fibre arises this way. Hence its root score is different. QED.

There are k possible score values in a k-component fibre. At least k
sidecar states are necessary to distinguish every member of that fibre;
the one score attains that local state count. This does not assert that
the displayed tensor-plus-score serialization uses fewer bits than an
adjacency mask.

### An explicit polynomial-time decoder, including rejection

Set r=0 and x_i=a_0i for i>0. Then

    a_ij=kappa_ij-x_i+x_j.                             (3)

For kappa_ij=1 the condition a_ij in {-1,1} is exactly x_i>=x_j;
for kappa_ij=-1 it is x_i<=x_j. The value 3 forces (x_i,x_j)=(1,-1),
and -3 forces (-1,1). Form a tournament on V\{0} by orienting i->j
when kappa_ij>0. Every feasible x is constant on each strong component
of this sign tournament, and is nonincreasing along its component chain.
It therefore comes from one of at most n prefix cuts: +1 on the prefix,
-1 on the suffix. Test these cuts against (3), the +/-3 constraints, and

    #{i:x_i=1}=s_T(0).

Every accepted candidate is literally a tournament with the requested
signature, and every realization appears. Theorem 2 says at most one
candidate passes the score test. This is a necessary-and-sufficient
realizability criterion, not just an invariant filter. The implementation
uses elementary transitive closure; SCC algorithms can reduce its cost.

The local alphabet alone is insufficient: at n=4 the three root curls
(3,3,3) are impossible. They force x_2=-1 from kappa_12=3 but x_2=1
from kappa_23=3. Thus actual arc feasibility provides the pigeonhole-type
constraint absent from the scalar transportation margins.

### A directed path is another exact sidecar

On the fixed path P=(0,1,...,n-1), all a_i,i+1=1. Starting with a_01=1,

    a_0,i+1=a_0i+1-kappa_i,i+1                         (4)

recovers the root star, and (3) recovers the remaining arcs. Hence root
curls alone are injective on the fixed-path tiling cube. Equivalently,
a nontrivial cut reversal cannot preserve a Hamiltonian path: the path
must cross the cut and its crossing arc is reversed.

## 3. Complete metagraph port acceptance

Keep an actual directed Hamiltonian path P of T. Let p_P be its signed
edge flow: +1 on directed path edges, -1 on their reverses, zero elsewhere.
The partner operation K_P retains the path arcs and reverses every other
arc, so

    a'= -a+2p_P,
    c'= -c+2 curl(p_P),
    s'(r)=n-1-s(r)+div(p_P)(r).                        (5)

Here div(p_P) is +1 at the start, -1 at the end, and zero at every internal
vertex. It follows that S_r(K_P(T)) is computable directly from S_r(T)
and the labelled path. The target signature is accepted only if the
decoder reconstructs its actual tournament. The directedness of P is a
separate input guard, checked on the decoded source.

**Corollary 3 (all-order iff port test).** Normalize two Hamiltonian ports
(T,P) and (U,Q) to the same position labels 0,...,n-1. They are complement
partners iff

    S_0(U along Q)=transport_P(S_0(T along P)).         (6)

Necessity follows from (5); sufficiency follows from Theorem 2. With the
path normalization, equality of root curls already suffices by (4).
Two paths in one isomorphism class have the same normalized tournament
iff an automorphism carries one path to the other. Thus (6) matches
Hamiltonian-path **orbits**, not merely vertex-score types.

For n>=3 and a complete inventory of actual normalized ports, the compatibility
relation (6) is a fixed-point-free involution. Each closed port has its
one actual open partner and conversely. Grouping by tournament class,
then by converse pairs, recovers the unique merged multigraph while
retaining each distinct port and its multiplicity. The colour is blue
iff the path reversal is an anti-automorphism of the decoded source;
otherwise black. This agrees with the grid reflection predicate and is
preserved by K_P. Labels may be changed simultaneously, but cannot be
discarded before (6) is checked.

This is an all-order **determination with complete port data**, not a proof
that node totals, even/odd fibre sizes, degrees, electrical potentials,
or a fixed low-order profile determine the metagraph. Any counterfeit
that changes the attachment of even one authenticated port fails (6).
The existence of the invariant formula does not compute the full class
inventory or solve black-core connectivity.

The restriction n>=3 in the pairing statement is essential: at n=2 the
free-edge set is empty and K is the identity. The transport and iff
reconstruction formulas remain valid there, but no fixed-point-free
closed/open pairing is asserted.

## 4. Why a cut can preserve tiling count and change bucket parity

For an ordinal sum of strong components, every Hamiltonian path traverses
the components in their linear order, and every automorphism preserves
each component set. Therefore

    H(T)=product H(C_i),   |Aut(T)|=product |Aut(C_i)|.

Cyclic rotations preserve both products and h=H/|Aut|. But self-converseness
and blue-port incidence can change.

For a concrete labelled n=5 witness, take C3 on {0,1,2} and singletons a,b.
The three arrays C3 -> a -> b, a -> b -> C3, and b -> C3 -> a have the same
signed curl. Each has H=3, |Aut|=3, and h=1. The middle-C3 array is
self-converse; the other two classes form a non-self-converse pair.
Consequently the same curl fibre meets an odd SC node fibre and an even
converse-merged node fibre. Identifying curl fibres with merged nodes is
impossible even though the unmerged tiling counts agree.

This supplies a precise version of the user's bucket question: a scalar
count does not remember where the linear order was cut. The missing cut
affects which isomorphism and converse bucket owns the port.

### Exact even/odd law on a cyclic component necklace

Now forget the vertex labels but retain the cyclic word w of **strong
component isomorphism types**, with their converse involution. Let rho
rotate the word one place, and let d be its least positive rotation
period. Exactly d distinct tournament classes occur among its linear
cuts. Let iota reverse the word and replace each component by its
converse. Unique ordered SCC decomposition proves both assertions and
the relation iota rho = rho^(-1) iota.

If iota(w) is not a rotation of w, this necklace and its distinct converse
necklace contribute d non-self-converse pairs and no SC class. Otherwise
write iota(w)=rho^a(w), with a taken modulo d. On the d linear-cut classes,

    iota(rho^t(w))=rho^(a-t)(w).

Consequently their SC classes are exactly the solutions of

    2t=a mod d.                                       (7)

There is precisely one when d is odd; when d is even there are two if a
is even and none if a is odd. If s is this SC count, the self-converse
necklace contributes (d+s)/2 converse-merged nodes: s singleton classes
with fibre h, and (d-s)/2 paired classes with fibre 2h. All have the same
unmerged h=product h(C_i), which is odd by the inherited tiling fibre law.
This is a structural source of odd/even buckets. It is scoped to one SCC
necklace, not a classification of all ports or blue-line multiplicities.

For self-converse component letters A,B, AAB illustrates one SC cut,
ABAB (least period 2) illustrates no SC cut, and AABB illustrates two.
These can be realized with A=K1 and B=C3. The n=5 witness above is the
AAB case. The checker also exhausts cyclic words through length 8 on
K1, a strong non-self-converse six-tournament T, and its converse, so the
nontrivial component-converse involution is included. Explicitly T has
upper-triangular positive arcs only 0->2 and 0->5; all its other arcs
point from larger labels to smaller labels. Its Hamiltonian cycle is
0->5->4->3->2->1->0, and the checker tests all 720 permutations to exclude
an anti-automorphism. Its mask in check.py is 18.

Audit repair: an initial draft used K1->C3 as an SCC letter. It is not
strong, so that realization was invalid even though the abstract dihedral
word test was valid. The strong chiral six-tournament replaces it; neither
the abstract law nor the necklace counts change.

## 5. Connected graphs: a real bridge with known prior art

Choose a fixed labelled root r. Every curl fibre has exactly one rotation
whose first strong component contains r. That representative is exactly
the tournament in the fibre from which r reaches every vertex. Hence

    signed-curl fibres <-> root-connected tournaments on K_n.            (8)

Let S(z) be the EGF for labelled strong tournaments, and
T(z)=sum_n 2^binom(n,2) z^n/n!, including the empty tournament. Strong
components form a linear sequence, so T=1/(1-S). A curl fibre is a cyclic
sequence of nonempty labelled strong components, giving EGF

    C(z)=sum_(k>=1) S(z)^k/k = log T(z).                (9)

All labelled simple graphs have the same EGF T, and decompose as sets of
connected components. Their connected EGF is also log T. Thus curl fibres
are equinumerous with connected labelled simple graphs, with counts

    n=1..6: 1,1,4,38,728,26704.

The general bridge is **known**, not a novel equinumerosity claim:
[Bernardi, *Tutte polynomial, subgraphs, orientations and sandpile model:
new connections via embeddings*](https://arxiv.org/abs/math/0612003),
Electronic Journal of Combinatorics 15 (2008), R109, gives an embedding-based
bijection between connected spanning subgraphs and root-connected
orientations, counted by T_G(1,2). Applying that cited construction to K_n
and composing with (8) provides a specified route to a bijection after
the root and embedding are chosen. The present checker does not implement
Bernardi's map, and no unlabelled or observer-free graph bijection is claimed.
[Backman's partial-orientation paper](https://arxiv.org/abs/1408.3962)
places cut-reversal quotients in the broader Tutte-polynomial theory.

## 6. Concurrent memory results and the exact scope of compression

The incoming [session33 memory
note](../session33_continuation_states_20261023/metagraph-memory/NOTE.md)
proves that for the sole operation K, node-observation memory needs only
(source node, target node), adding colour at mixed supports. Its
[full-flip addendum](../session33_continuation_states_20261023/metagraph-memory/ADDENDUM_FULL_FLIPS.md)
proves that allowing all individually labelled free-edge flips forces all
2^binom(n-1,2) tiling states. These are compatible with this packet:
(2) and (4) are lossless static coordinates, not smaller continuation
quotients. Formula (5) is an explicit update in those coordinates.

Source -> target: labelled tournament with marked root/path -> signed curl
plus one root score, or signed curl inside the fixed-path cube.
Preserved: all arcs, actual Hamiltonian ports, target, converse and colour.
Lost without the sidecar: precisely the cyclic SCC cut; after an unsigned
or sorted projection, much more is lost. Cheapest decisive test: run the
decoder and compare the transported target signature, not aggregate counts.

## 7. Exact reproduction and remaining targets

Run from the repository root:

    python -B 05-knowledge/results/metagraph_port_rigidity_20261023/check.py --out 05-knowledge/results/metagraph_port_rigidity_20261023/receipt.json
    python -O -B 05-knowledge/results/metagraph_port_rigidity_20261023/check.py

The checker imports no other producer, uses integer arrays throughout,
and explicit require/raise checks survive -O. The committed receipt records
the complete universes, checks, controls and source hash. It compares
brute-force curl fibres with both the SCC-rotation theorem and an independent
sign-SCC decoder; the fixed-path census checks affine transport, inverse,
closed/open opposition and colour preservation.

Next precise question: can a small portion of the root-curl tensor be
selected adaptively per ambiguous port to separate its surviving target
classes? The current theorem supplies an exact upper endpoint and a minimal
SCC-cut sidecar, not a uniform small selection bound. The bucket lane's
separate order-six scalar signatures can be tested against this complete
acceptance predicate as compressed certificates.
