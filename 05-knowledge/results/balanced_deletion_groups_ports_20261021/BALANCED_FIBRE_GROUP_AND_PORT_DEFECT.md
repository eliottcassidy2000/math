# The balanced deletion fibres are Klein-four orbits, but not metagraph fibres

2026-10-21, session20. Source checkout: `32579ff9c2f4`.

**PROVED:** literal insertion-toggle action; exact orbit/stabilizer types;
reflection semidirect product; an all-order two-port correction for path
reversal; sharp counts of correction sizes; obstruction to descending the
actual metagraph incidence through the undecorated deletion fibres.
**FINITE-EXACT:** independent controls in `check_balanced_group_ports.py`.
**NOT CLAIMED:** an arithmetic class-group identification follows from the
numbers 4,2,1; these results classify the unrestricted tournament metagraph;
path-set permutations are tournament automorphisms or metagraph moves.

## 0. Inheritance, hostile, and scope

The actual object and product decoder are inherited from
`05-knowledge/results/strong_blue_games_20261020f/double_return_counts.md`.
The pointed-set square, complement/swap involutions, literal centre-flanking pair deletion,
4/2/1 fibres, and scalar recurrences are inherited from
`05-knowledge/results/zero_heavy_ports_20261020g/microcosm.md`, Sections2–5.
Arbitrary interior thinning and the full endpoint-growth kernel are inherited
from `05-knowledge/results/balanced_fibre_endpoint_extensions_20261020/`
`BALANCED_DELETION_EXTENSIONS.md`, Sections1–4.

The relevant correction is the 2026-10-20 entry at the front of
`01-canon/MISTAKES.md`: an abstract permutation need not be a legal or black
Hamiltonian move. The newer marked-trijection correction retains the vertex
relabeling when literal path words move between isomorphic representatives.
The established homogeneous-twin hostile in the endpoint extension note
shows that even deleting one of a homogeneous pair need not take every
Hamiltonian path to a Hamiltonian path. We use the specific decoder below,
never an unqualified deletion principle. Reserved THM-4739/4740 are not used.

The least-used sidecar in this lane is the pair of surviving neighbours of
an inserted vertex, rather than the endpoints of the whole path. The board
is: pointed Boolean products, stabilizer subgroups, blue reflection,
conductor flags, and actual path-reversal incidence. The anchor is the
literal 4/2/1 connection; the niche is the forgotten cut-edge data; the
wildcard is the resulting faithful dihedral action.

## 1. Exact family and deletion map

For a,b>=2, T_(a,b) has vertices0,...,a+b. Every arc points upwards except
for a->0 and (a+b)->a. Put T_k=T_(k,k), with k>=2. Thus the central vertex is
k. Set A={0,...,k-1}, B={k+1,...,2k}, and sigma(v)=2k-v.

Let X_k=P({1,...,k-1}) disjoint-union {infinity}; write q=2^(k-1)+1.
Interpret an infinity coordinate as the full set A. For I,K in X_k set
J=sigma(A minus K). The unique Hamiltonian path with these coordinates is

    P_k(I,K) = I increasing, J increasing, k,
               (A minus I) increasing, (B minus J) increasing.       (1)

The decoder follows because removing the centre leaves a transitive
tournament. The portions before and after the centre must each increase.
Any A vertices after the centre must start with0, and any B vertices before
it must end with2k. These conditions give precisely the displayed states,
and they are sufficient for every joining arc. The infinity states are
not conventions that can be discarded: they are the exceptional choices
I=A and J=empty. In particular I=infinity iff the path starts at0, and
K=infinity iff it ends at2k.

Embed T_k in T_(k+1) by

    e(v)=v if v<k;  e(k)=k+1;  e(v)=v+2 if v>k.

The newly inserted vertices are L=k and R=k+2; the new centre is k+1.
Delete L and R and relabel increasingly. This is deletion of the two
vertices adjacent to the centre, **not deletion of the central vertex**.
It acts in (1) by forgetting the new bit k in each finite coordinate and
sending infinity to infinity. Denote the resulting HP map by pi_k.
The inherited decoder proves that it is an actual surjection of path sets.

## 2. New: the literal V4 action and all stabilizers

On X_(k+1), define t(S)=S symmetric-difference {k} for finite S, and
t(infinity)=infinity. On child paths, transport the coordinate actions

    tau_L(I,K)=(t(I),K),   tau_R(I,K)=(I,t(K)).              (2)

Both images are actual directed Hamiltonian paths by (1). Each generator
has order2, and the two commute. They act faithfully, since a finite/finite
state is moved by either nonidentity generator and their product. Let
G=<tau_L,tau_R>, identified with F2^2.

**Theorem1.** Every pi_k fibre is exactly one G-orbit. Its stabilizer is

    parent finite/finite:       {1};
    parent infinity/finite:     <tau_L>;
    parent finite/infinity:     <tau_R>;
    parent infinity/infinity:   G.                         (3)

**Proof.** Every finite coordinate has exactly the two lifts obtained by
choosing whether k belongs to it. An infinity coordinate has its single
infinity lift. The generators independently change exactly these choices,
so they act transitively on each fibre and fix exactly the coordinates
listed in (3). No other elements can fix a finite coordinate. QED.

Thus the orbit sizes4,2,1 are caused by two commuting binary insertions
with independent fixed boundaries. More precisely the complete G-set of
child paths is the disjoint union

    (q-1)^2 copies of G/{1},
    (q-1) copies of G/<tau_L>,
    (q-1) copies of G/<tau_R>,
    one copy of G/G.                                      (4)

This is a structure theorem, not just an orbit-size coincidence. In
particular both different size2 stabilizers really occur. A G-equivariant
isomorphism preserves point stabilizers because G is abelian. Therefore no
single G-set model using just one order-two subgroup for every size2 orbit
can represent (4), even after one global automorphism of G. Such an
automorphism permutes the two distinct subgroups, but cannot identify them.
Any comparison with a single arithmetic conductor flag must retain the
shore label and specify whether it compares kernels, quotients, or torsors.
It must not silently replace one type with another.

A natural marked section chooses the new bit0 in every finite coordinate,
and the unique lift at infinity. Relative to that section every fibre has
explicit quotient coordinates G/Stab. The section is invariant under shore
swap but not finite-set complementation. Without a chosen lift, the free
fibres are naturally G-torsors, not groups with an intrinsic identity.

For d reflected additions, the same proof gives the action F2^(2d): the d
new left bits and d new right bits toggle independently. Stabilizers are the
whole left d-space, the whole right d-space, zero, or the whole group.
The all-depth groups have exponent2; any arithmetic comparison retaining a
different group exponent at a later stage must stop or change its target.

## 3. Reflection, blue paths, and two different Klein-four groups

The actual anti-reversal is R(P)=sigma(reverse(P)). Formula(1) gives
R(I,K)=(K,I). Since T_k is rigid and sigma is its unique anti-automorphism,
a path is blue exactly when I=K. Rigidity follows from the score sequence:
only the two extreme adjacent pairs have equal scores, and their internal
arcs prevent either pair from being exchanged.

The inherited involution J complements each finite subset and fixes
infinity in both coordinates. It commutes with R, and <R,J> is an inherited
Klein-four group with a Pythagorean orbit interpretation. It is **not** the
new vertical group <tau_L,tau_R>. In particular

    R tau_L R = tau_R,    R tau_R R = tau_L.               (5)

The faithful group generated by R,tau_L,tau_R is
G semidirect-product C2, with the C2 swapping the two axes: the dihedral
group of order8. The action is faithful since G is faithful and R changes
some parent states whereas every element of G fixes all parent states.

Neither individual toggle preserves the blue diagonal: a finite diagonal
lift with equal new bits becomes off-diagonal after one toggle. In each
four-element fibre above a finite diagonal parent, two lifts are blue and
two black. A mixed fibre is entirely black. The infinity/infinity fibre is
one blue path. In contrast the simultaneous toggle tau_L tau_R does preserve
the diagonal. The inherited J also preserves blue but usually moves the
parent state, rather than remaining within a pi_k fibre.

This recovers the inherited recurrence from the actual G-set decomposition:

    H'=4H-4q+1,      q'=2q-1,
    Delta'=4Delta-2q+2,                                  (6)

where H=q^2 and Delta=H-q. A black parent has only black lifts, but a black
child can delete to a blue parent. Thus blue is not constant on all G-orbits,
and no blue-labelled quotient can discard the bits deciding the diagonal.

## 4. New: the exact cut-port correction for metagraph moves

Write E(P) for the set of unordered edges of a Hamiltonian path P, and
T^F for the tournament obtained by reversing the arcs on the unordered
edge set F. The repository's R move is T->T^E(P); its K move is the converse
of this target. Let D delete a prescribed set of vertices. Suppose deleting
those vertices from P gives a directed Hamiltonian path p of D(T).
Let C_D(P) be the set of edges of p whose two endpoints were separated in
P by at least one deleted vertex. Equivalently, each maximal nonempty block
of deleted vertices with a surviving vertex on each side contributes the
unordered pair of these two surviving neighbours. Deleted blocks at the
start or end contribute nothing.

**Theorem2 (exact defect identity).** With consistent relabeling,

    D(T^E(P)) = D(T)^(E(p) symmetric-difference C_D(P)).    (7)

The ordinary deletion/reversal square commutes if and only if C_D(P) is
empty. For deletion of r vertices, |C_D(P)|<=r. The identical correction
applies to K moves, since taking converse commutes with vertex deletion.

**Proof.** The only reversed edges visible after deletion are original
consecutive edges of P with both endpoints surviving. They are exactly the
edges of p except the newly created bypass edges C_D(P). Since C_D(P) is
contained in E(p), set subtraction equals the symmetric difference in(7).
An edge reversal changes the labelled tournament, so equality with the
uncorrected target holds exactly when the correction is empty. Each
contributing deleted block contains at least one of the r deleted vertices.
QED.

For our balanced deletion, |C|<=2. The sidecar is at most two **actual port
pairs**, not merely their count or the start/end of the entire path. It is
enough to restore the target of the move exactly. The full child path gives
these ports, but the fibre label, orbit size, and arithmetic group element
without the parent path/run boundaries do not.

More concretely the inserted L goes at the end of the pre-centre A run if
its bit is1, and at the end of the post-centre A run if its bit is0; infinity
forces the pre-centre choice. The inserted R goes at the start of the
pre-centre B run if its bit is0, and at the start of the post-centre B run
if its bit is1; infinity forces the post-centre choice. They can occupy the
same seam, so two inserted vertices need not create two correction edges.
Retaining the surviving predecessor and successor of each resulting block
is the precise insertion-port description.

There is a smaller exact interface than retaining the whole parent word.
Use START and END as formal endpoint markers and define the two parent seams

    pre = (last I if nonempty else START,
           first J if nonempty else centre),
    post = (last(A minus I) if nonempty else centre,
            first(B minus J) if nonempty else END).       (7a)

For this paragraph interpret an infinity coordinate as forcing its new bit
1. The pre seam receives an inserted vertex exactly when epsilon=1 or
eta=0; the post seam receives one exactly when epsilon=0 or eta=1. If both
inserted vertices occupy one seam they appear in the order L,R. Hence C is
exactly the set of occupied seams in(7a) having two genuine vertex endpoints.
This follows directly by inserting L at the end of an A run and R at the
start of a B run; every other old adjacency survives. Thus the two port
pairs together with the new bits reconstruct the correction. In particular
bit patterns00 and11 occupy both seams,01 occupies only post, and10 only
pre, subject to the forced infinity bits. This is an exact positive repair
of deleted-target reconstruction, not just a counterexample to the
undecorated version. It does not say that the reconstructed target is a
legal one-move labelled R or K neighbour of the parent tournament.

### A sharp finite hostile, already at T3 -> T2

Delete vertices2 and4 from T3 and relabel the surviving centre3 as2.
The two child states (empty,{1}) and ({2},{1}) are exchanged by tau_L.
Their literal paths are

    P00=(4,6,3,0,1,2,5),
    P10=(2,4,6,3,0,1,5).

Both are black, and both delete to p=(4,2,0,1,3). Their cut sidecars are
respectively {(1,3)} and empty. After reversing their respective child paths
and then deleting the same vertices, the score multisets are respectively

    (0,1,2,3,4),   (0,2,2,3,3).                          (8)

The first is transitive and the second is not. They are not isomorphic,
and not converse-isomorphic. Therefore even converse-merging does not make
this orbit-fibre quotient compatible with the move targets. This does not
prove that no other decorated or entirely different metagraph map exists;
it refutes the specific claim that the undecorated HP deletion descends the
actual move-incidence via the same vertex deletion.

There is an all-order single-fibre hostile as well. The unique lift of the
standard path is the standard increasing path. Both deleted vertices are
internal and separated by the retained centre, so its correction consists
of the two parent edges {(k-1,k),(k,k+1)}. Even the singleton fibre is not by
itself an incidence-preserving deletion.

## 5. New: the full all-order defect census without large enumeration

**Theorem3.** Put q=2^(k-1)+1. Among the (2q-1)^2 child Hamiltonian paths,
the numbers with |C|=0,1,2 are respectively

    N0=2q,
    N1=2q^2-4,
    N2=2q^2-6q+5.                                      (9)

**Proof.** Write s=q-1=2^(k-1). Each child coordinate has2s finite states
and one infinity state. Let e be the number of L,R occurring as endpoints
of the child path, and a be1 when L and R are consecutive,0 otherwise.
Counting internal deleted blocks gives |C|=2-e-a.

There are2s(s+1) paths with a=1: the two adjacency possibilities are
pre-centre L,R and post-centre L,R. In the first case the left coordinate
contains its new bit (s finite choices plus infinity), while the right
coordinate is finite and excludes its new bit (s choices). The second
case reverses these roles.

The endpoint incidences number6s+2 in total: L is first in2s+1 paths
(left coordinate exactly{L}); R is first in s paths (left coordinate empty
and right bit0); R is last in2s+1 paths (right coordinate exactly{L} in
common reflected coordinates); L is last in s paths (right coordinate
empty and left bit0). Exactly two paths have e=2, namely both coordinates
empty and both coordinates the new-bit singleton. Hence6s paths have at
least one deleted endpoint. There are exactly2s paths with both a=1 and
e>=1, obtained by the two adjacent blocks touching an endpoint; here e=1.

Consequently |C|=0 occurs in2+2s cases. By inclusion-exclusion, |C|=2
occurs in

    (2s+1)^2 - 6s - 2s(s+1) + 2s = 2s^2-2s+1

cases. The remaining count is2s^2+4s-2. Substitution s=q-1 gives(9).
QED.

Thus the proportion of child moves commuting literally with central
deletion is2q/(2q-1)^2, tending to0. This is a port-loss phenomenon, not a
failure of the exact product decoder or the Klein-four orbit statement.

## 6. Connection contracts and consequence

Path source -> pointed product: (1) is a bijection preserving every path
word. Counting only q^2 destroys the word and its actual edges. Keep the
coordinates and the tournament construction to recover them.

Child HP -> parent HP: pi_k is the actual centre-flanking pair deletion. It preserves
being a directed HP; its complete fibres are the G-orbits(3). It loses the
new bits and can lose blackness. Keep the boundary type, the bit pair, and
blue diagonal predicate when those matter.

Orbit types -> an arithmetic V4 carrier: the strongest automatic bridge is
an abstract marked group/quotient/torsor description, not an equality of
number-theoretic constructions. Both distinct stabilizer lines in(3), the
origin of a torsor, and the shore-swap action are mandatory sidecars. A
single distinguished arithmetic subgroup cannot represent both shores
under one globally fixed equivariance. A genuine arithmetic comparison
must specify the two characters or quotient maps and their provenance.

Metagraph incidence -> smaller incidence: naive deletion fails by(8).
Equation(7) is the exact repaired interface. Keep C_D(P) as actual labelled
port pairs; the bounded-size correction generalizes to any deletion for
which the path remains Hamiltonian. This is useful boundary information for
a metagraph classification, but it is not a complete classification.

## 7. Reproduction and limits

    python3 check_balanced_group_ports.py > check.out
    python3 -O check_balanced_group_ports.py > check.optimized.out
    cmp check.out check.optimized.out

The checker imports no repository producer. It compares the independent
literal HP DFS with the product decoder at source orders7,9,11; checks every
coded child state at source orders7,9,11,13,15; tests actual directedness,
blue anti-involutions, group composition, orbit equality, all stabilizers,
faithful eight-element dihedral action, the three independently calculated
cut sidecars (edge subtraction, deleted-block scan, and two-seam formula), and the literal corrected target tournament. The frozen
universe has5,709 source paths and91,344 group-composition checks. It retains
the explicit hostile and its different score spectra. All failures use
exceptions, so optimization cannot remove validation.

The all-order proofs do not depend on the finite controls. No literature
priority claim is made. Parent publication owns all repository writes.
