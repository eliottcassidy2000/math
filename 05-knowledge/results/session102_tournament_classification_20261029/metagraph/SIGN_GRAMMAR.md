# Native signs, switching obstructions, and SCC composition

Session102, 2026-10-29. **PROVED elementary statements; FINITE-EXACT controls.**
Source pin: `3b6e3ef82382da222fc59aa02914eea72d5a105b`. This is a local operation
and composition grammar, not a classification of all prime tournaments or the
merged metagraph. No literature-priority claim is made.

## 0. Inheritance and types

Closest mechanisms: THM-643, `01-canon/theorems/THM-643-blue-black-line-parity-structure.md`
(current correction governs); Session93 `metagraph/BLUE_QUOTIENT_FLATNESS.md` in
`05-knowledge/results/session93_reversal_lucas_clocks_20261029/`; the inherited
SCC center law in §3 of `05-knowledge/results/merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md`.
The last file's historical atom-classification status is not imported here.
**Provenance recovery:** §8 Theorem3 of
`05-knowledge/results/tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md`
already proves the cut-Cayley quotient and injectivity of its path generators
for n>=8, using the same degree/endpoint argument below. Our two-move
obstruction is its legal-word corollary, not a new discovery. This packet
refines the sharp n7 boundary by classifying the legal skeleton and colours,
and extends the degree obstruction to arbitrary native-word lengths.
See `CAYLEY_QUOTIENT.md` for the recovered theorem and fresh verification.
Hostiles: `C3 -> C3`, `C3 -> 1 -> C3`, and the inherited 81/27 decorated-C3
pair. Corrected near miss: unguarded sign covariance is not legality of a
native move. Least-used sidecar: the *actual edge set* of each named path.

Board: native path guards / cut switching / word parity / SCC reflection /
boundary endpoint ports / ordinary versus jointly decorated orbits.
Anchor is native-word obstruction; niche is exact SCC colour counting;
wildcard is the sharp bipartite two-path construction. The new connection is
from legal native words to parity subgraphs, preserving the final labelled
orientation and forgetting chronology unless intermediate guards are retained.

A tournament has vertices V and exactly one directed arc on each unordered
pair. Use a fixed vertex order for the sign gauge. Let J=E(K_V), and let E(P)
be the unordered edge set of an **ordered directed Hamilton path** P.
Writing orientations as bits x, a native move is

    K_P(x) = x XOR (J \ E(P)).                         (1)

Its guard is that P is directed in x. It remains directed in K_P(x). A native
presentation is blue iff its positional reversal rho_P is an anti-automorphism:

    A[p_i,p_j] = A[p_(n-1-j),p_(n-1-i)]  for i<j.       (2)

Otherwise it is black. This is the repository grid-symmetry colour, not a
colour assigned from a switching class or a sign of a permutation. Native
lines keep path presentations; merging to class support loses multiplicities.

## 1. Exact grammar and the guard obstruction

Vertex switching by S toggles the cut delta(S). Converse toggles J. Thus
switching and K_P commute as **unguarded bit maps**, and relabelling transports
K_P together with P. Converse transports the directed path to reverse(P):

    K_reverse(P)(T^op) = K_P(T)^op.                    (3)

But if P is directed in both T and its S-switch, no edge of E(P) crosses S.
Since the path spans V, S is empty or all V. Therefore:

**No nontrivial switch preserves a named Hamilton-path guard.** An algebraic
switching covariance cannot authorize the same native move after switching.
Converse is different: reverse(P) is the necessary transported guard.

**A single native move is never a switch when n>=3.** On three consecutive
positions of P, J\E(P) has exactly one edge, whereas every cut intersects a
triangle in zero or two edges. This is a labelled switching-orbit obstruction.
It does not exclude an isomorphic relabelling of a switched tournament.

**A single native move is never converse followed by a switch for n>=4.**
That would require E(P) to be a cut. A nonempty cut is complete bipartite, and
P_n is not complete bipartite for n>=4: if it has n-1 cut edges, the equation
s(n-s)=n-1 forces s=1 or n-1, giving a star, not a path. At n=3 it *is* a cut,
and K_P is converse plus switching the middle vertex. This is the sharp
small-order exception, independently pointed out in the session audit.

## 2. A sharp two-native-move switching threshold

For a **legal** word T_(j+1)=K_(P_j)(T_j), repeated XOR gives

    T_r = T_0 XOR ((r mod 2) J) XOR XOR_j E(P_j).       (4)

This is an exact endpoint formula, not a criterion that a prescribed path
word is legal. Legality must still be checked at each intermediate tournament.
For r=2 the change is E(P) symmetric-difference E(Q).

**Inherited injectivity corollary; refined sharp boundary.** For n>=8, two native moves cannot take a labelled tournament to
a nontrivial vertex switch of itself. The bound is sharp at n=7, even with
both moves black.

**Proof of obstruction.** If the endpoint difference is delta(S), each vertex
has degree at most four in E(P) symmetric-difference E(Q). A nontrivial cut has
degrees |S| and n-|S|, so n<=8. At n=8 both sides have size four. But any
endpoint of P has degree at most 1+2=3 in the symmetric difference, a
contradiction. This proof does not assume either line colour. ∎

**Sharp construction.** Partition V as L={0,1,2,3}, R={4,5,6}, and take

    P=(0,4,1,5,2,6,3),  Q=(1,6,0,5,3,4,2).

Their edge sets are disjoint and partition K_(4,3). Orient T along P and
opposite Q; orient all within-part pairs arbitrarily. P is directed in T;
K_P reverses every Q edge, so Q is directed in K_P(T). Formula (4) gives
K_Q(K_P(T))=switch_L(T). All 2^9=512 within-part orientations work. The
checker records their colours and an explicit black-black source matrix.
Thus the obstruction cannot be extended to n=7, including on black words.

**Complete extremal skeleton classification at n=7.** Conversely, any
such two-step nontrivial switching word must have cut sizes3 and4, by the
degree bound. The cut has12 edges, exactly the maximum total of two
Hamilton paths. Consequently their edges are disjoint and partition the
cut. Relabel P as(0,1,2,3,4,5,6). Its bipartition is then forced to be even
versus odd positions, and the complementary cut edges form precisely

    Q=(2,5,0,3,6,1,4), or reverse(Q).

Thus there is **one undirected two-path skeleton up to relabelling**, with
two directed choices for its second path. This is a complete classification
of this extremal operation pattern, not of all order-seven tournaments.
In particular every extremal word has equal-coloured moves, because the two
positional reversals coincide. The number of labelled ordered legal words
(T,P,Q) of this kind is7!·2·512=5,160,960;645,120 are blue-blue and4,515,840
are black-black. These counts do not quotient labels or merge metagraph lines.

The two paths have the **same positional reversal**, namely
(0 3)(1 2)(4 6), fixing5. Native movement preserves anti-reversal for P,
so both moves have the same colour. On the nine within-part pairs this
involution has six unordered-pair orbits (four on L and two on R); fixed
pairs exchange their vertices and impose no further restriction. Thus there
are exactly2^6=64 blue-blue fillings and512−64=448 black-black fillings.
The cross-part directions already obey the anti-reversal equation. This
explains the exact census rather than merely reporting it.

### Nine-bit syndrome and the actual skeleton action

In the displayed L/R labels, use free bits x_uv=A[u,v] for u<v on

    (01,02,03,12,13,23,45,46,56).

The blue conditions are the three independent linear equations

    sigma(x)=(x01+x23, x02+x13, x45+x56)=(0,0,0) in F2^3. (5a)

The three disjoint two-bit supports make the map surjective of rank3.
Each of its eight fibres has64 fillings: the zero fibre is blue-blue and
the seven nonzero fibres are black-black. The fixed-reversal pairs03,12,46
are completely free. In another reference-orientation gauge these same
conditions may be affine, but their rank and fibre sizes do not change.

The automorphism group of the **unordered pair of undirected path edge sets**
is exactly V4. Its elements are identity,

    r=(0 3)(1 2)(4 6), h=(0 1)(2 3)(4 6), rh=(0 2)(1 3).

Indeed an automorphism sends P to P or Q, then has only two endpoint choices,
so there are at most4; these four directly preserve the pair. If the two
undirected paths are individually distinguished, only identity and r remain.
The group is a symmetry of the skeleton, not generally of a filled tournament.

Raw r and h reverse all forced cross-part orientations of the source family.
To act on this fixed directed source family, explicitly compose their
relabellings with switch_L; identity and rh need no switch. All four preserve
L and R, so this correction obeys the group law and gives a genuine V4 action
on the512 fillings. It is a **switch-normalized action**, not an assertion
that its elements are tournament automorphisms or legal native moves.

Its action on sigma is **trivial**: each of the three paired supports is
preserved, possibly exchanging its entries and/or complementing both.
For instance r sends (x01,x23) to (1+x23,1+x01), while h sends it to
(1+x01,1+x23); in both cases their sum is unchanged. Directly the same holds
for the other two coordinates. Thus these realized skeleton operations do
not permute the seven nonzero syndrome values. There is no justified
GL3(2), Fano-plane, or transitivity promotion here. Independent choice of an
abstract basis for the syndrome space is not a realized tournament symmetry.

`check_extremal_syndrome.py` enumerates all7! vertex permutations to identify
this group, all512 fillings to verify the eight64-element fibres, and every
pair of group elements on every filling to check the corrected action law.
Its14,850 explicit checks agree under normal Python and `python3 -O`.

**General length obstruction, allowing converse.** Let r>=1 and suppose a
legal native word has endpoint `converse^epsilon(switch_S(T))`, where
 epsilon is0 or1. Put d=(r mod2) XOR epsilon. Formula (4) says that the path
parity graph is `delta(S) XOR dJ`.

* If d=0 and S is a nontrivial cut, then **n<=4r−1**.
* If d=1 (including empty/full S), then **n<=4r+1**.
* If d=0 and S is empty/full, no size bound follows: paths can cancel,
  as in the two-step involution K_P K_P. This exclusion is essential.

These are necessary, not claimed sharp for general r. Every degree of the
parity path-union is at most2r, and at an endpoint of any constituent path
it is at most2r−1. A nontrivial cut therefore has both part sizes at most2r,
giving n<=4r; equality requires degree2r everywhere, contradicting that
endpoint. The complement of a cut is the union of its two within-part
cliques. For two nonempty parts their sizes are at most2r+1, giving
n<=4r+2; equality likewise requires degree2r everywhere and fails.
For an empty part the graph is K_n, and the endpoint argument actually
gives the stronger n<=2r. This also implies the stated general bound.
For an actual switch (epsilon=0), this specializes to the even-r bound
4r−1 (nontrivial cut) and odd-r bound4r+1. For converse-switch targets the
parity cases exchange. No inequality here asserts existence or legality.

### Exact edge-parity guard (inherited mechanism)

The fixed-path complement-bit description of THM-643, and the edge-parity
calculation in §8 of the recovered `TRIJECTIONS_AND_TREE_GAUGES.md`, also give
an immediate necessary guard on every length bound above. This is an
application of inherited parity, not a new general parity theorem.
Write N=binom(n,2), beta=binom(n−1,2), and s=|S|. Each native step toggles
exactly beta pairs, so a word of length r reaching
`converse^epsilon(switch_S(T))` must satisfy

    r beta = epsilon N + s(n−s)  (mod2).              (5b)

Indeed the endpoint toggle is epsilon J XOR delta(S); cardinality modulo2
is additive under XOR. This proof applies to every legal word, regardless
of line colours. It is necessary, not sufficient, and applies to labelled
endpoints without an additional relabelling.

The four residues give a compact parity grammar:

| n mod4 | necessary equation mod2 |
|---|---|
|0|r=s|
|1|no extra restriction from this scalar equation|
|2|epsilon=s|
|3|r=epsilon|

In particular, for n=2 mod4, **a pure switch across an odd-size cut is
unreachable by any number of native moves** (epsilon=0). The cut reverses
an odd number of pairs, whereas every native move reverses an even number.
This is stronger than a finite word-length lower bound. Yet the switching
quotient is connected at these even orders, because taking that quotient
erases exactly this orientation-parity distinction. The order6 finite
component census is reported separately in `CAYLEY_QUOTIENT.md`; no claim
of all-order sufficiency or precisely two native components follows merely
from (5b).

## 3. SCC reflection recurrence (inherited center law made quantitative)

Write the **canonical strong-component chain** as T=T_1 -> ... -> T_k.
Every Hamilton path concatenates one path from each component; components
cannot be revisited. Automorphisms preserve every component set. Let

    a_i = H(T_i)/|Aut(T_i)|,
    b_i = number of blue ordinary path presentations of T_i.

Then the total presentation count is A=product_i a_i. The blue count B is zero
unless T_i is isomorphic to T_(k+1-i)^op for every outer pair. When all these
matches hold,

    B = (product_(i=1..floor(k/2)) a_i) * c,
    c = 1 if k is even, and b_((k+1)/2) if k is odd.    (5)

In particular black presentations number A-B. Equation (5) concerns line
*endpoints*, not the number of distinct target classes or support neighbours.

**Proof.** Represent an ordinary path orbit by its normalized adjacency word.
This is complete: equal normalized words give the unique path-coordinate
isomorphism. Reflection reverses the component order and converses each word.
A blue full word chooses one arbitrary word on each left component; its right
partner is uniquely forced. A center word must itself be fixed. A mismatch of
component types makes matching impossible. This proves (5), including every
zero case. No centralizer or automorphism enumeration is needed in the test.

This is a refinement/rederivation of the inherited center law, not a claim
that the SCC principle is new. Purity holds exactly when every outer a_i=1,
all outer converse matches hold, and the center (if any) is pure-blue.

| Chain | H | A | B | Black |
|---|---:|---:|---:|---:|
| C3 -> C3 | 9 | 1 | 1 | 0 |
| C3 -> 1 -> C3 | 9 | 1 | 1 | 0 |
| T5 -> T5 | 225 | 9 | 3 | 6 |
| C3 -> T5 -> C3 | 135 | 3 | 3 | 0 |

Here T5 is the regular cyclic tournament. The third row is a hostile to
“composing pure-blue atoms gives pure-blue”; the fourth identifies the exact
survivor: a nonrigid pure-blue atom is allowed *at the center*, not in an
independently varying paired position. This is positional recursion, not a
complete alphabet theorem for strong components.

## 4. Boundary-port defect of substitution

Let T=Q[T_1,...,T_k], and let P concatenate one directed Hamilton path P_i
from each child in a directed Hamilton order of Q, with k>=2. Native K_P(T)
is obtained from

    Q^op[K_(P_1)(T_1),...,K_(P_k)(T_k)]

by toggling exactly the k-1 bridge edges from the exit port of one child to
the entry port of the next. All other cross-child edges reverse; these bridges
are precisely the cross-child path edges and must stay forward. This is an
exact recursive construction with named entry/exit ports.

Every old nonsingleton child ceases to be a module: choose an incident bridge
and its external endpoint. That vertex sees the child uniformly in the
reversed-quotient substitution, except at the single bridge endpoint, so it
sees both directions after the toggle. Singleton children remain modules.
This statement is about the **old child partition**; new modules can arise,
and a strong quotient can have non-contiguous Hamilton paths outside this
construction. For an SCC chain all Hamilton paths are contiguous, so the
formula describes every native move from that chain.

For arbitrary paths through general substitutions, one needs full ordered
runs and their jointly named ports, not just one endpoint per child. This
aligns with Session92's `runs/RECOVERED_RUN_WITNESS.md` in
`05-knowledge/results/session92_joint_order_carriers_20261029/`: separately
matched child cursors need not determine the same native partner. The
81/27 same-cyclic-triples control remains a distinct warning: cyclic triples
and induced Hamilton counts do not determine automorphism transporters.
Neither the new sign grammar nor (5) erases that information requirement.

## 5. Reproduction and limits

Run `python3 check_sign_grammar.py > checks.json` and repeat with `python3 -O`;
outputs must agree. The script exhausts all 1,099 labelled tournaments on
orders1..5, every Hamilton path, and all switches modulo complementary sets;
it checks (1)--(3), blue invariance under native moves, and (5). It checks all
legal two-step words through order4, selected chains through order11, and all
512 within-part fillings of the sharp order7 construction. Runtime guards use
explicit exceptions and remain active under optimization.

The source matrices and path-edge toggles provide a second direct verification
route; the separate session auditor independently tests the sharp obstruction.
The count recurrence is checked using normalized word sets rather than assumed
automorphism orders. This guards quotient conventions but is not a new class
enumeration algorithm.

**Open boundaries:** unrestricted prime-tournament classification, an intrinsic
minimal state for arbitrary run-list substitution, all-order black-core
connectivity, and optimal switching word lengths beyond r=2. The lower bounds
are for fixed labelled switching endpoints; allowing arbitrary final relabelling
is a separate problem and is not covered.
