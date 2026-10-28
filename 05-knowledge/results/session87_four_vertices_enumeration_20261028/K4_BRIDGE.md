# Four vertices, six edges, and four port observations

Status: **PROVED** for the stated marked linear model; **FINITE-EXACT** field realization and exhaustive small controls. This is a coordinate bridge, not an identification of fields with tournaments or a tournament-class invariant.

## The inherited question

[Session86](../session86_relations_attachment_carriers_20261028/NORM_PORTS.md) identifies the nine norm-one elements of F64, under its binary alternating pairing, with the nonzero vectors in three mutually orthogonal binary two-dimensional subspaces. Write

\[
V=V_1\oplus V_2\oplus V_3,\qquad
U=(V_1\setminus0)\sqcup(V_2\setminus0)\sqcup(V_3\setminus0).
\]

Fixed block-local partner queries require six bits; adaptive partner queries require five; arbitrary linear queries require four. The four-query lower bound is simply that nine possible ports cannot be separated by three binary answers. The question here is whether the six coordinates and four answers have an actual K4 realization.

## The boundary-plus-parity map

Give K4 its six binary edge coordinates. Let B be its unsigned vertex-edge incidence matrix over F2. Define

\[
M(x)=((Bx)_0,(Bx)_1,(Bx)_2,\sum_e x_e).
\]

The omitted fourth boundary coordinate is determined by the other three. The map has rank four. Its kernel has dimension two and consists of zero and the three four-cycles of K4. Indeed Bx=0 describes the cycle space; its even-weight subspace consists precisely of these four vectors. Equivalently this is K4's bicycle space, the intersection of its cut and cycle spaces.

Partition the six edges into three pairs and identify the three V_i with those coordinate pairs. Then M separates U **if and only if all three pairs consist of adjacent edges**.

Proof: a collision between two U vectors is a nonzero kernel vector supported in at most two pairs. A four-cycle is the complement of an opposite pair. Thus it misses an entire coordinate pair exactly when that pair is opposite. Within a single pair the kernel has no nonzero vector. This proves both directions.

There are fifteen pair partitions in total and exactly eight successful ones. Inclusion-exclusion on the three possible opposite pairs gives 15−9+3−1=8. A successful partition has three distinct pair-centres. The omitted vertex is a root; the remaining three edges form a directed triangle when each is assigned to its pair-centre. Thus the eight choices are exactly four roots times two cyclic directions. Conversely, pair each root spoke rv with the outgoing triangle edge vv+.

Declaring the root a source identifies these eight scaffolds with the eight labelled tournaments 1→C3; declaring it a sink gives the other converse class. Literal full arc reversal also reverses the triangle, so it changes the wedge scaffold. The source/sink choice is additional data, not determined by the wedge partition. Each scaffold has stabilizer C3 under vertex relabeling.

## An explicit F64 realization

Use F64=F2[a]/(a^6+a+1), with integer masks encoding coefficients in the polynomial basis. Take ordered binary bases

\[
(1,58),\quad(6,26),\quad(11,20)
\]

for the three orthogonal planes. Map these basis elements respectively to the edge pairs

\[
(03,01),\quad(13,12),\quad(23,02).
\]

The nine ports are {1,58,59,6,26,28,11,20,31}. With B1(x,y)=Tr(xy^8), the four rows of M become the four field probes

\[
48,29,25,56.
\]

The kernel in field coordinates is {0,9,35,42}. It misses every difference of two distinct ports, so the four answers separate all nine. The field-to-edge map is linear on all 64 vectors, not merely a lookup table on the nine ports. It preserves B1 through the scaffold-dependent symplectic form J(x,y)=sum_(paired e,f)(x_e y_f+x_f y_e), not through ordinary edge dot product.

The producer and an independent checker verified the field arithmetic, all 64 coordinate vectors, all 15 pair partitions, and the 36 avoiding two-dimensional kernels among the 651 binary planes. The independent checker also tested all 1,536 pairs of a K4 orientation and a vertex permutation. Producer reproduction is given below by the co-located stdlib checker `check_connections.py` and its complete receipt `connections.json`. That code is independent of Session86's implementation; the separate auditor's code is in the coordinating packet's audit directory.

### All 36 kernels are realized, not newly invented

A linear response map separates U iff its kernel avoids every nonzero port difference. These differences are exactly the vectors supported on one or two planes. Thus every nonzero vector of a valid two-plane kernel has all three components nonzero, and each component projection is an isomorphism. Identifying the kernel through the first component gives

    {(x,Ax,Cx):x in V1},  A:V1→V2 and C:V1→V3 invertible.

There are |GL(2,2)|²=6²=36 choices. Changing the three block bases in the fixed successful wedge scaffold realizes all of them. This is elementary direct-sum linear algebra, already used in Session86, not a novelty claim for code classification. The checker exhausts all 216 ordered block-basis choices and recovers exactly those 36 kernels.

## Hostile controls and retained sidecars

1. **The natural opposite-edge pairing fails.** In the tested opposite-pair assignment, ports 28 and 31 have the same four-bit answer. Merely saying “six edges and four vertices” does not identify the successful map.
2. **The four-bit map does not classify tournaments.** Number edges 01,02,03,12,13,23 and use the transitive orientation as zero. Flip masks 0 and 45 have the same response but different isomorphism types. Across all 64 orientations, twelve output fibres each contain two transitive tournaments (H=1) and two strong tournaments (H=5). The other four fibres each contain four source-over-C3 or sink-under-C3 tournaments (H=3).
3. **Relabeling is affine in flip coordinates.** The vertex-boundary bits are changes in outdegree parity. Total flip parity changes by the permutation-sign cocycle when the reference transitive orientation is relabeled. A fixed orientation baseline must not be silently erased. In formulas x↦P_sigma x+c_sigma; P_sigma is an even permutation of the six unordered edge coordinates, but parity(c_sigma)=sign(sigma). These are distinct parity representations. With the stated transitive baseline, actual outdegree parity mask is 5 XOR boundary(x).
4. **This is not a multiplicative field quotient.** The nonzero proper kernel cannot be an ideal of F64. The faithful order-nine multiplication action on norm-one ports cannot be realized by K4 vertex relabelings, since S4 has no element of order nine. Concretely, 9 is in the displayed kernel but 2·9=18 has output 7, not zero. The multiplication of U can still be transported to its nine distinct outputs through the decoder as a separately defined, generally nonlinear operation; it is not multiplication of arbitrary quotient classes.
5. **The selected alphabet is not closed under addition.** The map is lossless on U but not on arbitrary sums of its elements. Full additive closure is all six-dimensional V and needs all six bits for exact reconstruction.

## What the connection buys

The six-versus-four gap has a genuine marked tetrahedral realization: a rooted cyclic scaffold selects a successful encoding, and boundary plus edge parity supplies four linear tests. The information saving comes from restricting the source to nine one-block states. It does not recover arbitrary joint attachment information for free, identify tournament classes, or replace the exact transporters required by the 81/27 hostile pair.

## Evidence pin and scope

Sources below were read at `b8672af00bcc6d440a9a00a6915325faf39fcdff`. The main proof was first checkpointed in this session at `cde5060a87990625f6b9cadaff27823b906e91fd`; this expanded version adds the recovered lineage and verified coding/arithmetic boundaries. No historical novelty claim is made.

## Source lineage, and what is actually being added

Closest proved mechanism: Session86 `NORM_PORTS.md`, sections 3–4, already
proves the orthogonal three-plane decomposition, costs **6 fixed / 5 adaptive /
4 unrestricted linear**, and the 36 avoiding kernels. None of these numbers
is a new claim here.

Canonical hostile: the K4 bicycle space is nonzero in characteristic two;
cut and cycle spaces are not complementary. This was already explicit in
THM-1405 and THM-1440. Corrected near miss: THM-1415's canonical switching
quotient is not the tournament isomorphism quotient. Least-used relevant
sidecar: a rooted cyclic scaffold partitions K4 edges into adjacent wedges,
rather than its intrinsic three opposite matchings.

The additional synthesis here is the exact criterion for when the bicycle
quotient separates Session86's port set, its eight scaffold realizations,
and the displayed field-to-edge transport. A bounded source search recovered
the prior bicycle and switching results; it does not establish publication
priority for the elementary coding/graph constructions.

### Source routes

All paths are repository-relative. Their statements, not chronological logs,
are the inherited evidence.

| Source | Inherited content and limit |
|---|---|
| `05-knowledge/results/session86_relations_attachment_carriers_20261028/NORM_PORTS.md` | Exact normalized U9 model, 3K3 pairing, 6/5/4 costs, 36 kernels. |
| `01-canon/theorems/THM-1405-star-quotient-is-the-cycle-space-transverse-to-isomorphism.md` | Cut/cycle quotient and its failure to descend to tournament isomorphism; bicycle caveat. |
| `01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md` | Canonical full cut switching; only two switching classes versus four ordinary classes at order four. |
| `01-canon/theorems/THM-1440-seidel-spectra-are-sine-and-the-odd-n-projection-is-one-parity-fact.md` | K_n bicycle dimension is n−2 at even n, hence two at n=4. |
| `01-canon/theorems/THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration.md` | S4 edge permutations are even; opposite matchings retain the S3 resolvent/sign sidecar. This is a different parity from orientation-bit parity. |
| `05-knowledge/results/fano_tournament_move_quotients_20261021/session23b-fano/README.md`, section 2 | Exact reduced-incidence Fano deletion `F7 \ omega = M(K4)`; metric/orientation data are not automatic. |
| `05-knowledge/results/duck_tournament_20260925.md`, section 4 | Four ordinary classes; cut-switching Pfaffian classes; marked root and gauge losses. |
| `01-canon/theorems/THM-4472-four-vertex-reading-of-3n-plus-minus-1.md` | Actual arithmetic four-vertex tournament map; converse is time reversal, not a free sheet swap. No arithmetic map from the present code is supplied. |
| `01-canon/theorems/THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy.md` | Integral rooted-tree incidence transport and graceful A/D discriminants; does not prove graceful existence. |
| `01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md` | Sum/difference incidence and the absolute-value collision sidecar. |
| `05-knowledge/results/graceful_square_exchange_20261019.md`, sections 3–5 | Four-point matching obstruction to graceful two-edge trades, and separate square-sum seam constraints. |

Live concept board: (i) restricted observations, (ii) block-distance kernels,
(iii) cut/cycle/bicycle codes, (iv) rooted cyclic edge scaffold, (v) Fano
syndrome, (vi) tournament affine gauge, (vii) integer arithmetic incidence.

## Fano deletion is the inherited syndrome, not a hidden tournament decoder

The reduced incidence columns of K4 are

    03↦1, 13↦2, 23↦4, 01↦3, 02↦5, 12↦6.

Adding omega↦7 gives all seven nonzero vectors of F2^3. Triangular faces
and opposite pairs completed with omega are exactly the XOR-zero triples.
This is the inherited Fano deletion map. The reduced boundary alone is its
three-bit syndrome; its kernel is the three-dimensional cycle code, not
the two-dimensional separating kernel. Adding edge parity restricts that
kernel to its even-weight bicycle subcode.

This distinguishes three dimensions often conflated:

* seven Fano columns have a four-dimensional relation kernel;
* deleting omega leaves six columns with a three-dimensional cycle kernel;
* adjoining a parity observation leaves a two-dimensional bicycle kernel
  and gives a four-dimensional **output**.

The Fano relation-space dimension four is not the assertion that four
vertex parities are independent. The 7→15 relation-object construction is
another operation and is not identified with this nine-state port set.

### The shortened codes, with their exact dimensions

Let H be the binary Hamming `[7,4,3]` code with parity-check columns all
seven nonzero vectors of F2^3, and let S=H-perp be the `[7,3,4]` simplex
code. **Shorten** a code at omega by first requiring that coordinate to be
zero and then deleting it. This is not puncturing, which only deletes.

* Shortening H gives the K4 cycle code `[6,3,3]`.
* Puncturing S gives the K4 cut code `[6,3,3]`.
* Shortening S gives the bicycle code `[6,2,4]`: its three nonzero words
  are the three four-cycles.
* The observation row space is bicycle-perp, a `[6,4,2]` code.

The cycle and cut weight enumerators are both `1+4z^3+3z^4`; the bicycle
one is `1+3z^4`; the observation row space has
`1+3z^2+8z^3+3z^4+z^6`. These classical identities are independently
checked here. They connect the old relation-object language to the new
observation language without treating a code and its dual as the same data.

The dimensions distinguish seven-column relation data (four free bits),
six-column cycle relations (three), bicycle ambiguity (two), and the four
output bits. No dimension coincidence by itself supplies a decoder.

## Graceful and square-sum incidence: a real shared map, a strict boundary

For integer vertex potentials, oriented differences are the transpose of
an oriented incidence matrix. Reducing modulo two makes + and − identical,
recovering the binary cut/gradient map. On a rooted tree the reduced
integral incidence is unimodular, so path integration really recovers all
vertex potentials (THM-2770). On K4 there are three cycle compatibility
conditions instead. This is an exact incidence connection, not a count
analogy.

Graceful labels additionally require distinct absolute edge differences
covering a specified interval. Square-sum edges require their actual integer
sums to be squares. Both predicates discard essential information upon
reduction modulo two. For example, sums 1+3=4 and 1+5=6 have identical
endpoint parity data but different square status. Binary incidence cannot
decide this predicate; magnitudes and the actual label assignment are needed.

The earlier four-vertex graceful trade obstruction is especially relevant.
For a<b<c<d, the three perfect matchings have length multisets

    {b−a,d−c}, {c−a,d−b}, {d−a,c−b}.

The first has smaller total than the other two, and the last has a larger
maximum than the middle. No two agree. Hence a degree-preserving two-edge
trade on these four vertices cannot preserve a graceful difference multiset.
The common K4/matching scaffold is real, but equality of binary syndromes
does not license that arithmetic trade. No graceful-tree, square-sum
Hamiltonicity, Collatz, or LRC consequence is claimed here.

## Connection contract and closure boundary

**Source:** marked F64, its norm-one U9, pairing B, and three symplectic bases.
**Target:** labelled K4 binary edge space with rooted cyclic wedge scaffold.
**Map:** phi sends bases `(1,58;6,26;11,20)` to edges
`(03,01;13,12;23,02)`, followed by vertex boundary and total edge parity.
**Preserved:** binary addition by phi, symplectic pairing through J, and
identity of U9 under h∘phi. **Destroyed by the quotient:** two additive
degrees of freedom on all E, tournament class/H, and general multiplicative
attachment information. **Required sidecars:** field marking, bases,
scaffold, and baseline if interpreting masks as orientations.

The graph quotient is injective only on the promised nine-state source,
not on its additive closure E. Multiplication does not descend to the
linear quotient: 9 is in its kernel, but 2·9=18 has response 7, not zero.
The U9 group operation can still be transported via the nine-state decoder
to its image as a separately tabulated, generally nonlinear operation.
That is not multiplication of arbitrary quotient classes.

**Cheapest decisive test:** the opposite-pair hostile 28 versus 31, followed
by the all-nine response table for the adjacent scaffold. **Strongest
survivor:** a complete and explicit classification of edge-pair scaffolds
for this graph decoder. No unmarked canonical bridge is asserted.

## Reproduction and finite scope

From this directory:

    python3 check_connections.py > connections.json
    python3 -O check_connections.py > connections-optimized.json
    cmp connections.json connections-optimized.json

Both runs pass and are byte-identical. The checker uses explicit exception
gates and only stdlib. It exhausts all 64 field elements; all 512 restricted
candidate subsets for adaptive costs; all normalized partner subsets through
size six; all 651 two-planes; all 15 pair partitions; all 216 block-basis
choices; all 64 tournaments and 24 relabellings; and all 16 graph-output
fibers. It independently verifies the transported pairing on all 4096
ordered field pairs. It does not test large-field generalizations, arithmetic
problem predicates, or general tournament generation. The coordinating independent audit has separately confirmed the field map,
pairing, all fifteen scaffolds, all 36 kernels, and affine-action identities.
Its own script and receipts live in the coordinating packet's audit directory;
this producer receipt is not a substitute for that independent implementation.

