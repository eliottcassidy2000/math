# One-bit incidence repair, full order-eight wiring, and dyadic switch tests

2026-10-23. **FINITE-EXACT:** all 6,880 ordinary tournament classes of order8,
all 2^21 fixed-path tilings, their aggregated ordinary-class ports and merged
coloured wiring. **PROVED, elementary:** the doubleton decoder, odd-cell
refinement rule, and electrical switching criterion below. **OPEN:** a
uniform short signature or efficient unrestricted enumeration algorithm.
These are repository-relative advances; no literature-priority claim is made.

## Inheritance and live board

The closest mechanism is [THM-4739, Kirchhoff/Hodge port
laws](../../../01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md),
refined by the [incidence/capacity packet](../metagraph_bucket_obstruction_20261023/README.md).
That packet separated 6,878 of 6,880 order-eight classes with J4 and all6,880
with a second refinement, but had not replayed their full port wiring.
The canonical hostile is its reflected order-six capacity switch. The
corrected near miss is promoting class separation into an unperformed port
census. The least-used sidecar is the correlation between two outgoing
incidences of the SAME vertex inside a repeated profile cell.

Anchor / niche / wildcard: minimal actual-port repair / electrical detection
of rewiring / recovery of the old order13 additive-energy collision.
Concept board: profile cells; common-source overlap; immutable endpoint
ports; binary gauges; dyadic energy; marked-source refinement.

The [signed-curl decoder](../metagraph_port_rigidity_20261023/PROOF.md)
already supplies an all-order lossless endpoint. This packet does not weaken
or replace it: it measures how little of that information is needed in the
currently ambiguous finite fibres.

## 1. Exactly one new bit on each surviving order-eight profile

For each vertex v retain its tag (s_v,c_v,f_v): outdegree, cyclic-triangle
incidence, and induced strong-four-tournament incidence. Set

    rho(v)=(tag(v), multiset of tag(w) over v->w),
    J4(T)=multiset of rho(v).

Rank the DISTINCT rho types lexicographically, starting at0. Even retaining
the complete directed quotient table

    E_ij = #{u->v : rho(u) has rank i, rho(v) has rank j}

does not separate the two known pairs. Each has seven types: six singletons
and one doubleton D. For singleton types i,j with vertices v_i,v_j, define

    W_D(i,j)=sum_{u in D} A[u,v_i] A[u,v_j].                    (1)

This is a two-incidence measurement retaining their common source. It is
invariant under relabelling; its observer types are intrinsic, not literal
vertex names.

| Ordinary classes | Canonical codes | D type | Chosen observers | W values |
|---|---|---:|---|---|
| 3490,3676 | 27467784,27467793 | 4 | 0,3 | 0,1 |
| 4239,5173 | 6427781,6427778 | 5 | 3,6 | 1,0 |

The second pair is the converse of the first. Both pairs have H=523 versus
527, but H is NOT an input to the repaired signature. The complete E table
is equal within each pair. Consequently every statistic formed solely from
E and the rho types, including arbitrary linear weighted type currents,
still misses the distinction.

Define J4* to equal J4 except that, on either of these two exceptional J4
values, append its displayed W bit. It separates all6,880 classes. This is
one bit per exceptional fibre, and zero on all other fibres. One bit is
minimal relative to either two-element surviving candidate set. No claim
of globally shortest serialization, minimum-cost computation, or uniform
all-order one-bit sufficiency is made.

## 2. The full wiring test, including an actual ambiguity witness

The path convention is 0->1->...->7. Free cells (i,j), j>=i+2, are ordered
lexicographically, with bit1 meaning i->j. Complement flips all21 free
bits; reflection sends (i,j) to (7-j,7-i). A closed source has 7->0.
If its endpoint scores are (a,b), its open target has scores (8-a,6-b).
Converse reflection sends (a,b) to (7-b,7-a).

The producer independently expands every directed Hamiltonian path of each
canonical representative. The resulting normalized masks cover every one
of the 2,097,152 presentations, with no overlap between ordinary classes;
each mask occurs exactly |Aut(T)| times. The Python consumer additionally
checks score-canonical representatives over403,193 permutations, whose
orbit masses sum to all2^28 labelled tournaments. Thus the universe is
complete, not a sample or just the self-converse subcube.

| Quantity | Exact value |
|---|---:|
| Directed paths enumerated in the6,880 representatives | 2,143,304 |
| Closed-to-open actual ports | 1,048,576 |
| Nonzero (source,target,endpoints,colour) records | 956,339 |
| Converse-merged nodes | 3,528 |
| Nonzero unordered node-pair/colour records | 438,374 |
| Blue tilings / literal blue lines | 4,096 / 2,048 |
| Literal black lines | 1,046,528 |

Given actual source rows, transported target J4, endpoint types, colour and
target capacities, J4 leaves370 rows and884 units of unresolved mass. Each
collision pair contributes185 rows of mass442. Each such row has two
available target columns of the same endpoint/colour type, with at least
two rows per residual component; singleton capacity peeling cannot proceed.
More decisively, the following real feasible switch proves nonuniqueness.
Keep blue unchanged. In the merged graph change literal black counts by

    (1763,3490): -2       (1763,3676): +2
    (3481,3490): +2       (3481,3676): -2.              (2)

The eight ordinary-class port edits, including their converse reflections,
are frozen in receipt.json. They preserve every source J4 row, target
endpoint capacity, converse pairing, nonnegativity, degree and black-twin
parity. They change the merged wiring. Twice all eight edits is also
feasible. This is an authenticated counterexample to THESE relaxed
constraints, not a claim that the counterfeit is geometrically realizable.

With J4* each target-signature fibre is a singleton. Every supplied source
row therefore has exactly one target column, with the native count. The
consumer verifies that these counts reproduce all956,339 actual port
records and the entire merged graph. The compressed output
[n8_merged_wiring.tsv.gz](n8_merged_wiring.tsv.gz) contains that full graph,
including colour, multiplicity and loops. Node names are the minimum
ordinary-class ID in each converse pair, not a new contiguous numbering.

Crucial scope: the source rows contain the signature of the ACTUAL
transported target. Computing it uses the actual source/path information;
the theorem does not recover that information from node margins, or claim
a speedup that ignores target-signature costs. Aggregated wiring is
determined; interchangeable literal paths inside an aggregated record are
not individually named by this signature.

## 3. An all-size decoder for one two-vertex ambiguity cell

Let a tournament carry a specified intrinsic partition with every cell a
singleton except D={u0,u1}. Retain E, and the common total outdegree s of
the two vertices in D. This last equality holds for a repeated rho cell.
For each singleton v_j put m_j=E_Dj in {0,1,2}.

If m_j is0 or2, its relation with D is fixed. If m_j=1, let x_j be the
binary choice of which D vertex beats v_j. Let M be these mixed observers,
with |M|=t. For i,j in M,

    W_D(i,j)=1 iff x_i=x_j,
    x_i+x_j=1-W_D(i,j) over F2.                        (3)

**Doubleton theorem.** E, s and the W bits on any spanning tree of M
reconstruct the full tournament up to swapping u0,u1, if realizable.
An arbitrary supplied connected observer graph is realizable as binary
row choices iff the right sides of(3) XOR to zero around every cycle.
The reconstructed rows must additionally pass their score and original
profile guards; parity alone does not certify those guards.

Proof: tree integration fixes all x after one initial choice. Changing
that choice swaps u0,u1. Entries between singleton vertices are given by
E. The nonmixed columns are fixed. Their contributions to the two D row
sums agree. The internal D arc contributes1 to exactly one row; equality
of total outdegrees therefore forces the external row sums to differ by
exactly1, and forces the lower external row to win the internal arc.
Reject unless the recovered common total is s and the retained profiles
match. This reconstructs every arc. The usual cycle-XOR criterion is both
necessary and sufficient for the binary potential. QED.

The degree condition sharpens the parity: t is odd and the two external
mixed row counts are (t-1)/2 and(t+1)/2. There are
binom(t,(t-1)/2) such assignments modulo swapping, before any additional
profile restrictions. For t=3 their three pair overlaps are one-hot: just
one of the three observer pairs agrees. In both order-eight collision
fibres a known unequal pair already excludes one possibility; one bit
chooses between the remaining two.

The t-1 tree bits give a uniform constructive upper bound for this
restricted family. They are minimal for UNRESTRICTED binary assignments
modulo swap, not necessarily for the degree-constrained family or a
particular known candidate list. Multiple unresolved cells introduce
further correlations; this theorem does not silently cover them.

## 4. Parity guides refinement, but a source mark is real information

**Odd-cell lemma.** Every cell of a stable directed neighbour-count
partition of a tournament has odd size. More locally, any even-sized cell
must split at the very next refinement.

Proof: if a cell C survives without splitting, all its vertices have the
same number d of outneighbors inside C. Summing those counts gives
|C|d=|C|(|C|-1)/2, hence d=(|C|-1)/2 is integral. QED.

This explains why a doubleton cannot persist indefinitely. It does not
imply that stable cells are singletons. The old
[March additive-energy experiment](../../../04-computation/H_determined_by_E.py)
supplies a useful hostile, recovered here and independently recounted.
On Z/13 let i->j iff j-i belongs to S, for

    S1={2,7,8,9,10,12},       H=3,683,797,
    S2={1,2,4,7,8,10},        H=3,704,857.

Both tournaments are vertex-transitive, with the same tag(6,21,112) at
every vertex. Their rho types, J4, complete type-incidence quotient, and
EVERY further unrooted neighbour-refinement round agree: the single
13-vertex cell is regular. Unequal H proves they are nonisomorphic. This
is a concrete boundary, not an inference from counting hardness.

Marking the actual source vertex0 in either example gives cell sizes

    [1,12] -> [1,6,6] -> [1,1,1,1,1,4,4] -> [1,...,1].

Three rounds separate all vertices. The source-marked discrete profiles,
with one further outgoing-colour read, recover every arc. This is an
example of using an already available port coordinate, not assigning an
arbitrary root to an unrooted class and forgetting its choice.

The old pair already had common additive energy118 and distinct cubic
sum-energies420/408. This packet's new use is the permanent refinement
obstruction and its marked-source repair, not rediscovery of the old H
collision. During recovery, Part3's historical energy display was repaired:
sum convolution has r(0)=0 for a skew set, while difference convolution has
r(0)=|S|. Their second moments agree, but their third moments need not.
The correction is logged in MISTAKES.md and checked on126 skew sets.

## 5. Electrical switching: an exact all-size detection criterion

Use the inherited literal multigraph convention G=L+2B, with blue L and
paired black B. Loops contribute twice to fibre degrees and zero to the
Laplacian. Suppose G is connected and write

    V=ker_F2 Lap(L) / <constant>,
    q2(x,y)=x^T Lap(G)y/2 mod2.                         (4)

The integrality and lift-independence of(4), and its identification with
the first dyadic Smith layer, are established in the incoming
[electrical packet](../session36_information_curl_contexts_20261023/metagraph/DYADIC-ELECTRICAL-METAGRAPH.md).
The known order-six counterfeit was detected by a single off-diagonal
entry. Here is its all-size criterion.

For four DISTINCT vertices a,b,c,d, perform the degree-preserving paired
black switch +ac+bd-ad-bc, assuming the resulting weights stay nonnegative
and the full graph remains connected. Set alpha(x)=x_a+x_b and
beta(x)=x_c+x_d on V. Direct expansion gives

    delta q2(x,y)=alpha(x)beta(y)+beta(x)alpha(y).        (5)

Thus delta q2 has rank2 exactly when alpha,beta are linearly independent;
otherwise it is zero. Equivalently it is invisible iff one of
e_a+e_b, e_c+e_d, e_a+e_b+e_c+e_d belongs to im_F2 Lap(L).
Indeed a functional vanishes on ker Lap(L) iff its coefficient vector
lies in that image; all three functionals annihilate constants. Formula(5)
is a two-dimensional alternating form when independent and zero otherwise.

In particular, every such unit switch among four blue-isolated vertices
is detected: their independent unit potentials survive in V. More
generally q2(e_i,e_j)=B_ij mod2 for distinct blue-isolated i,j, so it
recovers the entire paired-black off-diagonal parity matrix on those
vertices. Its diagonal alone misses EVERY degree-preserving switch(5).

The order-eight switch(2) has delta q2 rank2. Doubling it preserves the
ENTIRE q2 form and is still feasible in the actual port-capacity problem.
Consequently even this stronger electrical signature is insufficient at
order8. This is not merely a hypothetical large-multiplicity obstruction.
The incidence bit rejects both switches by restoring the actual target.

## 6. The full order-eight dyadic calculation and the next obstruction

The native full graph is connected. Its paired-black support has one
component of size3,525 and three isolated vertices. This confirms the
finite black-core connectivity check at order8, not an all-size theorem.
The blue graph has one176-node component and3,352 isolated nodes; the
connected component's binary bicycle dimension is11. Therefore

    dim V = 3,353-1+11 = 3,363.

On V, q2 has rank3,358 and radical dimension5. The first two layers alone
give3,358 critical-group factors with exact2-adic valuation1 and five with
valuation at least2.

To test the five survivors, root the full Laplacian at one blue vertex.
For x in rad(q2), choose an integer representative and solve

    Lap(G)y = -Lap(G)x/2 mod2.

This is solvable precisely because x is in the radical. Put v=x+2y;
then Lap(G)v is divisible by4. Define on rad(q2)

    q4(x,z)=v_x^T Lap(G)v_z/4 mod2.                    (6)

Changing a lift by2w changes(6) by w^T Lap(G)v_z/2, which is even.
Thus the form is well defined. Splitting off a unimodular block over Z_2
and repeating once identifies q2 and q4 with the reductions of the
successive divided symmetric blocks. Their ranks count Smith factors
with exact valuations1 and2, respectively. This also proves that the
quarter-energy radical counts factors divisible by8; no diagonal-only
shortcut is used.

The exact five-by-five q4 matrix is

    0 0 0 1 0
    0 0 1 1 0
    0 1 1 1 0
    1 1 1 1 1
    0 0 0 1 0

It has rank4. The native graph consequently has3,358 factors of exact
valuation1, four of exact valuation2, and ONE factor of valuation at
least3. Lifting that remaining direction to divisibility by8 gives the
one-by-one eighth-energy form[0], so its valuation is actually at least4.
For this lift, first cancel the obstruction against V using q2 and an
adjustment2w, then solve the remaining blue Laplacian equation and add4y.
The same unimodular-block argument identifies its zero form with the next
Smith layer. Independent all-minor controls verify this lifting procedure
on120 nonsingular symmetric integer matrices of orders1 through5.
This packet does not determine the last valuation or the odd
primary part of the critical group.

The doubled counterfeit has delta q4 rank2 on this SAME radical. Its
weights differ by a multiple of4, so the saved native lifts also remain
valid lifts for the counterfeit. This pins the lost information exactly
one layer deeper: mod2 misses any black pairing, q2 sees the unit switch,
q4 sees the doubled switch. Further even multiples can evade a fixed
layer. No fixed finite dyadic truncation is asserted to characterize all
possible wiring multiplicities.

## 7. Incoming work, connections, and next tests

The pulled [Session37 synthesis](../session37_statistics_parity_cuts_20261023/SYNTHESIS.md)
and its [independent audit](../session37_statistics_parity_cuts_20261023/statistics-audit/AUDIT.md)
distinguish exclusion of repeated vertices, exchange signs, and repeated
component occupation. Their exterior representation has the tournament
character and an exact even/odd component product, but does not supply a
converse- or port-compatible basis bijection. We therefore keep three
different parity predicates separate: odd profile-cell size, black-line
pairing, and exterior component signs. The shared useful move is retaining
the SAME source/component through a joint measurement; their objects are
not identified by their numerical parity alone.

| Connection | Map and preserved predicate | Lost information / needed sidecar | Decisive test |
|---|---|---|---|
| Incidence -> binary potential | Mixed two-row columns -> equations(3); preserves exact row correlations | Global row swap; score/profile guards still required | Cycle XOR and full decoded adjacency |
| Wiring -> electrical form | G=L+2B -> q2; retains tested switch parity | Even multiplicity digits and all Laplacian-invisible loops | Native doubled counterfeit; q4 |
| Symmetric class -> source-marked refinement | Retain actual port start; refine without moving it | Unrooted canonical choice is not supplied | Recovered order13 pair, all13 vertices discrete |
| Component statistics -> class census | Incoming character identity gives counts | Individual class and actual port map | Converse/parity mismatch in incoming order3 hostile |

Next precise target: build an adaptive port certificate which first uses
endpoint marks and neighbour counts, applies(3) only to the remaining
two-vertex fibres, and branches explicitly on odd regular cells that
survive. Its fallback is the already-proved signed-curl decoder. Measure
the number and cost of the retained incidences against actual target
canonicalization. The order13 hostile says to test source-marking before
adding another unrooted moment. A further electrical target is the one
remaining native factor divisible by16, not a claim that its valuation
will by itself establish all-size black-core connectivity.

## Reproduction

From the repository root:

    python3 -B 05-knowledge/results/incidence_wiring_20261023/check.py
    python3 -O -B 05-knowledge/results/incidence_wiring_20261023/check.py
    python3 -B 05-knowledge/results/incidence_wiring_20261023/controls.py
    python3 -O -B 05-knowledge/results/incidence_wiring_20261023/controls.py
    python3 -B 05-knowledge/results/incidence_wiring_20261023/dyadic_controls.py
    python3 -O -B 05-knowledge/results/incidence_wiring_20261023/dyadic_controls.py

The first command compiles wiring.cpp with clang++ and regenerates all raw
ledgers in a temporary directory. Explicit exceptions survive -O. Normal
mode checks the frozen receipts rather than overwriting them; --write is
an intentional regeneration option. --prefix may reuse raw producer
outputs but is not the fresh reproduction command. Source hashes, binary
ledger hashes, exact universes, counterfeits and dyadic ranks are frozen in
receipt.json. Controls exhaust all1,099 labelled tournaments through
order5, all54,253 tournament/partition pairs,2,502 equal-score doubleton
decodings and192 four-vertex electrical switch/blue-graph cases. Both
Hamiltonian counts at order13 are checked by two separate recurrences.
