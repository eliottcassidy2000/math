# Divisor ideals as marked graphs: a sharp edge budget and twisted complement

2026-10-25. **PROVED elementary**, with **FINITE-EXACT** controls for k=2..6.
This is a structural supplement, not a literature-priority claim. It supplies
an order and lattice representation with specified graph operations; it does
not identify arbitrary divisor arithmetic with arbitrary graph adjacency.

## Inheritance and the precise improvement

[Four-slot compression, §3](../collatz_four_slot_compression_20261004.md#3-the-strong-four-tournament-and-the-twelve-state-lattice)
already identifies the twelve divisors of p²qr with ideals of a two-element
chain plus two independent elements. It also records the failure of
intersection under the occupancy quotient of arbitrary subsets. This ideal
dictionary is inherited, not rediscovered here.

The [ten-state switched graph model](../reflection_context_20261020h/README.md#2-a-faithful-replacement-with-a-two-edge-switch)
retains both reflections but requires decoding for divisibility. The
[minimal twelve-state resolution](../divisor_reflection_resolution_20261021/REFLECTION_RESOLUTION.md)
retains ordinary graph complement and vertex reversal by splitting two
source fibres; its explicit graph table is not an order embedding.
The companion [graph bridge](graph_bridge.md) gives the exact distinction
between action embeddings, monotone injections, and order embeddings.

The construction below combines the inherited ideals with actual graph-edge
incidence. It preserves divisibility as edge inclusion and both lattice
operations. Its price is explicit: divisor complement is a **twisted
complement inside a marked host graph**, not ordinary complement in K_m.
The new small realization uses a square host on four vertices for k=2;
the number of Boolean edge coordinates is minimum for every k>=2.

Board: divisor order; prefix ideals; complement; prime swap; graph incidence;
quotient fibres. The canonical hostile is a proper divisor below its own
complement, and the least-used sidecar is the ordering of the repeated-prime
slots. The portfolio treats this as the structural supplement to the
session's reflection and seven-point questions.

## 1. The all-k marked host

For distinct primes p,q,r and k>=2 let L_k be all divisors of N=p^kqr.
Write d=(a,b,c)=p^a q^b r^c with 0<=a<=k and b,c in{0,1}. Let D_k omit
the unit and N. Use a marked star with hub h and leaves

    p_1,...,p_k,q,r.

Its edge universe E has k+2 edges, one incident with each named leaf. The
graph assigned to d has edge set

    I(d)={hp_i:1<=i<=a} union {hq:b=1} union {hr:c=1}.   (1)

All vertices remain present, including isolated ones. Thus this is a map to
labelled spanning subgraphs of one fixed host. Its image consists exactly
of the ideals of the slot poset p_1<...<p_k, with q,r incomparable to every
other slot. The inverse counts the p-edges and reads the q,r bits.

**Proposition I1 (PROVED).** The map I is an order embedding and a lattice
isomorphism onto its image:

    d|e iff I(d) subseteq I(e),
    I(gcd(d,e))=I(d) intersection I(e),
    I(lcm(d,e))=I(d) union I(e),
    |I(d)|=a+b+c.                                      (2)

The proof is coordinatewise minimum/maximum of prefix lengths and bits.
For D_k the first equivalence remains exact. Gcd and lcm are understood in
the full ambient lattice L_k: D_k itself is not closed. For example
gcd(p,q)=1 and lcm(p^kq,p^kr)=N. Intersection or union remains in the proper
image exactly when its result is neither the empty nor the full host.

The inherited faithful plane chart is directly readable from the edge set:

    f_k=2*number_of_p_edges-k*(q_bit+r_bit),
    g=q_bit-r_bit.                                     (3)

On D_k this pair is injective; on L_k its only collision is the omitted
bottom/top pair. For g=1 or-1, b,c are known and f_k determines a. For g=0,
the intervals f_k=2a (b=c=0) and f_k=2a-2k (b=c=1) meet only at0, at those
two endpoints. At k=2, f_k is twice the earlier coordinate a-b-c.

## 2. Exact operations, including the missing twist

Let K_E toggle all k+2 host edges. Let J be the vertex permutation
`p_i -> p_(k+1-i)`, fixing h,q,r. Let R swap the q,r leaves, fixing h and
all p-leaves. These are actual commuting vertex involutions on the host.
Their edge actions commute with K_E. The divisor operations are

    C(a,b,c)=(k-a,1-b,1-c),   Q(a,b,c)=(a,c,b).

**Proposition I2 (PROVED).** On every source, including the endpoints,

    I(Cd)=K_E J I(d),       I(Qd)=R I(d).                (4)

Indeed J sends the first a p-slots to the last a; complement in E then
retains precisely the first k-a. The q,r bits are toggled by K_E. R swaps
only their bits. Thus C,Q and CQ are represented by K_E J,R,K_E JR.

These preserve the image and reproduce the exact source fixed sets:

    Fix_D(C)=0,   |Fix_D(Q)|=2k,
    |Fix_D(CQ)|=2 if k even, and0 if k odd.              (5)

For Q one needs b=c and removes the two endpoints. For CQ one needs
a=k/2 and b=1-c. In chart (3), Q reflects g, CQ reflects f_k, and C negates
both coordinates. No source is split or discarded beyond the declared
unit/top deletion.

K_E alone generally leaves the prefix-ideal image: the complement of the
single first p-slot has the other p-slots without the first. The reversal
J repairs exactly this failure. It is therefore part of the operation,
not a relabelling that can be silently forgotten. The companion graph
bridge proves why bare graph complement cannot be retained by an order
embedding of D_k for any k>=2.

## 3. A square is enough for all ten divisors when k=2

Take four vertices0,1,2,3 and host E={01,12,23,03}, the edges of a square.
Use slot assignment

    p_1=12, p_2=03, q=01, r=23.

Here J=(01)(23) exchanges the p-slots and fixes the q,r edges, while
R=(03)(12) exchanges q,r and fixes both p-slots. They are reflections of
the square, and JR=(02)(13) is its half turn. Formula(1) gives:

| Divisor | Host edges |
|---|---|
|p|12|
|p²|12,03|
|q|01|
|r|23|
|pq|12,01|
|pr|12,23|
|qr|01,23|
|p²q|12,03,01|
|p²r|12,03,23|
|pqr|12,01,23|

This is an inclusion-faithful ten-state family. Its rank levels have
sizes3,4,3. These are **edge-count levels**, distinct from the chart's
horizontal g-levels, whose sizes also happen to be3,4,3. For example p
has rank1 and g=0, while p²q has rank3 and g=1. The two partitions are not
the same partition of the source.

The relative-complement qualification is visible here. The half turn JR
has no fixed edge in the square; it fixes the two omitted diagonals02,13
as unordered pairs. Thus K_E JR can fix the images of pq and pr without
contradiction. If the ambient universe were all six edges of K4, either
diagonal would require its bit to equal its opposite. Ordinary graph
complement cannot replace K_E in (4).

**Proposition I3 (PROVED; sharp coordinate budget).** Every order embedding
of D_k into subsets of an N-element set requires N>=k+2. The same bound
holds for L_k, and (1) attains it.

For the proper set, no image can be empty: it would be below every other
image and hence represent a least divisor, which D_k lacks. Similarly no
image can be the full coordinate set. The chain

    p < p² < ... < p^k < p^kq

has k strict inclusions. It starts with at least one coordinate and ends
with at most N-1, so N-1>=1+k. For L_k, its full maximal chain has k+2
strict inclusions and gives the same bound. At k=2 a simple graph on at
most three vertices has at most three edge coordinates, so four vertices
are necessary. The square realizes this minimum. No minimum vertex count
for general k is asserted for the star construction.

## 4. The occupancy quotient explains what a bare coordinate set forgets

Start instead with every subset S of the same k+2 slots, and map it to

    pi(S)=(|S intersection {p_1,...,p_k}|,q_bit,r_bit).  (6)

Two subsets have the same image exactly when a permutation of the p-slots
relates them. Thus L_k is the orbit poset of the Boolean lattice by S_k;
its order is the existence of comparable representatives. A fibre at
(a,b,c) has size binom(k,a). The map is order-preserving, rank-preserving,
and satisfies

    pi(K_E S)=C pi(S),    pi(R S)=Q pi(S).              (7)

Unlike the faithful section I, the quotient transports **bare** relative
complement. But it does not preserve intersection: singleton p_1 and
singleton p_2 have the same image p, while their intersection has image1.
Union fails by the same pair, whose union maps to p². The prefix section
satisfies pi I=id and restores both lattice operations, at the cost of
the twist J in (4).

For k=2, all sixteen subsets yield twelve divisors, with four fibres of
size2 and eight of size1, exactly as in the inherited Oct4 note. Removing
empty/full gives fourteen subset states over ten proper divisors, with
four double fibres and six single fibres. This natural Boolean cover is
different from the action-minimal twelve-state cover of the Oct21 packet.

There is an elementary coordinate-subspace reading: over any field with a
fixed basis indexed by the slots, S selects the span of those basis vectors.
Intersection, sum, and basis-coordinate complement reproduce the Boolean
operations. The prefix section and the p-basis permutation quotient then
give precisely (1) and (6). The basis and the chosen coordinate subspaces
are essential. This does not identify the full subspace lattice, assert
an unmarked linear equivalence, or supply a theory-independent F1 variety.
It is compatible with the distinction between a marked carrier and a
count-only model in [Session59](../session59_f1_geometry_and_information_20261025/SYNTHESIS.md).

The seven squarefree nonunit divisors at k=2 are the subfamily a<=1 with
rank counts3,3,1. Divisor complement C leaves this family: C(p)=pqr stays,
but C(q)=p²r does not. Any 3-1-3 rearrangement or replacement involution on
that family must therefore be specified separately; it is not the
restriction of (4).

## 5. Connection contract and validation

Source: exponent-labelled divisors with divisibility and the actual C,Q.
Target: labelled spanning subgraphs of a marked host with edge inclusion.
Map: the prefix-edge ideal (1). Preserved: source identity, divisibility,
ambient gcd/lcm, rank, faithful chart, and joint involution action via (4).
Lost on forgetting marks: prime identity, repeated-prime order, and the
host edge universe. Missing structure is not supplied by graph isomorphism.
The shortest hostile is complementing a first p-slot without J, or using
the omitted diagonal under full K4 complement.

[The checker](check_divisor_ideal_graph.py) exhausts k=2..6, all full-lattice
source pairs for gcd/lcm and the order equivalence, all source involutions,
all Boolean subsets for quotient fibres, and actual square vertex
permutations. It separately checks that rank and g partitions differ.
Normal and optimized modes must produce the same JSON receipt.

The [saved receipt](divisor_ideal_graph_receipt.json) records2160 source-pair
checks,100 source-action checks,496 quotient states, and12 square controls.
Normal and optimized stdout agree. Independent algebra review by the
session's graph referee confirmed (1)–(4), the square permutations, and
the minimum-coordinate proof; gcd/lcm boundary scope is explicit above.

    python 05-knowledge/results/reflection_bridge_20261025/check_divisor_ideal_graph.py
    python -O 05-knowledge/results/reflection_bridge_20261025/check_divisor_ideal_graph.py

The structural result closes the marked ideal-graph representation and
its minimum edge budget. It gives no equality between unmarked tournament
classes, Hamilton-path counts, or first-exit dynamical certificates.
