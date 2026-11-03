# Golden component words under actual path complements

**PROVED restricted all-size statements.** The Fibonacci word classification,
the nonstrong-to-strong complement lemma, and the general bounded-cutwidth
query compiler are inherited. Their combination below gives an explicit
seven-symbol source language and its literal complement image, not a closed
classification under arbitrary repeated re-marking.

## C1 What pairs and what does not

Fix the physical Hamilton path P=(0,1,...,n-1). Its tournament charts have
m=binom(n-1,2) free nonconsecutive edge bits. The actual complement K_P
retains path arcs and toggles every other arc, hence is a fixed-point-free
involution of charts for n>=3.

For an ordinary tournament class T, its chart fibre has size

    h(T)=H(T)/|Aut(T)|.

Indeed charts identify Hamilton words under the free action of Aut(T).
Both H and |Aut| are odd. The first is Redei's theorem; the second follows
because a tournament has no nonidentity involutory automorphism (a
transposed pair would reverse its own edge), followed by Cauchy's theorem.
Consequently all h(T) are odd. These facts and the parity conclusion are
inherited from
[THM-1420, no linear tournament invariants](../../../01-canon/theorems/THM-1420-no-linear-tournament-invariants-and-the-chirality-that-survives.md).

Quotient each complement pair to its two ordinary classes, retaining loops
and multiplicities. The resulting line multigraph has degree h(T), with a
loop counted twice. Handshaking makes the number of vertices even. This
recovers A000568(n) even for n>=3 without assuming that converse has no
fixed classes. The exact incidence convention is inherited from
[THM-646, line metagraph score skeleton](../../../01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md).

At n=4 the four classes are TT4, 1->C3, C3->1 and Q4. Their degrees are
1,1,1,5. Each of the first three has one line to Q4, and Q4 has one loop.
There is no perfect matching: any matching can cover only one of the
three leaves. The eight charts nevertheless have four disjoint complement
pairs. This is a minimal order for this obstruction; at n=3 the two classes
are joined by one line. An arbitrary pairing of an even finite set exists,
but does not by itself represent native tournament operations or symmetry.

The two black lines are converse twins; the transitive-to-Q4 line and the
Q4 loop are blue. Converse fixes TT4 and Q4 and swaps the remaining classes.
It is not the same involution as K_P. The SVG and PNG are generated from the
literal eight-chart receipt, not from an assumed symmetric drawing.

**Incoming constructive completion.** The independently authored
[odd-fibre cancellation lemma](../tournament_joint_parity_20261102/parity.md#32-odd-fibre-cancellation-lemma)
does pair ordinary classes, using additional re-marking steps. Order all
charts in each odd fibre, leave its least element unpaired, and pair the
remaining successive elements. Overlay this within-fibre matching with
the complement matching. Every remaining representative is an endpoint
of a finite alternating path to another representative. Those paths pair
all classes. The full ordering of each presentation fibre is required;
this is not a single-line involution of the class quotient.

In our bit convention, the n4 fibres are Q4={0,1,2,4,5}, TT4={7},
1->C3={3}, C3->1={6}. Within Q4 pair1 with2 and4 with5. The exact trails are

    0 --K-- 7;
    3 --K-- 4 --re-mark-- 5 --K-- 2 --re-mark-- 1 --K-- 6.

Thus TT4 pairs with Q4, and the two nonstrong converse classes pair with
each other. Our chart bit convention is the complement of the incoming
packet's convention, so its numerical mask choices must not be copied
without transport. The construction preserves vertex count but cannot
preserve H on every pair: at every n the transitive class is the unique
H=1 class. The incoming packet proves this obstruction and a refined
pairing preserving self-converseness. Both statements are inherited here.

## C2 The inherited golden family and its presentation sidecar

Let G consist of all ordered SCC words in the singleton, C3 and Q4.
Their respective vertex costs, Hamilton counts, automorphism orders and
chart multiplicities are

| Atom | Vertices | H | Aut order | Fixed-path charts | Blue charts |
|---|---:|---:|---:|---:|---:|
| Singleton |1|1|1|1|1|
| C3 |3|3|3|1|1|
| Q4 |4|5|1|5|3|

The ordinary class language and Fibonacci codec are proved in
[SCC arithmetic section2](../scc_basins_20261020d/arithmetic.md).
For a word with r singletons, a triangles and b Q4 atoms, the total vertex
count is r+3a+4b, H=3^a5^b, Aut order3^a, and chart count5^b. There are
(r+a+b)!/(r!a!b!) ordinary types with those multiplicities. Ordered SCC
decomposition gives all these statements uniquely; components cannot be
permuted across their dominance order by a tournament automorphism.

The exact scalar spectrum of this restricted family is therefore3^a5^b,
with minimum order3a+4b (and singleton padding at every larger order).
The forest criterion from
[CONFLICT_IMAGES section2](../session81_counterexample_families_20261028/CONFLICT_IMAGES.md)
is what makes this a complete structural stratum: no natural odd-cycle
conflict graph P3 or larger tree can supply a missing Fibonacci extension.

Replace each ordinary Q4 symbol by one of its five fixed-path charts.
Together with one singleton chart and one C3 chart, these form a seven-symbol
weighted alphabet. Its exact marked source count is

    M(x)=1/(1-x-x^3-5x^4).                              (C1)

This records an actual orientation sidecar. For a word with b Q4 atoms,
all5^b charts have distinct K_P outputs because K_P is an involution.
Thus any representation required to output the literal complemented
relation in the given path frame must distinguish those5^b inputs: at
least ceil(b*log2(5)) fixed-length binary bits in the worst case. This
lower bound concerns named, path-framed targets. No such claim for
unnamed target isotypes follows; multiple outputs can be isomorphic.

## C3 Exact image, colours and bounded width

Define the image family G^K_n as the path-marked tournaments whose actual
K_P complement belongs to G_n. This is a size-preserving bijection of
marked objects. Recognition is effective: apply K_P, find its SCCs in
their unique order, and require each to have size1,3 or4. Its inverse
also recovers all seven-symbol source data. Discarding P destroys this
specified inverse.

For n>=5, every member of G is nonstrong. Its path starts in its first SCC
and ends in its last, so the endpoint arc points forward. K_P reverses that
nonpath arc and closes the preserved Hamilton path into a directed cycle.
Every target is therefore strong and outside G. This is the inherited
[strong-core cut lemma, section2](../merged_metagraph_purity_and_strong_core_20261020/STRONG_CORE_COUNT_RECURSIONS.md),
specialized to G. In particular G and G^K are disjoint at n>=5; a golden
SCC recurrence alone is not closed under actual complements.

The involution given by converse plus path reversal reflects the ordered
SCC word. A fixed chart consists of arbitrarily paired outer atoms, with
an empty middle or one blue atom. Thus the blue source series is

    B(x)=(1+x+x^3+3x^4)/(1-x^2-x^6-5x^8).             (C2)

This is the restricted form of the inherited mirrored-component recurrence.
Since K_P commutes with this reflection, it preserves the blue/black label.
For every n>=5, exactly [x^n]M actual lines leave G, of which [x^n]B are
blue and [x^n](M-B) black. Every such line has just one endpoint in G,
so no factor of2 is required. These are marked line counts, not distinct
strong classes. For n=9 they are184 total,10 blue,174 black.

**Bounded-width consequence.** Supply every target with the reversed order
n-1,...,0. Its backward-arc cutwidth is at most4. At a cut between source
SCC blocks, all cross-block target arcs point forward in this new order
except the one retained adjacent path arc: width1. At a cut inside a
source block, no adjacent path arc crossing a different block boundary
crosses this cut. All backward crossing arcs are internal to that block,
whose size is at most4; at most2*2=4 arcs cross. This proves the bound.

The bound is sharp in this specified order, including strong targets:
one Q4 source chart has rows(2,4,9,3) as out-neighbor bitmasks; its K_P
image is the transitive chart(14,12,8,0), with reversed-order cutwidth4.
Appending a singleton to that Q4 source retains the same internal cut
and gives a strong five-vertex target, hence still width4. This is not a
lower bound on minimum cutwidth after optimizing the vertex order.

Apply the inherited
[DEFECT_TRANSFER D2 through D3](../overnight_rank_return_20261029/DEFECT_TRANSFER.md):
any fixed native k-port entry for these literal targets is computable in

    n * 2^O((k+4)log(k+6)) arithmetic operations.

The transfer must retain paired component endpoints, active vertices,
expired anonymous tails and completion guards. The seven-symbol grammar
supplies full target arcs; it is not itself a minimal carrier for arbitrary
joint weights or repeated re-marking. Cost is arithmetic, not a unit-cost
claim for arbitrarily large integer bit operations. Summing endpoint queries
can compute H; no scalar linear-time bound for all endpoint pairs is asserted.

## C4 Spectrum examples and the next obligation

All marked sources at n=5 give target counts9,11,13,15. At n=6 the set is
15,17,19,23,27,29,31,33. These are FINITE-EXACT, with full histograms in the
receipt. Thus the actual complement crosses from a3,5-only multiplicative
spectrum into prime lanes using a uniform source language. It does not
show that every required prime or7p occurs.

A useful next problem is to classify the H-values of G^K by a transfer that
composes its seven symbols and retains end conditions. Its fixed-width
structure gives a concrete route; its scalar value spectrum is still open.
Any proposed smaller transfer must distinguish the five Q4 phases whenever
the requested observer can see their different literal output relations.
