# Actual complements admit a small moving boundary; component pairing survives

2026-10-29, 12:00 UTC continuation. **PROVED all-size counting theorem for a
supplied ordered tournament; FINITE-EXACT controls and full twelve-vertex
three-run responses.** Same-agent review with independent counting and
elimination algorithms. No minimal general metagraph carrier, all-size
selected-path theorem, new A000568 coefficient or external priority is claimed.

## D1. Inheritance, sources and the operation

The closest mechanisms are the March interval-path bijection, the October
[future-boundary propositions 3--4](../tournament_future_boundary_20261020/README.md#5-a-marked-sparse-defect-grammar-for-the-tribonacci-family),
and [RUN_RECURRENCE R5](RUN_RECURRENCE.md). The exact, large recursively
closed literal-target carrier is already
[Session69 P2--P3](../session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md).
Its full complemented relation and ordered component endpoints cannot be
replaced by the scalar queries developed here.

The hostile is the twelve-vertex 81/27 pair. The corrected near miss is that
a degree law or two endpoint bags need not retain which endpoints belong to
the same path. The least-used sidecar is the October ordered boundary and
expired-coordinate mechanism, recovered alongside the capped-run carrier.
Anchor: actual complement queries. Niche: bounded crossing defects rather
than bounded total defects. Wildcard: a determinant connectivity test can
accept the wrong port matching. The live board is full relation / ordered
cut / paired components / endpoint expiry / integral counts / common phase.

Fix vertices 0,...,n-1 in a supplied order and write

    T = B_up xor D,

where B_up has i->j for i<j, and D contains the pairs reversed in T. At cut
t, put

    c_t = #{(u,v) in D: u<t<=v},        c=max_t c_t.       (D1)

This is the backward-arc cutwidth of this order. It is not the cutwidth of
the underlying complete graph, a minimum over all orders, a feedback-edge
count, or a maximum backward degree. Vertices and the pairwise arc predicate
are literal; no ties are introduced. The order is an orientation gauge.

The consumer is a single entry M_k[s,t]: spanning covers by k ordered,
nonempty directed paths, path i starting at s_i and ending at t_i. Starts
are distinct, as are ends; coincidences across the two lists are retained.
If s_i=t_j with i!=j, the entry is zero. This is an actual realization count,
not an arbitrary tensor of port weights.

## D2. The endpoint-expiry lemma and exact transfer

Process vertices in increasing order. The selected arcs on the processed
prefix must form a disjoint union of directed paths. Retain each component
as its ORDERED pair (head,tail), together with the set of completed port
pairs. A head is either its fixed source name i, or an actual vertex still
waiting for an incoming arc. A tail is either its fixed terminal name i,
an actual vertex, or one anonymous open-tail symbol.

Let a processed vertex be active when it meets an unprocessed defect edge.
There are at most c active vertices. The following expiry rules are exact:

1. An open head with no future defect neighbour is impossible: every future
   vertex is larger, so an incoming arc would have to be a backward defect.
2. An open tail with no future defect neighbour has an arc to EVERY future
   vertex. Its vertex name can be forgotten; its component association cannot.
3. A component with fixed source and terminal can be removed only when the
   two port names agree. Record that completed name.

To introduce v, choose at most one existing open tail as predecessor and
at most one existing open head as successor. An actual predecessor u is
legal exactly when (u,v) is not a defect; an anonymous tail is always legal.
A successor u is legal exactly when (u,v) is a defect. A designated source
has no predecessor; a designated terminal has no successor. Choosing neither
leaves the corresponding endpoint open. Join the selected components through
v, reject a choice using the SAME component on both sides, then apply the
expiry and matching rules. Equal resulting states have their integer counts
added. Initially the component list and completion set are empty; finally
every port must be complete and no open component may remain.

**Theorem (PROVED).** This transfer computes M_k[s,t] for every ordered
tournament, every n>=k>=1 and every valid pair of port lists.

**Proof.** Every cover restricts to a path forest on each prefix. Each arc
is introduced exactly once, when its larger endpoint appears. Its prefix
history therefore makes exactly the choices above. Conversely each accepted
history gives indegrees and outdegrees at most one, the prescribed sources
and terminals, no directed cycle, and no unfinished component. It is the
required cover. The only forgotten vertex names are tails whose entire
future adjacency row is constant. They have identical future choices, but
their retained component heads still distinguish the cycle and port-matching
tests. Thus expiry merges histories with the same continuation counts,
without changing their multiplicities. This proves both directions. QED.

Distinct components have distinct retained heads. An anonymous tail has no
independent multiplicity factor: its possible uses are counted once for
each component to which it belongs. This is why retaining only a number
of open tails would be insufficient.

## D3. A bound independent of the total size

At a cut every unfinished component has either one of k source names or
one of at most c active vertices as head. Hence there are at most k+c
components. A head's tail is absent, active, terminal-named or anonymous.
A sufficient bound on the number of states is

    N(k,c) = 2^k (k+c+2)^(k+c).                         (D2)

The bound deliberately includes impossible states. Each insertion has at
most (k+c+1)^2 choices, followed by component normalization. Consequently a
query uses

    n * 2^O((k+c) log(k+c+2))                           (D3)

integer arithmetic operations, and 2^O((k+c)log(k+c+2)) active count entries,
besides the input and its O(n+|D|) schedule. For fixed k,c the arithmetic
operation count is linear in n. This is not a bit-linear claim: integers
grow, and labels/schedules also occupy space. The order and defect list are
supplied; finding a best order is a separate problem. No optimality of D2
is asserted. Evaluating every entry separately gives the coarse polynomial
bound n^(2k+1) * 2^O((k+c)log(k+c+2)) for fixed k,c.

The original interval tournament S_n has width one in reversed vertex
order. Thus the earlier quadratic-in-n single-query algorithm can be
replaced by this transfer in that supplied frame. Its exact size 20,40,100
values from RUN_RECURRENCE agree; the measured peaks for these specific
queries are six states. This does not make six a universal minimum.

## D4. Every first interval complement has width at most two

Let P be any Hamilton path of S_n. Its inherited run composition consists
of ascending consecutive intervals, read from highest to lowest, with no
adjacent singleton runs. The actual complement is

    U=S_n^op xor E(P)=B_up xor (I xor E(P)).             (D4)

At each boundary between two consecutive label intervals, I xor E(P) has
the adjacent cut edge and the hinge joining the low interval's first vertex
to the high interval's last vertex. Inside a label interval, at most the
hinges to its two neighbours cross a cut. At an interval boundary, only its
own hinge and adjacent cut edge cross. Therefore c<=2. It is exactly two
if there is more than one run, and zero for the one-run path.

**Corollary (PROVED).** For fixed k, a specified joint-port count of EVERY
actual first complement of S_n is computable in O_k(n) arithmetic operations,
given its run composition and named ports. This includes the strong target
outside the interval family from RUN_RECURRENCE R6. The construction keeps
the literal target rather than applying the old interval formula to it.

For the composition (2,2,...,2), starts (0,1,2), ends (n-3,n-2,n-1):

| n | Exact three-run count | Peak states | Transfer contributions |
|---:|---:|---:|---:|
|20|38336592|35|1550|
|100|9755223338970122984743819506643436068412005776811746332|35|12510|
|1000|589 decimal digits; full integer in receipt|35|135810|

These are evaluations of the proved transfer. Independent direct checks
stop at the sizes specified in D8; no direct thousand-vertex census is claimed.

Repeated actual complement updates the SAME relation by

    (baseline direction,D) -> (opposite direction,D xor E(P)).  (D5)

Reversing the vertex order restores the increasing-baseline convention.
In that transported frame the width is at most c+c(P), where c(P) counts
path edges crossing a cut. No fixed width-two closure follows. Already
at four vertices, S4 with P=(2,3,0,1), then Q=(0,2,1,3), has second target
rows [12,9,2,4], strong, with decreasing-baseline defects
{(0,2),(0,3),(1,3)} and width three. This is a failure in the retained order,
not a lower bound on its best possible order or its unnamed carrier.

## D5. Pairing, cycles and Frobenius cannot be discarded

The actual S5 complement along P=(1,2,3,4,0) has rows [12,29,24,16,1].
Its query s=(1,0,2), t=(3,4,2) is zero. A deliberately damaged transfer
that independently sorts the head and tail bags returns one. It invented
a connection by changing the component pairing. This is a concrete failure
of that quotient, not a minimality theorem among all representations.
Here vertex 2 must be a singleton, so the requested remaining paths are
the edges 1->3 and 0->4. The latter is absent. The crossed edges 1->4 and
0->3 both exist, making the lost matching completely explicit.

Degree constraints alone also fail: C3->1, rows [10,12,9,0], has a degree
assignment consisting of the triangle cycle and singleton 3, but the
spanning one-path query s=t=(3) is zero. The same-component rejection is
therefore load-bearing.

The transfer is covariant under relabelling the source, its order and BOTH
port lists together. In Paley7, multiplication by two on exponent labels
is the inherited F8 Frobenius action. The receipt checks every simultaneous
three-port transport and records a nonzero query changed by transporting
only its starts. Common phase is retained; the orientation frame is not a
free reset. No new F64/F256 transport theorem is asserted.
Specifically s=(0,3,5), t=(1,6,2) has count one; replacing s alone by
(0,6,3) gives two.

## D6. Full three-run matrices of the 81/27 pair

**FINITE-EXACT.** The literal twelve-vertex sources in atom_precision.hostile
give the following full 1320-by-1320 native response matrices:

| Automorphism order | Three-run mass | Maximum entry | Rank over F3 |
|---:|---:|---:|---:|
|81|80614764|396|430|
|27|80614764|540|435|

The two-run ranks were 24 and 22, respectively. Thus the comparison reverses
at three runs. Equal mass and a numerical symmetry ordering cannot replace
the joint response. The maximum-entry difference is itself invariant under
port relabelling; the rank separator is not an artefact of the chosen names.

The C++ census visits every permutation with at most two forbidden joins,
and inserts every legal pair of cut positions. There are 20,338,857 such
permutations for each source. A third forbidden join makes a three-run
cover impossible, which proves the pruning rule. uint64 arithmetic is safe:
even the crude 12!*binom(11,2) bound is below 2^64.

An independent subset convolution proves the full mass, a mask/last/defect
DP proves the permutation count, and selected entries (including a maximum)
agree with both D2 and subset convolution. Two exact F3 eliminations agree:
NumPy column pivots and Python-integer two-bitplane row pivots. Every vertex
automorphism and both S3 generators preserve the full matrix in a common
frame. Binary matrix hashes, pivots and selected entries are recorded;
matrices are regenerated in a temporary directory rather than committed.

This is the first full three-run calculation for this pair in this packet.
It does not prove a full Smith form, minimality, or successful integral
complement recovery at twelve.

## D7. Algorithm manuscripts: connectivity is not prescribed linkage

The primary [Li, ICALP 2026, Lemma 20 and Theorem 21](https://drops.dagstuhl.de/storage/00lipics/lipics-vol374-icalp2026/html/LIPIcs.ICALP.2026.138/LIPIcs.ICALP.2026.138.html)
uses a determinant of a permutation Laplacian minor to recognize one cycle,
then builds derivative data for Hamiltonian-cycle counting. That published
method motivates an exact boundary test here; its complexity bound is not
imported into our algorithm.

Close a collection of k paths by the fixed cyclic splice t_i->s_(i+1).
If its actual endpoint matching is pi, the resulting permutation of path
indices is gamma*pi, with gamma the k-cycle. Exactly (k-1)! matchings pi
make it one cycle. The prescribed identity matching is unique only for
k<=2. At k=3, both pi=id and pi=gamma survive. In the transitive six-vertex
tournament with starts (0,1,2) and ends (3,4,5), each matching is realized
by three single-edge paths. The prescribed query is one, but the unmarked
one-cycle splice accepts two. Thus a connectivity determinant cannot replace
the port matching guard. This elementary counting argument is PROVED; the
receipt enumerates the port permutations through k=7.

The Big O Speedrun [Subset Sum entry](https://bigospeedrun.com/subset-sum)
leads to [Bills's pinned conditional manuscript, Lemma 2](https://raw.githubusercontent.com/yannouli/subset-sum-conditional-refinement/bf1f50f066f2c71a53293791f693f505a67a1540/paper.tex).
Its useful local mechanism is that an exactly zero certificate multiplier
allows a lower boundary to move without paying that row. The manuscript
explicitly keeps its upstream implementation and moment assumptions.
It also requires common phases in checksum and alias subtraction. Those
are source-checked claims, not an independently verified full Subset Sum
algorithm or a new exponent here. The transferable lesson is precise:
remove a coordinate only after its future effect vanishes. Our anonymous
tail rule supplies such a proof; our component pairing does not vanish.

## D8. Verification and next obligations

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/defect_transfer.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/defect_transfer.py
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/three_run_hostile.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/three_run_hostile.py
```

The last script requires a C++17 g++ compiler; it compiles the committed
source, checks the independently constructed source rows and regenerates
both matrices. All proof gates use exceptions, including under Python -O.

The frontier verifier checks 18,900 entries on 764 source/selected-target
occurrences covering the complete inherited 396 strong classes. It also
checks 37,908 queries over all 227 interval-path complements of sizes 3--9:
every three-port entry through five and twelve fixed-seed queries per target
thereafter. The twelve-vertex two-run mass/rank controls are replayed, with
64 nonzero entries per source independently checked. Long-query receipts,
paired-component erasure, cycle rejection and common Frobenius phase remain
separate controls. Finite tests support the supplied proofs, not all-size
classification by extrapolation.

Incoming 6c45d7f245 supplies the full Paley7 dyadic completion index 2^40 and
a ternary lower bound 3^8 with attainment open. It preserves joint clocks;
our boundary preserves joint components. Neither theorem can be obtained
by applying an independently valid operation to each marginal separately.
Those results are inherited, not rederived.

Next: optimize the sufficient boundary; extend the recurrence to output
literal complemented-target distributions, not only their port counts;
control width across selected repeated complements; and seek an all-size
integral reconstruction proof. The full relation (order,baseline,D) remains
available throughout. No economy for general strong-atom realization or
unnamed re-marking follows without those further obligations.

**Periodic continuation:** [PERIODIC_COMPLEMENT](PERIODIC_COMPLEMENT.md)
extracts a certified minimal order-18 recurrence from the 35-state transfer
for repeated two-vertex runs. Its exact ternary response has period 726 and
a degree-five F243 chart with order-121 multiplication and a necessary
nilpotent coordinate. A precise weighted bulk-observer quotient is minimal
at every modulus; its operation language is narrower than full metagraph
complement and re-marking.
