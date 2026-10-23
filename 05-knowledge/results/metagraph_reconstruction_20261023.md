# What fixes the tournament metagraph beyond its bucket counts

October 23, 2026. **PROVED:** the all-order reconstruction and component-cut
laws cited below. **FINITE-EXACT:** the small-order incidence determination
and counterfeit tests. **OPEN:** an all-order classification by small
invariants, black-core connectivity, and faster general A000568 enumeration.

The main conclusion is that the metagraph is fixed by **compatible actual
ports**, not by its node and line totals alone. This distinction now has
both a concrete counterfeit and an all-order rejection test. It also gives
a specific reason for part of the even/odd structure.

## Inheritance and live concepts

Closest proved mechanism: [THM-4739, Kirchhoff/Hodge and port
determination](../../01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md).
Its low-order profiles and endpoint constraints determine order5 but leave
48,000 transportation assignments at order6; only three project to the
actual merged wiring. Those are assignments, not 48,000 different graphs.
Canonical hostile: a balanced switch among the unresolved port buckets.
Corrected near miss: unsigned triangle counts are not signed labelled curl.
Least-used coordinate: the position of the cut in the ordered strong
components, restored by one vertex score or the marked Hamiltonian path.

The concept board is: bucket margins; incidence overlap; signed curl;
component cuts; two-height columns; allowed future operations. The anchor
is actual port acceptance, the niche is component-necklace parity, and the
wildcard is the connected-graph/curl-fibre correspondence.

## 1. The exact missing information in a triangle description

For a labelled tournament set a_ij=1 if i beats j and -1 otherwise, and
c(i,j,k)=a_ij+a_jk+a_ki. The [reconstruction proof](metagraph_port_rigidity_20261023/PROOF.md)
establishes:

* Equal signed triangle curls mean that the tournaments differ only by
  a cyclic rotation of their ordered strong components.
* One distinguished vertex's outdegree separates all these rotations.
* On a fixed directed Hamiltonian path, the path itself supplies the
  missing cut, so the signed root-triangle data already determines every arc.

The reason is restrictive: the difference of two equal-curl arc arrays is
a gradient, but its entries must lie in {-2,0,2}. Its potential therefore
has at most two levels. A nonconstant potential reverses exactly one
uniformly directed cut. Such cuts are precisely cuts between strong
components. This is an all-order argument, not a numerical pattern.

It also gives a polynomial-time feasibility decoder for a proposed curl
tensor and root score. Locally admissible curl values need not fit together:
the root tensor (3,3,3) would force one shared arc to have both signs.

For a marked Hamiltonian path P, write p_P for its signed path-edge flow.
The nonpath-complement operation is exactly

    a'=-a+2p_P,
    c'=-c+2 curl(p_P),
    s'(r)=n-1-s(r)+div(p_P)(r).

Thus a supplied source port has exactly one accepted target signature.
Normalizing two paths to the same position labels tests their actual
automorphism orbits. A counterfeit attachment fails this equality even
if its contribution to every coarse bucket total is correct. Given the
complete authenticated inventory, these pairings force the merged graph
up to the appropriate relabellings. This does not derive that inventory
from scalar margins or prove connectivity.

## 2. A literal two-slot capacity rule in the tiling triangle

The [column encoding](metagraph_two_height_20261023/README.md) makes the
obstruction more elementary. Along the fixed path define

    d_ij=c(i,i+1,j),     j>=i+2.

The array is realizable exactly when every column suffix satisfies

    sum_{k=i}^{j-2} d_kj in {j-i-1, j-i+1}.

After subtracting its baseline and dividing by two, each suffix has height
0 or 1: the single binary arc at its end. Two individually allowed upward
steps can exceed this capacity. At order4, d_13=d_03=3 would require
a_03=-3, an impossible tournament arc. This is a concrete version of the
user's proposed bucket obstruction.

Complement reflects the interior curls by d -> 2-d and the path-adjacent
curls by d -> 4-d. The difference is an exact boundary term. Reflection
of path positions commutes with complement. For n>=3, with m=(n-1)(n-2)/2,
f=floor((n-1)/2), and k=(m+f)/2, this yields 2^k blue tilings,
2^(k-1) fundamental blue lines, and (2^m-2^k)/4 fundamental black lines.
Parallel lines and loops retain their multiplicity. These totals alone
still do not specify which node ports meet.

## 3. Even and odd buckets come from a reflection of component cuts

Hamiltonian counts and automorphism orders multiply along the ordered
strong components. Consequently rotating their order preserves
h=H/|Aut|. But it need not preserve self-converseness.

For example, the three cuts of the necklace containing C3 and two
singletons all have h=1. The middle-C3 class is self-converse; the other
two cuts form a converse pair. The merged fibres therefore have sizes1
and2 even though the unmerged factors agree. The lost cut changes the
bucket, not its original scalar weight.

More generally, take a necklace of strong-component isomorphism types
with least rotation period d. If converse-reversal is rotation by a,
the self-converse cuts are exactly

    2t=a modulo d.

There is one when d is odd, and zero or two when d is even. Each fixed cut
has odd fibre h; each paired cut has even merged fibre2h. If the converse
necklace is distinct, all d cuts pair with that other necklace. This is
a complete parity law within each such necklace, not a classification of
all blue multiplicities or strong cores.

The same curl quotient has a known graph-theoretic connection: selecting
the unique cut from which a fixed root reaches everything identifies it
with root-connected orientations of K_n. Its EGF is log T(z), where
T(z)=sum 2^binom(n,2) z^n/n!, also the EGF for all labelled simple graphs.
Thus its counts equal connected-labelled-graph counts. The proof packet
identifies the existing embedding-dependent Bernardi bijection; it does
not claim a new canonical unlabelled bijection.

## 4. The small four-vertex structure repairs a genuine information loss

The [incidence packet](metagraph_bucket_obstruction_20261023/README.md)
constructs an order6 counterfeit preserving endpoint margins, existing
coarse target profiles, color and converse reflection. It changes merged
wiring. More detailed joint neighbor-incidence data removes the ambiguity.
This distinguishes a necessary degree/parity test from actual realization.

Let c3 count cyclic triangles, c_v count those through vertex v, and c_e
those through edge e. Then, for every tournament,

    D33 = binom(c3,2) - sum_v binom(c_v,2)
          + sum_e binom(c_e,2).

Here D33 counts unordered vertex-disjoint pairs of cyclic triangles. This
is inclusion-exclusion by intersection size0,1,2. The last term is exactly
the number of induced strong four-vertex tournaments: each such four-set
contains exactly one pair of cyclic triangles sharing an edge.

This makes a precise bridge to the incoming [context/port
work](session35_dimensions_dualities_20261023/recursions/CONTEXT-AND-PORTS.md).
Its cospectral six-tournament pair has equal individual cycle counts but
Hamiltonian counts29 and25. The discrepancy is exactly 4D33. The overlap
coordinate supplied by the strong four-tournament and the missing disjoint
coordinate are affine-equivalent once vertex incidences are retained.

These identities are all-order. The adequacy of particular short
signatures is finite and is recorded separately in the incidence packet;
it must not be extrapolated to every tournament order.

The finite progression is useful: joint score/triangle neighbor data plus
endpoint/color capacities forces all aggregated ordinary-class line
multiplicities at order6; adding each vertex's
strong-four incidence identifies all456 classes at order7 and forces its
aggregated line multiplicities. At order8, that first-round signature has two
collision pairs. The [full wiring and minimal-overlap
extension](incidence_wiring_20261023/README.md) now replays all2^21 fixed-path
ports: J4 leaves884 units of genuinely ambiguous port mass, while one
common-source overlap bit on each exceptional signature separates all6880
classes and fixes every aggregated port count. Even the full incidence table
between first-round profile types still misses these pairs. A complete second
round also works, but is more information than this finite repair needs.

## 5. Compression depends on which operations must survive

The incoming [session33 memory result](session33_continuation_states_20261023/metagraph-memory/NOTE.md)
shows that complement alone needs only the ordered source/target node pair,
with a color coordinate at mixed supports. Its [all-flips
extension](session33_continuation_states_20261023/metagraph-memory/ADDENDUM_FULL_FLIPS.md)
shows that allowing every individually labelled free-arc flip distinguishes
all2^m tilings. A singleton transitive node lets a suitable common flip word
separate any two states.

The present curl/height descriptions preserve all arcs, so they respect
that lower bound. Their advantage is a precise compatibility test and a
description of what goes wrong under forgetting. They do not promise a
smaller exact all-operations state space. The incoming spectral recurrence
work gives the complementary warning: a very small invariant can survive
infinitely many recursive contexts and still miss Hamiltonian counts.

### Final incoming work: compatibility and the next parity layer

The final pull brought [Session36's curl decoder](session36_information_curl_contexts_20261023/root/CURL_REALIZABILITY.md)
and [connected-graph derivation](session36_information_curl_contexts_20261023/root/CURL_CONNECTED_GRAPH_ENUMERATION.md).
They overlap the signed-curl/SCC mechanism developed here. These are
converging descriptions of one mechanism, not separate novelty claims.

Its [unsigned realization phases](session36_information_curl_contexts_20261023/context/REALIZATION_PHASES.md)
also sharpen the bucket question. Requiring positive cyclic triples gives
parity equations, but all other triples must simultaneously remain transitive.
The five-vertex demand of exactly cycles014 and023 passes positive parity
and fails the remaining constraints. Two cycles sharing an edge form S4;
two sharing just one vertex cannot be the only cyclic triples on their
five-vertex union. Positive incidence alone can miss the obstruction.

The incoming [dyadic electrical law](session36_information_curl_contexts_20261023/metagraph/DYADIC-ELECTRICAL-METAGRAPH.md)
writes the literal graph as G=Blue+2PairedBlack. Modulo2 its Laplacian sees
only blue, while the next layer retains the bilinear half-energy
q_G(u,v)=u^T Lap(G)v/2 modulo2 on blue-harmonic potentials.
Our [explicit bridge replay](metagraph_bucket_obstruction_20261023/dyadic_bridge.py)
tests the saved six-vertex counterfeit: its full Laplacian agrees modulo2
and its degrees agree, but q_G(e16,e12) changes from0 to1. Both named
vertices are blue-isolated, so these are admissible potentials. This
connects the incidence obstruction to a concrete electrical obstruction
one dyadic layer deeper. It detects this pinned switch; it is not a claim
that the half-energy form alone characterizes tournament wiring.

The [order-eight extension](incidence_wiring_20261023/README.md) proves the
all-size rank-two switch criterion and supplies a stronger hostile: doubling
its native feasible capacity switch preserves the entire half-energy form,
but changes the quarter-energy form. Its full wiring has3,358 critical-group
factors of exact2-adic valuation1, four of exact valuation2, and one of exact
valuation8, agreeing with the [independent incoming calculation](metagraph_order8_electrical_20261023/README.md).
The incoming [minimal menu](metagraph_minimal_signatures_20261023/README.md)
also removes the score and triangle seeds: v4 alone plus one round and its
guarded bit suffices at8. The extension recovers a historical order13 H collision that remains
invisible to every unrooted incidence-refinement round; marking the actual
source makes both examples discrete in three rounds. Stable cells must be
odd, but need not be singletons.

Next targets: choose only the triangle or incidence coordinates needed to
separate the surviving targets at each ambiguous port; certify their joint
capacities; and exploit the SCC-necklace description when assembling weak
sources from strong cores. Measure any enumeration speedup against its
canonicalization and target-signature costs. All-order short signatures,
black-core connectivity and universal Collatz coverage remain open.

The parallel [Collatz extension](weighted_tail_20261022/extension/README.md)
provides a concrete counterpart to the same preservation principle: it
keeps immutable source owners while aggregating synchronous values, checks
all intermediate cuts, and grounds new children without recycling credit.
This is a shared proof discipline, not an identification of the two systems.
