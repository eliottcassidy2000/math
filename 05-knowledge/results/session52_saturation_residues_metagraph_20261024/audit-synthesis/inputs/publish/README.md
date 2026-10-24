# Session 52: saturation patterns, residue filters, and metagraph carriers

24 October 2026. Session synthesis: the finite run correspondence and its boundary, exact
balanced metagraph port counts, marked regularization, and the MRB/CM residue
comparison. Every listed lane has an independent audit. See the [cross-thread synthesis](SYNTHESIS.md) for the exact common mechanisms, failed transfers, and remaining obligations.

The [pattern theorem](fortunate-square-pattern/THEOREM.md), with its
[independent audit](audit-fortunate-square/AUDIT.md), recovers

    3, 2, 3, 5, 1, 1, 2

as the run lengths of F_n=p_(n+1) versus F_n>p_(n+1), for 1<=n<=17.
Here F_n is the primorial-indexed Fortunate number and p_(n+1) is the next
indexed prime, not the number p_n+1. These are saturation bits, not the
numerical gap values. The equality-index sequence is already OEIS A035346.

For the square-sum graph G_N on 1..N, join distinct vertices whose sum is a
square. The same 17 bits occur on 10<=N<=26 for COMPONENTWISE traceability:
minimum vertex-disjoint path-cover number equals number of components.
Ordinary whole-graph Hamiltonian-path existence does not match at 10, 11, 12 because those
graphs are disconnected. The graph-side leading 3 and final 2 are clipped by
the selected window. This is an explicit interpretation, not a claim that
it was the only possible meaning of the original prompt.

The very next comparison fails: F_18=109>p_19=67, while G_27 has 35 Hamiltonian
paths up to reversal. An exact factor witness proves the smaller arithmetic
candidate composite. All graph deficits on 1..28 are certified, and all
nonzero ones in this range equal 1. No universal arithmetic/graph bridge or
Fortune-conjecture proof is inferred.

A genuine all-size result comes from the graph's birth operation. Adding a
vertex changes minimum path-cover number by at most 1 in either direction.
It drops by 1 exactly when a minimum old cover has two DISTINCT paths with
endpoints adjacent to the new vertex. If the new vertex meets r old
components, the path-cover excess above component count changes by

    change(path cover) -1 +r.

Thus a component merger and a path repair are distinct events even when they
produce the same Boolean run. The proofs retain the actual endpoints and
minimum-cover witnesses. Independent edge-subset computation, primality
certificates, normal/optimized replay, and hostile interpretation controls
are included in the audit packet.


## Balanced metagraph ports without a class atlas

The [balanced carrier theorem](balanced-metagraph/THEOREM.md) and its
[independent audit](audit-balanced-metagraph/AUDIT.md) refine the inherited
score partition function (THM-554/555), not claim that old product as new.
For even n>=4, a balanced source stays balanced under native path complement
exactly when its actual path starts at a high-score vertex and ends at a
low-score vertex. These internal ports biject, respecting C, with Eulerian
orientations of K_(n+1) with one fixed Hamilton cycle removed. Native C
becomes reversal of precisely those remaining edges.

At n=12, independent coefficient algorithms give 10,188,802,555,033 balanced
source ports and 1,606,912,987,280 internal C-lines, of which 3,401,824 are blue
and 1,606,909,585,456 are black. No target class atlas is used. Endpoint splits,
blue boundary counts, divisibility and rotation checks are retained. Original
metagraph C-lines are NOT silently quotiented again by reflection.

These are exact port/line counts, not A000568, loops or complete wiring.

## A marked regular carrier and its quotient loss

The [regularization packet](marked-regularization/THEOREM.md),
[independently audited](audit-marked-regularization/AUDIT.md), extracts
structure from the known Type-I construction of Hatem et al., Lemma 4.3.
A source arc becomes an edge-disjoint directed four-cycle in a larger regular
tournament. A source Hamiltonian path becomes a distinguished alternating
Hamilton cycle. Keeping the central vertex preserves the source class and
its automorphism group as the center stabilizer.

Forgetting that center is not generally faithful. The first ordinary
collision occurs at 4->9 and the first converse-merged collision at 5->11:
1+1+C3 and 1+C3+1 have the same unmarked regular image. The source complement
also differs from full target noncycle complement on 5(n-1) target edges.
The source frame and the actual operation cannot be discarded.

## MRB, mod 15, and the genuine CM bridge

The [MRB residue note](mrb-supersingular/RESEARCH.md) and
[independent audit](audit-mrb-supersingular/AUDIT.md) distinguish Ogg's 15-prime
set, reduction modulo 15, and supersingularity of one elliptic curve.
The m=8 MRB residue identity leads to j=8000 and E:y^2=x^3+4 x^2+2 x, a rational twist of
the previously studied Bolza elliptic quotient. Its CM splitting is governed
by mod 8. The prime 17 is Ogg but this E is ordinary; 37 is not Ogg but this E
is supersingular. No identification of the actual MRB defect follows.

For positive odd d put M(d)=sum_k>=1 (-1)^k[(dk)^(1/(dk))-1]. The exact unit
filter modulo 15 is M(1)-M(3)-M(5)+M(15), approximately 0.537584423937395341,
not 8/15 times the MRB constant. The mask retains period 30 because parity is
part of the state. Dropping parity instead produces divergence with leading
term(4/15)(log N)^2. All eight unit residues modulo 15 contain both Ogg and
non-Ogg prime witnesses; no residue-only classification is implied.

The old CM curve, score product and orbit-stabilizer mechanisms are explicitly
credited. New extractions retain their carrier, action, and witness data; no
Fortune conjecture, modular MRB-defect identity, planar-Jacobian conclusion,
or all-size native wiring theorem is claimed.


## Complete blue wiring at order ten

The [blue wiring certificate](balanced-blue-wiring/README.md) and
[independent audit](audit-balanced-blue-wiring/AUDIT.md) go beyond a scalar
port count. On all 229 self-converse almost-regular order-ten classes the
induced blue multigraph is connected and has no bridges or articulation
vertices: 11,452 edges, 24 loops,cycle rank 11,224. A455-edge spanning ear
certificate verifies two-connectivity. The inherited 20 positive-error
balanced nodes form a connected dominating set, with 90 core edges, 1,857
cross edges and 9,505 outside edges. Every other node directly meets this
core. The 18,795 blue boundary edges reach 786 distinct classes.

This is not the full black-and-blue metagraph or an all-size connectivity
proof. The four nonbalanced isolated nodes in the prior positive-error
study are not part of this balanced carrier.

## Recovered counting formulas and repaired interpretations

The [counting-lens recovery](counting-lenses/README.md) and
[audit](audit-counting-lenses/AUDIT.md) identify the old fast A000568 baseline,
correct self-converse partition filter, and minimal witnesses separating
single-arc spectral information from native complement-line incidence.
The [precision supplement](counting-precision/PRECISION_NOTES.md) and
[audit](audit-counting-precision/AUDIT.md) narrowly repair trace versus fixed
subspace, unsupported homological wording, and the necessary centering of
the cyclic-triangle eigenfunction. No count table or spectral value changes.

## Fifteen kernels, a real action, and a qualified CM connection

The [fifteen-kernel packet](fifteen-kernels/README.md), with
[independent audit](audit-fifteen-kernels/AUDIT.md), recovers the binary
symplectic carrier and constructs the full six-pentad outer-S6 transport.
The matching graph is KG(6, 2), but not under the ordinary untwisted action.
The cited Bolza kernel split 6+4+4+1 leads through its size-six orbit to the
j=8000 elliptic-square codomain. Its singleton is the loop. Edge-orbit counts
are not destination-vertex counts, and exceptional characteristics remain
explicit. This provides no Ogg-prime bijection or MRB-defect identity.

## All-size existence lead, not an all-size wiring proof

The [blue-port existence probe](blue-existence-probe/PROBE.md) and
[translation audit](audit-blue-existence/AUDIT.md) connect the regular-cycle
carrier to a Salvi-Zagaglia theorem reported in Farrugia's survey. The exact
carrier translation is proved and small orders checked. The original 1980
paper was not accessible, so the all-size consequence remains explicitly
CONDITIONAL pending primary-source verification. Even accepting it gives
internal blue-port existence, not connectivity or compatible global wiring.
