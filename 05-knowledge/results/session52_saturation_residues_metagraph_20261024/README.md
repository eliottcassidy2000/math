# Session52: saturation patterns, residue filters, and metagraph carriers

24 October 2026. Second checkpoint: the finite run correspondence and its boundary, exact
balanced metagraph port counts, marked regularization, and the MRB/CM residue
comparison. Every listed lane has an independent audit.

The [pattern theorem](fortunate-square-pattern/THEOREM.md), with its
[independent audit](audit-fortunate-square/AUDIT.md), recovers

    3,2,3,5,1,1,2

as the run lengths of F_n=p_(n+1) versus F_n>p_(n+1), for1<=n<=17.
Here F_n is the primorial-indexed Fortunate number and p_(n+1) is the next
indexed prime, not the number p_n+1. These are saturation bits, not the
numerical gap values. The equality-index sequence is already OEIS A035346.

For the square-sum graph G_N on1..N, join distinct vertices whose sum is a
square. The same17 bits occur on10<=N<=26 for COMPONENTWISE traceability:
minimum vertex-disjoint path-cover number equals number of components.
Ordinary whole-graph Hamiltonicity does not match at10,11,12 because those
graphs are disconnected. The graph-side leading3 and final2 are clipped by
the selected window. This is an explicit interpretation, not a claim that
it was the only possible meaning of the original prompt.

The very next comparison fails: F18=109>p19=67, while G27 has35 Hamiltonian
paths up to reversal. An exact factor witness proves the smaller arithmetic
candidate composite. All graph deficits on1..28 are certified, and all
nonzero ones in this range equal1. No universal arithmetic/graph bridge or
Fortune-conjecture proof is inferred.

A genuine all-size result comes from the graph's birth operation. Adding a
vertex changes minimum path-cover number by at most1 in either direction.
It drops by1 exactly when a minimum old cover has two DISTINCT paths with
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

At n=12, independent coefficient algorithms give10,188,802,555,033 balanced
source ports and1,606,912,987,280 internal C-lines, of which3,401,824 are blue
and1,606,909,585,456 are black. No target class atlas is used. Endpoint splits,
blue boundary counts, divisibility and rotation checks are retained. Original
metagraph C-lines are NOT silently quotiented again by reflection.

These are exact port/line counts, not A000568, loops or complete wiring.

## A marked regular carrier and its quotient loss

The [regularization packet](marked-regularization/THEOREM.md),
[independently audited](audit-marked-regularization/AUDIT.md), extracts
structure from the known Type-I construction of Hatem et al., Lemma4.3.
A source arc becomes an edge-disjoint directed four-cycle in a larger regular
tournament. A source Hamiltonian path becomes a distinguished alternating
Hamilton cycle. Keeping the central vertex preserves the source class and
its ordinary automorphism group.

Forgetting that center is not generally faithful. The first ordinary
collision occurs at4->9 and the first converse-merged collision at5->11:
1+1+C3 and1+C3+1 have the same unmarked regular image. The source complement
also differs from full target noncycle complement on5(n-1) target edges.
The source frame and the actual operation cannot be discarded.

## MRB, mod15, and the genuine CM bridge

The [MRB residue note](mrb-supersingular/RESEARCH.md) and
[independent audit](audit-mrb-supersingular/AUDIT.md) distinguish Ogg's15-prime
set, reduction modulo15, and supersingularity of one elliptic curve.
The m=8 MRB residue identity leads to j=8000 and E:y^2=x^3+4x^2+2x, a rational twist of
the previously studied Bolza elliptic quotient. Its CM splitting is governed
by mod8. The prime17 is Ogg but this E is ordinary;37 is not Ogg but this E
is supersingular. No identification of the actual MRB defect follows.

For positive odd d put M(d)=sum_k>=1 (-1)^k[(dk)^(1/(dk))-1]. The exact unit
filter modulo15 is M(1)-M(3)-M(5)+M(15), approximately0.537584423937395341,
not8/15 times the MRB constant. The mask retains period30 because parity is
part of the state. Dropping parity instead produces divergence with leading
term(4/15)(log N)^2. All eight unit residues modulo15 contain both Ogg and
non-Ogg prime witnesses; no residue-only classification is implied.

The old CM curve, score product and orbit-stabilizer mechanisms are explicitly
credited. New extractions retain their carrier, action, and witness data; no
Fortune conjecture, modular MRB-defect identity, planar-Jacobian conclusion,
or all-size native wiring theorem is claimed.
