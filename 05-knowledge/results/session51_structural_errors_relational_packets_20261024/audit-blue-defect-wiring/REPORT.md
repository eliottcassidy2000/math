# Independent audit of actual positive-error blue wiring

24 October 2026. **PASS, finite-exact within the stated induced slice.**

The frozen candidate is `session51-blue-defect-wiring`, manifest SHA256
`cd89cad144786c67f076c235409c90e3b4de21e0a5f81926139184068f1d676d`.
The inherited complete positive-error classification is frozen at manifest
`30b62a8fdf93065b19f071050d9fcfe5ad08c5f79cfcfbf30afcdae53babde5a`,
with its independent audit at
`7ca02dc5c0272039b552e7cd5538bbab04923a894e04c050e08fb69de721fc7f`.
Exact copied source hashes are in `SOURCE_PINS.json`. No frozen input was edited.

## Independent reconstruction

`independent.py` imports no producer code and does not use its anti-involution
certificates to generate paths. It enumerates **every directed Hamiltonian
path** on each of the 24 inherited rigid self-converse representatives using
ordinary adjacency DFS, obtaining 341504 paths. It applies the exact
anti-reflection predicate to each full normalized tournament. Precisely 4200
paths are blue, and all 24 per-node full-HP and blue counts agree with the
inherited certified counts.

The independent carrier consists of full 10-by-10 adjacency matrices, not the
producer's 36-bit word integers. Complement pairing is performed by literally
reversing every non-path arc while retaining the nine directed path arcs.
Only after reconstruction are matrices translated to the producer's external
word format to compare records. All 4200 actual words, all 90 internal edge
records with both actual source paths, and all 4020 boundary records agree
exactly. No target class is guessed from degrees or triangle counts.

Checks include the exact reflected matrix predicate, base-path preservation,
C being a free involution, R/C commutation, literal target normalization,
rigidity/injectivity of the entire actual presentation carrier, and every
per-node conservation equation. This independently recovers the following:

- 90 induced edges: four loops at nodes 9, 47, 74, 75, and 86 nonloops.
- 180 internal ports and 4020 boundary edges, so 4200 = 2(90) + 4020.
- 70 occupied unordered pairs; multiplicity histogram 1:52, 2:16, 3:2.
- Multiplicity-three pairs (28,74) and (31,65).
- One component of size 20 and four isolated nodes 0, 1, 2, 3.
- The isolated nodes are exactly the crossing-degree-one sector; the connected
  component is exactly the crossing-degree-two sector.
- Cycle rank 90 - 24 + 5 = 71, no bridges, and no articulation vertices.

Connectivity is recomputed independently with a disjoint-set algorithm.
Deleting each actual edge individually handles parallel edges correctly;
deleting each vertex handles isolated vertices and loops correctly. Loops
contribute two ports to the internal degree. The degree histogram is
0:4, 2:1, 4:1, 6:2, 7:1, 8:1, 9:5, 10:1, 11:4, 12:4.

## Carrier completeness and boundary scope

Exhaustiveness of the 24 positive-error self-converse classes is inherited
from the separately frozen and audited classification. Rigidity makes paths
bijective with their labelled base-path presentations: equal presentations
would yield a source automorphism. Distinct source classes cannot share a
presentation. Therefore the 4200-element carrier contains every blue
presentation of every positive-error self-converse class.

The literal complement of a blue presentation is blue and therefore
self-converse. If its target had positive error, the inherited complete
classification would put it in the carrier. A complement outside the carrier
thus has zero error. This is a completeness argument, not an inference from
an untested invariant. The present audit does **not** identify boundary target
classes, count their distinct nodes, assert their mutual adjacency, recover
black edges, or claim the full order-ten metagraph.

For an additional independent diagnostic, boundary matrices have stable-cell
patterns 1^10 on 3959 ports, 1^4+3^2 on 45, and 5^2 on 16. These are counts of
boundary ports, not distinct target classes. They do not change the scope.

## Endpoint obstruction and insufficiency

The Session40 actual-path degree identity is inherited, not newly claimed:

    d_K(v) = n-1-d_T(v) + 1_(v=s) - 1_(v=t).

Restoring the Hamiltonian-path arcs after taking the full converse cancels
incoming and outgoing corrections at each interior vertex and leaves only
the endpoint corrections. At n=10, positive self-converse target degrees must
be either 3^5,6^5 or 4^5,5^5. A crossing-degree-one source has interior target
degrees 3 or 6 but start target degree 7 or 4, so it cannot remain positive.
This proves the four induced isolated nodes without relying on lookup.

For crossing-degree-two sources, positive targets require start source degree
5 and end source degree 4. This is **necessary, not sufficient**. Independent
literal counts give 2037 eligible ports, of which only 180 are internal and
1857 are boundary. Another 1671 degree-two ports fail the endpoint filter.
The degree-one sources contribute 162 ports from degree 3 to 6 and 330 from
6 to 3, all boundary.

All 4200 literal source/partner pairs satisfy the full degree identity, its
squared-score consequence

    sum d_K^2 - sum d_T^2 = 2(d_T(t)-d_T(s)+1),

and the independently enumerated cyclic-triangle consequence

    c3(K)-c3(T) = d_T(s)-d_T(t)-1.

The candidate properly restricts the distinct-endpoint formulas to n>=2 and
its regular-source consequence to odd n>=3. Preserving scalar invariants does
not replace actual carrier membership. No general wiring claim follows.

## Replays, hostile controls, and artifacts

Normal and optimized Python produce byte-identical independent receipts.
`validate.py` makes actual copied directories under a temporary path, runs
both independent and producer implementations from an unrelated working
directory normally and under -O, and checks complete receipt equality. It
verifies the strict producer manifest both before and after relocation and
confirms the entire original producer inventory is unchanged.

Adversarial tests reject changed pinned input, C incorrectly made the
identity, reversed blue/nonblue membership, a changed candidate theorem,
loops counted only once, a wrong actual target path, a wrong boundary owner,
and the false claim that all 2037 endpoint-eligible ports are internal.
Source/semantic/script mutations are tested normally and under -O; all checks
use explicit exceptions, never removable Python assertions.

`normal.json` retains every actual edge and boundary word, all endpoint path
witnesses, the full adjacency and topology, per-node accounting, and source
pins. `validation.json` records frozen producer hashes and every replay/control.
No publication was performed.
