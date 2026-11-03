# Independent audit of the active-port weighted carrier

**PASS.** The all-order finite-carrier theorem in [ACTIVE_PORTS.md](ACTIVE_PORTS.md)
has the required boundary, exact transition and acceptance bijection. The
finite representation has unbounded integer weights; neither a finite set of
numerical states nor a solution of its target-output membership problem follows.

## Proof audit

**Boundary and clock.** Indegree never decreases. An inactive vertex can never
again be a selected root, so every later newborn beats it. Formula (1.1) is
literal on named vertices. The strong-tournament low-indegree bound applied
to a=d−1 yields at most 2d active vertices after each insertion. The initial
triangle at d=1 is a real exception; fixing its first root and starting the
matrices at the actual four-vertex child handles it correctly. The active
rank tuple preserves chronological port order and needs no hidden past arcs
for its degree update.

**Sufficiency of paired endpoints.** Each partial cover is a set of selected
arcs, with at most one selected incoming and outgoing arc per vertex and no
cycle. Existing selected arcs are never broken. A newborn can receive a
selected incoming arc only from the selected root, which must therefore be
a path end. It can send a selected outgoing arc to every component start
except that root. Forgotten starts all have the same future incidence.
Choosing both ends of the same component closes a cycle and is forbidden;
this is precisely why start/end pairing must survive.

An inactive path end can never receive a successor among future vertices.
Thus two such components cannot ever join, and both cannot survive in one
final Hamilton path. Discarding them is safe at every prefix, including the
moment a formerly active endpoint is forgotten. No other dead-state pruning
is claimed. At most one forgotten end and r distinct active ends imply at
most r+1 components. Endpoint labels, pairing, and multiplicity therefore
form a finite signature carrier.

Removing the newest vertex from a future Hamilton path recovers the unique
previous cover and the chosen incoming/outgoing components. Conversely every
accepted trajectory builds a spanning acyclic directed graph of degree at
most one in and out with one component. These are inverse constructions.
Consequently the transfer counts every actual Hamilton word once. Local
control, endpoint normalization and input letters retain their stated gauges.

**Finiteness and scope.** The recurrence for g_r,f_r correctly separates
absent ports, singleton paths, forgotten-start paths, the unique possible
forgotten-end path and paired active ports. A both-forgotten path accounts
for f_0=2. These are combinatorial upper bounds including empty signatures,
not ranks. Positive-support closure has no cancellation, so its finite
reachable coordinate count is meaningful. The 15/460/9672 values are the
producer's explicit closure receipts; this audit does not independently
recompute those full closures. Degree-control counts 1/6/24 are independently
recomputed. The cap-two/cap-three target conclusions remain the separately
audited bounded-history results, not consequences of finite dimension alone.

## Independent edge-set audit

`active_ports_audit.py` imports no producer or repository modules. Its cover
enumerator differs from both producer routines: it chooses a successor edge
or no successor for each vertex, rejects duplicate incoming edges, and rejects
directed cycles. This enumerates each unordered spanning path cover once,
without permutation/cut deduplication and without subset DP.

For every tested birth, it enumerates **all actual child covers**, deletes the
newborn from each, and groups these lifts by the actual parent cover. It then
checks that every two parent covers with the same endpoint signature have
identical complete distributions of child signatures. This is cover-by-cover
lumpability, not merely agreement of total H. It also proves directly in the
finite bank that a discarded two-forgotten-end parent is never rescued.
One-component edge covers are checked against literal Hamilton vertex words.

The universe is all actual cap-1/2/3 histories from the normalized four-vertex
seed through order six: 6,18,24 birth edges respectively. It includes 144 and
880 parent-cover observations with forgotten endpoints in the first two caps,
and multiplicity two occurs for cap two. A targeted cap-three extra edge
continues history (0,0,0) at root 0 from order six to seven, explicitly testing
forgotten vertices at that cap as well. This adds a missing hostile control,
not a target-spectrum census.

Reproduction:

    python3 active_ports_audit.py > active_ports_audit.json

Normal and optimized runs produce identical JSON. The all-order conclusion
rests on the paired-endpoint proof; the finite controls independently validate
the proposed sufficient state against actual named graphs.
