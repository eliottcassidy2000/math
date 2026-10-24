# Independent extension audit: one global Hamiltonian observer

Status: **PASS at the stated conditional all-size and finite-exact scopes.**
No producer correction is required. This is a separate extension packet;
`session47-audit-wiring-theory/` was left unchanged. No Git operation was used.

## Verdict and necessary scope

The directed-block compiler gives an exact bijection, not a polynomial-time
algorithm or a new implication-only/2-SAT representation. The source tournament
is a full fixed labelled input. Source completion, normalized resources,
electrical currents, and native helper supply remain separate obligations.

The semicomplete condition is sufficient without extra absolute clock pins,
not necessary and not sufficient after arbitrary global clocks are imposed.
The concrete witness disproves one global observer, not the incoming local
partner-host controller. Those distinctions are explicit in the producer note.

## Proof review

### Translating local partners to flags

For each covered pair, a partner retains its source orientation exactly when
the pair is consecutive in the relevant Hamiltonian path. Since the local
source orientations agree with the fixed full source S, equality of every
prescribed local partner restriction is exactly equality of those adjacency
flags. Reversing nonpath arcs introduces no additional pairwise ambiguity.

A required pair receives its directed orientation from S. Requiring both
adjacency and nonadjacency is immediately inconsistent. Two distinct required
successors or predecessors at one vertex are impossible in a directed path;
a required directed cycle is also impossible. Under indegree/outdegree bounds
one and acyclicity, the required arcs are disjoint directed paths, with isolated
vertices as singleton paths.

Every component is forced to occur as a contiguous oriented block: each
required arc makes its two endpoints adjacent in that order. Contracting the
blocks therefore preserves every possible global observer without ambiguity.
A block transition is permitted precisely when its last vertex points to the
next block's first vertex and their pair is not forbidden.

Conversely, expanding a Hamiltonian block ordering yields a source-Hamiltonian
path and preserves every required edge. A forbidden pair inside one block is
never made adjacent unless it was an internal required edge, whose conflict
was checked already. A forbidden pair between blocks can become adjacent only
at a transition, where it is explicitly excluded. Thus no negative flag is
lost. Expansion and contraction are inverse maps.

The endpoint transition carrier can contain both directions between a block
pair or neither. It cannot be assumed to be a tournament merely because S is
one. The independent checker includes explicit examples of each possibility.

### Semicomplete sufficiency

Insert a new vertex into an existing directed path immediately before the
first vertex to which it has an outgoing arc. If this position is not first,
its predecessor must have an arc into the new vertex, by semicompleteness and
the minimality of the chosen position. If there is no such position, every
existing path vertex points into the new vertex, so append it. Induction gives
a Hamiltonian path. Bidirectional pairs cause no problem.

This suffices for the unclocked contracted problem. The audit additionally
checks every insertion order on every labelled semicomplete digraph through
four vertices. A nonsemicomplete carrier can also have a Hamiltonian path, so
this is not a characterization.

### Absolute clocks and converse

If owner v has offset j in a forced block, a global clock requirement c fixes
the block start to c-j. Several pins in one block must agree. A block with m
owners can start only between 0 and |V|-m. Expansion starts each block at the
sum of preceding block lengths, making the cumulative-length constraints
necessary and sufficient. The producer may reject invalid starts through
exhaustion rather than early range testing; this does not change correctness.

An unclocked semicomplete carrier does not guarantee arbitrary prescribed
starts. The transitive three-vertex source with its last vertex pinned to
clock zero is a valid counterexample. The input's cell-local clocks must not
be silently interpreted as global absolute clocks.

Converse reverses every source arc, each forced block, and the block ordering.
The required/forbidden unordered pairs stay the same. A clock c transforms to
|V|-1-c; equivalently a block start s transforms to |V|-length(block)-s.
This proves the claimed converse feasibility correspondence, including pins.

## Actual shared-owner witness

The inherited JSON is authenticated by SHA-256

    5266665af6f013ff54a64a5b4b7ef83ac35a8b80c01507ef05cc9ed26d1b6880.

The audit independently checks that its 14-owner source matrix is a tournament,
that all four local paths are Hamiltonian in its induced cells, and that the
recorded global owner lists are exactly their local lists under the declared
embeddings. Here the JSON field `global_path` means a local path expressed in
global owner names, not a 14-owner Hamiltonian observer.

The authenticated paths are:

    core U:   (5,7,2,3,1,4,6,0)
    core V:   (5,7,6,3,1,4,2,0)
    helper U: (13,6,12,11,2,8,10,9)
    helper V: (13,2,12,9,11,6,8,10).

The exact allowed bit pairs in the source are (0,1) and (1,0). For the first,
required successor arcs include 2->3 and 2->12. For the second, they include
6->3 and 6->12. Either is already impossible for one directed Hamiltonian
path, even before adding forbidden pairs, global clocks, other cells, or
unspecified cross arcs. The obstruction remains valid for any host completion.

At the same time, the two cells' shared-pair flags agree on their shared owners
{2,6}; the audit verifies this directly. Hence the local partner restrictions
can glue while the global-observer condition fails. The incoming
`PARTNER_GLUING.md` explicitly distinguishes these predicates. The witness is
not a counterexample to that gluing theorem or its electrical controller.

## Independent executable coverage

`check_independent.py` imports no producer checker and uses explicit exceptions
rather than Python assertions. Its construction differs from the producer:
weak required-edge components are topologically ordered, block paths are built
by prefix extension, and results are compared with direct vertex permutations.

Both normal and optimized runs pass **494,223 checks**, with identical results
apart from elapsed time and the optimization flag:

- All four-valued pair-flag assignments (unspecified, required, forbidden,
  contradictory) on every labelled tournament of orders 1..4: 262,665 total
  instances, including all 262,144 order-four instances and 215,786 contradictory
  instances across the stated orders.
- 103,743 selected clock configurations, including two simultaneous owner pins,
  incompatible clocks, out-of-range clocks, and reflected-clock converse.
- 109,606 converse checks in total, counting clocked and unclocked cases.
- 512 additional deterministic seeded order-five/six instances with two pins.
- All 760 semicomplete digraphs of orders 1..4 under all 17,665 insertion orders.
- Explicit neither-direction, both-directions, required-cycle, unnecessary
  semicompleteness, and pin-obstruction controls.
- Direct authentication of the four native local paths, both allowed states,
  both successor branches, and compatible shared-pair flags.

The producer's normal and optimized reruns each reproduce **145,809 checks**:
46,656 noncontradictory four-owner flag instances, 93,312 single-pin instances,
and its additional converse/native controls. These are separate from the
independent checks above.

Six frozen input files are recorded in `INPUT_PINS.json` and copied into
`inputs/`. The four producer packet artifacts and the inherited source JSON
and gluing note were unchanged at the end of this audit. The prior wiring
audit manifest was rechecked without modification.

## Reproduction

From this directory:

    python check_independent.py
    python -O check_independent.py
    python verify_manifest.py
    python -O verify_manifest.py

The default checker writes its receipt only to standard output. Supply an
explicit `--output` to save a new receipt, noting that replacing a frozen timed
receipt changes its hash. The manifest verifier checks input hashes, artifact
hashes, and semantic equality of both normal/optimized result pairs.
