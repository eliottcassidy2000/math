# Current integration: path reconstruction and operation-specific information

2026-10-26. **PROVED results are attributed to their linked proofs;
FINITE-EXACT controls retain their stated universes.** This is the current
integration of the [growing-sieve/Paley packet](../growing_sieve_path_20261026/README.md)
with the concurrent guard, attachment and fresh re-marking work. The earlier
packet's mathematical bytes and historical source bindings remain frozen.

## Concurrent reconstruction proofs and a sharper decoder

Three continuations converge on the same theorem: the union U(T) of all
Hamilton-path arcs determines every finite labelled tournament T. They use
the same CITED Busch separator theorem, not three unrelated classification
principles:

- [Growing-sieve P1–P4](../growing_sieve_path_20261026/paths/RESULT.md) proves
  vertex-deletion reachability by retaining cycle edges via prepended and
  appended paths, then uses the separator to recover missing orientations.
- [Session70](../session70_attachment_guard_alignment_20261026/session70-guards/REPORT.md)
  first obtains two visible length-two detours through nonseparating middle
  vertices, and derives the vertex-deletion reachability property from them.
- [Guard reconstruction G1–G2](../guard_reconstruction_20261026/README.md)
  uses those middle vertices directly, giving an O(n³) decoder from U.

For the last decoder, first restore cross-component directions from the
ordered strong components of U and retain every existing arc. For a
**missing pair x,y inside one strong component**, set

    c_U(x,y) = #{w in that component : x→w→y in U}.

Exactly one of c_U(x,y), c_U(y,x) is at least two; orient in that direction.
The classical separator puts at least two middle components between the
tail and head. A vertex in either is nonseparating, making its incident
arcs visible. Every reverse path has to pass through the single separator,
so the reverse two-step count is at most one. This proves the rule on the
promised support image, using the original support throughout.

The restrictions are essential. In the transitive triple the missing long
arc has only one two-step route; component order handles it. In Paley7 all
arcs are visible, and an actual arc can have forward count1 and reverse
count2; applying the rule to existing arcs would reverse a correct arc.
This supplies an exact connection to the Paley branch cover: its objects
lie at the no-missing-arcs boundary of the same reconstruction algorithm.
O(n³) measures decoding from supplied U, not the cost of acquiring the full
Hamilton language or an oracle-query bound.

The shared all-order conclusion is the exact `2^binom(n−1,2)` count of
valid fixed-path states under full named legality guards. It is recorded
once as a closure, with independent proof routes. The full support and its
isomorphism type still preserve H and Aut; this transports the existing
7/21 exclusions but does not enumerate the support image.

## What the incoming fresh-re-marking audit changes

The corrected [HYP-9416](../../hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md)
now proves that the **whole sandwich tournament** e*h*e determines every
strong joint two-port state for n≥3. Its added degree-one entrance and exit
are structurally recognizable. This part is no longer an open minimality
question. The scalar H(e*h*e) is a different observation and does compress.
The [independent audit](../../../04-computation/experiments/joint_port_bisimulation_20261026/audit_hyp9416/REPORT.md)
also distinguishes fixed-port successor multisets from normalized
distributions and from full named guards. MISTAKE-676 records the repairs.

The remaining all-order problems include separation by fixed-port dynamics,
port-forgotten atom-response classification and the asymptotic compression
ratio. The older packet's references to open fresh questions must be read
with this repaired partition of proved and open clauses; not every clause
of HYP-9416 is open. Its finite counts and degree-reconstruction proof are
inherited from that audit, not independently recomputed here.

## A stronger unifying account

The [method note](../growing_sieve_path_20261026/METHOD.md) asks whether a
target is constant on observation fibres, then tests the allowed operations.
The incoming [Session70 synthesis](../session70_attachment_guard_alignment_20261026/session70-root/SYNTHESIS.md)
adds useful contrasting outcomes:

| Specified experiment | Information retained | Scope |
|---|---|---|
|Full named tournament guards|Entire labelled presentation|All orders; reconstructible even when arcs are invisible|
|Unrestricted fresh connections between r disjoint nontrivial paths|Pairing of their 2r distinct named endpoints|`(2r−1)!!` states; fixed arithmetic hosts need not realize every distinguishing context|
|Specified graceful edge-overwrite validity guards|Color-to-component permutation|k! states in that family; realization signs can be forgotten|
|Conic branch transport|Pole plus determinant bit|16 oriented Paley objects; no recovery of arbitrary arithmetic quality|
|Twin primality on a bounded interval|Prime tests through a height-dependent cutoff, with self-hit exceptions|Exact finite recurrence; no infinitude conclusion|

The first and last rows explain why a universal fixed-size quotient is the
wrong target for some protocols; the graceful row shows that exact forgetting
is possible for others. The endpoint-pairing row gives the precise reason
ports must remain joint. Session70's Boolean moment hierarchy likewise
distinguishes k bits for one alignment from 2^k linear coordinates for an
arbitrary law on alignments. Those are different targets and costs.

The cohesive theory is therefore operation-relative recoverability, with
three practical repairs: structural inference, a composable finite extension,
or a growing certified horizon. Repeated small numbers suggest experiments;
the successful connections are the maps, context distinctions and hostile
examples above. No reduction among twin-prime infinitude, abc, RH, P versus
NP or black-core connectivity is claimed.

## Version and verification boundary

The source record binds the incoming guard/Session70 material at
`c2039c32f9` and the fresh audit at `9f2a9722ef`. A separate finite check tests
the sharper direct decoder on every labelled tournament through order6,
using the original packet's independently implemented subset-DP support.
Its normal and optimized receipts are equal; it is a new test of the sharper
rule, not another independent implementation of Hamilton support.

Reproduce the merge control from the repository root with
`python3 -I -B 05-knowledge/results/sieve_path_merge_20261026/verify.py`.
The [review](REVIEW.md) audits this synthesis and the preserved decoder
control; the manifest binds both source versions and all local files.

The old packet's default verifier continues to replay its historical source
and correction records. Its `--live` option intentionally rejects maintained
files that have since advanced. The merge preserves the new reservations,
audit corrections and all three reconstruction references while retaining
the conditional all-k Mersenne-fibre scope. No RESERVED theorem is promoted.
