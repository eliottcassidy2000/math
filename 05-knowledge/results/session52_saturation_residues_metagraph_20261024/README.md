# Session52: saturation patterns, residue filters, and metagraph carriers

24 October 2026. First checkpoint: a precise reconstruction of the supplied
run word, its boundary, and a general graph-birth theorem. Other mathematical
lanes are being independently audited for subsequent checkpoints.

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
