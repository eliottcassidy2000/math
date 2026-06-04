---
id: HYP-2188
status: SUPPORTED computational speedup and search diagnostic; not a proof of u(22)=60
source: user-2026-06-03; codex-2026-06-03-S617
tags: [unit-distance, n22, Moser-carrier, beam-search, frontier-gain, counting-speedup, tournament-analysis]
---

# HYP-2188: unit-distance counts speed up by retained Moser unit-shell carriers and frontier-gain extension ledgers

Small unit-distance candidate counting should not recount all pairwise
distances after every child extension when the search is already inside a fixed
unit-shell carrier.

For the rank-4 Moser coordinate carrier extracted in S614, the directed unit
shell has `18` vectors, or `9` antipodal unit directions. If `S` is a connected
cluster and `q` is a frontier point, then the child edge count is

```text
E(S union {q}) = E(S) + gain(q),
gain(q) = |{u in U : q+u in S}|.
```

Thus one frontier-gain table per parent state replaces per-child pairwise
unit-distance recounts. The quotient preserved here is the unit-edge increment
observable on the fixed carrier. It deliberately forgets global planar
embeddability beyond the carrier, totally-unfaithful obstruction labels, and
automorphism-class data; those must be reattached if this is to become a proof
route.

## Evidence

S617 implements this carrier-local count in
`04-computation/unit_distance_moser_beam_speedups_s617.py`. With beam width
`1200` and target `n=22`, the run recovers the dense Moser-carrier lane

```text
n:       2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
best E:  1  3  5  7  9 12 14 18 20 23 27 30 33 37 41 43 46 50 54 57 60
```

At `n=22`, the ledger records `228629` unique children, `242190` frontier
evaluations, best edge count `60`, span `5`, and an edge-count-check speedup of
about `211.4x` relative to naive child-local recounting. The speedup number is
only for the unit-edge count itself; canonicalization and beam ordering remain
real runtime costs.

The same run adds the dense-core extension diagnostic suggested by HYP-2176.
Inside the retained beam, `57`-edge `21`-cores accept only gain-`3` frontier
extensions, while `56`-edge cores reach gain `4`; both lanes top out at `60`,
not `61`.

```text
core_edges=56: gains {1: 162157, 2: 51932, 3: 7673, 4: 1156}, max child E=60
core_edges=57: gains {1: 13848,  2: 4776,  3: 648},  max child E=60
```

This is not an upper bound for planar `u(22)`. It is a beam-local obstruction
ledger: the currently retained Moser-carrier search recovers the known
`60`-edge lower-bound lane and shows exactly why it did not stumble into `61`.

Rebase integration: incoming opus S599w-x/HYP-2187 independently applies the
same state-local frontier-gain principle to LRC via survival bitmasks on
`Z/(2n-1)`. That does not change the unit-distance evidence, but it makes the
abstraction sharper: after choosing the right carrier, update the local
observable by a frontier ledger rather than recomputing the whole child state.

## Tournament Analysis

S617 challenges the default vertex assumption before using Tournament Analysis.
Useful vertex sets include points, frontier candidates, unit directions,
`21`-cores, obstruction filters, and proof routes. The script uses counting
tricks/proof routes as vertices, because the user asked for speedups and the
preserved predicate is "how quickly and faithfully does this route count or
certify dense unit-distance extensions?"

The resulting route tournament is transitive. It ranks:

1. `21-core extension ledger`;
2. `totally-unfaithful obstruction library`;
3. `frontier-gain incremental count`;
4. `Moser unit-shell carrier`;
5. `bitset adjacency popcount window`;
6. `raw F-free graph enumeration`;
7. `triangular-lattice-only beam`;
8. `naive pairwise distance recount`.

The lesson is that the fastest raw counter is not automatically the best proof
route. The extension ledger and obstruction library are higher-value because
they preserve the side information that HYP-2176 identified as decisive.

## Next targets

1. Add automorphism canonicalization for the Moser unit shell so equivalent
   child states do not compete in the beam.
2. Replace set lookups with bounded-window bitsets and popcounts after choosing
   a finite Moser box.
3. Turn the gain ledger into a `21`-core extension solver: enumerate candidate
   gain-`5` and gain-`4` extensions against exact or known dense cores.
4. Compile the totally-unfaithful subgraphs from Alexeev-Tikhonov into a
   reusable obstruction library and run it before expensive embedding checks.
5. Compare the Moser-carrier lower-bound lane with CM/class-field asymptotic
   carrier data from HYP-2184 and the incoming HYP-2186 equidecomposability
   packet: both are retained unit-shell/angle carriers, but they use very
   different side channels.
6. Port the state-local frontier-gain pattern back and forth with incoming
   S599w-x/HYP-2187: unit-distance uses unit-shell neighbor gains; LRC uses
   survival bitmasks on the `2n-1` shell.

## Artifacts

- `04-computation/unit_distance_moser_beam_speedups_s617.py`
- `05-knowledge/results/unit_distance_moser_beam_speedups_s617.out`
- `07-reflections/unit-distance-moser-beam-speedups-s617.md`
