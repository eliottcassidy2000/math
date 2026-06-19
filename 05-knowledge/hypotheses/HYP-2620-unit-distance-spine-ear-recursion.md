# HYP-2620: Unit-Distance Spines Persist by Endpoint-Compatible Ear Recursion

**Status:** SUPPORTED by S16 computation; OPEN for all true optimal planar
unit-distance configurations.

## Claim

The intrinsic version of the user's point-to-tournament question is still
traceability of the unit-distance graph `G_1(P)`, but S16 adds a stronger
recursive invariant:

> A dense unit-distance witness should have an endpoint-compatible ear: a
> vertex `v` such that `G_1(P\v)` has a Hamiltonian path and `v` is adjacent
> to at least one possible endpoint of such a path.

If endpoint-compatible ears persist, then an all-unit Hamiltonian path is not
a lucky after-the-fact labeling.  It is an inductive certificate section.  A
genuine geometric flop would need to destroy this recursion, not merely make a
fixed lexicographic point-order tournament lose an all-unit directed path.

## Evidence

`04-computation/unit_distance_spine_ear_recursion_codex_s16.py` computes
Hamiltonian endpoint masks for each unit graph and audits every vertex
deletion.  Output is stored at
`05-knowledge/results/unit_distance_spine_ear_recursion_codex_s16.out`.

The run covers:

- the S625 compact triangular carrier through `n=22`;
- the S622/S626 Moser beam through exact rows `n<=14`, exact `n=21`
  with `57` unit edges, and the `n=22` Moser `60`-edge lane.

Main findings:

```text
graph-level unit-spine failures: none in either carrier
fixed lexicographic unit-flip flop: first at n=7 in both carriers
endpoint-ear failures: none in either carrier
nonunit complement HP: true at n=6, false at n=7, true again from n=8
```

The Moser rows are the strongest evidence.  In the width-`1200` beam:

```text
n=14: 33 edges, every vertex endpoint-compatible, 1 exact-edge deletion ear
n=21: 57 edges, every vertex endpoint-compatible, 1 exact-edge deletion ear
n=22: 60 edges, every vertex endpoint-compatible, 1 exact-edge deletion ear
```

Here "exact-edge deletion ear" means deleting that vertex leaves the known
edge target for the previous row (`33 -> 30`, `57 -> 54`, `60 -> 57`).
The same rows have no deletion separator with more than one component in this
beam witness, so the usual traceability obstructions are absent.

The triangular carrier data is useful for the gauge/flop pattern, but after
`n=8` it is a compact beam carrier rather than an exact Harborth optimum.  Its
role here is as a persistence lab, not an exact-optimal certificate.

## Answer to the Flop Question

No graph-level unit Hamiltonian-path flop was observed in the tested carriers.
The first observed flop remains the fixed-order tournament artifact at `n=7`:
the lexicographic unit-flip tournament cannot carry an all-unit directed
Hamiltonian path, while the undirected unit graph has many unit Hamiltonian
paths and every vertex is still endpoint-compatible.

Thus the current best reading is:

```text
geometric unit-spine flop: not observed
fixed tiling-order / lex gauge flop: n=7
nonunit complement exception: n=7 compact hexagon center
```

The "mandatory Hamiltonian path" should therefore be read with its retained
tie-order witness.  If the tie path is chosen from the unit spine, either flip
convention can expose an all-unit mandatory path whenever `G_1(P)` is
traceable.

## What a Genuine Flop Would Need

A true intrinsic flop needs more than a bad coordinate order.  It would likely
require one of the following:

- no endpoint-compatible deletion ears;
- a cut vertex separating three or more path-obligatory branches;
- multiple frontier ears whose endpoint requirements are mutually incompatible;
- a dense core extension where every removable vertex misses all smaller-spine
  endpoints;
- a non-lattice true optimum at or beyond the `n=22` frontier whose extra edge
  comes from a spine-breaking attachment rather than a spine-preserving ear.

This makes the next proof target concrete: classify `21`-core extensions by
whether the degree-`4/5` new vertex is endpoint-compatible with some unit
spine of the core.

## Tournament Analysis

S16 deliberately challenges the assumption that the tournament vertices must
be points.  It considers points, unit edges, nonunit edges, deletion ears,
endpoint masks, direction shells, core-extension obligations, and S/U/L
impurity words.

The proof-lens tournament uses the pairwise observable:

```text
which lens preserves the unit-spine predicate while supporting recursion?
```

The score order is transitive:

```text
unit graph traceability
> endpoint-ear recursion
> spine-order tournament
> S/U/L impurity word
> nonunit complement HP
> lexicographic flip gauge
> raw edge count
```

It has score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, no directed
`3`-cycles, singleton SCCs, and one Hamiltonian path.  This says the present
comparison is a clean proof-obligation hierarchy, not yet a nonlinear
tournament obstruction.

## Next Problems

1. Run the endpoint-ear audit on the five stored exact `n=21`, `57`-edge
   graph6 cores from HYP-2176/HYP-2188, not just the Moser beam
   representative.
2. Add endpoint-ear compatibility to the `n=22` core-extension ledger.
3. Build the S/U/L impurity-word scanner promised by HYP-2249/T746 for exact
   and candidate `n=21/22` point sets.
4. Prove a sufficient condition: high edge density plus no multi-branch
   separator forces at least one endpoint-compatible ear.
5. Search specifically for spine-breaking attachments rather than higher edge
   counts alone.
