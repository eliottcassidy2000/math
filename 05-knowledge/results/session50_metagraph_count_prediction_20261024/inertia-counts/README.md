# Native observer holonomy and exact class certificates

This packet proves a finite completeness certificate for owner-labelled
reobservations, rather than assuming an all-size local-connectivity theorem.
Read `THEOREM.md` for all quantifiers and inherited boundaries.

Main results:

- A normalized reobservation component with `b` words and closed-walk
  owner group `G` satisfies `b|G|<=H`. Equality certifies both the entire
  isomorphism fibre and the entire source automorphism group.
- Complete unbounded-span single-vertex-relocation graphs saturate on every
  component through order seven, recovering `1,1,2,4,12,56,456` ordinary
  tournament classes without an isomorphism canonicalizer in generation.
- Partial graphs give the rigorous sandwich
  `sum_B b_B|G_B|/H_B <= A000568(n) <= number of components`.
- Fixed contiguous-window observer moves fail all-size completeness.
  The source `C3[TT_k]` even has a frozen observer with nontrivial source
  isotropy, so its local holonomy misses an actual symmetry correction.
- A strictly increasing full-reobservation rescue always finishes when
  local moves fail. It need not be efficient or confined to small windows.
- The native cycle-type loop-weighted count agrees coefficientwise with
  inherited Burnside. That identity is not claimed as a new formula.

Reproduce with Python 3.10+ and only its standard library:

```sh
python3 check.py --output fresh-normal.json > fresh-normal.out
python3 -O check.py --output fresh-optimized.json > fresh-optimized.out
cmp fresh-normal.json normal.json
cmp fresh-optimized.json optimized.json
cmp fresh-normal.out normal.out
cmp fresh-optimized.out optimized.out
```

For the actual manifest-only, empty-directory reproduction:

```sh
python3 verify_packet.py
```

The verifier hashes all manifest entries, copies only those entries into a
fresh temporary directory, and reruns both modes there. Normal and optimized
JSON and stdout must each match the frozen outputs byte-for-byte. There are
no parent-directory reads, third-party dependencies, network calls, or Git
operations. Historical source Markdown is archived verbatim as inert `.txt`
with source paths, commit, sizes, and byte hashes in `sources/SOURCE_PINS.json`.

The generation universe is every fixed-path word through order seven.
Independent controls enumerate actual Hamiltonian orders and source
automorphisms of every discovered representative; every labelled tournament
through order five additionally gets a full-permutation canonicalization
and cycle-index fixed-point check. Explicit window hostiles go through
order eight, and four cyclic-chain inflation controls reach order twelve.
The latter are structured examples, NOT a complete order-twelve census.

All-size connectivity of the unrestricted-span relocation graph and a
compressed predictive cross-order recurrence remain open in this packet.
