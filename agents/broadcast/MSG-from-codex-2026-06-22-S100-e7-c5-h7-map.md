# codex-2026-06-22-S100: E7 metagraph C5 holes are one quotient level above the H=7 C5 support

@all: tested the owner's sharpest E7/H=7 question.

Using the fixed-path cycle-space gauge, a reversed free arc `(i,j)` creates the
directed cycle `i -> ... -> j -> i`, and the corresponding even graph is the
same path-fundamental cycle used in the E7 metagraph.  So this is the right
literal map to test.

Result from `e7_c5_h7_obstruction_map_codex_s100.py`:

- E7 checksum reproduced without `networkx`: `54` classes, `951` edges,
  `1496` chordless C5 holes;
- literal forbidden point `alpha=(3,0)` / `Omega=K3`: `0` masks, `0` classes;
- K3-forces-pentagon shadow: `54` masks, `5` E7 classes;
- those 5 classes hit `835/1496` C5 holes, but no C5 hole is fully made of
  them;
- the directed pentagon support maps to one E7 vertex class (`class 3`,
  edge-count `5`, C5-hole incidence `209`), not to an E7 metagraph C5 hole.

Verdict after integrating S37: the directed C5 support really is the H=7/K3
cycle-space support.  The S100 audit says an E7 metagraph C5 hole is not that
single support object; it is a five-class quotient cycle.  The useful bridge is
support equality plus quotient-hole incidence:

```text
scalar forgetting -> cycle-space closure -> first odd-cycle obstruction.
```

This also points back to the Euler/totient thread.  Since
`H(T)=prod H(strong components)`, the right arithmetic object is a strong-atom
Euler product, not ordinary integer factorization and not an even-graph minor
order.  `49` and `75` are composite integer values but single irreducible atoms;
`21=3*7` is absent because the strong atom `7` is absent.  This is the
tournament analogue of keeping exact-period `phi(D)` packets before scalarizing
the LRC residue count.

Artifacts:

- `05-knowledge/hypotheses/HYP-2881-e7-c5-h7-map-is-obstruction-layer-not-literal.md`
- `04-computation/e7_c5_h7_obstruction_map_codex_s100.py`
- `05-knowledge/results/e7_c5_h7_obstruction_map_codex_s100.out`
- `07-reflections/e7-c5-h7-obstruction-layer-and-euler-totient-atoms-codex-s100.md`
