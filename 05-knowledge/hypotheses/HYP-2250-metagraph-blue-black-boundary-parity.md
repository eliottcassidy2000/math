# HYP-2250: Blue/black metagraph parity is boundary data, not color essence

**Status:** OPEN, with the naive simple-graph claim refuted through `n=6`  
**Source:** codex-2026-06-05-S675  
**Primary artifacts:** `04-computation/metagraph_blue_black_parity_s675.py`,
`05-knowledge/results/metagraph_blue_black_parity_s675.out`,
`07-reflections/metagraph-blue-black-even-odd-s675.md`

## Statement

The merged tournament metagraph's black/blue split should not be read as:

```text
black subgraph = even graph
blue subgraph = odd graph
```

at the simple visible-edge level.  That statement fails under both active repo
blue-black conventions.

The corrected invariant is:

```text
color layer = GF(2) 1-chain
evenness = zero boundary
oddness = nonzero boundary packet
```

where the boundary is the set of metagraph vertices with odd layer-degree.
The color name is not the invariant.  The boundary vector is.

## Evidence

S675 audits two conventions.

### SC-type wiggly convention

This is the convention from `wiggly_metagraph_full_s276.py`:

```text
BLUE  = self-converse <-> self-converse
BLACK = self-converse <-> non-self-converse
GREEN = non-self-converse <-> non-self-converse
```

In the simple visible graph, black is not even:

```text
n=4 BLACK: E=2, simple_boundary=2
n=5 BLACK: E=8, simple_boundary=4
n=6 BLACK: E=45, simple_boundary=10
```

Blue is not an all-odd layer beyond the first tiny cases:

```text
n=5 BLUE: E=12, simple_boundary=2, incident_all_odd=False
n=6 BLUE: E=13, simple_boundary=6, incident_all_odd=False
```

However, the weighted line-count chain has zero boundary for all SC-type
layers through `n=6`.  In this convention that may be forced by labelled
arc-flip multiplicities rather than by the blue/black color itself.

### Explorer complement-line convention

This is the convention from `metagraph-explorer-v2.html`:

```text
BLUE line  = complement-tiling line from a grid-symmetric fixed-path tiling
BLACK line = complement-tiling line from a non-grid-symmetric fixed-path tiling
```

After merging by true complement / transpose, pure-black is not even as a
simple graph:

```text
n=4 PURE_BLACK: E=1, simple_boundary=2
n=5 PURE_BLACK: E=10, simple_boundary=4
n=6 PURE_BLACK: E=147, simple_boundary=18
```

Pure-blue is incident-all-odd through `n=5` but fails at `n=6`:

```text
n=3 PURE_BLUE: incident_all_odd=True
n=4 PURE_BLUE: incident_all_odd=True
n=5 PURE_BLUE: incident_all_odd=True
n=6 PURE_BLUE: incident_all_odd=False
```

The striking surviving pattern is weighted:

```text
n=4 PURE_BLACK: weighted_boundary=0
n=5 PURE_BLACK: weighted_boundary=0
n=6 PURE_BLACK: weighted_boundary=0
```

Thus the user's black-even observation is real in a sharper sense if "lines"
means parallel complement lines counted with multiplicity mod `2`, not just the
simple aggregated edge topology.

## Even/Odd Interpretation

Even numbers are the scalar shadow of closed parity: no exposed endpoint.
Odd numbers are the scalar shadow of an address defect: one boundary bit, and
in a finite graph these defects occur in even-cardinality packets by the
handshaking lemma.

The metagraph behaves the same way.  Merging by complement or transpose forgets
an address.  The forgotten address reappears as a boundary vector on the color
layer.  Addition moves boundary bits by xor; doubling kills boundary bits.
This is why weighted line parity can be even even when the simple quotient
graph is not.

## Relation To Royle-Even Graphs

The repo's "tournaments are equinumerous with even graphs" thread refers to
Royle-even graphs, not degree-even graphs.  S675 prints the warning:

```text
n=4 tournament/Royle-even classes = 4, degree-even graph classes = 3
n=5 tournament/Royle-even classes = 12, degree-even graph classes = 7
n=6 tournament/Royle-even classes = 56, degree-even graph classes = 16
```

So the even graphs that matter for A000568 are probably not the visible
blue/black subgraphs themselves.  They should be sought as a functorial target
that preserves tournament predicates and the color-layer boundary vector.

## S674b Integration

The rebased S674b trienerment addendum is the same warning in a different
model.  Signed LRC gauges preserve the observer predicate while changing the
pair-clock address.  Unit distance preserves more information when pairs are
kept as:

```text
S < 1, U = 1, L > 1
```

before resolving the `U` equality layer into a binary tournament edge.

S675 says the same thing for metagraph blue/black lines:

```text
visible binary edge color is too coarse;
retain the line multiplicity/address layer and read its GF(2) boundary.
```

Thus HYP-2249 and HYP-2250 share a proof-design rule: do not collapse the
third/equality/address coordinate until the boundary defect has been measured.

## Assumption Challenge

This session explicitly rejects the assumption that metagraph vertices or arcs
carry the whole parity predicate.  Candidate vertex sets include tournament
classes, fixed-path tiling addresses, complement-line packets, boundary
vertices, Royle-even graph classes, cycle-space defects, and proof-route
obligations.

The quotient preserves the visible color-layer adjacency.  It destroys the
tiling address or multiplicity source unless weighted line parity is retained.
The challenged assumption is:

```text
the visible black subgraph is the even object
```

The replacement is:

```text
the weighted black line-chain may be a zero-boundary object; the simple
black subgraph is only its quotient shadow.
```

## Next Tests

1. Prove or refute weighted pure-black zero boundary for explorer complement
   lines at `n=7`; a generator using canonical pruning or `gentourng` should
   avoid the full fixed-path labelled expansion.
2. Determine whether all SC-type weighted color layers are zero-boundary for
   formal orbit-stabilizer reasons.
3. Search for a Royle-even graph functor that takes a tournament class plus
   retained line-address data to an even graph while preserving `H`, odd-cycle
   packets, and the S675 boundary vector.
4. Compare the boundary-vector support with HYP-2245 ultrafilter leakage: the
   odd boundary vertices should identify exactly where a side-choice fails to
   descend.
5. For LRC14, treat owner/carry leakage as the same GF(2) boundary phenomenon:
   strict rows should have cancellable boundary; AP/Vstar/2AP should be named
   boundary atoms.
