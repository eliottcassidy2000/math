# Fixed Hamiltonian frames, smaller graphs, and a semi-inverse for H + 2E

**Chan–Pak remains open.** This packet provides an exact ordered-graph
translation and a constructive inverse on a particular graph geometry,
with an unbounded two-parameter family of adjacent physical counts.

- **[Ordered graph chart](FRAMED_SIMPLE_GRAPH.md).** Fixing the
  Hamiltonian spine leaves one free tile for each unordered pair
  of gaps, giving a simple graph on one fewer vertex. The chart
  carries an order. Reflection has the expected converse meaning,
  but arbitrary graph relabeling is not a metagraph identification.
  The outlet edit is a last-labelled leaf when the root isn't the
  final vertex. Undirected leaf paths count full paths at that root;
  the tournament reads even prefixes and weighted complements.
  A different order-averaged relationship uses fence-poset weights
  for every mixed adjacency word, not just the two ordinary
  Hamilton path counts of a graph and its complement.
- **[Port spectra and controls](ONE_EDGE_PORTS.md).** Even one
  backward tile gives exact mixed powers of 2 and 3 as the root
  moves. Two labelings of the same P4 both have one ordinary
  undirected Hamilton path but give strong rooted totals 37 and 31.
- **[Complement-star inverse](COMPLEMENT_STAR_INVERSE.md).** The
  positive undirected Hamilton count is d(d+1)(m−2)!/2 on the
  smaller graph. Positions of the missing star edges directly
  encode E in binary; the actual target is 1+4E minus an exact
  pair-separation energy. A bounded inverse decodes a target
  and a chosen graph shape. At a final root, the re-framed outlet
  instead adds one isolated backward gap before that root.
- **[Two-scale physical adjacent pairs](TWO_PARAMETER_BAND.md).**
  The two high words 8A−5 and 8A−3 give the same leading coefficient
  with opposite low-bit adjustments. At A=2^(t+2), t>=0, k>=4,
  their exact counts are
  (8A−5)(8A−3) 2^(k−3) + {3−A,5−A},
  and both have the same entire unlabelled graph on the smaller
  gap set. They are strong. The first pair is 1565 and 1567,
  and their larger family has an exact power-of-two target decoder.
  They do not preserve the full continuation state, and the
  intervals between their centres need not be covered.
- **[Neményi's semi-inverse analogy](SEMI_INVERSE.md).** Choose
  part of the geometry, then distinguish equilibrium from
  compatibility; equilibrium alone is not sufficient. The Cartesian
  elasticity control exhibits the passage from free equilibrium
  functions to five compatible constants, under its assumptions.
  In the directed setting, arc-flow balance admits subtours.
  A subset/last-state lift tests compatibility; the actual outlet
  readout needs the joint split count, not the endpoint marginals.
  The continuum theorem supplies the design principle, not a
  five-parameter bound for arbitrary graphs.

The next constructive problem is to find a second independently
selectable readout or additional graph geometries with complementary
absolute bands. The current complement-star inverse has every
fixed odd-compatible residue and still provably omits a positive
fraction of odds, so more local tests alone will not suffice.

## Replay

Run `replay_frame.py`, `replay_ports.py`, and
`replay_inverse.py` with Python 3. The companion receipts
include literal subset-DP controls for both directed tournaments
and the smaller undirected charts.
