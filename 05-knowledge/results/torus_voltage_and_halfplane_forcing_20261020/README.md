# Geometric forcing and winding audit

- [Proof packet](geometric_forcing_bridges.md): G1 half-plane equality; G2 triangular lattice → K7 torus; G3 voltage-sensitive traces; G4 carousel/Paley oriented lifts.
- [Primary-source scope audit](source_audit.md): current ProperColoring definition is standard; no Lean build or complete analytic proof audit performed.
- [Source manifest](source_manifest.json): pinned OpenAI source commit, blob hashes, and inspected declaration locations.
- [Exact checker](check_geometric_forcing.py) and [receipt](exact_checks.json): 160 checks, normal/optimized runs identical.

- [Independent audit](audit/AUDIT.md): 272 exact checks, matched primary-source blob hashes, and explicit proof-status boundaries.

- [Classic lattice nonduplication](CLASSIC_LATTICE_NONDUPLICATION.md): a self-contained finite-boundary proof for arbitrary pieces, with exact square and D6-invariant hexagonal bounds.

The packet proves scoped bridge/obstruction statements and exact controls. It does not prove a new Euclidean six-colour lower bound or claim novelty.
