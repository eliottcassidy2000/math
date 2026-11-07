# Session 210: coherent permutohedra and Grassmannian path cells

The full [proof](PROOF.md) establishes:

- an exact, order-preserving bijection between tournaments and orientations
  of the permutohedron that agree on parallel edge classes; Hamilton paths
  are exactly sinks;
- a product law for the sinks of every face and its incidence polynomial;
- a sharp boundary: generic convex-linear objectives only give one sink,
  whereas arbitrary acyclic orientations admit every count through n!/2;
- a Gaussian-binomial refinement of the classical one-reversed-arc family,
  with an explicit decorated-path bijection to all finite-field subspaces;
- local feedback minima with quadratic gap from the global optimum.

The missing counts 7 and 21 are inherited facts. Their translation to this
carrier is not a new exclusion proof, and global spectrum completeness is
not established. The finite-field lift applies to the specified one-arc
family, with explicit coordinates and decorations, not arbitrary tournaments.
The classical scalar formula is credited separately from the refinement.

Reproduce from this directory:

    python3 verify.py > /tmp/session210-permutohedron.json
    cmp RESULT.json /tmp/session210-permutohedron.json
    python3 -O verify.py > /tmp/session210-permutohedron-optimized.json
    cmp RESULT.json /tmp/session210-permutohedron-optimized.json

All tests are deterministic and use only Python's standard library.
SOURCES.json pins inherited repository claims. REVIEW.md, when supplied,
records independent mathematical audit rather than an author self-review.
