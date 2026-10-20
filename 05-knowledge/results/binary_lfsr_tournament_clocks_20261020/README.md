# Binary LFSR and tournament source clocks

A self-contained research packet, 20 October 2026. Status: elementary proofs plus exact finite controls; historical novelty has not been established.

- [Binary primary clocks](BINARY_LFSR_PRIMARY_CLOCKS.md): complete companion graph formulas, exact dyadic-profile classification, a sharp degree-18 information-loss example, and compatible supernatural clock limits.
- [Tournament source-adjunction odometer](TOURNAMENT_SOURCE_ADJUNCTION_ODOMETER.md): a native-operation bridge from the repository's source-free/Hasse module to a completely classified binary automaton. Its completed dual has every nonzero orbit closure explicitly conjugate to addition on Z₂.
- [Exact checker](verify_lfsr_tournament_clocks.py), [machine-readable report](verification.json), [text report](verification.out).

Reproduction: `python3 verify_lfsr_tournament_clocks.py`.

The generic LFSR cycle ingredients are classical. The strongest new-to-this-repository deliverables are the sharp information-loss boundary and the actual source-adjunction Pascal embedding. Matching finite counts or supernatural numbers does not identify modules, whole phase spaces, coloring geometries, or physical dynamical systems.
