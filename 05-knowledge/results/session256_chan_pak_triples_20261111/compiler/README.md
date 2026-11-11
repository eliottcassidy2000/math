# Session256: scalar-fibre compiler

- `PROOF.md`: positive nine-coordinate graph-address machine; terminating
  all-height Pareto optimizer; synchronized primitive-triple driver;
  rational word series; positive-polynomial substitution extension.
- `engine.py`: standalone implementation with graph reconstruction and
  subset-DP controls. Run `python3 -B engine.py`.
- `results.json`: matrices, statistics and explicit minimizing graphs.

New all-height result is exact optimization in fixed modular/regular
fibres, not a universal height bound or Chan–Pak proof. Normal and `-O`
runs produced identical result files; assertions are finite normal-mode
implementation controls, not the termination proof.
- `morphic_optimizer.py` / `morphic_results.json`: actual positive32-coordinate
  substitution compiler, all-height modular minima, regular B-occurrence
  constraint, and exact membership/inverse examples. Run
  `python3 -B morphic_optimizer.py`.

The morphic image misses13 and15 at every depth; their modulo16 fibre minima
are3181 and1615. This narrows the usefulness of the restricted Berggren lift
without implying any obstruction in the five-core grammar.
- `morphic_mod16_certificate.json`: complete633-vector semigroup closure
  certificate (about140KB), with every reversed witness word.
- `verify_morphic_certificate.py`: lightweight standalone closure/replay
  check, without optimizer execution. Checks633 words and1899 transitions.
- `morphic_certificate_verified.json`: exact minima, chronological words,
  triples and graph rows produced by the verifier.

The all-depth morphic exclusions are FINITE-EXACT COMPLETE SEMIGROUP
CERTIFICATES, supported by full positive domination closure, not depth caps.
