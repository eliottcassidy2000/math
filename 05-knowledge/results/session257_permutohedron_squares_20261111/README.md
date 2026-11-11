# Square-sum ports, permutohedra and filtered arithmetic

**Session257, 11 November 2026. COMPLETE research packet.**

This session recovers the user's original square-sum coloring, valuation cycle
and Monster prime shells, then develops actual geometric and operation maps.
Inherited results are inputs, not targets of a broad re-audit. Chan–Pak remains
**OPEN**. General square-sum existence from order25 and cycles from32 are
already cited in the repository; this session studies richer structure.

The [inheritance atlas](history/README.md) resolves the original references:

- The three edge colors are birth-role classes in Q25. They are neither its
  three old components nor its three first square-sum labels.
- The cycle4→2→5→3→7→13=4 comes from valuations of |M|/|B| and an explicit
  4↔13 identification. It is not a proved arithmetic recurrence.
- The prime sets are B minus Fi22 and M minus B, with inherited modular-level
  and Sylow-fixed-point meanings.
- The recent tournament permutohedron is a coherently oriented braid graph:
  Hamiltonian paths are sinks. Earlier transitive-tournament geometry alone
  does not capture its full role.

[Concept board](BOARD.md). The completed construction notes and exact experiments are below. Starting main:
`6512e16e0729401eab11508348d491992bff3cf3`.

## Completed mechanisms (checkpoint 2)

- [Singleton braid chambers](chambers/PROOF.md): square-sum Hamiltonian words
  are exactly singleton fibers of the graphical-zonotope quotient.
- [Weighted path-cover atlas](permutohedron/PROOF.md): exact block geometry,
  cut refinement and surgery; Q25's elementary chart is completely frozen
  even though two-edge surgery connects all ten paths.
- [Two-seam Gram and primitive decoder](gram/PROOF.md): positive sharp probe
  bank, four-seam negative boundary, and a triple-preserving tournament map.
- [Fixed faces and valuation normalization](filtration/PROOF.md): prime-shell
  port regimes and an exact360-edge permutohedron contraction.

[Read the synthesis](SYNTHESIS.md) and [next questions](OPEN_QUESTIONS.md).
Run `python3 -B reproduce_all.py` and
`python3 -B reproduce_all.py --optimized` for the four exact engines.
Source provenance and publication hashes are recorded in
[SOURCE_PINS.json](SOURCE_PINS.json) and [MANIFEST.json](MANIFEST.json).
