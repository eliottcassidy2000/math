# Dynamic repair cores

The main result is an all-context gap-two repair on five core vertices and a gap-four window (three consecutive odds) on six core vertices, valid under arbitrary singleton binary prefix/suffix padding. `PROOF.md` provides physical definitions, exact states, all-order proof, small-universe minimality scope, and the reason repeated cores multiply rather than fill the gaps.

- `probe.py`: portable subset path-state and ordinal-composition routines.
- `verify.py`, `CHECK.json`, `CHECK.optimized.json`: separate literal permutation/core and whole-graph checks, no assertions.
- `search.cpp`, `SEARCH.txt`, `SEARCH.stats`: exhaustive all-labeled-tournament census through six vertices, preserving the entire interface except T.
- `search.py`, `UNIVERSAL_WITNESS.json`: early discovery route for the minimal five-vertex pair, not required by verification.

No all-target coverage, arbitrarily long scalar interval, or literature-priority claim. No repository files are imported.
