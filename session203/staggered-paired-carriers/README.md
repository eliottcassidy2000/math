# Paired staggered carriers and a full-state grammar boundary

**PROVED in the stated scopes; F(8,137) and the eventual threshold 23 remain OPEN.**

- Four explicit pairs of rooted tournaments retain simultaneous staggered residue certificates under every common O/Z word, with starts at most 17 and shifts at most 5.
- An explicit unbounded start-two family lies off every possible O/Z/J last-letter face. No finite physical seed bank is complete for all nine-coordinate states, even with start cap 2. This does not rule out completeness for the smaller pair (B,s).
- The separate sibling packet `../carrier-menu-density/` proves that the selected four carrier orbits cannot cover an odd tail.

Read PROOF.md. Exact dependency provenance is in dependencies/PROVENANCE.json. The reviewed original and the status-only promotion diff are retained. The independent audit is in audit/REVIEW.md and preserves all files named by its manifest.

Reproduce the producer:

    g++ -std=c++17 -O3 verify.cpp -o verify
    ./verify > replay.json
    cmp replay.json NORMAL.json

The producer checks 504 literal carrier graphs and 11 interior-state family graphs. The independent audit additionally checks 120 literal four-core readouts, 855 one-vertex transfers, and 1140 direct cover interpretations. No completeness of the discovery bank is required.
