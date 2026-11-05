# A dyadic residual replacement

This is a restricted constructive continuation of the balanced-source normal
form. It characterizes an entire physical source fiber and gives a replacement
that preserves vertex count, Hamilton-path count and the C1 condition, while
changing a graph with no balanced root into one with a balanced root.

The resulting scalar values already have the simple sink-tail realization on
the replacement side. The new content is the exact physical entry condition
and the same-order replacement of the residual graph, not new scalar attainment
or a proof of cofinite coverage.

Run from this directory:

    python3 verify.py > replay-normal.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-normal.json replay-optimized.json
    cmp replay-normal.json VALIDATION_NORMAL.json
    sha256sum -c MANIFEST.sha256

The verifier uses only Python's standard library and local frozen files. It
checks the complete capped dyadic profile fiber, canonical small graph inverses,
literal Hamilton counts and root conditions for fifteen graph pairs, and the
symbolic balanced-input identities. Independent review is delivered alongside
this packet. The proof is all-order, but explicitly restricted to its stated
fiber. Universal odd coverage remains open.
