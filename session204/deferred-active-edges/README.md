# Exact registers with deferred active edges

This packet supplies an exact finite nonnegative path-forest state representation.
At a boundary, edges between active ports have not yet been included in the
counted digraph. They are processed when a port is forgotten or at completion.
This allows different intermediate completed readouts from one partial state.

It does **not** enlarge the final tournament class of ordinary register
constructions: each completed word is realized by choosing its final active-edge
orientations in advance. It proves neither arithmetic coverage nor an eventual
odd Hamilton-path threshold.

Boundary width W and peak introduction width W+1 are distinct. For ordinary H,
p=1 bulk-start component suffices; rooted L and B additionally require p=2 and
retaining the designated root through completion.

## Exact standalone reproduction

    python3 verify.py > /tmp/deferred-normal.json
    python3 -O verify.py > /tmp/deferred-optimized.json
    cmp NORMAL.json /tmp/deferred-normal.json
    cmp OPTIMIZED.json /tmp/deferred-optimized.json

Only Python's standard library is needed. Every sparse transfer column is
regenerated, with ordinary arbitrary-precision Python integers. No cached matrix,
seed bank, NumPy installation, external sibling file, or optimization-sensitive
assertion is required.

The verifier checks 65,712 sparse columns, 1,464 literal transitions, 292 terminal
completions, and 48 larger permuted transitions. The proof establishes the
all-order transfer semantics and final-class equivalence; the finite controls
check the implementation. Earlier exploratory digit searches are not included
or used as dependencies.
