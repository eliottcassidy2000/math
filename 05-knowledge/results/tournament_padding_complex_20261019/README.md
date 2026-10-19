# Exact tournament padding complex

Read [the proof](PADDING_COMPLEX.md), [independent proof and computation review](review/INDEPENDENT_REVIEW.md), and [algebra audit](REVIEW_ALGEBRA.md). Deletion summation d and universal-source evaluation h satisfy d²=0 and dh+hd=identity over F₂, giving exact kernel/image equality and explicit primitives. Embedding padding and source-deletion realize the operators. Tournament class counts determine all operator ranks.

This describes the entire function space, with exact graph embedding identities. It does not assert that every primitive is a single tree or classify all good trees.

Run from this directory:

    python3 -B verify_padding_complex.py
    python3 -B -O verify_padding_complex.py
    python3 -B review/independent_controls.py
    python3 -B -O review/independent_controls.py

The producer's 35 checks include tournament classes through order five and direct ranks through d₄. The independent 8,832 checks include every oriented graph through order four and padded paths through order five. Each normal/optimized pair matches exactly.
