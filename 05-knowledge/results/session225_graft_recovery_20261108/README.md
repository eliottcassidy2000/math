# Session225: reconstructed full-rank rooted graft

**PROVED scoped construction. Chan-Pak completeness remains OPEN.**

This recovery packet gives a self-contained proof of an exact four-coordinate
rooted tournament operation, an explicit invertible four-vertex gadget,
and an infinite strong family with Hamilton-path counts

    1, 37, 1113, 31349, 886465, 25060261, ...

The initial value 1 belongs to the two-vertex seed. Every later graph is strong;
the jth iterate has 4j+2 vertices.

- [Complete proof and definitions](PROOF.md).
- [Portable fresh verifier](verify.py), [normal receipt](CHECK.json),
  [optimized receipt](CHECK.optimized.json).
- [Independent review and verification scope](AUDIT.md), separately written
  [independent enumeration](independent.py), and [its receipt](INDEPENDENT.txt).

The gadget's determinant is -8. Because physical rooted states span the
four-coordinate space, no one fixed universal O/Z/J word has the same state
action. This says nothing about minimum gadget size, novelty of individual
counts, input-dependent words, unrestricted coverage, or literature priority.

## Recovery provenance

The original Session225 draft files were unavailable during recovery.
This is a newly reconstructed proof and newly generated verification packet,
not a byte-for-byte restoration of those files. Retained mathematical statements
and exact gadget data were rechecked; the accompanying 5625-pair receipt was
generated afresh without the original source bank or audit files.

The four-coordinate framework and named root operations are inherited research
objects. This packet records a recovered construction and its rigorous scope,
not a priority claim for the surrounding framework. Later inverse refinements,
other Session225 lanes, and their missing artifacts are not included.

## Replay

    python3 verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

The verifier uses only the Python standard library, no network or external data,
and explicit checks that remain active under optimization.

Independent small enumeration:

    python3 independent.py > independent.replay.txt
    cmp independent.replay.txt INDEPENDENT.txt
