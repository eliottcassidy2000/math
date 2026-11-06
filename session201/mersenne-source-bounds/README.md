# Mersenne targets require growing source coordinates

**PROVED, independently audited in this scope.** If a representation of
n=2^k-1 by one of the twelve actual-source templates has both g coordinates
at most G, then n<=18630 G^8. No sufficiently large Mersenne target can be
represented while keeping both source coordinates bounded by a fixed G.

This does **not** show that any target is missing when G may grow. The uniform
twelve-template occupancy lemma and the global spectrum problem remain OPEN.

- [Proof](MERSENNE.md)
- [Independent audit](audit/AUDIT.md)
- [Status-only promotion](PROMOTION.diff)
- [Frozen arithmetic menu](dependencies/TWELVE.json)
- [Reviewed physical compiler and product bound](dependencies/CORE_PROOF.md)

Reproduce from this directory:

    python3 verify_mersenne.py
    python3 -O verify_mersenne.py
    python3 audit/check.py
    python3 -O audit/check.py

Author controls check 189,165 form/source pairs. The independent controls
check 13,536,024 pairs and deliberately test zero-constant non-Mersenne cases.
The finite controls support the proof; they are not a coverage census.
G=1 already represents 255, and G=3 represents 65535, so the result must not
be read as excluding all Mersenne values at small G.

The original audited source hashes are in audit/AUDITED-HASHES.json.
MERSENNE.md differs only by its reviewed status. The independent checker's
input path is changed only to the bundled, byte-identical MENUS.json.
Unaudited exploration, inverse experiments, and selected-target scans are
excluded from this release.
