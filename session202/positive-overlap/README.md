# A free edge-order sector and exact positive overlap arithmetic

**PROVED, independently audited in the stated physical scope.**

The explicit tournament E(u) has exactly u Hamilton paths placing its marked
arc's tail before its head, for every integer u>=1. Its opposite sector is
correlated with u by the inherited function g. For u>=2 these are strong
one-Hamilton-cycle graphs of order 2+bit_length(u-1).

Two literal shared-edge constructions implement the sector changes u->2u
and u->2u-1. Their total counts are exactly 2u+g(u) and
3u+g(u)-(u mod2)-2. The triangle construction is genuinely nonproduct: its
ratio to the product of constituent counts approaches one half along an
explicit family.

This is not yet total-count surjectivity. The family misses13. The proposed
physical constructors with outputs6u-1 and6u+1 remain OPEN. The proof shows
that finding both for every u>=2 would imply positive existence for every
odd count other than7 and21, using ordered triangle products and the actual
E(39) witness with63 paths. No eventual threshold is established here, and
the separate global exclusions of7 and21 are not reproved.

Read EDGE_SECTOR.md for the construction, counting proofs and open target.
The independent report is audit/AUDIT.md. It also retains the exact three-line
scope clarification made before approval. The initial and revised audited
proofs are preserved in originals/; PROMOTION.json records the status-only
promotion. Frozen dependencies identify inheritance and comparison sources.
No exploratory graft probe is part of this reviewed release.

Reproduce with the Python standard library from this directory:

    python3 verify_sector.py > /tmp/free-sector-author-normal.json
    python3 -O verify_sector.py > /tmp/free-sector-author-optimized.json
    cmp /tmp/free-sector-author-normal.json SECTOR_NORMAL.json
    cmp /tmp/free-sector-author-normal.json /tmp/free-sector-author-optimized.json
    python3 audit/check.py > /tmp/free-sector-audit-normal.json
    python3 -O audit/check.py > /tmp/free-sector-audit-optimized.json
    cmp /tmp/free-sector-audit-normal.json audit/NORMAL.json
    cmp /tmp/free-sector-audit-normal.json /tmp/free-sector-audit-optimized.json
    sha256sum -c RELEASE_FILES.sha256

The author and independent checks each cover129 direct indices and256 literal
extensions through order11. The independent checker also verifies deletion
counts, endpoint exceptions, the u=1 boundary and the nonproduct arithmetic
family through k=100. Its conditional factorization test is arithmetic under
the stated OPEN hypothesis, not a claim that missing extractors were built.
