# Physical counts inside congruence intervals

This Session 199 packet has independent proof and code review in audit/.
The main new result forces compatible joint one-apex residues inside every
aligned binary interval of an explicit fixed length. For each fixed modulus,
it gives a physical balanced C1 count below any sufficiently large odd target,
congruent to that target, with a square-root error bound and logarithmic graph
order. This is not exact coverage: the constructed family has positive-density
gaps, as explicitly recorded in PROOF.md.

Reproduce from this directory with Python 3, without third-party packages:

    python3 verify.py > /tmp/congruence-primary-normal.json
    python3 -O verify.py > /tmp/congruence-primary-optimized.json
    cmp /tmp/congruence-primary-normal.json /tmp/congruence-primary-optimized.json
    python3 audit/independent.py > /tmp/congruence-audit-normal.json
    python3 -O audit/independent.py > /tmp/congruence-audit-optimized.json
    cmp /tmp/congruence-audit-normal.json /tmp/congruence-audit-optimized.json
    sha256sum -c MANIFEST.sha256

The primary proof and verifier were frozen before independent review. The
reviewer's separate implementation reconstructs the physical source compiler
and the modular correction without importing primary functions. See
DEPENDENCIES.md for the exact inherited sources and the historical connection.
