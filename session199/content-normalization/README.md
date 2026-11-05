# Primitive-normalization hostile

An actual seven-vertex source with profile `(39,9,39,9)` has the nonphysical
primitive quotient `(13,3,13,3)`. The same failure persists after any number
of sink extensions. Complete enumeration certifies that U=39 is the first
cover-height failure, with a unique failing profile there.

This rules out bare gcd division as a physical height-reduction procedure.
It does not give a new scalar spectrum gap or refute arbitrary replacement.

Reproduce:

    python3 verify.py > replay-normal.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-normal.json replay-optimized.json
    cmp replay-normal.json VALIDATION_NORMAL.json
    python3 audit/independent.py
    python3 -O audit/independent.py

The independent audit counts the witness by literal path permutations,
separately from the author's subset dynamic program and response recurrence.
Frozen proof dependencies and their receipt hashes are included.
