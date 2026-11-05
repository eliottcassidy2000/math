# Field-invariant physical source compiler

Read `PROOF.md` for the mathematical claims and limits. The packet contains
explicit quadratic equations encoding actual source words, logarithmic word
bounds, characteristic-zero integrality, and complete C1/C2 target recognition.
It does not prove universal odd coverage or give a single fixed Diophantine
formula for all unbounded targets.

Reproduce the finite checks with Python 3 and its standard library:

    python3 verify.py > replay-normal.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-normal.json replay-optimized.json

The supplied receipts are `VALIDATION_NORMAL.json` and
`VALIDATION_OPTIMIZED.json`. The checker uses explicit exceptions, not assert
statements. It independently counts literal Hamilton paths and cycles for
small source-pair constructions, checks all source words through length seven,
and checks exact length-capped inverses against complete height-capped banks.
The field selector tests include characteristics two and three, where omitting
orthogonality would be unsound.

The positive example decodes target 101 from physical source words OO and
ZJOJ and directly verifies the resulting nine-vertex C1 tournament.
Frozen mathematical inputs are under `dependencies/`. Their status and scope
must remain distinct from the new claims; finite checks alone do not prove
the inherited all-order classifications.
