# Unequal covered cyclic modules

PROOF.md proves a three-query separation theorem under six cross-profile
Turán inequalities, with no equal-profile assumption. It covers every
cyclic blowup with three blocks of orders 2 or 3. It also derives the exact
unconditional unequal-profile three-query formula and gives a sign boundary
showing why three separate log-concavity hypotheses cannot replace the
cross conditions in that proof.

The mixed seven-vertex example retains the full named polynomial contrast,
including all cubic and quartic coefficients, and all 64 block-uniform
queries. Each polynomial coefficient is independently checked through
literal two-port counts on the full coefficient-determining query set.

Reproduce with standard Python, no packages:

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json

No universal log-concavity, full-response injectivity, or unrestricted
module-sequence claim is made.
