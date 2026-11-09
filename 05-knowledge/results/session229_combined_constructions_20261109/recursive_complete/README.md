# Recursive-template extension

Current portable deliverable: PROOF.md, verify.py, CHECK.json,
CHECK.optimized.json, ternary.py, TERNARY.json, TERNARY.optimized.json.

- Sharp size bound and complete compiler for arbitrary three-child cyclic substitution.
- Literal three-moment construction separating the former binary language.
- Complete finite images at4095:577 ->620 ->660 after full ternary reuse and an inherited H19 seed.
- General cycle-certified finite-template logarithmic node bound and compiler.
- Structural small-target barrier:19 without the seed,29 after it.

All-odd Chan–Pak completeness remains OPEN. The core substitution/Newton laws and H19 module are inherited, with explicit sources in PROOF.md. Computation is exact but grammar-relative.

Reproduce:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B ternary.py > ternary.replay.json
    cmp ternary.replay.json TERNARY.json

Add -O for optimized replay; expected output is unchanged.

checkpoint1/ is the immutable earlier proof+17-escape snapshot. inherited.py,
explore.py,new.json and progress logs in the producer directory are exploratory
provenance, not dependencies of the portable deliverable.
