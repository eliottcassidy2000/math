# Physical adjoint and square-root inverse

An exact symmetric unimodular pairing transfers every actual source word
by reversal and exchange of O with Z. The transferred source is physical,
and the number of root changes is preserved. Every represented C1 target n
has a pair presentation with both cover heights at most floor(2 sqrt(n-1)).
The same cutoff gives an exact inverse for balanced counts with a prescribed
root-change budget. Universal odd-target coverage remains open.

Reproduce in this directory:

    python3 verify.py > replay-normal.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-normal.json replay-optimized.json
    cmp replay-normal.json VALIDATION_NORMAL.json
    python3 audit/independent.py
    python3 -O audit/independent.py
    python3 verify.py --target 1000001 --budget 3
    python3 verify.py --target 879001 --budget 3
    python3 verify.py --target 27781 --budget 3

The last target is rejected only in the balanced class with at most three
precursor root changes. Rejection in that class is not a tournament gap.

The independent audit uses separately pinned positive-coordinate matrices,
a complete forward bank, and literal permutation counts on 32 graphs. It
imports no candidate code. Both programs use only the Python standard library.
Frozen dependency texts are included; they supply source grammar completeness,
C1 reconstruction, and the four-block language. The proof and verifier hashes
in audit/INDEPENDENT.json identify the reviewed versions. The manifest records
all release files except itself.
