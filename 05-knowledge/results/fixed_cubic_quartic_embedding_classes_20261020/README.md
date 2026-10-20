# Finite embedding classes and fixed cubic/quartic cylinders

20 October 2026. **Independently reviewed PASS.**

Every fixed polynomial curve algebra with degree semigroup <3,4> has at most
two plane-embedding classes. After normalizing
h=y³+a y and k=y⁴+b y²+c y, the only candidates are
(h,k) and (k,h²+(3c/2)h). The second candidate generates the algebra
exactly when c≠0 and 2a=3b.

The separate `EXACT_EMBEDDING_COUNT_APPENDIX.md` proves the exact count:
two classes when c≠0 and 2a=3b, one otherwise.

The finite-class theorem restores the
fixed-cylinder period obstruction without falsely asserting uniqueness.
The period obstruction handles normalization collisions; the chain rule
handles ramification. Together they exclude entire generating Keller pairs
on every such fixed cylinder.

This does not exclude all moving cubic/quartic Keller sources.

Read `FINITE_EMBEDDING_CLASSES_AND_FIXED_CYLINDERS.md`.
The independent review is `review/AUDIT.md`: 79 core producer checks,
15 appendix checks, and 70 independent controls, with matching normal
and optimized outputs.
Run from this directory:

    python3 verify_finite_embedding_classes.py
    python3 -O verify_finite_embedding_classes.py

Appendix controls:

    python3 verify_exact_embedding_count.py
    python3 -O verify_exact_embedding_count.py
