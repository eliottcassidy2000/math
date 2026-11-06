# Physical section of the 6/7 signed response

Proof candidate for an intrinsic physical inverse and target-preserving
rebase on the actual slice E=2(U−D), d=3(U−D). The signed response6U−7D
on that slice is exactly an unsigned source path count. An explicit
five-letter source extension gives the section; the converse is proved
by forcing physical predecessor types rather than filling lost coordinates.

The section has two disjoint branches and an inverse cover-height reduction
by at least a factor six. It does not prove that every signed fiber meets
the slice. Universal odd-target coverage remains open.

Reproduce serially:

    python3 verify.py > VALIDATION-normal.json
    python3 -O verify.py > VALIDATION-optimized.json
    cmp VALIDATION-normal.json VALIDATION-optimized.json

The default complete source bank has U<=20,000. Literal graph controls count
all four response coordinates, root-start paths, and Hamilton cycles.
