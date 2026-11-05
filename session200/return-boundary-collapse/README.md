# Boundary returns collapse after two blocks

Independently audited all-order obstruction to a proposed joint-fiber
amplification: any two nonempty words returning the physical triangle boundary
to its own scalar ray force the synchronized auxiliary boundary onto that ray.
Later choices change only the common scale. Every resulting Mersenne target
is composite, so this restricted scheme cannot solve prime targets.

One return can escape: the nine-vertex C1 construction with139 Hamilton paths
is checked literally. This is not an obstruction to the unrestricted source
inverse or to universal odd-target coverage, which remains open.

Run serially:

    python3 verify_return_collapse.py > check-normal.json
    python3 -O verify_return_collapse.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json RETURN-VALIDATION-normal.json

Independent proof review, separate matrix/permutation checker, exact replays,
and status-only promotion verification are in audit/. Frozen receipts retain
their pre-promotion scope wording. Dependencies are preserved in dependencies/.
