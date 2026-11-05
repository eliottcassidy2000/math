# Four-block physical coverage probe

Read PROOF.md for the new bounded-root-change construction and its exact
remaining coverage lemma. The proof does not close the {7,21} spectrum.

Run serially:

    python3 verify.py > VALIDATION_NORMAL.json
    python3 -O verify.py > VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json

The default exact cap is1,000,001 and can use several GiB of memory.
`--cap 50001` is a shorter complete probe. The small independent positive-
coordinate closure, symbolic block controls and literal graph controls run
in either case.

GAPS_THREE.json lists all73 gaps in the tested at-most-three-J class. The
largest is27,781. This stabilized finite gap set motivates the open all-order
four-block lemma; stabilization is not a proof of it.
