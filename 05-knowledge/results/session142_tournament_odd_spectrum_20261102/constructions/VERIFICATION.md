# Verification

Commands (all passed):

    python3 -B verify_binary.py
    python3 -B -O verify_binary.py

For the relocated modes, only `verify_binary.py` was copied into a fresh directory, and both commands were run from `/`. All four receipts have SHA-256

    f20e98dbafd4c04a419bcd064e722fad73f1745425afa5f2b819d7fef7bb9731

The checker solves all lemma sign cases algebraically, then separately controls 199 eligible formal seeds in coordinate box0..6 and 1,352 valid binary-history edges through depth8. These finite controls are not the all-length proof. Literal subset-DP counts independently verify two actual H5 source graphs, all listed S/R continuations, and the named H-neutral S1/reset gadgets. There are no asserts, floating-point calculations, external packages, graph-atlas assumptions, or arbitrary numerical source promotion.

Independent auditor `/root/audit_session142_gap_proofs` confirmed the graph laws, the algebraic trap, the strengthened at-most-one-output-per-depth theorem, and the tail formulas. Its audit artifact is separate from this producer packet. The unrestricted {7,21}-only gap conjecture is still open.
