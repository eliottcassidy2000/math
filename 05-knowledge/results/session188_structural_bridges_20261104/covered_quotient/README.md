# Unequal transitive cyclic blocks

`PROOF.md` gives an exact named two-port separator for quotient converse of
C3[TT_a,TT_b,TT_c], with a,b,c>=2 not all equal. It proves a general cyclic-run
start identity and an elementary Stirling-ratio positivity lemma. Equal sizes
require a different argument; arbitrary child tournaments are not claimed.

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The checker is self-contained and standard-library only. It validates start
recovery by actual port attachments, the cyclic-run identity on arbitrary
small children, Stirling minors, literal separating batteries, and a precise
arbitrary-child ratio counterexample. Explicit checks remain active under -O.

Publication allowlist: PROOF.md, README.md, verify.py, receipt.json,
receipt_optimized.json, SHA256SUMS. Other files in this scratch directory are
exploratory and are not part of the result packet.
