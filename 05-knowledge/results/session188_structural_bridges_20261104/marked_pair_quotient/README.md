# Marked-pair quotient separation

PROOF.md gives an exact grouped-quartic criterion for quotient converse with
unchanged arbitrary child interiors. A pair child is enough when deleting its
quotient vertex leaves an ordered cut. Strong prime W5-type examples are
included; the regular five-cycle quotient remains outside this corollary.

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The checker is standalone and standard-library only. It retains all sixteen
literal finite-difference queries for each grouped quartic. It checks fifteen
seeded general grouped-coefficient cases, 121 arbitrary small ordered cuts,
and three explicit strong-prime controls, including a ten-name coefficient224.
The proof dependency is the adjacent response/PROOF.md, whose hash is pinned.

Publication files: PROOF.md, README.md, verify.py, receipt.json,
receipt_optimized.json, SHA256SUMS. Independent review passed in the adjacent marked_pair_quotient_audit packet.
