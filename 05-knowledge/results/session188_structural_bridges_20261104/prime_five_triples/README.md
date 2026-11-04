# Prime-five quotients with arbitrary triple children

PROOF.md and the complete retained certificate establish a finite fifteen-name
class: every strong prime five-vertex quotient inflated by arbitrary named
three-vertex tournaments is separated from its quotient converse. This is not
an all-order theorem. All cubic/quartic effects at the actual query are retained.

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The verifier requires g++ and the adjacent marked_pair_quotient/verify.py helper.
It builds the independent physical subset-DP control in a temporary directory;
no installation is performed. It checks all 1,024 quotients, all 8 triple children,
all 96 representative/profile patterns, and eighteen physical 17-name counts.
Both per-pattern query values and 32-coefficient profile polynomials are retained.

Publication allowlist: PROOF.md, README.md, verify.py, literal.cpp, receipt.json,
receipt_optimized.json, SHA256SUMS. Exploratory probe.py/probe_receipt.json and the
local compiled binary are not part of the packet. Independent review passed in the adjacent prime_five_triples_audit packet.
