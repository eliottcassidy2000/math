# Short physical slope transport

Status: PROVED and independently audited. The frozen proof snapshot retains its pre-audit status; `audit/AUDIT.md` records completed review of its exact SHA256. This is a separate reviewed addendum to [signed fiber repairs](../signed-fiber-repair/README.md).

- A fixed physical block O^ell Z^j transports every reduced interior slope with odd denominator to an even numerator, preserving the denominator and adding at most 2 ceil(log2 b) vertices.
- Every strong unit profile with 5U-9D=24 is exactly q.OOOZ for a strong unit profile with 8U_q-9D_q=29. The predecessor profile is unique.
- For the selected C2 left source, the entire repaired family has n=79232j+30133 and physical membership requires g(18j+7)=16j+3. Universal positive-target coverage remains OPEN.

Reproduce from this directory:

    python3 verify.py > /tmp/block-transport-normal.json
    python3 -O verify.py > /tmp/block-transport-optimized.json
    cmp /tmp/block-transport-normal.json /tmp/block-transport-optimized.json
    python3 audit/check_block_transport.py
    python3 -O audit/check_block_transport.py

The author verifier tests 1053 slopes through denominator 101, full-coordinate identities, literal small graphs, a declared bounded equivalence check, and the large exact positive witness. Independent checks cover 16382 slopes through denominator 401, including 8187 odd numerators, and 100000 fiber controls. The proofs establish the all-order statements.

The frozen helper module is imported only for source operations, graph counting, and the g recurrence. Its main routine is not executed by this verifier. No external package is needed.
