# Independent Session213 proof and computation audit

`AUDIT.md` records the proof review and its boundaries. `APPROVED_PACKETS.json`
pins every approved research file by SHA-256. `REPLAY_RECEIPT.json` records
normal and optimized replays of the authors' frozen verifiers.

The checkers here were written without importing an author's research code.
The C++ endpoint classifier uses degrees and a transitive remainder, whereas
the author's classifier uses strong components and a chain reconstruction.

From this directory, reproduce the independent checks:

```bash
python3 independent_counts.py > /tmp/s213-independent.json
cmp independent-counts.json /tmp/s213-independent.json
python3 audit_comb.py > /tmp/s213-comb.json
cmp comb-normal.json /tmp/s213-comb.json
python3 audit_covers.py > /tmp/s213-covers.json
cmp covers-normal.json /tmp/s213-covers.json
python3 audit_deletions.py > /tmp/s213-deletions.json
cmp deletion-normal.json /tmp/s213-deletions.json
python3 audit_support.py ../cycle-to-path/RESULT.json > /tmp/s213-support.json
cmp support-independent.json /tmp/s213-support.json
g++ -std=c++17 -O2 -Wall -Wextra -pedantic audit_endpoint_equality.cpp -o /tmp/s213-endpoint
/tmp/s213-endpoint ../../session212/cycle-levels/classes7.rows ../../session212/cycle-levels/classes8.rows > /tmp/s213-endpoint.jsonl
cmp endpoint-equality-audit.jsonl /tmp/s213-endpoint.jsonl
```

Repeat Python checks with `python3 -O`; the first four have separately retained
byte-identical optimized receipts. The C++ check was also rebuilt at `-O0`
and reproduced its complete receipt. It enumerates every labeled tournament
through order six without external data. Its supplementary order-seven/eight
class universe is inherited from Session212, with the exact input hashes in
`SOURCE_VALIDATION.json`; it is not a new graph-isomorphism census.

Every proof's all-order conclusion is justified by the written argument,
not by extrapolation from these finite checks. The unrestricted tournament
Hamiltonian-path spectrum conjecture remains open.
