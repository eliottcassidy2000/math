# Independent review of the serial low-cycle obstructions

Verdict: **PASS** for the frozen three-file packet below. No correction is required.

## Reviewed snapshot

- `session197-lowcycle/SERIAL_OBSTRUCTIONS.md`: SHA256 `42a4c6c0bd5aabb1b47e3b1a84359a79fdc3d5bdc90ab71306a07ac93396ad40`.
- `session197-lowcycle/verify_serial_obstructions.py`: SHA256 `2b684e8b6565cfcecd9ce18a9dd1e46112be1b32da86b457ab60c5be38720489`.
- `session197-lowcycle/serial_obstructions_certificate.json`: SHA256 `ac6191c6849b5f2eca64af15523c993af14331709a96119e10e814cfa3e263a4`.

The verifier's named dependencies `verify.py`, `verify_rooted.py`, and `unique_start.py` were read/copied only as required by this audit. The source packet was not edited.

## Proof checks

1. The cycle-extension argument preserves the cyclic order of the named original vertices, including the two-new-vertices insertion when every outside vertex is uniform toward the current cycle. Strong connectivity guarantees the required cross edge. This proves injectivity and the heredity of the componentwise cycle bound.
2. Every cycle crossing a block cut must use its unique backward spine edge. In a spanning cycle the edges through the singleton spine are consecutive. The residual traversal is precisely one A-path starting at the designated root and one B-path ending at its root, yielding exactly four cycles.
3. The minimality argument covers every strong induced subset. Missing an intermediate singleton or a needed endpoint root leaves a forward cut. If both endpoint blocks remain and the induced subset is proper, one selected endpoint has at most three vertices, and its rooted path response is at most one. The other response is at most two. The order-two exceptions cannot be strong. Applying this to each component handles arbitrary weak proper induced subsets.
4. The matrix readout `(5,-1,9) S^k (5,1,9)^T` agrees with the stated rooted data and one-seam operation. `S^3=S^2+S+I` proves the tribonacci law, rather than fitting it to sampled counts.
5. The two endpoint triangles and the disjoint consecutive spine triples prove the stated feedback-vertex lower bound. No equality claim about the feedback number is being made.
6. Distinct orders plus induced minimality give an infinite antichain and rule out a finite forbidden-induced list for this class. This does not imply any obstruction to a finite-state description in a separately specified ordering.

## Independent computation

- Replayed the frozen checker in an isolated copy, both normally and with `python -O`. Outputs agreed byte for byte, and the regenerated certificate had the exact frozen SHA256 above.
- Independently constructed every graph by literal adjacency rules, without using the source `glue` function.
- Used predecessor-recursion subset dynamic programming, separate from the forward anchored dynamic program in the source, to count cycles on every subset of `F_0,...,F_7` (orders 8–15).
- All 65,272 proper subsets had at most two Hamiltonian cycles; each full graph had four. Since every component vertex set is itself among those proper subsets, this also checks the componentwise conclusion.
- Independently counted H as `105,187,349,641,1177,2167,3985,7329`, agreeing with the exact recurrence.

The check is `independent_check.py`; detailed results are `INDEPENDENT_RESULTS.json`.

## Additional exact-span-two observation

The author improved the preliminary span-three ordering, and this improvement was independently checked. Put A's local vertices in order `[3,0,1,2]`, then the singleton spine, then B's local vertices `[2,1,0,3]`. Both endpoint blocks have maximum backward distance two. The designated A root is last and B root first, so all reversed spine seams become adjacent; every other cross edge remains forward. This is a direct all-k span-two upper bound, independently checked through `k=100`.

The bound is exact. Every strong tournament with a span-one ordering is the reversed-adjacency chain `R_n`, by the audited ribbon classification, and has exactly one Hamiltonian cycle: all adjacent backward arcs are forced and the forward first-to-last arc closes the cycle. But every `F_k` has four cycles. Hence its minimum backward span is exactly two.

The author added this exact-span-two conclusion and a self-contained span-one exclusion to the final frozen packet. The amended proof and checker were reviewed and the full normal/optimized replay again reproduced the new certificate hash listed above. This addendum is included in the PASS verdict.
