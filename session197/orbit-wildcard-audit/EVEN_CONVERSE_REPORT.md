# Independent even-converse congruence and spectrum audit

**PASS.** The all-order congruence `H=1 modulo4` for even-order self-converse tournaments is correct. The complete arbitrary-order self-converse spectrum through 73 is exactly the 25-value set stated in the frozen packet. The omissions, including the prime targets 53 and 59 and the induction seed 63, are genuine limitations of this construction class. No source correction is needed.

## Frozen source snapshot

All inputs are in `session197-orbit-wildcard/`:

- `EVEN_CONVERSE.md`: `0fc0e7ca35a1d48626ffd1b49251b64dcc0564f66594bfd0e70307287c4c4cbb`.
- `verify_even_converse.py`: `ee8c15a3a27a2b5f1540d17eb04f9e71084fc3fec612659903c53bed2a413e77`.
- `even-converse-normal.json`: `dd1e6716e9fa29d7cf2764c013c4ac75b2531cfac4098f36da6e479d50edb4df`.
- `even-converse-witness.json`: `08c74a5182999514685de72edd9c776970aa49324fa1e20816c93f6ceb5cd2cc`.

The earlier mass and through-41 spectrum packets remain unchanged and have separate audit reports.

## Mathematical review

The inherited odd-cycle collection formula gives `H=1+2c modulo4`, where `c` counts individual directed odd cycles modulo cyclic rotation. Applying an anti-involution to vertices and reversing cyclic order is a well-defined involution on those directed cycles. Consequently `c` has the parity of its individually invariant cycles.

An even-order anti-involution has no fixed vertex: it has at most one, because two fixed vertices would contradict arc reversal, and its fixed-point count is even. No odd-sized vertex set can then be invariant. Thus no directed odd cycle is individually invariant, proving the stated congruence. At odd order an invariant odd cycle must contain the unique fixed vertex, but need not be every cycle through it. The proof does not confuse these two assertions. The supplied order-three and order-four controls correctly rule out dropping the even-order hypothesis or strengthening the result to `1 modulo8`.

For the complete support, the already audited SCC identity writes every self-converse count as `d^2 h`, with `h` a strong self-converse count or one. The Busch bound was independently checked against the official paper during the preceding spectrum audit. The exact inequality `390625>389017` bounds every relevant middle strong component by order eight. This does not bound an arbitrarily padded weak tournament.

The finite normalized universe is complete by the anti-involution existence and fixed-Hamilton-word argument from the mass theorem. Closing its support under squares `1,9,25,49` creates no extra values below 74. Allowing the factor `d=7` at this upper-bound stage is harmless: it enlarges the possible source factors, without asserting that an unrestricted count-seven graph exists. Every value retained in the resulting set has its own explicit self-converse witness, so no spurious attainment is inferred from that enlargement.

The shorter odd-order reduction for the omitted `3 modulo4` targets is also sound. For 63, the only relevant nontrivial square factor is nine, which would require a central self-converse count seven; the complete bounded strong bank excludes it. Thus paired components and singleton padding cannot restore that seed.

## Independent exact verification

`check_even_converse.py` imports no author verifier. It independently generates the normalized pair-orbit universe and counts all **4696** presentations through order eight using vertex-cover inclusion-exclusion over directed walks, rather than either source endpoint recurrence. The full support is

    {1,3,5,9,11,13,15,17,25,27,29,31,33,37,41,
     45,49,51,57,61,65,67,69,71,73}.

All **4165 even-order presentations** satisfy the congruence. For all **600 presentations through order seven**, the audit independently enumerates every directed odd cycle, checks closure under the reflected-cycle involution, identifies its fixed cycles, and verifies both congruence expressions. Invariant cycles at odd order contain the fixed vertex; at even order there are none. These are direct finite controls of the all-order cycle proof, not its replacement.

The checker verifies every emitted positive witness and anti-involution. It reconstructs all six new unrestricted two-apex graphs from their physical words, checks the literal rows, and independently obtains their Hamilton counts by both walk inclusion-exclusion and full vertex permutations. Exact converse-isomorphism tests additionally certify the displayed escape graphs are non-self-converse.

The independent normal and optimized receipts agree exactly. Source replay uses an audit-only copy, so generated witness files do not touch the frozen source directory; its normal and optimized receipts and witness bytes reproduce the frozen packet.

## Reproduction and scope

Run `python check_even_converse.py` and `python -O check_even_converse.py` in this audit directory. Its helper dependencies are the earlier independent graph, mass and spectrum audit checkers.

- Independent checker SHA-256: `59b115b2a774ebdd1661b9324d1866a0080cf33c78aca4a96029853c5efbaaf1`.
- Independent receipt SHA-256: `6c0137183ac9e14efc0cbf6d2e3efa466175e931751cec134e78bc8e0bdf661a`.

The certified result excludes these values only from the self-converse class. It neither introduces new unrestricted spectrum gaps nor proves universal prime or `7p` coverage.
