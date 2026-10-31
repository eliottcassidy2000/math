# Independent audit of the minimal threshold42 extension

**PASS.** A separate verifier, `audit_independent.py`, imports no producer module.
It independently decodes the old base data, checks all eleven new small bases,
replays all 196 universal certificates, and constructs its own recursive pairs
and streamed cycles. Normal and optimized runs produce byte-identical receipts.

- 22,149,425 checks;
- exactly 2,040 base pairs, indices42 through2081;
- 98 paired and98 terminal-cycle states;
- eight independent cycle controls, including a control consuming B51,
  both sides of the deeper recursion cutoff, and N=5,000,082;
- B2082 prefix regression checked;
- three hostile small-base mutations rejected.

The frozen tables do not algebraically require input index53. It was the lower
end of the previously available certified input supply, not a hypothesis of the
transition identities. At42 all blocks and tails are nonempty and required
prefix/endpoints exist. Affine ownership, forced member choice, endpoint square
identities, and paired position-shift equations depend only on the weak B
invariant and input parity. Thus the universal proof is valid at42.

The changed cutoff2081 is essential: larger weak bases do not certify P1. The
independent implementation uses bases only up to2081 and recurses from2082,
proving the stronger prefix invariant with S_(l+1)=49*S_l+24, S0=42.
The unchanged direct-cycle consumer then proves the unscaled P_l statement
for every N>=S_l and l>=1. This is still one explicit nested ray, not arbitrary
prescribed-path embedding.

Run independently (Python standard library):

```sh
python audit_independent.py --extension-dir . \
  --core-dir ../session124_recursive_square_prefix_20261031 \
  --base-dir ../square_sum_B_base_interval_20261019
python -O audit_independent.py --extension-dir . \
  --core-dir ../session124_recursive_square_prefix_20261031 \
  --base-dir ../square_sum_B_base_interval_20261019
```

`independent.normal.json` and `independent.optimized.json` are the corresponding
receipts. Their shared cycle hashes agree with the separate producer replay.
No claim about impossibility of smaller seed architectures is needed here.
