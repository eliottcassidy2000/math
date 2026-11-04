# Independent prime-flow extension audit

See [AUDIT.md](AUDIT.md) for the interval proof, primitive-prime exceptions,
exact arithmetic examples, propagation soundness, and scope review.

    python3 verify_independent.py > receipt.json
    python3 -O verify_independent.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The checker imports no producer module. It reads the sibling
`primeflow_extension` packet to verify replay equality and pin the audited
proof and checker. Complete modular periods independently establish the
capped row domains used in the arithmetic counterexamples.
