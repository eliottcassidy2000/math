# Tail-scale certificate obstruction

New Session22 proof and finite native witness subset. A literal obstruction refutes the proposed depth+6 quadratic-tail contraction; it does not refute uniform153N coverage or prove failure at every later cut. Conditional sufficient scale-cone inequalities, the exact source-domain recursion, and the equivalence of uniform old-source completion with entry retain their full scope.

- [Proof and producer evidence](session22-tail-scaling/README.md)
- [Independent review](session22-tail-review/REVIEW.md)
- [Required published base residual packet](../routefree_suffix_gate_20261021/README.md)

The independent audit authenticates 4,096 supplied pending witnesses through 454,656 literal native edges with full clocks and mod15/mod35 checks. Its 1,396,405 checks do not constitute a complete new-annulus sweep. Unlisted new sources and full new-scale energy remain uncomputed here.

## Reproduce

From this packet directory inside the repository:

```sh
python3 session22-tail-scaling/check.py --repo-root ../../.. --base-dir ../routefree_suffix_gate_20261021/session22-residual --base-receipt ../routefree_suffix_gate_20261021/session22-residual/receipt.json --verify session22-tail-scaling/receipt.json
python3 -O session22-tail-scaling/check.py --repo-root ../../.. --base-dir ../routefree_suffix_gate_20261021/session22-residual --base-receipt ../routefree_suffix_gate_20261021/session22-residual/receipt.json --verify session22-tail-scaling/receipt.json
python3 session22-tail-review/audit.py --repo-root ../../.. --packet session22-tail-scaling --base-dir ../routefree_suffix_gate_20261021/session22-residual --verify session22-tail-review/receipt.json
python3 -O session22-tail-review/audit.py --repo-root ../../.. --packet session22-tail-scaling --base-dir ../routefree_suffix_gate_20261021/session22-residual --verify session22-tail-review/receipt.json
```

All selected producer/review files are unchanged. The original producer manifest's pending-review wording is historical; the included final independent review supersedes it. Paths are explicit command arguments, not source changes. Both checkers authenticate dependency pins. Baseline repository commit: fe180d8ea4356bd31c7c3258a23982b0fd879d4a. The top manifest pins every included file except itself.
