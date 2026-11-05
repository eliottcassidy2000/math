# Height-bounded local lifting for physical sources

The main result is an exact target-dependent quadratic witness system for
strong C1 Hamilton-path counts. Its finite-field criterion uses one prime
p>4n and keeps explicit ordinary-height bounds. A complete CRT hostile at
n=21 explains why different prime witnesses cannot be merged as one physical
source word, even when the product modulus exceeds every possible count.

- `PROOF.md`: all-order construction, lifting theorem, scope, and hostiles.
- `verify.py`: explicit polynomial-circuit checks and independent modular
  trace closure; no disabled Python assertions.
- `VALIDATION_NORMAL.json`, `VALIDATION_OPTIMIZED.json`: identical receipts.
- `dependencies/PHYSICAL_ADJOINT.md`: frozen Session 199 physical C1 interface.

Reproduce:

    python3 verify.py > /tmp/hasse-source-normal.json
    python3 -O verify.py > /tmp/hasse-source-optimized.json
    cmp /tmp/hasse-source-normal.json /tmp/hasse-source-optimized.json

The proof is conditional on the named inherited physical grammar and C1
classification. Finite controls do not establish all-target coverage. No
unrestricted Hamilton-spectrum gap, number-field integrality definition,
or external-paper correctness conclusion is claimed.
