# Mahler safety, binary entropy, and Chan–Pak height

Session253. New result: every fixed finite core bank with BOTH no11 flanks has only O(X^log₂φ log X) outputs<=X. Requiring both finite Mahler-safe flanks improves the exponent to log₂(3/2). The proof combines the new sharp coercivity n<=3g(n) with the positive compiler and a hyperbola count. It assumes no rational witness selector.

The natural asymmetric repair was tested and refuted:52477 is the first target through65535 lacking any five-bank certificate with one Mahler-safe flank. Its four exact certificates expose either11 on both sides or11 on one side and10101 on the other. One-no11-flank coverage survives through65535. All one-safe successes in that range already use a safe side<=5.

- `PROOF.md`: proof, full scope, inheritance, bridge contract and live survivor.
- `probe.py`, `results.json`: new two-sided safety restrictions through4095 and coercivity controls.
- `asymmetric.py`, `asymmetric.json`: new one-sided and bounded-safe-side restrictions through65535.
- `hostile.py`, `hostile.json`: complete52477 inverse and independent graph DP.

Run each Python file from its directory or any working directory. It locates the inherited Session251 compiler relative to the repository. No theorem ID reserved. Chan–Pak, Mahler and LRC(14) remain OPEN.
