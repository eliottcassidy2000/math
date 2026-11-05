# Modular universality of physical one-apex gadgets

For every even modulus, the gadget index L and Hamilton count g(L) realize every joint residue pair with g odd. The construction gives arbitrarily large witnesses. For modulus 2^k, a witness uses O(k²) vertices.

Consequently strong one-apex tournaments with exactly one Hamilton cycle attain infinitely many counts in every parity-compatible congruence class, and each fixed rooted low-cycle prefix hits every odd dyadic residue.

These are residue results. They do not prove exact-target coverage, cofiniteness, or infinitely many prime or7p counts.

- [Proof](PROOF.md)
- [Standalone verifier](verify.py), with normal and optimized receipts
- [Independent review](audit/REPORT.md) and standalone checker
- [Inherited arithmetic law](dependencies/ARITHMETIC_INVERSE.md)
- [Physical gadget compiler proof](dependencies/ROOTED.md)

Run `python3 verify.py` and `python3 -O verify.py`; run `audit/check.py` similarly. Both checkers use only the Python standard library, have no data dependencies and write only stdout. Dependency proofs are included unchanged; session-relative references describe their original packets.
