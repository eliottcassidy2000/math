# Complete congruence obstructions for a fixed prefix

For a fixed physical rooted prefix (h,e,R), let d=gcd(h,R). An integer passes every congruence obstruction for its final-gadget image exactly when it is odd and is congruent to0 or−e modulo d. The proof gives the exact residue image for every modulus and constructive physical witnesses.

This is local admissibility, not exact-target sufficiency. The same parity branch must be used across the prime factors of a composite modulus.

- [Proof](PROOF.md)
- [Standalone verifier](verify.py), with both receipts
- [Independent physical/BFS audit](audit/REPORT.md) and standalone checker
- [Joint modular universality](dependencies/MODULAR_GADGET_UNIVERSALITY.md)
- [Inherited arithmetic law](dependencies/ARITHMETIC_INVERSE.md)
- [Physical source grammar](dependencies/ROOTED.md)

Run `python3 verify.py` and its `-O` form; run `audit/check.py` similarly. Neither checker needs external data or writes files. The mathematical dependency proofs are included unchanged.
