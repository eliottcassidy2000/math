# Minimum-cover factorization of clone amplitude

The leading coefficient of a physical transitive-clone family factors exactly as a_(nu+1)=nu! L_C R_C for every minimum crossing cover. The two positive factors count explicitly defined endpoint-constrained half covers. This gives a direct all-order Hamilton-word bijection, including zero matching rank and empty sides.

- [Proof](PROOF.md).
- [Independent review](REVIEW.md): PASS, no correction.
- [Standalone C++ verifier](verify.cpp) and [Python runner](verify.py).
- [Independent literal permutation/slot verifier](independent/literal.py).

Run `python3 verify.py` and `python3 -O verify.py` in this directory; C++17 and Python3 are sufficient. The source census checks 9,774 roots and 17,002 minimum covers. Run `python3 independent/literal.py` for the separate complete labelled order-five check. No external data bank or code dependency is needed.

All audited files are preserved byte for byte. Scratch-style source paths in the frozen proof/review identify their original session packets; the corresponding files in this directory have identical hashes. The separately cited full polynomial lift is an optional stronger companion, not a dependency of the elementary scalar proof or either verifier. This checkpoint makes no arbitrary-factor realization, unit-offset congruence, or full-spectrum claim.
