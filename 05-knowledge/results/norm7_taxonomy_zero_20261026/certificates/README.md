# Definedness and finite comparison certificates

**Status: PROVED elementary statements + CITED comparison hypotheses + FINITE-EXACT computations; independent audit PASS after three bounded repairs.** The [independent certificate audit](../taxonomy/AUDIT_CERTIFICATES.md) records the proofs, primary-source checks and frozen receipts.

- [Zero, bottom, quantifiers and the13-prime TT7 witness](ZERO_AND_FINITE_CERTIFICATES.md).
- [Galois comparison gates and explicit failures when data are forgotten](GALOIS_COMPARISON_GATES.md).
- [Frozen13-row integer input](tt7_rank13_input.json), [exact checker](check_certificates.py), [results](results.json).
- Independent audits: [quadratic fibre atlas](AUDIT_QUADRATIC_FIBRE_ATLAS.md) and [odd-prime conductor support](AUDIT_CONDUCTOR_SUPPORT.md).

The concrete arithmetic output is a nonzero13×13 determinant at primes≤103, certifying independence of the four exact raw TT7 sums modulo the nine inherited functions. This recovers and compresses the existing rank observation; it does not identify the proposed four motives. The second concrete output is a rank-eight Artin-representation counterexample to a dimension-free noncubic comparison shortcut: fourteen explicit primes pass, while prime43 separates the representations.

Run from this directory:

```sh
gp -q replay_modular_rows.gp > modular_rows.replay.out
cmp modular_rows.out modular_rows.replay.out
python3 check_certificates.py --output results.replay.json
python3 -O check_certificates.py --output results.optimized.json
cmp results.replay.json results.optimized.json
```

The checker uses Python's standard library; fresh modular coefficient replay uses PARI/GP. It does not import inherited producer code. The proof notes explicitly retain candidate-membership promises, rank, coefficient fields, residual hypotheses, semisimplification and arithmetic losses.
