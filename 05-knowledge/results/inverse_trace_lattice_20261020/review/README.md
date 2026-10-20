# Independent inverse-trace-lattice review

The uniform mathematical audit is in [AUDIT.md](AUDIT.md).

The independently written controls use multiplication matrices for perturbed
non-binomial packets, exact numerator/denominator Newton arithmetic against a
known root, exhaustive small split-lattice allocations, and the stronger
full-rank unimodular hostile.

Reproduce from this directory:

```bash
python3 independent_controls.py > independent-results.json
python3 -O independent_controls.py > independent-results.optimized.json
cmp independent-results.json independent-results.optimized.json
```

Both modes pass **505 checks** with byte-identical JSON. Python 3 and SymPy
are required. The script itself performs no file writes. Run with
`--progress` for an optional diagnostic on standard error.

The independent controls supplement the proof audit. They do not prove
unconditional trace regularity, construct a Keller map, or resolve JC(2).

The final `REVIEW_RECEIPT.json` records the exact reviewed producer hashes,
the normal/optimized producer replay, and the independent review artifacts.
