# Independent cusp-arithmetic review

**Verdict: PASS.** See `independent_cusp_operations_review.md` for the
proof audit and its precise scope.

## Reproduce

Run with Python 3 and SymPy:

    python independent_cusp_operations_controls.py
    python -O independent_cusp_operations_controls.py

The script has no input-file, network, or filesystem-write dependency.
It was checked with Python 3.12.14 and SymPy 1.14.0. Its Decimal interval
arithmetic and exact-rational arithmetic use Python's standard library.
The two captured outputs are byte-identical. Normal and optimized replays
of the author's mpmath-based verifier are also included.

The interval controls independently enclose pi by
pi=16 arctan(1/5)−4 arctan(1/239), using consecutive alternating partial
sums. All Decimal arithmetic endpoints are rounded outward. Correctly
rounded square roots are enlarged by one representable value at each
endpoint. The C_2 enclosure is an inherited certified input, not a
numerical estimate fitted by this verifier. Original x-variable iterates
are propagated without evaluating an escape-threshold table.

`source_fingerprints.json` identifies the audited author packet and its
dependencies. `MANIFEST.sha256` covers this review packet, excluding
itself. The scripts and receipts are finite controls, not a calibration of
the existential guard or a proof of prime abundance.
