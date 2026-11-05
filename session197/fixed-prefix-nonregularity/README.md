# Fixed-prefix synchronous nonregularity

The complete binary graph of the exact final-gadget count for any fixed rooted low-cycle prefix is not synchronous regular, even after a finite initial cutoff. The proof isolates a balanced-run language and transfers it through regular fixed-coefficient arithmetic.

This is not a scalar-spectrum omission, a nonregular-range theorem, or an obstruction to a suitable selected-witness subrelation. Range-only and cofinite-coverage approaches remain open.

- [Proof](PROOF.md)
- [Standalone verifier](verify.py), with normal and optimized receipts
- [Independent review](audit/REPORT.md) and standalone checker
- [Inherited arithmetic law](dependencies/ARITHMETIC_INVERSE.md)
- [Physical gadget compiler proof](dependencies/ROOTED.md)

The independent review identified a projection-padding gap in the first draft. The final proof explicitly trims forced common zeros by a regular left quotient before canonical normalization. The frozen review binds this corrected snapshot.

Run `python3 verify.py` and `python3 -O verify.py`; run the independent `audit/check.py` similarly. Neither checker writes files. All sources use only the Python standard library. The dependency proofs are included unchanged; their session-relative references describe their original packets.
