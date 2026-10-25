# Session 55: exact CRT collision wiring

Main result: fixed-cell equality partitions compose by intersection under CRT. Their minimal-mask antichains are sufficient for exact non-Parker detection and constructive witnesses. The least counterexample to classification by local Boolean Parker flags is 35=5·7. Complementary three-block partitions prove a larger infinite family: every product of coprime odd integers at least 5 is non-Parker.

Read THEOREM.md for proofs and boundaries. Reproduce with:

    python verify.py --limit 160 --output census.json

`census.json` contains the prime-power local antichains, all 159 modulus results, and reconstructible witnesses with roots. Two full-modulus enumerators and CRT composition must agree exactly on minimal masks. `SHA256SUMS` pins the main packet files.

Additional result: CLASSIFICATION.md reduces all quotient rings to 203 CRT candidates using affine Lo Shu constructions at prime powers, the orthogonal AP construction, and the separately proved prime-field classification. `python classify_candidates.py` reproduces classification_candidates.json: 136 Parker moduli n≥2, maximum3216; odd maximum201. Independent odd-ring review reproduced all 203 classifications; see `session55-audit-odd-rings`. Add `--verify` to either reproducer to compare the recomputed data with the saved certificate rather than overwrite it.
