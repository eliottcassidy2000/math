# Uniform equivariance obstruction

20 October 2026. **HYP-9172(a) PROVED for every n≥5; independent audit PASS.**

There is no S_n-equivariant map inducing a bijection of isomorphism types from tournaments to sign-even graphs, or from tournament switching classes to sign-even Euler graphs, for any n≥5. Arbitrary bijections of unlabelled types remain possible.

The [proof](UNIFORM_EQUIVARIANCE_OBSTRUCTION.md) uses only primes 3 and 5. A signed Walsh determinant forces every hypothetical matched stabilizer to preserve all odd-prime valuations. Overlapping 3-cycles would create a forbidden double transposition; disjoint 3-cycles force a triple-block partition with insufficient 5-part. The empty graph requires both full ambient Sylows, yielding a contradiction. No additive-prime conjecture or group classification is used.

The [independent audit](review/AUDIT.md) checks the general reversal quotient, all invariant edge subspaces through n=4, and both principal settings through n=5. The producer's independent finite controls extend through n=6. The all-order argument is a proof, not an extrapolation.

The packet also proves an all-order joint-prime obstruction and supplies an explicit Hall-deficient family. HYP-9172(b)'s branch integrality and pairing question are handled separately.

Reproduce: `python3 -B exact_controls.py --max-n 6`; independent checks: `python3 -B review/independent_controls.py` and `python3 -B review/group_controls.py`. The original proof and audit bytes are preserved; staging paths in the audit identify the reviewed originals.
