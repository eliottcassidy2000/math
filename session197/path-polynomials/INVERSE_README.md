# Exact arithmetic inverses for a physical module language

These three packets extend the earlier arbitrary-module theorem. Each proof and verifier is frozen and independently reviewed. Run `python3 audit_inverse.py` to replay all three verifiers in normal and optimized Python and compare every byte with the frozen receipts.

- `RATIO_LATTICE.md`: exact cancellation of the binary mode is equivalent to integer prime-valuation lattice membership. Smith normal form gives a physical ordered-join witness or an all-order obstruction for a specified finite bank.
- `TARGET_RAY_INVERSE.md`: prescribed amplitude and offset add bounded nonnegative-copy and parity constraints. The complete finite algorithm preserves both prescribed coefficients and nonempty modules.
- `TWO_GENERATOR_CLASSIFICATION.md`: for bank K1/H5, one singleton leaf and any nonempty transitive leaf, cancellation occurs only at transitive length four. The complete family is 15360*4^t*1024^s*3^m-65*25^s; every ray has exact common divisor five and no prime or 7p output.

Independent reports:

- `session197-bridges/ratio-lattice-audit/REPORT.md`, SHA256 `0c2d6d4f2566c93105dee79abd68ea3a78e5646c8ff72a502ecc1e5d1353c8e3`.
- `session197-target-ray-audit/REPORT.md`, SHA256 `3de256278edf4b20a79b895b6aff5947299f4ce239203d18a0b377b76f5c3566`.
- `session197-bridges/two-generator-audit/REPORT.md`, SHA256 `0d8e66dae27e917bbc9a61ee0edfddb0fcbb1f6610eb20efec90825c48f3c15c`.

These are scoped positive and negative results for specified physical construction languages. They do not establish the full tournament spectrum. The separately frozen primitive-cancellation construction is not part of this manifest.
