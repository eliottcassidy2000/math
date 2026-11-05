# Independent dyadic fiber and height audit

Status: PASS. Frozen proof SHA256 `6fe7ce5a86b2e819e9bba72038e8d6335556869d8c1b0166a532686f9f521bd9`; verifier SHA256 `c7d3ff81d5eea7190eb857109c4535a9df19e4c73750bdef69d0d1bd94eadce5`.

The exact O/Z polynomial and division by two were checked. Its odd quadratic support is exactly the square of the path on word positions. The generalized-binomial bit argument is valid for negative coefficients/values by finite formal-series truncation. At critical length the full Boolean monomial can arise only from disjoint quadratic terms in every top-degree bit component; the matching partition multiplicity is odd. The first-vertex matching split proves the Fibonacci recurrence, and the critical Fibonacci index is a power of two, never divisible by three. The empty-product/empty-cube case r=1 is explicitly valid.

The theorem counts indexed words and distinct labelled constructions, not exact scalar fibers or graph isomorphism classes. Different O/Z words change a named root/newborn edge, and later moves preserve that edge. The proof keeps this scope.

The same-pair height criterion is sound: oddness upgrades congruence to lcm(m,2), and the strict product bound excludes both neighboring representatives. The primitive balanced family v=w=(2^a+1,1,2^a+1,1), with its actual zero-J word t.Z^a.O^a.J, attains product n+m and shows the inequality cannot be weakened. Its balanced/total-J bookkeeping matches the inherited adjoint.

Independent `check_dyadic.py` checks four physical prefix profiles and both closures at critical and next depth for r=1..4 (64 fiber checks), exact mixed coefficients through block length8, matching parity through r=13, and the sharpness family through a=8. Normal and optimized receipts match. The full author verifier also passes normally and optimized with identical receipts. No unresolved mathematical issue found; all-target coverage is expressly open.
