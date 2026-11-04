# Independent review of the ordered-cover ambient reduction

**PASS.** The exact reduction and its scope check. The proof does not assert coordinatewise domination for all tournament profiles.

An ordered h-cover places both ports in the same slot in h ways, or in distinct slots in h(h-1) ways. The other slots are nonempty old paths. Prefix and suffix old segments of a port-containing slot may be empty. Reading slots in their specified order gives the global order of each child's runs, so the common profile functional is correct without extra factorials.

The common-slot term is h F E^(h-1). The distinct-slot term is h(h-1) Fx Fy E^(h-2); it is absent for h=1. The old-only slot series E is symmetric. Swapping the first two child variables reverses the triangle and exchanges the two independent single-port factors. Hence the distinct-slot contrast is antisymmetric and its six-permutation sum vanishes. Combining this with the previously audited two-port numerator gives exactly h/3 times L_q(K E0^(h-1)/(1-XYZ)^(h+1)). The power of the denominator and the factor h are correct.

If all these coordinates have the proposed sign, nonnegative ambient run coefficients and the positive one-run coefficient do give strict transfer. The condition is not proved for every tournament profile. The artificial sequence with q2=0 and later positive entries is explicitly unrealizable as a tournament profile, and its positive h2 coordinate only blocks a universal nonnegative-sequence proof. It is not an actual ambient collision.

The independent checker imports only the separately written covered_quotient primitives. It verifies the separated-slot numerator antisymmetry by integer polynomial arithmetic, recomputes every ordered-cover coordinate directly for TT1, TT2, TT3 and C3 children, and confirms the artificial h2 value 146. No covered_ambient or covered_response producer code is imported. Normal and optimized receipts agree, and the exact reviewed source hashes are pinned.
