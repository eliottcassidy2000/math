# Independent review of the order-ten arbitrary-ambient certificate

**PASS.** Frozen producer proof SHA256 `e2ce7ca477a2eec99b16e87f5bafad62826040e153c7b3fbba6544ff693fc3b1`; verifier `faba9a4dfdcc00dbb28fe28dc55914ab59307fbade91b5a36a60bb57f13ceade`. The independent audit adds only order ten. The earlier range is reviewed separately in `equal_nine_ambient_audit`.

## 1. Exact finite Hamilton inputs and the split bound

The producer's exhaustive endpoint recursion covers all 64 labelled order-four and 32,768 labelled order-six tournaments. Its integer storage is safe: every path count is at most 720, and all histogram bins and universe counts fit the declared types. The encoding is explicitly lexicographic with bit one meaning the lower name points to the higher name.

The independent audit uses a different counting direction and no tournament path dynamic program. For every permutation, its five (or three) successive directed edges are fixed; it adds one to every tournament completion of the remaining unordered edges. Thus each tournament receives exactly one contribution for each of its directed Hamilton words. The complete histograms, maxima and first attaining codes agree with the producer: maxima 5 and 45, attained first at codes 4 and 408. A hash of the complete count vector in the common code order is retained for each universe.

For an order-ten Hamilton word, the first-six-name set and its directed prefix and suffix uniquely recover the word. Ignoring the connecting edge only enlarges the image. Hence

    H <= sum over six-name subsets of H(prefix set) H(complement)
      <= binom(10,6) 45 5 = 47,250.

No independence of the two induced tournaments is assumed, and 47,250 is not called the exact ten-name maximum.

The inspected [OEIS A038375 entry](https://oeis.org/A038375) agrees with the numerical provenance. Its order-ten value is not used in this proof. The review did not inspect the books cited there and does not attribute a first proof to them.

## 2. Complete continuous-domain certificate

The actual order-ten odd-cycle-cover profile is `t10 + w2 t2 + w4 t4 + w6 t6 + w8 t8`, with nonnegative weights and Hamilton count `1+sum(w)`. The split bound therefore gives exactly the continuous simplex budget 47,249. The unit top coefficient and nonnegative slack are both essential; arbitrary unbounded positive profile cones are not covered.

The independent script rederives the full physical two-port run-word kernel through exponent ten, rather than taking the producer numerator as an input. It multiplies that series by the independently enumerated ordinary cyclic path-word series for each additional cover slot. Transitive profiles are computed using inclusion-exclusion for ordered surjections rather than the producer Stirling recurrence. Finally it evaluates the symmetric trilinear functional at simplex vertex profiles; this independently performs the homogeneous cubic substitution, including the multiplicities of repeated simplex indices.

All 1,085 coefficients are explicitly compared with the complete source list, including zeros, in all 31 coordinates `h=2,...,32`. Every coefficient is nonpositive. The exponent support vanishes afterwards, so no higher coordinate is omitted. The audit retains each coordinate's complete-list hash, extrema and count; the complete numerical lists remain in the producer receipts. The independent normal and optimized receipts agree byte-for-byte, and both producer modes are independently replayed and match their frozen receipts.

## 3. Strictness and the physical outer lift

The inherited scalar identity gives `D1=-P(q)<0` for every nonzero nonnegative finite profile. Combining it with the checked nonpositive higher coordinates gives componentwise domination of the two port-augmented cores, strict at the Hamilton coordinate.

For an arbitrary fixed rooted exterior, both ports receive the original root's incidence to every exterior name. This is the precise physical assignment that makes the augmented core a common module. The exact module-run substitution has nonnegative coefficients and first coefficient `H(Q)>0`, so the actual named-response contrast is at most `-H(Q)P(q)<0`. The factorization is a linear combination of full ordered-cover counts; no invalid multiplicativity of the child profile functional is used. Target interiors may change only while preserving the common complete profile on the same physical child sets.

This is an exact arbitrary-fixed-ambient result at order ten, and with the prior packets through ten. It does not assert a result at order eleven, an unequal-profile theorem, or invariance under changes involving exterior names. The earlier positive coefficient on a larger continuous domain is compatible with this smaller valid domain and was never an exact response collision.
