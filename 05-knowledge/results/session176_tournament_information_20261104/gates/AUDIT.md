# Independent audit: shared gates and exact memory

The reviewer checked the arguments and independently reran the supplied checker in normal and optimized Python modes on 2026-11-04. Both outputs exactly match the committed certificate files.

- Deleting the apex yields two increasing parts. Endpoint terms and the two possibilities for their extreme vertices give the stated quadratic polynomial, including the adjacent-pair correction.
- Its mixed derivative is the corresponding quadratic coefficient, negative for every pair. This establishes a complete interaction graph, not a count of distinct gates.
- Across an ordered cut, the geometric coefficient factors as the displayed outer product plus one adjacent-entry correction. The interior 2-by-2 minor has determinant -4, establishing rank exactly two over Q; end cuts have rank one.
- The simultaneous register update uses at most six binary arithmetic operations with constants 0, 1, 2. Growing operand bit length is explicitly distinguished from Boolean gate count.
- The two suffix responses differ by twice the prefix's binary integer. Hence every length-k prefix has a different residual response, requiring at least 2^k deterministic states, or k bits, under the stated final-output streaming model. Sending the entire prefix attains the one-way bound.

Finite rerun coverage: 511 literal tournaments, 1,023 derivative base assignments, 4,950 cut factorizations, and 8,191 prefixes. Checks use explicit exceptions and survive `python3 -O`.

The all-length statements follow from the proofs, not the finite audit. The packet does not prove a general Boolean circuit lower bound, novelty priority, or a lower bound for arbitrary nonlinear cut encodings. No mathematical correction was required.
