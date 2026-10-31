# Prime-coordinate covariance and the numerical-label boundary

**PROVED elementary; FINITE-EXACT controls.** No novelty claim for the standard prime-exponent encoding of FRACTRAN.

## Covariance theorem

Let P=(p1,...,pr) be distinct primes containing every prime dividing the numerators and denominators of an ordered reduced FRACTRAN list. Let Q=(q1,...,qr) be any other distinct primes. On the free multiplicative monoids define E(product p_i^e_i)=product q_i^e_i. Replace every fraction a/b by E(a)/E(b), preserving order.

Then E is a bijection of these monoids, and conjugates the entire partial computation: an instruction is eligible exactly when its image is eligible; the first eligible index, halting, instruction word, and exponent-vector trajectory are preserved. Proof: b|N iff every denominator exponent is at most the corresponding exponent of N. E preserves these coordinate inequalities and addition/subtraction of valuations. Apply induction to successive states. The same result extends to all positive integers after extending the finite prime bijection to a permutation of all primes.

Therefore a statement intrinsic to this abstract register machine must be invariant under prime renaming. The numerical sizes of prime labels, ordinary additive identities between encoded integers, triangular membership, and Fibonacci membership are extra structure. They may be meaningful for Conway's chosen encoding, but are not consequences of the prime-register transition structure alone.

## Exact hostile and priority sidecar

Map the first ten primes (2,3,5,7,11,13,17,19,23,29) to (31,37,41,43,47,53,59,61,67,71). The instruction numerator55, both T10 and F10 in the conventional indexing F0=0,F1=1, becomes41*47=1927. This is neither triangular (8*1927+1 is not a square) nor Fibonacci. Every instruction trace remains conjugate.

The product of an instruction list is also insufficient to describe computation. At seed2, Conway's list starts with15/2, producing15. The reversed list starts with55, producing110. Their total rational products are equal. A block-product map must retain eligibility and original priority if it is claimed to simulate computation.

## Reproduction

Run `python3 check_encoding.py` and `python3 -O check_encoding.py`. The independent integer-level simulation verifies25,500 paired transitions over all seeds<=500 supported on the first ten primes, at most100 transitions each. This finite test supplements the all-state proof above. It does not claim all seeds<=500 have the required prime support, or that all unbounded trajectories were enumerated. Failure checks remain active under optimization.

Primary input: J.H. Conway, *FRACTRAN: A Simple Universal Programming Language for Arithmetic* (1987), printedpp4 and9, [chapter scan](https://gwern.net/doc/cs/computable/1987-conway.pdf). The exact ordered14-fraction list is embedded in the checker. The theorem uses the definition, not the prime-generating correctness theorem.
