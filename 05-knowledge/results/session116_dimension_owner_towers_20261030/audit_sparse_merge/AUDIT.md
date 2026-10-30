# Independent incoming sparse-tower audit (C1–C5 only)

**PASS with a distribution-generalization wording caveat.** Source is `05-knowledge/results/cayley_dickson_integral_20261030/PROOFS.md`, pinned at incoming main `7f70dcd6012439b7cf75fb43f26850f3807126b9`. This is independent verification and source integration, not a novelty claim. No C6+ claims or automaton minimality are audited.

## Proof review

* **C1:** Left multiplication by a basis element is integral orthogonal. Multiplication by its transpose converts the two-support operator into `aI+bP`; the nontrivial XOR translation pairs coordinates. Each nonzero 2-by-2 block has entry gcd `g` and determinant `a²−epsilon b²`, proving the displayed cyclic cokernel decomposition, including singular and zero-coefficient cases. Global Smith factors need valuation reordering. Counting cyclic kernels proves the formula for every modulus, not only primes. This does not transport named right-hand sides without transformations.
* **C2:** The displayed doubling matrices are correct under the stated Baez convention. For distinct nonzero lower indices the right/left column-sign discrepancy flips exactly the four quaternionic coordinates. This yields the reflected `h−4−2k` count. The proof uses basis alternativity, not alternativity for arbitrary sums. Same-half and exceptional cross-half cases are consistent; recursion has one smaller input.
* **C3:** The three disjoint pair classes establish the polynomial recurrence. The support induction and constant term formula hold starting at dimension 8; differentiation cancels inherited first moments. The extrema are for two-support elements, not unrestricted elements.
* **C4:** The pair `(1,10)` has four positive blocks in dimension 16 in this convention, and embedding doubles this number. For any odd prime p, coefficients `(1,p−1)` make the positive determinants zero modulo p and negative determinants 2 modulo p; positive blocks retain rank one. Both real determinants are nonzero. This proves all odd primes with no distribution conjecture.
* **C5:** Expansion of the composition defect has only distinct-coordinate mixed rectangles; their coefficients are 0 or ±4. Independence/centering makes distinct rectangle monomials orthogonal. The rectangle count equals the C3 first moment. Unit-length Rademacher inputs therefore give exactly the stated mean and variance, and Chebyshev proves convergence in probability. This is not uniform control over the second input for a fixed multiplier.

**Wording caution:** The extension to general centered variance-one independent coordinates is valid for `D(x,y)=||xy||²−||x||²||y||²` and its second moment. It does not replace the fixed-length formula by `E(||XY||²−1)²` for arbitrary coordinate distributions. For instance normalized independent Gaussian inputs in dimension 8 have no composition defect but `Var(||XY||²)=4/8+4/8²=9/16`, not zero. No fourth moments are needed for D, but may be needed for product norm variance. The original random-sign claim is unaffected.

## Independent finite checks

`check.py` builds the Baez sign table directly from `(a,b)(c,d)=(ac−d conjugate(b), conjugate(a)d+cb)`. It imports no existing checker or parent tower witness. All unordered basis pairs and their XOR rectangles are enumerated at dimensions 2 through 128. The direct block counts agree with the named-pair recurrence, polynomial recurrence, support, zero count, first moment, and rectangle count. At dimension 128: 8,128 pairs, 1,186 zero-count pairs, 234,360 violated rectangles, variance 29,295/2,097,152.

Direct Gaussian elimination modulo 3,5,7,11 at dimensions 16,32,64 gives ranks 12,24,48 respectively for the C4 family. Another 336 exact exhaustive two-variable kernel counts cover positive/negative determinants, zero coefficients, singular blocks, nonprimitive coefficients, and moduli 1 through 12.

Both `python check.py` and `python -O check.py` pass with byte-identical JSON outputs; gates use explicit exceptions. `receipt.json` records exact counts and checker SHA-256. These finite checks corroborate the reviewed proofs; they do not replace the all-dimension arguments.

## Interface with Session116

The parent dense dimension-16 classification uses a different doubling convention; this audit deliberately does not transport named coordinates or witnesses between conventions. The two packets complement one another at the level of their explicitly stated results: arbitrary real multipliers at dimension 16 versus all-dimension integral two-support blocks and pair counts. Sparse worst-case cancellation failure can proliferate while two independent dense random-sign inputs show a shrinking average norm defect. Neither statement implies an arbitrary-dimension dense classification or operator-norm convergence.
