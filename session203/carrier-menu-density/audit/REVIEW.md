# Independent audit of the four-carrier density bound

6 November 2026. **PASS, restricted to the four specified paired O/Z orbits.** The addendum proves those orbits do not contain an odd tail. It says nothing comparable about all physical sources, other seed banks, J continuations, or F(8,137).

The reviewed proof has SHA-256 b0ebc37202e641c3eb0f55245b54dda8c61624489a7ac5a58dfeb402215732d1. Its producer checker has hash a4f54844e1a273d79f30b85b1a403f9af3283f2d6f3d256442c06a3e0e52419a. The input rows are the four actual coupled rows independently audited in the main staggered-allorder packet. This addendum does not infer coverage from them.

For a seed (X0,e0,V0,d), let A=V0-d>0. Both row transitions preserve V_k>=2^k A+d. The O step increases the spare d term; the Z step subtracts exactly d. Initially e0<=d, and every later e is either 0 or d. X is nondecreasing. Thus a word of length n whose last letter is Z increases X on that last step by at least 2^(n-1)A and has X>=X0+2^(n-1)A.

Deleting only the terminal O run preserves X. Therefore every represented base has a normalized representative which is empty or ends in Z. There are exactly 2^(n-1) such words of length n. The last-Z inequality bounds their possible lengths for X<=M. Summing the geometric series gives the safe per-seed bound 1+2M/A. Word collisions can only reduce the number of represented bases, so no injectivity assumption is used.

The four A values are 1012,898,1582,1046. Exact rational arithmetic gives

    sum 1/A = 693609731/187977134884.

Summing the four word bounds and dividing by M gives upper integer density at most 693609731/93988567442. Since all bases are odd, dividing instead by the number of odds at most M gives upper relative density at most 693609731/46994283721, strictly less than 3/200. An odd tail would have relative density one, so these particular four orbits cannot contain one. The additive constant four disappears in the limit.

The independent checker verifies the exact fractions, all 262,140 normalized last-Z words through length 16, and the closed form V_n=2^n V0-d times the binary subtraction sum for all 524,284 prefix nodes. It checks the geometric count at threshold and empty-word boundaries. Normal and optimized receipts match. The producer's normal and optimized replays also agree with its frozen receipt. These finite controls check the implementation; the induction and geometric count prove the all-order bound.

Reproduce with:

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

This result is separately scoped from the positive carrier identities and the full-state grammar obstruction. It does not claim that a finite bank including other seeds or a richer operation alphabet has the same density. Publication and status promotion require a separate final package check.
