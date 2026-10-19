# Focused correction to the THM-4716 intercept index

Status: PROVED local repair, with FINITE-EXACT regression. This is not a whole-theorem audit or a status promotion.

## Exact scope and mechanism

For integers q >= 2, 1 <= p < q, gcd(p,q)=1 and 0 <= t <= q-2, let

    f_t(k) = floor((kp+t)/q).

Its unique changed position when t increases to t+1 is the representative

    k1 = p^(-1)(q-1-t) mod q,  1 <= k1 <= q-1.

Indeed the increment occurs exactly when kp+t is congruent to q-1 modulo q. Inverting multiplication by the unit p gives the displayed index. The index cannot be zero because q-1-t is nonzero modulo q. Taking successive differences of f_t proves that the mechanical word changes from 01 to 10 at positions (k1-1,k1), without wrapping. This repairs the exact index while preserving the adjacent-swap conclusion.

The published proof omitted the modular inverse and wrote k1 = q-1-t mod q. The first hostile in lexicographic (q,p,t) order is (3,2,0): the old/new words are 011 and 101; actual k1=1, while the formula gives 2. At (5,2,0), the words 00101 and 01001 give actual k1=2, not 4. The original formula is correct in the control lane p=1.

## Verification

Run from this directory:

    python3 -B verify_intercept_index.py
    python3 -B -O verify_intercept_index.py

Both runs check 52,807 transitions with 2 <= q <= 64. The direct floor-increment location and adjacent word-change positions are computed independently of the modular-inverse expression. There are 48,994 failures of the original expression. Explicit checks remain enabled under optimization. Matching outputs are retained as verification.out and verification_optimized.out.

The changed theorem paragraph is recorded in PATCH_PROPOSAL.json, with its baseline and replacement hashes. No other theorem claim, scope, dependency, or status is changed. Other claims, including boundary-case formulations, cycle consequences, literature attribution, and Mandelbrot identifications, were not re-audited here.
