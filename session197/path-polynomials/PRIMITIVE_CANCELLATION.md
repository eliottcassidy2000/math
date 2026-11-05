# A primitive cancellation language containing prime and 7p counts

Session 197. **PROVED physical construction and complete fixed-bank cancellation classification.** Unlike the two-generator language in `TWO_GENERATOR_CLASSIFICATION.md`, this language has no fixed prime divisor. It includes an explicitly certified prime count and an explicitly certified 7p count. No infinitude of primes or complete spectrum claim is made.

## 1. Two concrete module profiles

Let K=K1 and let the tournaments X,Y have outgoing-row encodings

    X=[54,500,152,243,64,340,261,113,157],
    Y=[32,17,3,7,13,30].

Thus X has nine vertices and Y has six. Their exact labelled-bin polynomials are

    B_X(x)=x^9+40x^7+354x^5+648x^3+104x,
    B_Y(x)=x^6+12x^4+24x^2.

In particular

    h_X=1147=31*37, b_X=22352=2^4*11*127,
    h_Y=37,         b_Y=352=2^5*11,

where h=H and b=B(2). The checker derives these polynomials from actual Hamilton-path counts on every vertex subset, followed by the exact set-partition/labelled-bin transform. They are not virtual input profiles.

Use the seven-vertex quotient `[14,32,2,22,103,13,47]` with nonempty modules

    (Tr_m, K, L, K, Tr_7, K, M), m>=1.                    (1)

The two fixed singleton modules displayed between the other modules are distinct physical vertices. Allow L,M to be any nonempty ordered joins of copies of the bank K,X,Y.

## 2. Complete cancellation classification for this bank

For the fixed leaves P=K and R=Tr_7, the ratio criterion from `RATIO_LATTICE.md` says that the binary mode vanishes exactly when

    rho_L/rho_M=127/62, rho_Z=B_Z(2)/H(Z).                (2)

Indeed d_P=1,e_P=2,d_R=127,e_R=-124. The bank ratios are

    rho_K=2,
    rho_X=2^4*11*127/(31*37),
    rho_Y=2^5*11/37.

Let z_K,z_X,z_Y be the copy counts in L minus those in M. The 127-adic valuation of (2) forces z_X=1. The 37-adic valuation then forces z_Y=-1. The 2-adic valuation gives z_K=0. The remaining prime valuations, at 11 and 31, agree automatically. Thus the unique difference vector is

    (z_K,z_X,z_Y)=(0,1,-1).

All nonnegative allocations are therefore exactly

    L: t copies of K, s+1 copies of X, j copies of Y,
    M: t copies of K, s copies of X, j+1 copies of Y,
    t,s,j>=0.                                            (3)

Conversely every allocation (3) has ratio (2), so it cancels. Both modules are automatically nonempty. Any ordering of each module's factors has the same responses; a canonical construction lists K factors, then X factors, then Y factors. This is an exact all-order classification, not a bounded search over copy counts.

## 3. The actual ternary rays and their content

The audited arbitrary-module formula in `UNIT_OFFSET_RIGIDITY.md`, specialized to the leaves in (1), is

    H_m=A*3^m+B*2^m+C,
    A=254*b_L*b_M,
    B=254*b_M*h_L-124*b_L*h_M,
    C=-125*h_L*h_M.

For (3), B=0 and

    A=A_0*4^t*(22352^2)^s*(352^2)^j,
    C=C_0*(1147^2)^s*(37^2)^j,
    A_0=1998447616=2^10*11^2*127^2,
    C_0=-5304875=-5^3*31*37^2.                           (4)

This gives physical strong tournaments of order

    m+25+2t+18s+12j.                                    (5)

To verify strongness, the quotient has directed spanning cycle 0->2->1->5->3->4->6->0. Every block is nonempty, so the blowup is strong even when L or M is an ordered join.

All prime factors of A lie in {2,11,127}, while all prime factors of C lie in {5,31,37}. Thus gcd(A,C)=1, C is odd, and 3 does not divide C. Since H_2-3H_1=-2C,

    gcd(H_1,H_2)=gcd(3A+C,2C)=1.

Consequently every ray in this entire cancellation language has exact common divisor 1 over positive m. This is stronger than merely failing to exhibit an obvious common divisor. It does not guarantee a prime in every ray.

There is also an exact 7-divisibility targeting rule. Modulo 7, A_0=4, C_0=5, 22352^2=1147^2=1, and 352^2=37^2=4. Therefore

    7 divides H_m  iff  m+4t=4 modulo 6.                 (6)

Here 3 has multiplicative order 6 modulo 7, and the factors 4^j cancel. The rule does not depend on s or j. It targets divisibility only; primality of H_m/7 still needs a separate certificate.

## 4. Explicit certified prime and 7p witnesses

Take t=s=j=0, so L=X and M=Y. The exact all-order count is

    H_m=1998447616*3^m-5304875.

At m=1, the constructed 26-vertex graph has

    H_1=5990037973, a prime.

At m=4, the constructed 29-vertex graph has

    H_4=161868952021=7*23124136003,

and 23124136003 is prime. Complete integer trial division certifies the two prime assertions: the respective square-root cutoffs are 77395 and 152066, and none of the integers from 2 through the relevant cutoff divides the candidate. The verifier performs every test, checks the exact floor-square-root inequalities and records a deterministic digest of all nonzero remainders.

Both full adjacency-row lists appear in the receipt, with vertex blocks specified in (1). Their large Hamilton-path counts are derived from the proved arbitrary-module identity and the independently reproducible small module profiles. We do not claim to have run a 26- or 29-vertex Hamilton-path subset DP. The verifier checks both full graphs are tournaments, checks all block orientations, their orders and strongness, and checks the root crossing graph has its explicit two-vertex cover and a two-edge matching.

The bank contains a nine-vertex tournament whose path count 1147 is itself composite. These examples are therefore an arithmetic construction of new rays from known physical modules, not an inverse reduction of every prime or 7p target to smaller prime inputs.

## 5. Reproduction and precise scope

Run `python3 verify_primitive_cancellation.py` and `python3 -O verify_primitive_cancellation.py`. The deterministic verifier recomputes the two complete physical input profiles, checks the displayed coefficients and factorization, reconstructs both large witnesses, completes the two primality proofs, and tests the allocation and residue formulas on finite parameter boxes. The proof above supplies the unbounded classification and content statements; no sampling or probabilistic primality test is used by the release verifier.

Dependencies:

- `UNIT_OFFSET_RIGIDITY.md`, SHA256 `3c116086bfd697eaeaf5fcbf46308f1e31cdf653eb2a651053f01bce46ce374b`.
- `RATIO_LATTICE.md`, SHA256 `a6f9c0b09d715865d01bbd62d95f2261c59e1f4a7eea2ed05cc09ac90aff60a3`.
- `verify_unit_offset_rigidity.py`, SHA256 `756ea2c09f32d4830a05ebf92bee079812b3b013d9526b173a04ca2f70bbab0f`, with its named exact graph/counting dependencies.

The finite exploratory search that found X,Y is not a premise of any theorem here and is not needed to reproduce the construction. The scope is precisely the quotient, fixed leaves and bank in Section 1. We do not assert infinitely many prime or 7p values, arbitrary prime targeting, or full {7,21} spectrum coverage.
