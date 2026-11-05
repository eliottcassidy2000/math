# Exact count targeting inside the primitive cancellation language

Session 197. **PROVED finite constructive inverse for the complete language in `PRIMITIVE_CANCELLATION.md`.** The input is a target count N, rather than a prescribed pair of exponential coefficients. The output is every parameter recipe in that language, or an exact statement of absence from that language.

## 1. The bounded inverse

Retain the constants

    A0=1998447616, D0=5304875,
    bx=22352, by=352, hx=1147, hy=37,
    H0=3*A0-D0=5990037973.

The physical family has counts

    N=A0*4^t*bx^(2s)*by^(2j)*3^m-D0*hx^(2s)*hy^(2j),
    t,s,j>=0, m>=1.                                      (1)

If N is not an integer or N<H0, it is absent. Otherwise, for every pair s,j>=0 satisfying

    bx^(2s)*by^(2j) <= N/H0,                             (2)

form

    a=A0*bx^(2s)*by^(2j), d=D0*hx^(2s)*hy^(2j).

If a does not divide N+d, reject this pair. Otherwise put Q=(N+d)/a. Remove all factors of 2 and 3 from the positive integer Q. Accept exactly when the remaining factor is 1, the 2-adic exponent is even, and the 3-adic exponent is positive. The unique remaining parameters for this pair are

    t=v_2(Q)/2, m=v_3(Q).                                (3)

Every accepted tuple gives the actual graph recipe in `PRIMITIVE_CANCELLATION.md`. There is no empty-module boundary: the left module contains at least one X and the right at least one Y. Keep all accepted tuples if all recipes are wanted; distinct tuples are not assumed to have distinct counts.

## 2. Why this terminates and is complete

Let

    W=4^t*bx^(2s)*by^(2j)*3^(m-1).

Since bx>hx, by>hy and all exponents are nonnegative, W>=hx^(2s)*hy^(2j). Thus every count (1) satisfies

    N=3*A0*W-D0*hx^(2s)*hy^(2j) >= H0*W.                 (4)

In particular (2) is necessary. Both bases bx^2 and by^2 exceed 1, so (2) contains only finitely many pairs. Integer multiplication and comparison enumerate them without logarithmic rounding or factorization of N. For each pair, division and extraction of 2- and 3-adic exponents are finite and exact.

For necessity, an actual recipe obeys (2) by (4), so its pair is reached; dividing (1) by a then gives Q=4^t*3^m, and the algorithm recovers exactly t,m. For sufficiency, a successful pair and (3) satisfy (1) identically, and the already-proved physical construction realizes that count. If no pair succeeds, no graph in the specified language has count N.

This is a global absence certificate for this infinite language, not merely a failed search up to a chosen graph order. It is not an absence theorem for other tournaments.

## 3. Explicit positive and impossible targets

The inverse returns:

- N=5990037973: (t,s,j,m)=(0,0,0,1), a 26-vertex prime-count witness.
- N=161868952021: (t,s,j,m)=(0,0,0,4), a 29-vertex witness with count 7*23124136003.

Both primality statements are independently reproducible in `PRIMITIVE_CANCELLATION.md`.

The prime N=5990037979 is absent from this language. Here N/H0<2, so (2) forces s=j=0 and (4) also forces t=0,m=1. That would give H0, not H0+6. The checker certifies N's primality by complete trial division through floor(sqrt(N))=77395. Thus even this primitive prime-producing language misses a prime immediately above its first output. Content 1 and some successful prime targets do not imply full prime coverage.

## 4. Reproduction and dependencies

Run `python3 verify_primitive_target_inverse.py` and `python3 -O verify_primitive_target_inverse.py`. The checker verifies both positive targets, the missing prime, and every recipe in a finite parameter box against the exact inverse. It also builds a separately enumerated complete count table below 10^13 using bounds from (4), and compares exact inverse membership against that table on the members, their close neighbours and deterministic additional targets. The finite checks are controls; Section 2 proves all-order completeness.

The proof depends on `PRIMITIVE_CANCELLATION.md`, SHA256 `6c577db04d939196b04cbec9ba322f7f5abff7ae95e2f0fa19f7365c3acfd6ea`. The checker imports `verify_primitive_cancellation.py`, SHA256 `535d68803ccce6f1253e7ebf27c797bd9097665c2edc6e72cf897994c72ac0ce`, only for the exact constants and complete trial-division certificate routine. Its coefficient formula, physical realization and bank classification are inherited explicitly, not reproved by target enumeration.
