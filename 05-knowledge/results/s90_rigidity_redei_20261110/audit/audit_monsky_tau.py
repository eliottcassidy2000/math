#!/usr/bin/env python3
"""THM-5125 Remark 1 / README section 3 item 2 gloss the cubic x^2 + x + 1 as "the action of an element of order 3
on the 2-torsion of the curve".  Neither the packet nor the cited note proves a statement of that kind.  What CAN be
checked at matrix level (cited note bsd_selmer_growth_20261013a.md, eqs. (4)-(5)):

   for n = p_1 ... p_t with every p_i = 3 mod 8, the kernel of the Monsky matrix [[A + I, I], [I, A]] is stable
   under the order-3 map  tau(x, y) = (y, x + y)  of the coefficient pairs, and tau acts on the x-coordinate as
   A + I (= omega^2 for the F_4-structure "omega acts by A").

This script verifies that on random sets of primes = 3 mod 8, and shows as a control that tau-stability fails for
general odd squarefree n (so it is a feature of the class 3 mod 8, not of 2-descent in general).  It does NOT verify
that tau is induced by an automorphism of E_n[2] in Monsky's coordinates (not checked; E_n itself has no
automorphism of order 3, Aut(E_n) = mu_4).
"""
import random


def jacobi(a, n):
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
    return [i for i in range(N + 1) if s[i]]


def monsky(primes):
    """rows of the 2t x 2t Monsky matrix as (left t bits | right t bits), cited note eq. (4)."""
    t = len(primes)
    A = [0] * t
    for i in range(t):
        for j in range(t):
            if i != j and jacobi(primes[j], primes[i]) == -1:
                A[i] |= 1 << j
        if bin(A[i]).count("1") & 1:
            A[i] |= 1 << i
    D2 = [1 if jacobi(2, p) == -1 else 0 for p in primes]
    Dm2 = [1 if jacobi(p - 2, p) == -1 else 0 for p in primes]
    top = [((A[i] ^ (D2[i] << i)), (D2[i] << i)) for i in range(t)]
    bot = [((D2[i] << i), (A[i] ^ (Dm2[i] << i))) for i in range(t)]
    return top + bot, A


def kernel(rows, t):
    ker = []
    for x in range(1 << t):
        for y in range(1 << t):
            if all(((bin(r[0] & x).count("1") + bin(r[1] & y).count("1")) & 1) == 0 for r in rows):
                ker.append((x, y))
    return ker


def apply(A, x):
    out = 0
    for i, row in enumerate(A):
        if bin(row & x).count("1") & 1:
            out |= 1 << i
    return out


def main():
    random.seed(5125)
    P = [p for p in primes_upto(3000) if p > 2]
    P38 = [p for p in P if p % 8 == 3]
    ok = True
    stable = 0
    acts = 0
    trials = 0
    dims = {}
    for _ in range(300):
        t = random.randint(1, 6)
        ps = sorted(random.sample(P38, t))
        rows, A = monsky(ps)
        ker = set(kernel(rows, t))
        trials += 1
        d = len(ker).bit_length() - 1
        dims[d] = dims.get(d, 0) + 1
        if all((y, x ^ y) in ker for (x, y) in ker):
            stable += 1
        if all(y == apply(A, x) ^ x for (x, y) in ker):
            acts += 1
    print(f"primes = 3 mod 8: {trials} random n (1..6 prime factors): kernel tau-stable in {stable}; y = (A + I) x on the "
          f"kernel in {acts}; kernel dimensions {dict(sorted(dims.items()))}")
    ok &= stable == trials and acts == trials and all(d % 2 == 0 for d in dims)
    bad = 0
    tot = 0
    for _ in range(300):
        t = random.randint(2, 5)
        ps = sorted(random.sample(P, t))
        rows, A = monsky(ps)
        ker = set(kernel(rows, t))
        if len(ker) > 1:
            tot += 1
            if not all((y, x ^ y) in ker for (x, y) in ker):
                bad += 1
    print(f"control, arbitrary odd primes: among {tot} random n with nonzero Monsky kernel, tau-stability fails for {bad}")
    print("AUDIT MONSKY TAU:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
