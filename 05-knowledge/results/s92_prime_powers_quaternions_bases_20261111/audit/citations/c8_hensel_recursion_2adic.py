#!/usr/bin/env python3
"""Auditor 'citations' - an elementary proof-by-recursion of THM-5135 (3)(b) at the prime 2, modulo EVERY 2^k.

Let Q be an integer-valued quadratic form in n variables whose polar form B(x, y) = Q(x+y) - Q(x) - Q(y) is unimodular
at p (det B a p-adic unit); p any prime (for p = 2 this forces n even and Q 'even unimodular').
N_k(c) = #{x mod p^k : Q(x) = c mod p^k}.  Split x into primitive (x != 0 mod p) and x = p x'.
  * primitive: each primitive solution modulo p^k has exactly p^(n-1) lifts modulo p^(k+1) in each class above c
    (the map y -> B(x, y) mod p is onto), so  P_(k+1)(c') = p^(n-1) P_k(c' mod p^k);
  * imprimitive: Q(p x') = p^2 Q(x'), so I_k(c) = 0 unless p^2 | c, and I_k(p^2 c'') = p^n N_(k-2)(c'')  (k >= 3),
    I_2(0) = p^n, I_1(0) = 1.
Hence all N_k are determined by the class counts modulo p alone.  For p = 2 and an even unimodular form these are
(Z, 2^n - Z) with Z = 2^(n-1) + 2^(n/2-1) (Arf invariant 0) or 2^(n-1) - 2^(n/2-1) (Arf invariant 1).
E8 and the hyperbolic form of rank 4 both have Z = 136: SAME class counts modulo every 2^k.

This script (a) runs the recursion from Z = 136 and compares it with the brute-force class counts of the E8 form
modulo 2, 4, 8, 16 (my enumeration, c2_local_class_counts*.out) and with the formula of THM-5135 (1), r = 4, for
k <= 14; (b) does the same for the rank-2 seam (Z = 10) against the table square, k <= 10; (c) runs the odd-p
version for four squares modulo 3^k, 5^k, 7^k against the rank-2 formula.
"""
from fractions import Fraction


def recursion_counts(p, n, counts_mod_p, K):
    """N_k as lists for k = 1..K from the class counts modulo p."""
    N = {0: [1], 1: list(counts_mod_p)}
    P = {1: [counts_mod_p[0] - 1] + list(counts_mod_p[1:])}
    for k in range(2, K + 1):
        q = p ** k
        Pk = [p ** (n - 1) * P[k - 1][c % (p ** (k - 1))] for c in range(q)]
        P[k] = Pk
        Nk = list(Pk)
        if k == 2:
            Nk[0] += p ** n
        else:
            for c2 in range(p ** (k - 2)):
                Nk[(p * p * c2) % q] += p ** n * N[k - 2][c2]
        N[k] = Nk
    return N


def vp(c, p, k):
    if c % (p ** k) == 0:
        return k
    v = 0
    while c % p == 0:
        c //= p
        v += 1
    return v


def rank_law(p, k, r, c):
    v = vp(c, p, k)
    if v < k:
        return (1 - Fraction(1, p ** r)) * sum(Fraction(1, p ** ((r - 1) * i)) for i in range(v + 1))
    return (1 - Fraction(1, p ** r)) * sum(Fraction(1, p ** ((r - 1) * i)) for i in range(k)) + Fraction(p ** k, p ** (r * k))


def main():
    ok = True
    # (a) E8 / rank-4 seam at 2
    brute = {
        1: [136, 120],
        2: [17536, 15360, 17280, 15360],
        3: [2246656, 1966080, 2211840, 1966080, 2242560, 1966080, 2211840, 1966080],
    }
    N = recursion_counts(2, 8, [136, 120], 14)
    for k in (1, 2, 3):
        same = N[k] == brute[k]
        print("(a) recursion from Z = 136, modulo %d: equals the enumerated E8 class counts: %s" % (2 ** k, same))
        ok &= same
    bad = 0
    for k in range(1, 15):
        q = 2 ** k
        for c in range(q):
            if Fraction(N[k][c] * q, q ** 8) != rank_law(2, k, 4, c):
                bad += 1
    print("    recursion = (1 - 2^-4) sigma_3(2^v)/2^(3v) law of THM-5135 (1), r = 4, modulo 2^k for k = 1..14: failures %d" % bad)
    ok &= bad == 0
    print("    modulo 16 from the recursion:", N[4])
    # the other 2-adic class (Arf 1, H^3 + A) for contrast
    N1 = recursion_counts(2, 8, [120, 136], 4)
    print("    contrast, Z = 120 (H^3 + A): modulo 8:", N1[3])

    # (b) rank-2 seam at 2: Z = 10 of 16
    N2 = recursion_counts(2, 4, [10, 6], 10)
    bad = 0
    for k in range(1, 11):
        q = 2 ** k
        for c in range(q):
            if Fraction(N2[k][c] * q, q ** 4) != rank_law(2, k, 2, c):
                bad += 1
    print("(b) rank-2 seam at 2 (Z = 10): recursion = (1 - 2^-2) sigma(2^v)/2^v law modulo 2^k, k = 1..10: failures %d" % bad)
    ok &= bad == 0

    # (c) odd p: four squares, counts modulo p by brute force, then the recursion against the rank-2 law
    for p, K in ((3, 8), (5, 6), (7, 5), (11, 4)):
        cnt = [0] * p
        for a in range(p):
            for b in range(p):
                for c in range(p):
                    for d in range(p):
                        cnt[(a * a + b * b + c * c + d * d) % p] += 1
        Np = recursion_counts(p, 4, cnt, K)
        bad = 0
        for k in range(1, K + 1):
            q = p ** k
            for c in range(q):
                if Fraction(Np[k][c] * q, q ** 4) != rank_law(p, k, 2, c):
                    bad += 1
        print("(c) four squares modulo %d^k, k = 1..%d: recursion from the counts modulo p = rank-2 law: failures %d" % (p, K, bad))
        ok &= bad == 0
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
