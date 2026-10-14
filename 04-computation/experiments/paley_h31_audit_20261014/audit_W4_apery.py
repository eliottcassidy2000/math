#!/usr/bin/env python3
"""AUDIT of HYP-9274 (B): Sigma*(W_4) = q * W(q),  W(q) = sum_{a,b,c,d distinct nonzero} chi(abcd(c-a)(d-a)(c-b)(d-b)),
claimed W(q)/(q-1) = -a_q - 4 with f = eta(2z)^4 eta(4z)^4.

Independent of the author's scripts:
 (1) a_n of f by Euler's pentagonal theorem for prod(1 - x^n), then powers/convolutions (author: in-place products);
 (2) a_p cross-checked by Beukers' supercongruence A((p-1)/2) = a(p) mod p^2 (Apery numbers
     A(n) = sum_k C(n,k)^2 C(n+k,k)^2; proved by Ahlgren-Ono 2000), which pins a(p) for p >= 17 since |a(p)| < p^2/2;
 (3) W(q) by brute force from the definition (own loops) at small primes and at q = 113, 127, 131;
 (4) W(q) by the exact O(q^2) reduction (referee's derivation):
        W = (q-1) [ sum_{l != 0,1} chi(l) g(l)^2 + (q - 3) ],   g(l) = sum_t chi(t(t-1)(t-l)),
     for all primes q < QMAX;
 (5) the Greene / Ahlgren-Ono link:  q^3 4F3(1) (Greene normalisation) = sum_{u != 0} chi(u) g(u)^2,
     checked against the triple sum chi(-1) sum_{x,y,z} chi(xyz(1-x)(1-y)(1-z)(1-xyz)) and against -a(p) - p;
 (6) bonus: the W_4 sum over GF(27) against the Frobenius trace alpha^3 + beta^3 (a_3 = -4).
Usage: python audit_W4_apery.py [QMAX=1000]
"""
import sys, math
import numpy as np
import numba as nb

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def eta_coeffs(N):
    """a_n, n <= N, of x * E(x^2)^4 E(x^4)^4 with E(x) = prod (1 - x^n) = sum_k (-1)^k x^(k(3k-1)/2) (pentagonal)."""
    E = [0] * (N + 1)
    k = 0
    while True:
        done = True
        for kk in (k, -k) if k else (0,):
            e = kk * (3 * kk - 1) // 2
            if e <= N:
                E[e] += (-1) ** (abs(kk) % 2)
                done = False
        if done and k > 0:
            break
        k += 1
    E = np.array(E, dtype=object)

    def mul(A, B):
        C = [0] * (N + 1)
        for i in range(N + 1):
            if A[i] == 0:
                continue
            ai = A[i]
            for j in range(N + 1 - i):
                if B[j]:
                    C[i + j] += ai * B[j]
        return np.array(C, dtype=object)

    def subst(A, m):            # A(x^m)
        C = [0] * (N + 1)
        for i in range(N + 1):
            if i * m <= N:
                C[i * m] = A[i]
        return np.array(C, dtype=object)
    E2 = subst(E, 2)
    E4 = subst(E, 4)
    P = E2 * 0
    P[0] = 1
    for _ in range(4):
        P = mul(P, E2)
    for _ in range(4):
        P = mul(P, E4)
    a = [0] * (N + 2)
    for i in range(N + 1):
        if i + 1 <= N + 1:
            a[i + 1] = int(P[i])
    return a


def apery_mod(n, m):
    s = 0
    for k in range(n + 1):
        s += (math.comb(n, k) ** 2) * (math.comb(n + k, k) ** 2)
    return s % m


def chi_table(q):
    chi = np.full(q, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, q):
        chi[(x * x) % q] = 1
    return chi


@nb.njit(cache=True)
def W_brute(q, chi):
    tot = 0
    for a in range(1, q):
        for b in range(1, q):
            if b == a:
                continue
            for c in range(1, q):
                if c == a or c == b:
                    continue
                for d in range(1, q):
                    if d == a or d == b or d == c:
                        continue
                    tot += chi[a] * chi[b] * chi[c] * chi[d] * chi[(c - a) % q] * chi[(d - a) % q] \
                        * chi[(c - b) % q] * chi[(d - b) % q]
    return tot


@nb.njit(cache=True)
def g_table(q, chi):
    g = np.zeros(q, dtype=np.int64)
    for l in range(q):
        s = 0
        for t in range(q):
            s += chi[(t * ((t - 1 + q) % q)) % q] * chi[(t - l + q) % q]
        g[l] = s
    return g


@nb.njit(cache=True)
def triple_sum(q, chi):
    tot = 0
    for x in range(q):
        cx = chi[x] * chi[(1 - x + q) % q]
        if cx == 0:
            continue
        for y in range(q):
            cy = chi[y] * chi[(1 - y + q) % q]
            if cy == 0:
                continue
            xy = (x * y) % q
            for z in range(q):
                tot += cx * cy * chi[z] * chi[(1 - z + q) % q] * chi[(1 + q * q - xy * z) % q]
    return tot


def gf27():
    """GF(27) = F_3[y]/(y^3 + 2y^2 + 1) (referee's own model), elements a0 + 3a1 + 9a2"""
    def dg(a):
        return [a % 3, (a // 3) % 3, (a // 9) % 3]

    def num(d):
        return d[0] + 3 * d[1] + 9 * d[2]

    def mul(a, b):
        A, B = dg(a), dg(b)
        pr = [0] * 5
        for i in range(3):
            for j in range(3):
                pr[i + j] += A[i] * B[j]
        for k in (4, 3):
            c = pr[k]
            pr[k] = 0
            pr[k - 1] += c
            pr[k - 3] += 2 * c
        return num([x % 3 for x in pr[:3]])

    def sub(a, b):
        A, B = dg(a), dg(b)
        return num([(A[i] - B[i]) % 3 for i in range(3)])
    sq = {mul(x, x) for x in range(1, 27)}
    chi = {0: 0}
    for x in range(1, 27):
        chi[x] = 1 if x in sq else -1
    return mul, sub, chi


def main():
    N = max(QMAX, 140) + 5
    a = eta_coeffs(N)
    print('a_p (pentagonal construction), p < 60:', {p: a[p] for p in primes_upto(60)})
    # (2) Beukers' supercongruence check
    bad = []
    for p in primes_upto(min(N, 1200)):
        if p < 17:
            continue
        r = apery_mod((p - 1) // 2, p * p)
        # lift to the symmetric residue
        if r > p * p // 2:
            r -= p * p
        if r != a[p]:
            bad.append(p)
    print(f'(2) Beukers A((p-1)/2) = a(p) mod p^2 with |a(p)| < p^2/2 pins a(p) for 17 <= p <= {min(N, 1200)}: '
          f'{"all agree" if not bad else "MISMATCH at " + str(bad)}')
    # (3) brute force at small primes and 113, 127, 131
    brute_q = [p for p in primes_upto(70) if p >= 3] + [113, 127, 131]
    allok = True
    for q in brute_q:
        chi = chi_table(q)
        w = int(W_brute(q, chi))
        ok = (w == -(q - 1) * (a[q] + 4))
        allok &= ok
        print(f'(3) q = {q:3d} ({q % 4} mod 4): brute W = {w:9d};  -(q-1)(a_q+4) = {-(q - 1) * (a[q] + 4):9d};  {"OK" if ok else "FAIL"}')
    print('(3) brute-force agreement:', allok)
    # (4) O(q^2) reduction for all primes < QMAX; (5) Greene/Ahlgren-Ono link
    fails4, fails5 = [], []
    maxratio = (0, 0)
    for q in primes_upto(QMAX):
        if q < 3:
            continue
        chi = chi_table(q)
        g = g_table(q, chi)
        tw = int(sum(int(chi[l]) * int(g[l]) ** 2 for l in range(2, q)))     # l != 0, 1
        W = (q - 1) * (tw + q - 3)
        if W != -(q - 1) * (a[q] + 4):
            fails4.append(q)
        G = int(sum(int(chi[u]) * int(g[u]) ** 2 for u in range(1, q)))      # = q^3 4F3(1) (Greene), claimed
        if G != -a[q] - q:
            fails5.append(q)
        r = abs(a[q]) / (2 * q ** 1.5)
        if r > maxratio[0]:
            maxratio = (r, q)
    print(f'(4) O(q^2) reduction, all primes 3 <= q < {QMAX}: W/(q-1) = -a_q - 4 fails at {fails4 if fails4 else "none"}')
    print(f'(5) sum_(u != 0) chi(u) g(u)^2 = -a(p) - p (Ahlgren-Ono, Greene normalisation) fails at {fails5 if fails5 else "none"}')
    print(f'    max |a_p|/(2 p^1.5) over p < {QMAX}: {maxratio[0]:.4f} at p = {maxratio[1]}')
    for q in [p for p in primes_upto(100) if p >= 3]:
        chi = chi_table(q)
        ts = int(triple_sum(q, chi)) * int(chi[q - 1])
        g = g_table(q, chi)
        G = int(sum(int(chi[u]) * int(g[u]) ** 2 for u in range(1, q)))
        if not (ts == G == -a[q] - q):
            print(f'(5) TRIPLE-SUM MISMATCH at q = {q}: chi(-1) triple = {ts}, twisted moment = {G}, -a-p = {-a[q] - q}')
            break
    else:
        print('(5) chi(-1) sum_{x,y,z} chi(xyz(1-x)(1-y)(1-z)(1-xyz)) = sum_(u!=0) chi(u) g(u)^2 = -a(p) - p for all primes 3 <= p < 100')
    # (6) GF(27)
    mul, sub, chi = gf27()
    tot = 0
    for A in range(1, 27):
        for B in range(1, 27):
            if B == A:
                continue
            for C in range(1, 27):
                if C in (A, B):
                    continue
                cc = chi[A] * chi[B] * chi[C] * chi[sub(C, A)] * chi[sub(C, B)]
                for D in range(1, 27):
                    if D in (A, B, C):
                        continue
                    tot += cc * chi[D] * chi[sub(D, A)] * chi[sub(D, B)]
    a3 = a[3]
    # alpha + beta = a3, alpha beta = 27 ; alpha^3 + beta^3 = a3^3 - 3*27*a3
    tr27 = a3 ** 3 - 3 * 27 * a3
    print(f'(6) GF(27): W = {tot}, W/26 = {tot / 26}; -(alpha^3+beta^3) - 4 = {-tr27 - 4} (alpha^3+beta^3 = {tr27}); '
          f'-a_27 - 4 = {-a[27] - 4} (a_27 = {a[27]})')


if __name__ == '__main__':
    main()
