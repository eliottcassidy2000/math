#!/usr/bin/env python3
"""Test: W(q)/(q-1) = -a_q(f) - 4 for f = eta(2z)^4 eta(4z)^4, the weight-4 level-8 newform (Apery / zeta(3)).

W(q) is the wheel-W_4 character sum of paley_W4_charsum_20261014.py (hub 0, distinct nonzero a, b, c, d):
  W(q) = sum chi(a b c d) chi((c-a)(d-a)(c-b)(d-b)).
f = q * prod_n (1 - q^(2n))^4 (1 - q^(4n))^4 = sum a_n q^n is expanded exactly here (no tables used).
Also printed: the Deligne bound |a_p| <= 2 p^(3/2).  Session opus-2026-10-14-S34.
"""
import numpy as np
import numba as nb

NMAX = 140


def eta_product_coeffs(N):
    # P(x) = prod_{n>=1} (1 - x^(2n))^4 (1 - x^(4n))^4 up to x^N; f = x * P
    P = [0] * (N + 1)
    P[0] = 1
    for n in range(1, N + 1):
        for step in (2 * n, 4 * n):
            if step > N:
                continue
            for _ in range(4):
                for i in range(N, step - 1, -1):
                    P[i] -= P[i - step]
    a = [0] * (N + 2)
    for i in range(N + 1):
        if i + 1 <= N + 1:
            a[i + 1] = P[i]
    return a


@nb.njit(cache=True)
def W(q, chi):
    tot = 0
    for a in range(1, q):
        ca = chi[a]
        for b in range(1, q):
            if b == a:
                continue
            cab = ca * chi[b]
            for c in range(1, q):
                if c == a or c == b:
                    continue
                cabc = cab * chi[c] * chi[(c - a) % q] * chi[(c - b) % q]
                for d in range(1, q):
                    if d == a or d == b or d == c:
                        continue
                    tot += cabc * chi[d] * chi[(d - a) % q] * chi[(d - b) % q]
    return tot


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(2, n + 1) if s[i]]


a = eta_product_coeffs(NMAX)
print('f = eta(2z)^4 eta(4z)^4 coefficients a_p (computed):', {p: a[p] for p in primes_upto(60)})
ok_all = True
for q in primes_upto(110):
    if q < 5:
        continue
    chi = np.zeros(q, dtype=np.int64)
    sq = set((x * x) % q for x in range(1, q))
    for x in range(1, q):
        chi[x] = 1 if x in sq else -1
    w = int(W(q, chi))
    lhs = w // (q - 1) if w % (q - 1) == 0 else None
    ok = (lhs is not None) and lhs == -a[q] - 4
    ok_all &= ok
    print(f'   q = {q:3d} ({q % 4} mod 4): W/(q-1) = {lhs},  -a_q - 4 = {-a[q] - 4},  match: {ok};  |a_q|/(2 q^1.5) = {abs(a[q]) / (2 * q ** 1.5):.3f}')
print('ALL MATCH' if ok_all else 'MISMATCH')
