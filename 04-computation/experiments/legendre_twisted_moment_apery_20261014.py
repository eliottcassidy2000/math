"""Twisted second moment of the Legendre family vs the Beukers-supercongruence form f = eta(2z)^4 eta(4z)^4.

M(q) = sum_l chi(l) g(l)^2 with g(l) = sum_c chi(c(c-1)(c-l)); claim M(q) = -(q + a_q(f)) at every odd prime
(Ahlgren-Ono Lemma 2.2 + Theorem 6; S34 audit).  Scaling one rim vertex gives U(W_4) = q(q-1) M(q).
Session opus-2026-10-14-S35.
"""
import numpy as np
N = 160
P = [0] * (N + 1); P[0] = 1
for n in range(1, N + 1):
    for step in (2 * n, 4 * n):
        if step > N: continue
        for _ in range(4):
            for i in range(N, step - 1, -1):
                P[i] -= P[i - step]
a = [0] * (N + 2)
for i in range(N + 1): a[i + 1] = P[i]

def primes(n):
    return [p for p in range(3, n) if all(p % d for d in range(2, int(p ** 0.5) + 1))]

print(' q  q%4   M(q)=sum_l chi(l) a_l^2   -(q + a_q)   M2(q)=sum_l a_l^2   q^2-2q-3')
for q in primes(110):
    chi = np.zeros(q, dtype=np.int64)
    sq = set((x * x) % q for x in range(1, q))
    for x in range(1, q): chi[x] = 1 if x in sq else -1
    xs = np.arange(q)
    base = chi[(xs * (xs - 1)) % q]
    al = np.array([int(np.dot(base, chi[(xs - l) % q])) for l in range(q)])   # sum_c chi(c(c-1)(c-l))
    M = int(np.dot(chi, al * al))
    M2 = int(np.sum(al * al))
    print(f'{q:3d}  {q % 4}   {M:10d}              {-(q + a[q]):8d}     {M2:8d}        {q * q - 2 * q - 3}')
