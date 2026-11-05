#!/usr/bin/env python3
"""gf_identity.py (klein-2026-11-04): the gluing count.  E(s) = sum over ordered pairs k + l = s, k, l >= 2 of F(2k-3) F(2l-3)
counts the gluings of a uniquely Hamiltonian tournament rooted at a sole entry (order k+1; F(2(k+1)-5) = F(2k-3) rooted
classes, THM-4868 (F)) with one rooted at a sole exit (order l+1); D_n = E(n-1) (vertex junction) + E(n-2) (arc junction).
Checks D_n = (n-2) F(2n-9) - L(2n-10) and the generating function x^5 (1+x)(1-x)^2 / (1-3x+x^2)^2 for 6 <= n < 60."""
F = [0, 1]
for i in range(130): F.append(F[-1] + F[-2])
L = [2, 1]
for i in range(130): L.append(L[-1] + L[-2])
def E(s): return sum(F[2 * k - 3] * F[2 * (s - k) - 3] for k in range(2, s - 1))
N = 60
ok = all(E(n - 1) + E(n - 2) == (n - 2) * F[2 * n - 9] - L[2 * n - 10] for n in range(6, N))
num = {5: 1, 6: -1, 7: -1, 8: 1}; den = [1, -6, 11, -6, 1]; ser = [0] * N
for i in range(N): ser[i] = num.get(i, 0) - sum(den[j] * ser[i - j] for j in range(1, min(i, 4) + 1))
ok2 = all(ser[n] == E(n - 1) + E(n - 2) for n in range(6, N))
print('E(n-1) + E(n-2) = (n-2)F(2n-9) - L(2n-10) for 6 <= n < %d: %s' % (N, ok))
print('series of x^5(1+x)(1-x)^2/(1-3x+x^2)^2 equals E(n-1)+E(n-2) for 6 <= n < %d: %s' % (N, ok2))
print('D_6..D_16 =', [E(n - 1) + E(n - 2) for n in range(6, 17)])
# the gluing counts by (|S2|, |PR|) at orders 6-12 (vertex junction k + l = n - 1; arc junction k + l = n - 2)
for n in range(6, 13):
    print('order %d: vertex junction (k,l) -> F(2k-3)F(2l-3): %s ; arc junction: %s' % (n,
          {(k, n - 1 - k): F[2 * k - 3] * F[2 * (n - 1 - k) - 3] for k in range(2, n - 2)},
          {(k, n - 2 - k): F[2 * k - 3] * F[2 * (n - 2 - k) - 3] for k in range(2, n - 3)}))
