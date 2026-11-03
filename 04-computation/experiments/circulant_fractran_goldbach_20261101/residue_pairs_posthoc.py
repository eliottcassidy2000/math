#!/usr/bin/env python3
"""mac-mini-2026-11-02-triangles: POST-HOC (observed in quadratic_triangle_n4.out, then proved; not pre-registered).

For even m >= 6 let P(m) be the least prime = +-1 (mod m). Let R(m) be the set of unordered residue pairs {r, m - r}
(r odd, 3 <= r <= m - 3) such that both r and m - r occur as residues mod m of odd primes p < P(m), and let G(m) be
the set of Goldbach partitions {p, m - p} of m into odd primes.
Proposition (elementary). If m - 1 or m + 1 is prime, then R(m) = G(m). In general G(m) is contained in R(m).
Proof: if P(m) = m - 1 or m + 1, the odd primes below P(m) are exactly the odd primes <= m - 3 (m - 1 is composite when
P(m) = m + 1), whose residues are themselves; so a pair of residues summing to m is a pair of primes summing to m.
In general every odd prime < m lies below P(m) (P(m) >= m - 1), so G(m) is contained in R(m); extra pairs come
only from primes in (m, P(m)), whose residues wrap around.
This script checks both statements for all even m in [6, 20000] and counts the m with R(m) strictly larger.
"""
import sys


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return s


def main():
    N = 20000
    S = sieve(400000)
    eq_ok = sub_ok = True
    strict = []
    adj = 0
    for m in range(6, N + 1, 2):
        k, P = 1, None
        while P is None:
            for c in (k * m - 1, k * m + 1):
                if S[c]:
                    P = c
                    break
            k += 1
        res = {p % m for p in range(3, P) if S[p]}
        R = {tuple(sorted((r, m - r))) for r in res if r % 2 == 1 and 3 <= r <= m - 3 and (m - r) in res}
        G = {tuple(sorted((p, m - p))) for p in range(3, m - 2) if S[p] and S[m - p]}
        if not G <= R:
            sub_ok = False
        if S[m - 1] or S[m + 1]:
            adj += 1
            if R != G:
                eq_ok = False
        elif R != G:
            strict.append(m)
    total = (N - 6) // 2 + 1
    print(f'  m with m-1 or m+1 prime: {adj} of {total}; R(m) = G(m) for all of them: {eq_ok}')
    print(f'  G(m) contained in R(m) for every even m in [6, {N}]: {sub_ok}')
    print(f'  m (neither neighbour prime) with R(m) strictly larger than G(m): {len(strict)}; first {strict[:10]}')
    ok = eq_ok and sub_ok
    print('POST-HOC PROPOSITION CONFIRMED' if ok else 'POST-HOC CHECK FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
