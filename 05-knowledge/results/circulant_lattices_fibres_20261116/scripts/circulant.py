#!/usr/bin/env python3
"""circulant.py (H-1, H-4): circulant tournaments of prime orders p = 3, 5, 7, 11, 13 (every connection set S with
S disjoint from -S, |S| = (p-1)/2): h, P_2, P_3 by DP and subset convolution; the coefficients a_m of P_q (q^p + sum a_m
q^m, m = p (mod 2)) solved from P_1..P_{(p-1)/2}; the congruences p | P_q, a_m ≡ 0 (1 < m < p), a_1 ≡ -1 (mod p); the
orbit numbers H/p; distinct counts per order; the Paley tournaments."""
import itertools, sys
from fractions import Fraction
exec(open('nearjoin_check.py').read().split('# (1)')[0])
def ham_table(T):
    n = len(T); N = 1 << n; dp, Hm = ham_all(T); return Hm
def conv(f, g, n):
    N = 1 << n; out = [0] * N
    for W in range(N):
        s = 0; sub = W
        while True:
            s += f[sub] * g[W ^ sub]
            if sub == 0: break
            sub = (sub - 1) & W
        out[W] = s
    return out
def circ(p, S): return [[1 if (j - i) % p in S else 0 for j in range(p)] for i in range(p)]
def solve_int(rows, rhs):
    # Gaussian elimination over Fractions
    n = len(rows); M = [[Fraction(x) for x in r] + [Fraction(b)] for r, b in zip(rows, rhs)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[piv] = M[piv], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0: M[r] = [a - M[r][c] * b for a, b in zip(M[r], M[c])]
    return [M[r][n] for r in range(n)]
fails = 0; checks = 0
for p in (3, 5, 7, 11, 13):
    pairs = [(d, p - d) for d in range(1, (p + 1) // 2)]
    counts = {}; ms = list(range(p - 2, 0, -2)); k = len(ms)
    for choice in itertools.product([0, 1], repeat=len(pairs)):
        S = frozenset(pr[c] for pr, c in zip(pairs, choice)); T = circ(p, S); Hm = ham_table(T); N = 1 << p; full = N - 1
        h = Hm[full]; F = Hm[:]; Pq = [h]
        for q in range(2, max(3, k) + 1):
            F = conv(F, Hm, p); Pq.append(F[full])
        P2, P3 = Pq[1], Pq[2]
        checks += 3; fails += (h % p != 0) + (P2 % p != 0) + (P3 % p != 0)
        # coefficients: P_q - q^p = sum_m a_m q^m for q = 1..k
        rows = [[q ** m for m in ms] for q in range(1, k + 1)]; rhs = [Pq[q - 1] - q ** p for q in range(1, k + 1)]
        a = solve_int(rows, rhs) if k else []
        assert all(x.denominator == 1 for x in a); a = [int(x) for x in a]
        for m, am in zip(ms, a):
            checks += 1
            if m == 1: fails += (am + 1) % p != 0
            else: fails += am % p != 0
        counts.setdefault(h, []).append((sorted(S), P2, P3, a))
    qr = frozenset((x * x) % p for x in range(1, p)) if p % 4 == 3 else None
    print(f"p = {p}: {2 ** len(pairs)} connection sets, distinct counts {sorted(counts)}; orbit numbers H/p {sorted(set(h // p for h in counts))}; "
          f"Paley (quadratic residues) count {next((h for h, L in counts.items() for (S, P2, P3, a) in L if qr is not None and frozenset(S) == qr), None)}")
    for h in sorted(counts)[:3]:
        S, P2, P3, a = counts[h][0]; print(f"      h = {h}: S = {S}, (P_2, P_3) = ({P2}, {P3}), coefficients a_m for m = {ms}: {a}")
print(f"(H-1) congruence checks: {checks}, failures: {fails}")
