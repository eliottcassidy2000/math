#!/usr/bin/env python3
"""Wave 38: the Z[sqrt -2] / mod-8 bridges between cyclic joins, Paley tournaments and the Collatz tables.

(1) Diagonal joins mod 16 / 32 (refining N = 3 mod 8): for a block A of order n with H = h, Q = Pi_2, beta = (Q - h) + h sqrt(-2),
    N(beta) = H(cyc(A, A', *)) satisfies
        odd n:  N = 3 h^2 + 8 HC(A)        (mod 16),   HC = number of directed Hamiltonian cycles,
        even n: N = 3 h^2 + 16 c_{n-1}(A)  (mod 32),   c_{n-1} = number of directed (n-1)-cycles,
    from Q = 4 HC (mod 8) at odd n and Q = 8 c_{n-1} (mod 16) at even n (chanpak_pair_line_census PAIR_LINE_STRUCTURE s.4).
(2) Paley primes p = 3 mod 4, p < LIMIT: equivalence of  p = 3 (mod 8),  (2/p) = -1 (halving reverses the Paley tournament),
    -1 in <2> mod p (the 3x+1 / 3x-1 tables alternate, THM-4792), and p = x^2 + 2y^2.
(3) Hasse's density 17/24 of primes with ord_p(2) even (= alternating primes of THM-4792), counted to LIMIT.
(4) 2-adic: the subgroup <3> mod 2^j equals the unit norms {x^2 + 2y^2 : x odd} mod 2^j = {u = 1, 3 mod 8}, j <= 16.
(5) Paley tournaments P_3, P_7, P_11, P_19: H, Q, beta, the mod-8 data, and |Aut| = p(p-1)/2 dividing H and Q.
(6) Diagonal-join primes: every prime value of (h + P2)^2 + 2 h^2 over blocks of order <= 10; the census is COMPLETE below
    (2^11 - 1)^2 + 2 = 4,190,211, since a block of order >= 11 has P2 >= 2^11 - 2.
Usage: python3 w38_norm_bridges.py HQ_BINARY PAIRS_DIR
"""
import random
import subprocess
import sys
from math import isqrt

HQ, PD = sys.argv[1], sys.argv[2]
rng = random.Random(9792)


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]], s


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def hq_batch(tours):
    inp = "".join(f"{len(t)} " + " ".join(map(str, t)) + "\n" for t in tours)
    out = subprocess.run([HQ], input=inp.encode(), capture_output=True, check=True).stdout.decode().split("\n")
    return [tuple(map(int, l.split())) for l in out if l.strip()]


def cycles_through_all(out, verts):
    """number of directed Hamiltonian cycles of the subtournament on verts (list)."""
    k = len(verts)
    if k < 3:
        return 0
    idx = {v: i for i, v in enumerate(verts)}
    adj = [[(out[verts[i]] >> verts[j]) & 1 for j in range(k)] for i in range(k)]
    N = 1 << k
    dp = [[0] * k for _ in range(N)]
    dp[1][0] = 1
    for S in range(1, N):
        if not S & 1:
            continue
        for v in range(k):
            c = dp[S][v]
            if not c:
                continue
            for u in range(k):
                if not (S >> u) & 1 and adj[v][u]:
                    dp[S | 1 << u][u] += c
    return sum(dp[N - 1][v] for v in range(1, k) if adj[v][0])


# (1) diagonal-join congruences
fails = 0
tests = 0
tours = []
for _ in range(300):
    n = rng.randint(3, 10)
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    tours.append(out)
for out, (h, q) in zip(tours, hq_batch(tours)):
    n = len(out)
    Nb = (q - h) ** 2 + 2 * h * h
    if n % 2:
        HC = cycles_through_all(out, list(range(n)))
        ok = (q - 4 * HC) % 8 == 0 and (Nb - 3 * h * h - 8 * HC) % 16 == 0
    else:
        c = sum(cycles_through_all(out, [v for v in range(n) if v != w]) for w in range(n))
        ok = (q - 8 * c) % 16 == 0 and (Nb - 3 * h * h - 16 * c) % 32 == 0
    tests += 1
    fails += not ok
print(f"(1) diagonal-join law N = 3h^2 + 8 HC (mod 16) [odd n], 3h^2 + 16 c_(n-1) (mod 32) [even n]: "
      f"{tests} random tournaments of orders 3-10, {fails} failures")

LIMIT = 10 ** 7
P, sieve = primes_upto(LIMIT)

# (2) Paley-prime equivalence
bad = 0
cnt = 0
for p in P:
    if p % 4 != 3 or p < 7:
        continue
    cnt += 1
    a = p % 8 == 3
    b = pow(2, (p - 1) // 2, p) == p - 1          # (2/p) = -1
    m = p - 1
    while m % 2 == 0:
        m //= 2
    c = pow(2, m, p) != 1                          # ord_p(2) even  <=>  -1 in <2>
    if p < 10 ** 6:
        d = any(isqrt(p - 2 * y * y) ** 2 == p - 2 * y * y for y in range(1, isqrt(p // 2) + 1))
    else:
        d = a
    if not (a == b == c == d):
        bad += 1
print(f"(2) Paley primes 7 <= p < 10^7 ({cnt} primes): p = 3 mod 8 <=> (2/p) = -1 <=> -1 in <2> <=> p = x^2 + 2y^2 "
      f"(the last tested below 10^6): {bad} exceptions")

# (3) Hasse's 17/24
ev = 0
for p in P[1:]:
    m = p - 1
    while m % 2 == 0:
        m //= 2
    ev += pow(2, m, p) != 1
print(f"(3) odd primes < 10^7 with ord_p(2) even: {ev} of {len(P) - 1} = {ev / (len(P) - 1):.5f} (Hasse: 17/24 = {17 / 24:.5f})")

# (4) 2-adic norm group
ok = True
for j in range(3, 17):
    M = 1 << j
    g3 = set()
    x = 1
    for _ in range(M):
        g3.add(x)
        x = 3 * x % M
    norms = {(x * x + 2 * y * y) % M for x in range(1, M, 2) for y in range(0, M)} if j <= 10 else None
    target = {u for u in range(1, M, 2) if u % 8 in (1, 3)}
    if g3 != target or (norms is not None and norms != target):
        ok = False
print(f"(4) <3> mod 2^j = unit norms x^2 + 2y^2 mod 2^j = {{u = 1, 3 mod 8}} for 3 <= j <= 16 (norms for j <= 10): {ok}")

# (5) Paley tournaments
rows = []
tours = []
for p in (3, 7, 11, 19):
    QR = {(x * x) % p for x in range(1, p)}
    tours.append([sum(1 << j for j in range(p) if (j - i) % p in QR) for i in range(p)])
for p, (h, q) in zip((3, 7, 11, 19), hq_batch(tours)):
    aut = p * (p - 1) // 2
    x, y = q - h, h
    Nb = x * x + 2 * y * y
    rows.append(f"  P_{p}: H = {h}, Q = {q}, H/|Aut| = {h // aut}{'' if h % aut == 0 else ' (!)'}, Q/|Aut| = {q / aut:g}, "
                f"beta = {x} + {y} sqrt(-2), N(beta)/|Aut|^2 = {Nb / aut ** 2:g}, H mod 8 = {h % 8}, Q mod 8 = {q % 8}")
print("(5) Paley tournaments (|Aut| = p(p-1)/2 acts freely on paths and on path pairs):")
print("\n".join(rows))

# (6) diagonal-join primes
BOUND = (2 ** 11 - 1) ** 2 + 2
vals = set()
small = [(1, 0), (1, 2), (1, 6), (3, 6)]  # orders 1, 2, 3
for h, p2 in small:
    vals.add((h + p2) ** 2 + 2 * h * h)
for n in range(4, 11):
    for line in open(f"{PD}/pairs{n}.txt"):
        h, p2 = map(int, line.split())
        vals.add((h + p2) ** 2 + 2 * h * h)
assert all(v % 8 == 3 for v in vals)
pv = sorted(v for v in vals if is_prime(v))
target = [p for p in P if p % 8 == 3 and p < BOUND]
hit = [p for p in pv if p < BOUND]
print(f"(6) diagonal-join values over blocks of order <= 10: {len(vals)} distinct, all = 3 mod 8; {len(pv)} prime")
print(f"    complete range p < {BOUND}: {len(hit)} of the {len(target)} primes = 3 mod 8 are diagonal joins "
      f"({len(hit) / len(target):.4f}); the first realised: {hit[:12]}")
miss = [p for p in target if p not in set(hit)]
print(f"    first primes = 3 mod 8 that are NOT diagonal joins: {miss[:15]}")
for lo, hi in ((0, 10 ** 4), (10 ** 4, 10 ** 5), (10 ** 5, 10 ** 6), (10 ** 6, BOUND)):
    t = sum(1 for p in target if lo <= p < hi)
    hh = sum(1 for p in hit if lo <= p < hi)
    print(f"    [{lo}, {hi}): {hh} / {t} = {hh / max(t, 1):.4f}")
