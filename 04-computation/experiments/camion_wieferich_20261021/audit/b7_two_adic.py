# Audit b7 (independent code): THM-4750 statement 6 (the prime 2).
#  (A) EXACT s(n) for n <= NEX (Python big integers, Hurwitz inverse W = T^(*-1); multiplication by T(k) = 2^C(k,2) is a shift),
#      f(n) = s(n)/n! in Z_(2):  e(n) = v_2(s(n)) - v_2(n!).
#      Checks (a) e(n) = 0 iff n != 2 (mod 3) (n >= 1, n != 2), i.e. the parity of the number of strong classes;
#             (c) e(n) = 2 + v_2(n - 3250) for every n = 2 (mod 3), 5 <= n <= NEX (|n - 3250| < 4096, so this is exact);
#             and determines n* mod 2^k from the exact data alone;
#             (b) the quasi-period f(n + 3*2^j) = f(n) mod 2^(j+2) (1 <= n <= NEX - 3*2^j), and the extra claim
#                 "mod 2^(j+3) on the class n = 1 (mod 3)".
#  (B) f(n) mod 2^K for huge n by Kitamasa (x^n mod the characteristic polynomial), a different algorithm from the session's
#      companion-matrix powers; n* determined to 2^NB bits, the law v_2(f(n)) = 2 + v_2(n - n*) tested at random huge n,
#      comparison with the session's 64-bit value, and a test of the session's own rational reconstruction
#      765490658/1598160473 at the higher precision.
# Usage: python3 b7_two_adic.py [NEX]
import sys, random, time
from math import comb
NEX = int(sys.argv[1]) if len(sys.argv) > 1 else 990
t0 = time.time()
W = [1] + [0] * NEX
binrow = [1]
for n in range(1, NEX + 1):
    binrow = [1] + [binrow[i - 1] + binrow[i] for i in range(1, n)] + [1]
    acc = 0
    for k in range(1, n + 1):
        acc += (binrow[k] * W[n - k]) << (k * (k - 1) // 2)
    W[n] = -acc
s = [0] + [-W[n] for n in range(1, NEX + 1)]
del W
print(f"exact s(n), n <= {NEX}: {time.time() - t0:.1f}s; s(1..8) = {s[1:9]}")
def v2(x): return (x & -x).bit_length() - 1 if x else None
fv = [0] * (NEX + 1)
for n in range(1, NEX + 1): fv[n] = fv[n - 1] + v2(n)
e = {n: v2(s[n]) - fv[n] for n in range(1, NEX + 1) if s[n]}
# (a)
bad_a = [n for n in e if (e[n] == 0) != (n % 3 != 2)]
print(f"(a) parity law e(n) = 0 iff n != 2 mod 3, 1 <= n <= {NEX} (n = 2 has s = 0): failures {bad_a}")
print("    e(1..30) =", [e.get(n) for n in range(1, 31)])
# number of strong isomorphism classes for small n, from the groupoid count f(n) = s(n)/n! = sum 1/|Aut| and the
# OEIS values of A051337 (n = 1..12) -> parity check only (the OEIS numbers are quoted, not recomputed)
A051337 = [1, 0, 1, 1, 6, 35, 353, 6008, 178133, 9355949, 884464590, 152310149735]
print("    A051337 parity vs (n mod 3 != 2):", all((A051337[n - 1] % 2 == 1) == (n % 3 != 2) for n in range(1, 13)))
# (c)
cls2 = [n for n in range(5, NEX + 1) if n % 3 == 2]
bad_c = [(n, e[n], 2 + v2(n - 3250)) for n in cls2 if e[n] != 2 + v2(n - 3250)]
print(f"(c) e(n) = 2 + v_2(n - 3250) for all {len(cls2)} n = 2 mod 3, 5 <= n <= {NEX}: failures {bad_c[:10]}")
# n* mod 2^k from exact data alone: largest k such that a unique residue r mod 2^k is consistent
r, k = 0, 0
while True:
    c = [x for x in (r, r + (1 << k)) if all((e[n] >= 3 + k) == ((n - x) % (1 << (k + 1)) == 0) for n in cls2 if e[n] >= 2 + k)]
    if len(c) != 1: break
    r, k = c[0], k + 1
print(f"    exact data determine n* = {r} mod 2^{k}")
# (b) quasi-period, exact
K2 = 64
MOD = 1 << K2
F = [0] + [None] * NEX
oddf = 1
for n in range(1, NEX + 1):
    i = n
    while i % 2 == 0: i //= 2
    oddf = oddf * i % MOD
    F[n] = ((s[n] >> fv[n]) * pow(oddf, -1, MOD)) % MOD if s[n] else 0
for j in range(1, 9):
    P = 3 * 2 ** j
    rng = range(1, NEX - P + 1)
    fails = [n for n in rng if (F[n + P] - F[n]) % (1 << (j + 2))]
    extra = [n for n in rng if n % 3 == 1 and (F[n + P] - F[n]) % (1 << (j + 3))]
    extra_other = [n for n in rng if n % 3 != 1 and (F[n + P] - F[n]) % (1 << (j + 3)) == 0]
    print(f"(b) j={j}: f(n+{P}) = f(n) mod 2^{j+2} fails at n = {fails[:8]} (n <= {NEX - P}); "
          f"mod 2^{j+3} on n = 1 mod 3 fails at {len(extra)} n {extra[:6]}; "
          f"(n != 1 mod 3 that happen to hold mod 2^{j+3}: {len(extra_other)})")
print(f"[part A {time.time() - t0:.1f}s]")

# ---------------- (B) huge n via Kitamasa ----------------
K = 160; M = 1 << K
def tau(n):
    ex = n * (n - 1) // 2; odd = 1
    for i in range(2, n + 1):
        while i % 2 == 0: i //= 2; ex -= 1
        odd = odd * i % M
    return (pow(odd, -1, M) << ex) % M if ex < K else 0
d = 2
while tau(d + 1) != 0: d += 1
Q = [tau(j) for j in range(d + 1)]
assert Q[0] == 1 and all(tau(j) == 0 for j in range(d + 1, d + 8))
# g = 1/Q: g_n = -sum_{j=1}^{d} Q[j] g_{n-j}
g = [1]
for n in range(1, 2 * d + 2): g.append(-sum(Q[j] * g[n - j] for j in range(1, min(d, n) + 1)) % M)
# characteristic polynomial R(x) = x^d + Q1 x^(d-1) + ... + Qd ; reduce x^n mod R
def polmulmod(a, b):
    prod = [0] * (2 * d - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b): prod[i + j] += ai * bj
    for t in range(2 * d - 2, d - 1, -1):        # x^t = -sum_j Q[j] x^(t-j)
        c = prod[t] % M
        if c:
            for j in range(1, d + 1): prod[t - j] -= c * Q[j]
        prod[t] = 0
    return [x % M for x in prod[:d]]
def g_at(n):
    if n < len(g): return g[n]
    res = [1] + [0] * (d - 1); base = [0, 1] + [0] * (d - 2); e_ = n
    while e_:
        if e_ & 1: res = polmulmod(res, base)
        base = polmulmod(base, base); e_ >>= 1
    return sum(res[i] * g[i] for i in range(d)) % M
f = lambda n: (-g_at(n)) % M
# sanity vs the exact values
assert all(f(n) % (1 << 64) == F[n] for n in range(1, 200)), "Kitamasa disagrees with exact s(n)/n!"
print(f"(B) K = {K} bits, recurrence order d = {d}; Kitamasa values agree with exact f(n) for n < 200 (mod 2^64)")
def v2m(x): x %= M; return v2(x) if x else K
rnd = random.Random(20261021)
NB = 140
rr, kk = 2, 0
while kk < NB:
    hits = []
    for cand in (rr, rr + (1 << kk)):
        ok = True
        for trial in range(2):
            # n = cand mod 2^(kk+1), n = 2 mod 3, n huge and random
            base_n = cand + (1 << (kk + 1)) * rnd.randrange(1 << 40)
            while base_n % 3 != 2: base_n += 1 << (kk + 1)
            if v2m(f(base_n)) < 3 + kk: ok = False; break
        if ok: hits.append(cand)
    if len(hits) != 1: print("    ambiguous at k =", kk, hits); break
    rr = hits[0]; kk += 1
nstar = rr
print(f"    n* mod 2^{kk} = {nstar}")
print(f"    n* mod 4096 = {nstar % 4096};  n* mod 2^64 = {nstar % (1 << 64)}  (session: 16221783116972227762) "
      f"{'AGREE' if nstar % (1 << 64) == 16221783116972227762 else 'DISAGREE'}")
# the law at random huge n (not of the special form used to find the bits)
viol = 0; tested = 0
for trial in range(300):
    n = rnd.randrange(1 << 60, 1 << 100)
    while n % 3 != 2: n += 1
    if rnd.random() < 0.5:                       # force deep agreement with n* for half of the tests
        j = rnd.randrange(1, 120); n = nstar % (1 << j) + (1 << j) * rnd.randrange(1 << 50) + (1 << j) * 0
        while n % 3 != 2: n += 1 << j
    vt = v2((n - nstar) % (1 << kk)) if (n - nstar) % (1 << kk) else kk
    if vt + 2 >= K - 2: continue
    tested += 1
    if v2m(f(n)) != 2 + vt: viol += 1
print(f"    law v_2(f(n)) = 2 + v_2(n - n*) at {tested} random huge n = 2 mod 3: violations {viol}")
# other classes: v_2(f(n)) = 0 for n != 2 mod 3 (parity law), at random huge n
v0 = sum(1 for t in range(100) for n in [rnd.randrange(1 << 80) * 3 + rnd.choice((0, 1))] if v2m(f(n)) != 0)
print(f"    v_2(f(n)) = 0 at 100 random huge n = 0,1 mod 3: violations {v0}")
# the session's rational reconstruction candidate
a, b = 765490658, 1598160473
for bits in (64, 80, 100, NB):
    m = 1 << bits
    print(f"    candidate {a}/{b} = n* mod 2^{bits}? {(a * pow(b, -1, m) - nstar) % m == 0}")
def ratrec(r, m):
    from math import isqrt
    a0, a1, b0, b1 = m, r % m, 0, 1; bound = isqrt(m // 2)
    while a1 > bound:
        q = a0 // a1; a0, a1 = a1, a0 - q * a1; b0, b1 = b1, b0 - q * b1
    return a1, b1
for bits in (64, 96, 128, NB):
    m = 1 << bits; A_, B_ = ratrec(nstar, m)
    good = B_ % 2 == 1 and abs(B_) <= (m // 2) ** 0.5 and (A_ * pow(B_, -1, m) - nstar) % m == 0
    print(f"    rational reconstruction at 2^{bits}: {A_}/{B_} (|num| ~ 2^{A_.bit_length()}, |den| ~ 2^{abs(B_).bit_length()}), "
          f"valid {good}")
A128, B128 = ratrec(nstar, 1 << 128)
print(f"    the 2^128 reconstruction {A128}/{B128} matches n* mod 2^{NB}? {(A128 * pow(B128, -1, 1 << NB) - nstar) % (1 << NB) == 0}")
print(f"[total {time.time() - t0:.1f}s]")
