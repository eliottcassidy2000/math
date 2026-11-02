#!/usr/bin/env python3
"""prime_path_mod36.py (klein-2026-11-01) -- the mod-36 partner law behind the 2-matching obstruction.
If p = 1 (mod 12) and q = 11 (mod 12) are adjacent then p + q = 0 (mod 12), so the square p + q is divisible by 36 and
q = -p (mod 36): the three classes 1, 13, 25 (mod 36) of primes = 1 (mod 12) are adjacent only to the classes 35, 23, 11
(mod 36) respectively (and to the prime 3); likewise the classes 7, 19, 31 (mod 36) of primes = 7 (mod 12) only to 29,
17, 5 (mod 36) (and to the prime 2).  Each restricted class a therefore needs, in a Hamiltonian path, 2 c_a - e_a path
edges inside the capacity 2 c_a' + 2 - e' of its partner class (the +2 for the shared prime 3 or 2), so c_a <= c_a' + 2.
Sharper (segment condition): the path edges inside X_a = class a u class a' u {special prime} form a forest with
at least 2 c_a - e_a edges on c_a + c_a' + 1 vertices, so c_a - c_a' <= e_a, where e_a is the number of global
endpoints in class a; with |B| = |A| + 1 both endpoints lie in B (e_a = 0 for the A-classes 1, 13, 25), with |A| = |B|
one endpoint lies on each side, with |A| = |B| + 1 both in A.  This script (i) verifies the partner law on all edges
below 10^5, (ii) prints the six pair counts at N = 615,731, (iii) scans every prime bound N < BOUND for the mod-4
balance and the six weak conditions c_a <= c_a' + 2, and (iv) scans the segment conditions, reporting the first
passing bound and the status of every bound in [615,731, 617,059].  usage: python3 prime_path_mod36.py [BOUND=5000000]"""
import math, sys, bisect
BOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 5_000_000
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
PAIRS = {1: 35, 13: 23, 25: 11, 7: 29, 19: 17, 31: 5}
# (i) brute force below 10^5
P = primes_below(100000); Ps = set(P); top = P[-1]; viol = 0; edges = 0; by36 = {}
for p in P:
    s = math.isqrt(p + 2)
    while s * s <= p + top:
        q = s * s - p
        if q > p and q in Ps:
            edges += 1
            for a, b in ((p, q), (q, p)):
                if a > 3 and a % 36 in PAIRS and b > 3 and b % 36 != PAIRS[a % 36]: viol += 1
                if a > 3 and a % 36 in PAIRS and b > 3: by36[(a % 36, b % 36)] = by36.get((a % 36, b % 36), 0) + 1
        s += 1
print('(i) edges below 10^5: %d; violations of the mod-36 partner law: %d; observed partner pairs: %s' % (edges, viol, sorted(by36)))
assert viol == 0
def counts(Pn):
    c = {r: 0 for r in range(36)}
    for p in Pn:
        if p > 3: c[p % 36] += 1
    return c
P = primes_below(BOUND)
k = bisect.bisect_left(P, 615731); c = counts(P[:k])
print('(ii) N = 615731: pairs (c_a, c_partner, c_a - c_partner):', {a: (c[a], c[b], c[a] - c[b]) for a, b in PAIRS.items()})
# (iii) scan
c = {r: 0 for r in range(36)}; passing = []
for i, p in enumerate(P[:-1]):
    if p > 3: c[p % 36] += 1
    N = P[i + 1]
    if N <= 100: continue
    A = 1 + sum(c[r] for r in range(36) if r % 4 == 1); B = 1 + sum(c[r] for r in range(36) if r % 4 == 3)
    if abs(A - B) > 1: continue
    if all(c[a] <= c[b] + 2 for a, b in PAIRS.items()): passing.append(N)
print('(iii) bounds N in (100, %d) passing mod-4 balance and all six mod-36 conditions: %d' % (BOUND, len(passing)))
if passing:
    win = []; start = passing[0]; prev = passing[0]
    for N in passing[1:]:
        if N - prev > 20000: win.append((start, prev)); start = N
        prev = N
    win.append((start, prev)); print('      windows:', win[:20])
    # the worst pair inside the second window
    for N in passing[:1] + passing[-1:]:
        kk = bisect.bisect_left(P, N); cc = counts(P[:kk])
        print('      N = %d: pair margins c_partner + 2 - c_a =' % N, {a: cc[b] + 2 - cc[a] for a, b in PAIRS.items()})

# (iv) segment conditions: c_a - c_a' <= e_a, endpoints distributed by the side balance; the surplus pieces on one side
#      cannot exceed the endpoints available on that side
def seg_ok(c):
    A = 1 + sum(c[r] for r in range(36) if r % 4 == 1); B = 1 + sum(c[r] for r in range(36) if r % 4 == 3)
    d = A - B
    if abs(d) > 1: return False, d
    endsA = {-1: 0, 0: 1, 1: 2}[d]; endsB = 2 - endsA
    surplusA = [c[a] - c[b] for a, b in PAIRS.items() if a % 4 == 1]
    surplusB = [c[a] - c[b] for a, b in PAIRS.items() if a % 4 == 3]
    needA = sum(max(0, x) for x in surplusA); needB = sum(max(0, x) for x in surplusB)
    return (max(surplusA) <= endsA and max(surplusB) <= endsB and needA <= endsA and needB <= endsB), d
c = {r: 0 for r in range(36)}; first = None; passing4 = []; status = []
for i, p in enumerate(P[:-1]):
    if p > 3: c[p % 36] += 1
    N = P[i + 1]
    if N <= 100: continue
    ok, d = seg_ok(c)
    if ok:
        passing4.append(N)
        if first is None: first = (N, i + 1, d, {a: c[a] - c[b] for a, b in PAIRS.items()})
    if 615731 <= N <= 617059 and (ok or N in (615731, 617053, 617059)):
        status.append((N, i + 1, 'PASS' if ok else 'fail', d, c[25] - c[11]))
print('(iv) segment conditions: first passing bound N = %d (%d primes, |A|-|B| = %+d, surpluses c_a - c_partner %s)' % first)
print('      bounds in [615731, 617059] passing, and the three landmarks (N, primes, status, |A|-|B|, c_25 - c_11):', status)
print('      passing bounds below %d: %d; windows:' % (BOUND, len(passing4)), end=' ')
win = []; start = passing4[0]; prev = passing4[0]
for N in passing4[1:]:
    if N - prev > 20000: win.append((start, prev)); start = N
    prev = N
win.append((start, prev)); print(win[:20])
