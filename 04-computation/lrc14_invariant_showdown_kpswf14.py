#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) INVARIANT SHOWDOWN.

Side-by-side comparison of ALL candidate tournaments on the discriminating row set, with a
FULL iso-fingerprint (score-seq, c3, c5, H=#Ham-paths). Goal: find the cleanest magnitude-aware
tournament whose iso class separates TIGHT {AP,GW} from the apex-twins / loose rows.

The apex twins (AP, 12->26, 12->96) are the acid test: identical under any periodic winding.

Definitions of each tournament are printed. Exact Fraction arithmetic. n=13.
c5 = number of directed 5-cycles (full enumeration over 5-subsets, counts cyclic 5-perms).
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
from functools import reduce

# ---------- iso fingerprint ----------
def H_paths(A):
    n = len(A); size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c == 0: continue
            Av = A[v]
            for w in range(n):
                if mask & (1 << w): continue
                if Av[w]: dp[mask | (1 << w)][w] += c
    return sum(dp[size - 1])

def scores(A):
    n = len(A); return tuple(sorted(sum(A[i]) for i in range(n)))

def c3(A):
    n = len(A); cnt = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]: cnt += 1
        if A[a][c] and A[c][b] and A[b][a]: cnt += 1
    return cnt

def c5(A):
    """number of directed 5-cycles. For each 5-subset, count cyclic orderings that form a
    directed cycle. There are (5-1)!/2 = 12 undirected cyclic perms; each may be a directed
    5-cycle. We count directed Hamiltonian cycles within the induced 5-tournament."""
    n = len(A); cnt = 0
    for sub in combinations(range(n), 5):
        s = sub
        # fix s[0], count directed cycles through all 5
        for perm in permutations(s[1:]):
            cyc = (s[0],) + perm
            if all(A[cyc[k]][cyc[(k+1) % 5]] for k in range(5)):
                cnt += 1
    return cnt // 1  # each directed cycle counted once (start fixed at s[0], but direction free)

def is_tournament(A):
    n = len(A)
    return all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(i+1, n))

def fp(A):
    return dict(scores=scores(A), c3=c3(A), c5=c5(A), H=H_paths(A), tourn=is_tournament(A))

# ---------- candidate tournaments ----------
def winding(S, t):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            rel = (F(S[i]-S[j]) * t) % 1
            if 0 < rel < F(1, 2): A[i][j] = 1
            elif rel > F(1, 2): A[i][j] = 0
            else: A[i][j] = 1 if S[i] < S[j] else 0
    return A

def divis_lt_half(S):
    k = len(S); raw = [[2*(S[i] % S[j]) < S[j] for j in range(k)] for i in range(k)]
    A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            a, b = raw[i][j], raw[j][i]
            if a and not b: A[i][j] = 1
            elif b and not a: A[j][i] = 1
            else: A[i][j] = 1 if S[i] < S[j] else 0; A[j][i] = 1 - A[i][j]
    return A

def floor_odd(S):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            ai = (S[i] // S[j]) % 2 == 1; aj = (S[j] // S[i]) % 2 == 1
            if ai and not aj: A[i][j] = 1
            elif aj and not ai: A[j][i] = 1
            else: A[i][j] = 1 if S[i] < S[j] else 0; A[j][i] = 1 - A[i][j]
    return A

def reciprocal_stack(S):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            votes = 0; tot = 0; di = S[i] - S[j]
            for s in S:
                rel = (F(di, s)) % 1
                if rel == 0 or rel == F(1, 2): continue
                tot += 1
                if rel < F(1, 2): votes += 1
            if 2*votes > tot: A[i][j] = 1
            elif 2*votes < tot: A[j][i] = 1
            else: A[i][j] = 1 if S[i] < S[j] else 0; A[j][i] = 1 - A[i][j]
    return A

def cf_parity(S):
    def cl(a, b):
        n = 0
        while b: a, b = b, a % b; n += 1
        return n
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            L = cl(S[i], S[j])
            small_beats = (L % 2 == 1)
            if small_beats: A[i][j] = 1 if S[i] < S[j] else 0
            else: A[i][j] = 1 if S[i] > S[j] else 0
            A[j][i] = 1 - A[i][j]
    return A

CANDS = [
    ("winding T(1/14) apex [periodic]", lambda S: winding(S, F(1, 14))),
    ("winding T(3/41) Farey [periodic]", lambda S: winding(S, F(3, 41))),
    ("divis  s_i%s_j<s_j/2  [magnitude]", divis_lt_half),
    ("floor(s_i/s_j) odd    [magnitude]", floor_odd),
    ("reciprocal-stack 1/s  [magnitude]", reciprocal_stack),
    ("CF-depth parity       [magnitude]", cf_parity),
]

ROWS = [
    ("AP        TIGHT", list(range(1, 14))),
    ("GW        TIGHT", list(range(1, 12)) + [13, 24]),
    ("12->26    loose", list(range(1, 12)) + [13, 26]),
    ("12->36    loose", list(range(1, 12)) + [13, 36]),
    ("12->96    loose", list(range(1, 12)) + [13, 96]),
]
TIGHTNAMES = {"AP        TIGHT", "GW        TIGHT"}

for cname, fn in CANDS:
    print("="*92)
    print(cname)
    print("-"*92)
    tight_H = set(); loose_H = set()
    tight_fp = set(); loose_fp = set()
    for rname, S in ROWS:
        f = fn(S); d = fp(f)
        sig = (d['scores'], d['c3'], d['c5'], d['H'])
        print(f"  {rname}: c3={d['c3']:4d} c5={d['c5']:5d} H={d['H']:>10d}  tourn={d['tourn']}")
        if rname in TIGHTNAMES:
            tight_H.add(d['H']); tight_fp.add(sig)
        else:
            loose_H.add(d['H']); loose_fp.add(sig)
    print(f"  --> H separates 5-row tight/loose: {tight_H.isdisjoint(loose_H)}   "
          f"full-fp separates: {tight_fp.isdisjoint(loose_fp)}")
    # the acid test: do AP, 12->26, 12->96 get DISTINCT iso classes?
    twins = {}
    for rname, S in ROWS:
        if rname.split()[0] in ("AP", "12->26", "12->96"):
            d = fp(fn(S)); twins[rname.split()[0]] = (d['scores'], d['c3'], d['c5'], d['H'])
    distinct = len(set(twins.values())) == len(twins)
    print(f"  --> apex-twins AP/12->26/12->96 get DISTINCT iso classes: {distinct}")
    print()
