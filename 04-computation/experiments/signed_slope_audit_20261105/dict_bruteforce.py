#!/usr/bin/env python3
"""Independent brute-force check of THM-4882 (1), the dictionary.

T(b): apex x = vertex 0; chain v_1..v_m = vertices 1..m with v_i -> v_j for i < j;
x -> v_i iff b_i = 1, else v_i -> x.  Profile in Session 198's conventions:
  D = H(T), E = #Hamiltonian paths ending at x, U = #ordered covers (P1, x P2) with P1 possibly empty,
  d = H(T - x), s = #Hamiltonian paths starting at x.
Claims: (D,E,U,d) = (H, 1 - b_m, 1 + bin(b), 1); H = g(1 + bin(b)) with Session 198's g;
appending bit 1 acts as O(D,E,U,d) = (D,0,2U,d); appending bit 0 acts as Z = (D-E+U, d, 2U-d, d).
All words of length <= LMAX are enumerated (both leading bits).  Path counts by subset DP, cross-checked
against permutation enumeration for order <= 7.
"""
import itertools, sys
from functools import lru_cache

LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7

def build(bits):
    m = len(bits); n = m + 1
    adj = [[False]*n for _ in range(n)]
    for i in range(1, n):
        for j in range(i+1, n):
            adj[i][j] = True
    for i, b in enumerate(bits, start=1):
        if b: adj[0][i] = True
        else: adj[i][0] = True
    return adj

def ham_counts(adj, verts):
    """Return dict: (end[v]) and (start[v]) counts of Hamiltonian paths of the induced subtournament on verts, plus total."""
    verts = list(verts); n = len(verts); idx = {v:i for i,v in enumerate(verts)}
    if n == 0: return {'H': 1, 'end': {}, 'start': {}}
    full = (1 << n) - 1
    end = [[0]*n for _ in range(1 << n)]
    for i in range(n): end[1 << i][i] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S >> v) & 1 or end[S][v] == 0: continue
            for w in range(n):
                if (S >> w) & 1: continue
                if adj[verts[v]][verts[w]]:
                    end[S | (1 << w)][w] += end[S][v]
    # starts: reverse
    start = [[0]*n for _ in range(1 << n)]
    for i in range(n): start[1 << i][i] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S >> v) & 1 or start[S][v] == 0: continue
            for w in range(n):
                if (S >> w) & 1: continue
                if adj[verts[w]][verts[v]]:
                    start[S | (1 << w)][w] += start[S][v]
    H = sum(end[full])
    return {'H': H, 'end': {verts[v]: end[full][v] for v in range(n)}, 'start': {verts[v]: start[full][v] for v in range(n)}}

def ham_perm(adj, verts):
    verts = list(verts)
    return sum(all(adj[p[i]][p[i+1]] for i in range(len(p)-1)) for p in itertools.permutations(verts))

def profile(bits):
    adj = build(bits); n = len(bits) + 1; V = list(range(n))
    hc = ham_counts(adj, V)
    D = hc['H']; E = hc['end'][0]; s = hc['start'][0]
    d = ham_counts(adj, V[1:])['H']
    # U: ordered covers (P1, x P2), P1 a Hamiltonian path of T[X], X subset of V - x (possibly empty), P2 starts at x on V - X
    U = 0
    others = V[1:]
    for r in range(len(others) + 1):
        for X in itertools.combinations(others, r):
            h1 = ham_counts(adj, X)['H'] if X else 1
            rest = [0] + [v for v in others if v not in X]
            h2 = ham_counts(adj, rest)['start'][0]
            U += h1 * h2
    return D, E, U, d, s

@lru_cache(None)
def g(n):
    if n == 1: return 1
    if n % 2 == 0: return g(n // 2)
    m = (n + 1) // 2
    return g(m) + m - (m % 2)

def binval(bits):
    v = 0
    for b in bits: v = 2 * v + b
    return v

def O(p): D,E,U,d = p; return (D, 0, 2*U, d)
def Z(p): D,E,U,d = p; return (D - E + U, d, 2*U - d, d)

bad = 0; total = 0; permchecked = 0
prof = {}
for m in range(1, LMAX + 1):
    for bits in itertools.product([0, 1], repeat=m):
        D, E, U, d, s = profile(bits)
        prof[bits] = (D, E, U, d, s)
        total += 1
        v = binval(bits)
        ok = (E == 1 - bits[-1]) and (U == 1 + v) and (d == 1) and (D == g(1 + v)) and (s == bits[0]) and (D % 2 == 1)
        if m + 1 <= 7:
            hp = ham_perm(build(bits), range(m + 1)); permchecked += 1
            ok = ok and (hp == D)
        if not ok:
            bad += 1; print("DICTIONARY FAIL", bits, (D, E, U, d, s), "expected", (g(1+v), 1-bits[-1], 1+v, 1), "s expected", bits[0])
# transitions
tbad = 0; tt = 0
for bits, p in prof.items():
    if len(bits) >= LMAX: continue
    D, E, U, d, s = p
    q1 = prof[bits + (1,)]; q0 = prof[bits + (0,)]
    tt += 2
    if q1[:4] != O((D, E, U, d)): tbad += 1; print("O FAIL", bits, p, q1)
    if q0[:4] != Z((D, E, U, d)): tbad += 1; print("Z FAIL", bits, p, q0)
print(f"words of length 1..{LMAX}: {total} checked ({permchecked} also by permutation enumeration); dictionary failures: {bad}")
print(f"transitions checked: {tt}; O/Z failures: {tbad}")
# seed
print("seed word (1):", prof[(1,)], " Session 198 seed t = (1,0,2,1)")
# the orbit statement: for odd n, point (g(n), 1, n, 1) is the word of value n-1 (leading bit 1), ends with Z (bit 0)
