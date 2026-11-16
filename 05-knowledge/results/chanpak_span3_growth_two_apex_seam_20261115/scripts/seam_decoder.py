#!/usr/bin/env python3
"""seam_decoder.py (klein-2026-11-15): the signed-seam product decoder of chanpak_signed_seam_bank_20261116/PRODUCT.md,
run in reverse.  Family: tape TT_k -> TT_t -> TT_r (ordinal sum of transitive blocks), sentinel s with tape -> s, head x
with s -> x, x beats every tape vertex except the marked ones V -> x.  Left marks: bit i of b (0 <= b < 2^(k-1)) marks
left position k-1-i (backward; position 0 never marked).  Right marks: bit p of a (0 <= a < 2^(r-1)) marks right
position p (forward; the last position r-1 never marked).  Claimed exact count (PRODUCT.md (*)):
    H = P * 2^t + C,  P = (2^k - b)(2^r + a),  C = 1 + rev_k(b) - rev_r(a) - J(b) - J(a),
    J(z) = sum over set bits i > j of z of (2^(i-j-1) + [i == j+1]).
Modes:
  verify kmax rmax tmax   : brute-force Hamiltonian-path counts of the literal tournaments against (*) for all masks
  decode N tmin tmax kmax rmax : every (k, r, t, a, b) with H = N
  windows kmax rmax       : for each (k, r) and each product P, the set of remainders C over admissible factor pairs,
                            and the longest run of consecutive odd values P*2^t + C (independent of t)."""
import sys
from collections import defaultdict
def rev(z, w): return int(format(z, f'0{w}b')[::-1], 2) if w > 0 else 0
def J(z):
    bits = [i for i in range(z.bit_length()) if z >> i & 1]; s = 0
    for x in range(len(bits)):
        for y in range(x):
            i, j = bits[x], bits[y]; s += 2 ** (i - j - 1) + (1 if i == j + 1 else 0)
    return s
def C_of(k, r, a, b): return 1 + rev(b, k) - rev(a, r) - J(b) - J(a)
def P_of(k, r, a, b): return (2 ** k - b) * (2 ** r + a)
def build(k, t, r, a, b):
    """adjacency (out-sets) of the literal tournament; vertices: tape 0..k+t+r-1, sentinel s, head x"""
    n = k + t + r; S = n; X = n + 1; N = n + 2
    out = [set() for _ in range(N)]
    for i in range(n):
        for j in range(i + 1, n): out[i].add(j)          # transitive tape (ordinal sum of transitive blocks)
        out[i].add(S)                                      # tape -> sentinel
    out[S].add(X)                                          # sentinel -> head
    marked = set()
    for i in range(k - 1):
        if b >> i & 1: marked.add(k - 1 - i)               # left, backward
    for p in range(r - 1):
        if a >> p & 1: marked.add(k + t + p)               # right, forward
    for u in range(n):
        if u in marked: out[u].add(X)
        else: out[X].add(u)
    return out, N
def ham_paths(out, n):
    dp = [dict() for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        if not row: continue
        for v, c in row.items():
            for w in out[v]:
                if not (S >> w) & 1:
                    T = S | (1 << w); dp[T][w] = dp[T].get(w, 0) + c
    return sum(dp[(1 << n) - 1].values())
def admissible(k, r):
    """(a, b, d, e) with d = 2^k - b in (2^(k-1), 2^k], e = 2^r + a in [2^r, 3*2^(r-1))"""
    for b in range(2 ** (k - 1)):
        for a in range(2 ** (r - 1)):
            yield a, b, 2 ** k - b, 2 ** r + a
mode = sys.argv[1]
if mode == 'verify':
    kmax, rmax, tmax = map(int, sys.argv[2:5]); bad = 0; tot = 0
    for k in range(1, kmax + 1):
        for r in range(1, rmax + 1):
            for t in range(1, tmax + 1):
                if k + t + r + 2 > 16: continue
                for a, b, d, e in admissible(k, r):
                    out, N = build(k, t, r, a, b); H = ham_paths(out, N); pred = P_of(k, r, a, b) * 2 ** t + C_of(k, r, a, b)
                    tot += 1
                    if H != pred: bad += 1; print("MISMATCH", k, t, r, a, b, H, pred)
    print(f"verified {tot} literal tournaments; mismatches {bad}")
    out, N = build(5, 1, 4, 5, 13); print("control k=5,r=4,t=1 (a,b)=(5,13): H =", ham_paths(out, N), "formula", P_of(5, 4, 5, 13) * 2 + C_of(5, 4, 5, 13), "(packet: 801)")
    out, N = build(5, 1, 4, 3, 11); print("control (3,11): H =", ham_paths(out, N), "formula", P_of(5, 4, 3, 11) * 2 + C_of(5, 4, 3, 11), "(packet: 803)")
elif mode == 'decode':
    N = int(sys.argv[2]); tmin, tmax, kmax, rmax = map(int, sys.argv[3:7]); hits = []
    for k in range(2, kmax + 1):
        for r in range(2, rmax + 1):
            Cmax = 2 ** (k + 1) + 2 ** (r + 1)
            for t in range(tmin, tmax + 1):
                lo = (N - Cmax) >> t; hi = (N + Cmax) >> t
                for a, b, d, e in admissible(k, r):
                    P = d * e
                    if P < lo or P > hi: continue
                    if P * 2 ** t + C_of(k, r, a, b) == N: hits.append((k, r, t, a, b, d, e, C_of(k, r, a, b)))
    print(f"N = {N}: {len(hits)} realisations (k, r, t, a, b, d, e, C):")
    for h in hits[:40]: print("  ", h)
elif mode == 'windows':
    kmax, rmax = map(int, sys.argv[2:4]); best = (0, None)
    for k in range(2, kmax + 1):
        for r in range(2, rmax + 1):
            byP = defaultdict(set)
            for a, b, d, e in admissible(k, r): byP[d * e].add(C_of(k, r, a, b))
            for P, Cs in byP.items():
                cs = sorted(Cs); run = 1; bestrun = 1; start = cs[0]; bs = cs[0]
                for i in range(1, len(cs)):
                    if cs[i] == cs[i - 1] + 2: run += 1
                    else: run = 1; start = cs[i]
                    if run > bestrun: bestrun = run; bs = start
                if bestrun > best[0]: best = (bestrun, (k, r, P, bs, len(Cs)))
            longest = max(len(v) for v in byP.values())
            print(f"k={k} r={r}: products {len(byP)}, max factor pairs per product {longest}, best consecutive-odd run {best[0]} at {best[1]}")
    print("overall best run:", best)
