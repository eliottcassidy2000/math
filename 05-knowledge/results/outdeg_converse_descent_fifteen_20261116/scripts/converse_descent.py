#!/usr/bin/env python3
"""converse_descent.py: (1) replay of the digit-0 characterisation on all tournaments of orders 3-8 (every vertex);
(2) for strong D of orders 3-7: h(D), P_D(2) and whether P_D(2)/h(D) is a power of two; (3) for strong S of orders 4-8 with a
digit-0 vertex z: the terminal component D of S - z (size 1 = near-join proper), and the identity H(S) = h(A)h(D) + B_t0(A)P_D(2);
(4) the dead-end descent: delete digit-0 vertices until none remain; max/min number of steps; the sterile endpoints."""
import subprocess
exec(open('nearjoin_check.py').read().split('# (1)')[0])
from collections import Counter
def digit0_vertices(S):
    n = len(S); dp, H = ham_all(S); N = 1 << n; res = []
    for u in range(n):
        st, _ = starts(S, u); B = sum(st[m] * H[(N-1) & ~m] for m in range(1, N) if (m >> u) & 1)
        if B == H[N-1]: res.append(u)
    return res
def P2(D):
    n = len(D); dp, H = ham_all(D); N = 1 << n
    return sum(H[S] * H[(N-1) & ~S] for S in range(N)), H[N-1]
def term_component_of_base(S, z):
    rest = [v for v in range(len(S)) if v != z]; B = sub(S, rest)
    # terminal component: vertices from which nothing outside their SCC is reachable
    n = len(B)
    def reach(s):
        r = {s}; f = [s]
        while f:
            a = f.pop()
            for b in range(n):
                if B[a][b] and b not in r: r.add(b); f.append(b)
        return r
    term = [v for v in range(n) if reach(v) == scc_of(B, v)]
    comp = scc_of(B, term[0]); return B, sorted(comp), rest
# (1)
bad = 0; tot = 0
for n in range(3, 9):
    for S in classes(n):
        dp, H = ham_all(S); N = 1 << n
        for u in range(n):
            outs = [x for x in range(n) if S[u][x]]
            st, _ = starts(S, u); B = sum(st[m] * H[(N-1) & ~m] for m in range(1, N) if (m >> u) & 1)
            isd0 = (B == H[N-1])
            if len(outs) == 0: pred = True
            elif len(outs) == 1:
                t0 = outs[0]; rest = [v for v in range(n) if v != u]; Sr = sub(S, rest); dpr, Hr = ham_all(Sr)
                pred = (dpr[(1 << (n-1)) - 1][rest.index(t0)] == 0)
            else: pred = False
            tot += 1
            if pred != isd0: bad += 1
print(f"(1) digit-0 characterisation (sink, or sole exit with dead-end port) on all tournaments of orders 3-8: {tot} vertices, {bad} mismatches")
# (2)
pow2 = 0; notpow2 = []; cnt = 0
for n in range(3, 8):
    for D in classes(n):
        if not strong(D): continue
        p, h = P2(D); cnt += 1
        if p % h == 0 and (p // h) & (p // h - 1) == 0: pow2 += 1
        elif len(notpow2) < 6: notpow2.append((n, h, p))
print(f"(2) strong D of orders 3-7: {cnt}; P_D(2)/h(D) a power of two in {pow2}; first non-examples (n, h, P): {notpow2}")
# (3)
single = 0; multi = 0; idbad = 0; multi_counts = Counter()
for n in range(4, 9):
    for S in classes(n):
        if not strong(S): continue
        for z in digit0_vertices(S):
            B, comp, rest = term_component_of_base(S, z)
            t0 = rest.index([x for x in range(n) if S[z][x]][0])
            A_idx = [v for v in range(len(B)) if v not in comp]; A = sub(B, A_idx); D = sub(B, comp)
            pD, hD = P2(D); hA = ham_all(A)[1][-1] if A_idx else 1
            BtA = 0
            if A_idx:
                ta = A_idx.index(t0); dpA, HA = ham_all(A); NA = 1 << len(A); stA, _ = starts(A, ta)
                BtA = sum(stA[m] * HA[(NA-1) & ~m] for m in range(1, NA) if (m >> ta) & 1)
            HS = ham_all(S)[1][-1]
            if HS != hA * hD + BtA * pD: idbad += 1
            if len(comp) == 1: single += 1
            else: multi += 1; multi_counts[(hD, pD)] += 1
print(f"(3) strong digit-0 pairs (S, z), orders 4-8: single-vertex terminal component (near-join proper): {single}; strong terminal D (|D| >= 3): {multi}; identity H = h(A)h(D) + B_t0(A) P_D(2) failures: {idbad}; (h(D), P_D(2)) of the right sources: {dict(multi_counts)}")
# (4)
steps_max = Counter(); sterile_counts = Counter(); sterile_strong = 0; strong_total = 0
import functools
def key(adj): return tuple(tuple(r) for r in adj)
@functools.lru_cache(None)
def descent(k):
    adj = [list(r) for r in k]; n = len(adj)
    d0 = digit0_vertices(adj)
    if not d0 or n == 1: return 0, 0
    lo = hi = None
    for z in d0:
        rest = [v for v in range(n) if v != z]; a, b = descent(key(sub(adj, rest)))
        lo = a + 1 if lo is None else min(lo, a + 1); hi = b + 1 if hi is None else max(hi, b + 1)
    return lo, hi
for n in range(3, 9):
    for S in classes(n):
        if not strong(S): continue
        strong_total += 1
        d0 = digit0_vertices(S)
        if not d0:
            sterile_strong += 1; sterile_counts[ham_all(S)[1][-1]] += 1; continue
        lo, hi = descent(key(S)); steps_max[hi] += 1
print(f"(4) strong tournaments orders 3-8: {strong_total}; sterile (no digit-0 vertex): {sterile_strong}; max descent steps distribution over the rest: {dict(sorted(steps_max.items()))}")
print("    counts of sterile strong tournaments (first 30):", sorted(sterile_counts)[:30])
