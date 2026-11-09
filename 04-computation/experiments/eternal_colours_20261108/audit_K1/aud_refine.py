#!/usr/bin/env python3
"""audit K1: independent recomputation of the note's section-5 refinement table (written from scratch).
For every value h >= 3 with least order n(h) <= 9, the birth fibre is the set of classes of order n(h) with H = h.
A refinement keeps the classes extremal in a statistic; the owner's dichotomy holds on the refined fibre if it is one
self-converse class, or two non-self-converse classes that are converses of each other.
Class fields (H, c3, SC, |Aut|) come from THM-4922's census, re-verified class by class by aud_census.c; the converse
pairing of every two-class refined fibre is re-checked here by a direct isomorphism search between the converse of one
class and the other (it does not use the census's ccode field).
usage: aud_refine.py DIR_WITH_cls_n.txt
"""
import sys, os, collections
D = sys.argv[1]

def decode(code, n):
    out = [0] * n; b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> b & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
            b += 1
    return out

def isomorphic(A, B, n):
    da = [bin(x).count("1") for x in A]; db = [bin(x).count("1") for x in B]
    if sorted(da) != sorted(db): return False
    p = [0] * n; used = [False] * n
    def bt(k):
        if k == n: return True
        for w in range(n):
            if used[w] or db[w] != da[k]: continue
            if all(((A[j] >> k) & 1) == ((B[p[j]] >> w) & 1) for j in range(k)):
                p[k] = w; used[w] = True
                if bt(k + 1): return True
                used[w] = False
        return False
    return bt(0)

def converse(A, n):
    full = (1 << n) - 1
    return [(full & ~A[v]) & ~(1 << v) for v in range(n)]

fib = {}
seen = set()
for n in range(1, 10):
    cur = collections.defaultdict(list)
    for line in open(os.path.join(D, f"cls_{n}.txt")):
        f = line.split("|")[0].split()
        code, H, c3, SC, aut = int(f[0]), int(f[1]), int(f[2]), int(f[4]), int(f[5])
        if H not in seen: cur[H].append((code, n, c3, SC, aut))
    for H, L in cur.items(): fib[H] = L
    seen |= set(cur)
vals = sorted(h for h in fib if h >= 3)
conv_cache = {}
def holds(L):
    if len(L) == 1: return L[0][3] == 1
    if len(L) == 2:
        a, b = L
        if a[3] or b[3]: return False
        key = (a[0], b[0])
        if key not in conv_cache:
            n = a[1]
            conv_cache[key] = isomorphic(converse(decode(a[0], n), n), decode(b[0], n), n)
        return conv_cache[key]
    return False
def ext(L, i, mx): v = (max if mx else min)(c[i] for c in L); return [c for c in L if c[i] == v]
rules = [("none", lambda L: L),
         ("max |Aut|", lambda L: ext(L, 4, True)),
         ("min |Aut|", lambda L: ext(L, 4, False)),
         ("min c3", lambda L: ext(L, 2, False)),
         ("max c3", lambda L: ext(L, 2, True)),
         ("max |Aut| then max c3", lambda L: ext(ext(L, 4, True), 2, True)),
         ("max |Aut| then min c3", lambda L: ext(ext(L, 4, True), 2, False))]
print("values with birth order <= 9 (h >= 3):", len(vals))
for name, r in rules:
    print(f"{name:24s}: dichotomy holds for {sum(holds(r(fib[h])) for h in vals)}")
