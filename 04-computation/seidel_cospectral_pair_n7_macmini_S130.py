#!/usr/bin/env python3
"""
The first Seidel-cospectral pair of switching classes: n = 7  (mac-mini-2026-07-20-S130)
========================================================================================
THM-1440(B) found 11 distinct skew-Seidel spectra against A049313(7) = 12 switching
classes, so the spectrum is a COMPLETE invariant of the switching class for n <= 6 and
first fails at n = 7 with exactly ONE cospectral pair.  This isolates that pair.

Method: work from the 456 iso-class representatives (not all 2^21 tournaments).  Two iso
classes are switching-related iff some switching of one is isomorphic to the other, so
456 reps x 64 switchings = 29184 canonicalizations, batched over the 5040 permutations.
"""
import numpy as np
from itertools import permutations

n = 7
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
idx = {p: k for k, p in enumerate(pairs)}
E = len(pairs); pow2 = (1 << np.arange(E, dtype=np.int64))

PERM = []
for p in permutations(range(n)):
    src = np.empty(E, dtype=np.int64); fl = np.zeros(E, dtype=np.uint8)
    for e, (i, j) in enumerate(pairs):
        a, b = p[i], p[j]
        t = idx[(min(a, b), max(a, b))]
        src[t] = e; fl[t] = 1 if a > b else 0
    PERM.append((src, fl))

def codes_to_A(codes):
    c = np.asarray(codes, dtype=np.int64).reshape(-1, 1)
    return ((c >> np.arange(E, dtype=np.int64)) & 1).astype(np.uint8)

def canon(codes):
    A = codes_to_A(codes); best = None
    for src, fl in PERM:
        c = (A[:, src] ^ fl) @ pow2
        best = c if best is None else np.minimum(best, c)
    return best

def skew(code):
    S = np.zeros((n, n))
    for e, (i, j) in enumerate(pairs):
        if code >> e & 1: S[j, i], S[i, j] = 1.0, -1.0
        else:             S[i, j], S[j, i] = 1.0, -1.0
    return S

# ---- iso-class reps by vertex extension --------------------------------------
reps = {0}
for k in range(2, n + 1):
    pk = [(i, j) for i in range(k) for j in range(i + 1, k)]
    ik = {p: t for t, p in enumerate(pk)}
    op = [(i, j) for i in range(k - 1) for j in range(i + 1, k - 1)]
    cand = []
    for r in reps:
        base = 0
        for e, (i, j) in enumerate(op):
            if r >> e & 1: base |= 1 << ik[(i, j)]
        for mask in range(1 << (k - 1)):
            v = base
            for b in range(k - 1):
                if mask >> b & 1: v |= 1 << ik[(b, k - 1)]
            cand.append(v)
    if k < n:
        Ek = len(pk); p2 = (1 << np.arange(Ek, dtype=np.int64))
        Ak = ((np.array(cand, dtype=np.int64)[:, None] >> np.arange(Ek)) & 1).astype(np.uint8)
        best = None
        for p in permutations(range(k)):
            src = np.empty(Ek, dtype=np.int64); fl = np.zeros(Ek, dtype=np.uint8)
            for e, (i, j) in enumerate(pk):
                a, b = p[i], p[j]
                t = ik[(min(a, b), max(a, b))]
                src[t] = e; fl[t] = 1 if a > b else 0
            c = (Ak[:, src] ^ fl) @ p2
            best = c if best is None else np.minimum(best, c)
        reps = set(int(x) for x in best.tolist())
    else:
        reps = set(int(x) for x in canon(cand).tolist())
reps = sorted(reps)
print(f"n={n}: iso classes = {len(reps)}")

# ---- switching masks (X subset of {0..n-2}, kills the X ~ X^c redundancy) -----
switches = []
for m in range(1 << (n - 1)):
    b = 0
    for e, (i, j) in enumerate(pairs):
        bi = (m >> i & 1) if i < n - 1 else 0
        bj = (m >> j & 1) if j < n - 1 else 0
        if bi ^ bj: b |= 1 << e
    switches.append(b)
print(f"switchings = {len(switches)}")

# ---- union-find over iso classes, joined by switching ------------------------
pos = {r: t for t, r in enumerate(reps)}
parent = list(range(len(reps)))
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[max(ra, rb)] = min(ra, rb)

flat = [r ^ s for r in reps for s in switches]
canon_flat = canon(flat).tolist()
for t, r in enumerate(reps):
    for u in range(len(switches)):
        union(pos[r], pos[canon_flat[t * len(switches) + u]])

groups = {}
for r in reps: groups.setdefault(find(pos[r]), []).append(r)
print(f"switching classes up to isomorphism = {len(groups)}   (A049313(7) = 12)")

# ---- spectra -----------------------------------------------------------------
spec_of = {}
for g, members in groups.items():
    s = tuple(np.round(np.sort(np.linalg.eigvalsh(1j * skew(members[0])).real), 6))
    spec_of.setdefault(s, []).append((g, members))
print(f"distinct spectra = {len(spec_of)}")
print()
for s, gs in spec_of.items():
    if len(gs) > 1:
        print("*** COSPECTRAL PAIR OF SWITCHING CLASSES (n=7) ***")
        print(f"  shared spectrum of iS: {[round(float(v),5) for v in s]}")
        for g, members in gs:
            rep = members[0]; S = skew(rep)
            w = (S > 0).sum(axis=1)
            c3 = n*(n-1)*(n-2)//6 - sum(int(x)*(int(x)-1)//2 for x in w)
            Hs = []
            for m in members[:200]:
                Sm = skew(m); wm = (Sm > 0)
                dp = [[0]*n for _ in range(1 << n)]
                for v in range(n): dp[1 << v][v] = 1
                for Ssub in range(1 << n):
                    for v in range(n):
                        c = dp[Ssub][v]
                        if not c or not (Ssub >> v & 1): continue
                        for u in range(n):
                            if Ssub >> u & 1 or not wm[v][u]: continue
                            dp[Ssub | 1 << u][u] += c
                Hs.append(sum(dp[(1 << n) - 1]))
            print(f"  class: {len(members)} iso classes, rep code {rep}, "
                  f"scores {sorted(w.tolist())}, 3-cycles {c3}, H values {sorted(set(Hs))}")
