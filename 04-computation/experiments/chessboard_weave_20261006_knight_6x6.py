#!/usr/bin/env python3
"""6x6 closed knight's tours vs rings and scaffolds; edge counts c(e) and their
parity (chessboard-weave, 2026-10-06).

Compiles and runs chessboard_weave_20261006_knight_tours.c (exhaustive DFS,
positive controls 5x6=8, 3x10=16, 3x12=176, 5x8=44202, 6x6=9862) and
chessboard_weave_20261006_knight_hp.c (all directed open tours by
(start, second, end); positive control 6x6 = 6637920).
Run: python3 chessboard_weave_20261006_knight_6x6.py
"""
import os, subprocess, tempfile
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="knight6x6_")
def build(name):
    exe = os.path.join(TMP, name)
    subprocess.run(["cc", "-O2", "-o", exe, os.path.join(HERE, f"chessboard_weave_20261006_knight_{name}.c")], check=True)
    return exe
TOURS, HP = build("tours"), build("hp")

N = 6
def ring(s, n=N):
    c = (n - 1) / 2
    return int(max(abs(s[0] - c), abs(s[1] - c)) - 0.5)
SQ = [(i, j) for i in range(N) for j in range(N)]
RG = [ring(s) for s in SQ]
COL = [(i + j) % 2 for i, j in SQ]
def alg(k):
    return "abcdef"[SQ[k][1]] + str(SQ[k][0] + 1)
E = set()
for k, (i, j) in enumerate(SQ):
    for di, dj in [(1, 2), (2, 1), (-1, 2), (-2, 1)]:
        a, b = i + di, j + dj
        if 0 <= a < N and 0 <= b < N:
            E.add(tuple(sorted((k, a * N + b))))
E = sorted(E)
DEG = Counter(x for e in E for x in e)
assert len(E) == 80

# positive controls
for r, c, want in [(5, 6, 8), (3, 10, 16), (3, 12, 176), (5, 8, 44202)]:
    o = subprocess.run([TOURS, str(r), str(c)], capture_output=True, text=True, check=True).stdout.strip()
    print(o, " (known:", want, ")"); assert o.endswith(str(want))
tf = os.path.join(TMP, "t66.txt")
o = subprocess.run([TOURS, "6", "6", tf], capture_output=True, text=True, check=True).stdout.strip()
print(o, " (known: 9862)")
tours = [list(map(int, l.split())) for l in open(tf)]
assert len(tours) == 9862 and all(sorted(t) == list(range(36)) for t in tours)
Es = set(E)
def moves(t):
    return [tuple(sorted((t[k], t[(k + 1) % 36]))) for k in range(36)]
assert all(all(m in Es for m in moves(t)) for t in tours)
assert len({frozenset(moves(t)) for t in tours}) == 9862   # distinct as undirected cycles

# D4 on squares
def d4maps():
    fs = [lambda i, j: (i, j), lambda i, j: (j, N-1-i), lambda i, j: (N-1-i, N-1-j), lambda i, j: (N-1-j, i),
          lambda i, j: (i, N-1-j), lambda i, j: (N-1-i, j), lambda i, j: (j, i), lambda i, j: (N-1-j, N-1-i)]
    return [[f(*SQ[k])[0] * N + f(*SQ[k])[1] for k in range(36)] for f in fs]
G = d4maps()
tourset = {frozenset(moves(t)) for t in tours}
print("tour set D4-invariant:", all(frozenset(tuple(sorted((g[a], g[b]))) for a, b in T) in tourset
                                    for g in G for T in list(tourset)[:500]), "(checked on 500 tours x 8 maps)")

# ---- c(e)
c = Counter(m for t in tours for m in moves(t))
print("\n== c(e) = number of undirected closed tours through edge e ==")
print("sum c(e) =", sum(c.values()), "= 36 x 9862 =", 36 * 9862)
orbits, seen = [], set()
for e in E:
    if e in seen:
        continue
    orb = {tuple(sorted((g[e[0]], g[e[1]]))) for g in G}
    seen |= orb
    orbits.append(sorted(orb))
print("D4-invariance of c:", all(len({c[x] for x in orb}) == 1 for orb in orbits), f"({len(orbits)} edge orbits)")
print(f"{'rep edge':10s} {'orbit':>5s} {'rings':>6s} {'degs':>6s} {'c(e)':>6s} parity")
for orb in sorted(orbits, key=lambda o: (tuple(sorted((RG[o[0][0]], RG[o[0][1]]))), -c[o[0]])):
    a, b = orb[0]
    print(f"{alg(a)+'-'+alg(b):10s} {len(orb):5d} {str(tuple(sorted((RG[a], RG[b])))):>6s} "
          f"{str(tuple(sorted((DEG[a], DEG[b])))):>6s} {c[orb[0]]:6d} {'odd' if c[orb[0]] % 2 else 'even'}")
odd_edges = [e for e in E if c[e] % 2]
print("edges with odd c(e):", len(odd_edges), "of 80;  per ring pair:",
      dict(sorted(Counter(tuple(sorted((RG[a], RG[b]))) for a, b in odd_edges).items())))
# parity of c(e) at each square: sum_{e ni v} c(e) = 2 * 9862 (every tour uses 2 edges at v)
print("check sum_{e at v} c(e) = 2*9862 at every square:", all(sum(c[e] for e in E if v in e) == 2 * 9862 for v in range(36)))
# odd-c subgraph degrees
oddeg = Counter(x for e in odd_edges for x in e)
print("odd-c subgraph: degree of each square (row 6 at top):")
for i in reversed(range(N)):
    print("   ", " ".join(f"{oddeg[i*N+j]}" for j in range(N)))

# ---- ring statistics
print("\n== ring statistics over the 9862 tours ==")
def pairs(t):
    return Counter(tuple(sorted((RG[a], RG[b]))) for a, b in moves(t))
P = [pairs(t) for t in tours]
keys = sorted({k for p in P for k in p} | {tuple(sorted((RG[a], RG[b]))) for a, b in E})
print("ring pairs with knight edges:", keys, " edge counts:",
      dict(sorted(Counter(tuple(sorted((RG[a], RG[b]))) for a, b in E).items())))
for k in keys:
    vals = Counter(p[k] for p in P)
    print(f"  m{k[0]}{k[1]}: min {min(vals)} max {max(vals)}  distribution {dict(sorted(vals.items()))}")
W = Counter(sum(p[k] for k in keys if k[0] == k[1]) for p in P)
print("within-ring moves W=m11+m22: distribution", dict(sorted(W.items())))
blocks = Counter(sum(1 for k in range(36) if RG[t[k]] != RG[t[(k + 1) % 36]]) for t in tours)
print("ring blocks: distribution", dict(sorted(blocks.items())), " (= 36 - W:",
      all(sum(1 for k in range(36) if RG[t[k]] != RG[t[(k+1) % 36]]) == 36 - sum(p[k] for k in keys if k[0] == k[1])
          for t, p in zip(tours, P)), ")")
d2 = Counter(sum(p[k] for k in keys if k[1] - k[0] == 2) for p in P)
print("|ring change|=2 moves (= m02): distribution", dict(sorted(d2.items())))
vec = Counter(tuple(p[k] for k in keys) for p in P)
print(f"distinct ring-transition vectors {tuple('m%d%d' % k for k in keys)}: {len(vec)}")
for v, n_ in sorted(vec.items()):
    print("   ", v, n_)
print("identities (PROVED by half-edge counting) hold on all tours:",
      all(p[(0, 1)] + p[(0, 2)] == 8 and 2 * p[(1, 1)] + p[(0, 1)] + p[(1, 2)] == 24
          and 2 * p[(2, 2)] + p[(0, 2)] + p[(1, 2)] == 40 and p[(2, 2)] == 12 + p[(1, 1)] - p[(0, 2)] for p in P))
# count-admissible vectors: m01+m02=8, m12=24-2m11-m01, m22=4+m01+m11, each within its edge supply
supply = Counter(tuple(sorted((RG[a], RG[b]))) for a, b in E)
adm = []
for m01 in range(9):
    for m11 in range(supply[(1, 1)] + 1):
        v = (m01, 8 - m01, m11, 24 - 2 * m11 - m01, 4 + m01 + m11)
        if all(0 <= x <= supply[k] for x, k in zip(v, keys)):
            adm.append(v)
print(f"count-admissible vectors (half-edge identities + edge supplies): {len(adm)}; realised {len(vec)};"
      f" admissible but NOT realised: {[v for v in adm if v not in vec]}")
# rim within moves (ring 2 inner matching): each tour uses how many
# ring-2 within edges are the 8 K2's next to corners; check forcedness pattern
rim_within = [e for e in E if RG[e[0]] == RG[e[1]] == 2]
print("ring-2 within edges:", [alg(a) + "-" + alg(b) for a, b in rim_within],
      " c(e):", [c[e] for e in rim_within])
inner_c4 = [e for e in E if RG[e[0]] == RG[e[1]] == 1]
print("ring-1 within edges (two C4):", [alg(a) + "-" + alg(b) for a, b in inner_c4], " c(e):", [c[e] for e in inner_c4])

# ---- cross-check c(e) against open-tour table: c(uv) = #Ham paths u->v
hf = os.path.join(TMP, "hp66.txt")
o = subprocess.run([HP, "6", "6", hf], capture_output=True, text=True, check=True).stdout.strip()
print("\n" + o, " (known 6637920)")
Cxyz = defaultdict(int)
for l in open(hf):
    x, y, z, n_ = map(int, l.split())
    Cxyz[(x, y, z)] = n_
HPuv = defaultdict(int)
for (x, y, z), n_ in Cxyz.items():
    HPuv[(x, z)] += n_
print("c(uv) == #directed Ham paths u->v for every knight edge:", all(HPuv[(a, b)] == c[(a, b)] == HPuv[(b, a)] for a, b in E))
