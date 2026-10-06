#!/usr/bin/env python3
"""Owner's 8x8 board, decoded exactly (chessboard-weave session, 2026-10-06).

Objects: rings (Chebyshev shells about the centre), the two ferz scaffolds
(colour classes under diagonal adjacency), sliding lines, the knight graph.
Prints every count the owner's prompt states plus the cross-structure.
"""
from itertools import product
from collections import Counter, defaultdict

N = 8
SQ = [(i, j) for i in range(N) for j in range(N)]

def ring(s, n=N):
    # Chebyshev shell index about the centre ((n-1)/2,(n-1)/2); 0 = innermost
    i, j = s
    c = (n - 1) / 2
    return int(max(abs(i - c), abs(j - c)) - 0.5)

def colour(s):
    return (s[0] + s[1]) % 2

def on_board(s, n=N):
    return 0 <= s[0] < n and 0 <= s[1] < n

def leaper_edges(a, b, n=N):
    E = set()
    for s in [(i, j) for i in range(n) for j in range(n)]:
        for da, db in {(a, b), (b, a)}:
            for sa, sb in product((1, -1), repeat=2):
                t = (s[0] + sa * da, s[1] + sb * db)
                if on_board(t, n) and t != s:
                    E.add(frozenset((s, t)))
    return E

def components(V, E):
    adj = defaultdict(set)
    for e in E:
        u, v = tuple(e)
        adj[u].add(v); adj[v].add(u)
    seen, comps = set(), []
    for v in V:
        if v in seen:
            continue
        stack, comp = [v], []
        seen.add(v)
        while stack:
            x = stack.pop(); comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); stack.append(y)
        comps.append(comp)
    return comps, adj

print("== 1. rings (Chebyshev shells) ==")
rs = Counter(ring(s) for s in SQ)
print("ring sizes:", [rs[k] for k in range(4)], "(owner: 4,12,20,28)")
W = leaper_edges(0, 1)          # wazir = one rook step
# is each ring a cycle in the wazir graph?
for k in range(4):
    Vk = [s for s in SQ if ring(s) == k]
    Ek = [e for e in W if all(ring(x) == k for x in e)]
    comps, adj = components(Vk, Ek)
    degs = Counter(len(adj[v]) for v in Vk)
    print(f"  ring {k}: |V|={len(Vk)} |E|={len(Ek)} components={len(comps)} degrees={dict(degs)}"
          f" -> cycle C_{len(Vk)}: {len(comps)==1 and set(degs)=={2} and len(Ek)==len(Vk)}")
spokes = Counter(tuple(sorted(ring(x) for x in e)) for e in W if len({ring(x) for x in e}) == 2)
print("  wazir edges between rings (spokes):", dict(spokes), "total wazir edges", len(W))

print("\n== 2. ferz scaffolds (diagonal one-step) ==")
F = leaper_edges(1, 1)
comps, adj = components(SQ, F)
print("components:", len(comps), "sizes", [len(c) for c in comps])
for c in comps:
    col = {colour(s) for s in c}
    anti = Counter(s[0] + s[1] for s in c)    # anti-diagonals i+j = const
    diag = Counter(s[0] - s[1] for s in c)    # diagonals i-j = const
    print(f"  colour {col}: anti-diagonal lengths {[anti[k] for k in sorted(anti)]}; "
          f"diagonal lengths {[diag[k] for k in sorted(diag)]}")
# odd-length lines partition the board; so do even-length lines
lines = []
for s0 in range(2 * N - 1):
    lines.append(("anti", s0, [s for s in SQ if s[0] + s[1] == s0]))
for d0 in range(-(N - 1), N):
    lines.append(("diag", d0, [s for s in SQ if s[0] - s[1] == d0]))
odd = [L for L in lines if len(L[2]) % 2 == 1]
even = [L for L in lines if len(L[2]) % 2 == 0]
cover_odd = Counter(s for L in odd for s in L[2])
cover_even = Counter(s for L in even for s in L[2])
print("  odd-length diagonals:", len(odd), "lines; every square covered exactly once:",
      all(cover_odd[s] == 1 for s in SQ))
print("  even-length diagonals:", len(even), "lines; every square covered exactly once:",
      all(cover_even[s] == 1 for s in SQ))
print("  odd line lengths:", sorted(len(L[2]) for L in odd))
print("  ring k size 4(2k+1) = 4 x (odd line length 2k+1):",
      [rs[k] == 4 * (2 * k + 1) == sum(1 for L in odd if len(L[2]) == 2 * k + 1) * (2 * k + 1) for k in range(4)])

print("\n== 3. knight graph: weave of scaffolds and rings ==")
K = leaper_edges(1, 2)
deg = Counter()
for e in K:
    for x in e:
        deg[x] += 1
print("knight edges:", len(K), " degree histogram:", dict(sorted(Counter(deg.values()).items())))
print("every knight edge joins the two colours (scaffolds):", all(colour(a) != colour(b) for a, b in map(tuple, K)))
RR = Counter(tuple(sorted(ring(x) for x in e)) for e in K)
print("knight edges by ring pair:", dict(sorted(RR.items())))
M = [[0] * 4 for _ in range(4)]
for e in K:
    a, b = tuple(e)
    M[ring(a)][ring(b)] += 1
    M[ring(b)][ring(a)] += 1
print("ring-incidence matrix (half-edges, row r = edges leaving ring r):")
for r in range(4):
    print("   ", M[r], " row sum", sum(M[r]), " ring size", rs[r])
# ring change distribution of a knight move
dr = Counter(abs(ring(a) - ring(b)) for a, b in map(tuple, K))
print("|ring change| of knight edges:", dict(sorted(dr.items())))

print("\n== 4. D4 orbits by ring ==")
def d4(s):
    i, j = s; n = N - 1
    imgs = [(i, j), (j, n - i), (n - i, n - j), (n - j, i), (j, i), (n - i, j), (i, n - j), (n - j, n - i)]
    return min(imgs)
orb = defaultdict(set)
for s in SQ:
    orb[d4(s)].add(s)
byring = Counter(ring(r) for r in orb)
print("orbits:", len(orb), " per ring:", [byring[k] for k in range(4)], " orbit sizes:",
      sorted(len(v) for v in orb.values()))

print("\n== 5. sliding lines (rook/bishop/queen on the empty board) ==")
def slide_moves(s, dirs):
    out = 0
    for d in dirs:
        t = (s[0] + d[0], s[1] + d[1])
        while on_board(t):
            out += 1
            t = (t[0] + d[0], t[1] + d[1])
    return out
R4 = [(1, 0), (-1, 0), (0, 1), (0, -1)]
B4 = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
for name, dirs in (("rook", R4), ("bishop", B4), ("queen", R4 + B4)):
    tot = sum(slide_moves(s, dirs) for s in SQ)
    byr = defaultdict(set)
    for s in SQ:
        byr[ring(s)].add(slide_moves(s, dirs))
    print(f"  {name}: total moves {tot} (edges {tot//2}); mobility by ring {dict(byr)}")
print("  bishop mobility = 7 + 2*ring (ring = 0..3)?",
      all(slide_moves(s, B4) == 13 - 2 * ring(s) for s in SQ), "(bishop: 13,11,9,7 by ring 0..3)")
print("  king mobility by ring:", {k: sorted({slide_moves(s, []) + sum(1 for d in R4 + B4 if on_board((s[0]+d[0], s[1]+d[1]))) for s in SQ if ring(s) == k}) for k in range(4)})
