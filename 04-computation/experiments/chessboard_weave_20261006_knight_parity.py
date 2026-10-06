#!/usr/bin/env python3
"""Parity of Hamiltonian path / cycle counts in rectangular knight graphs
(Redei/Thomason-type tests; chessboard-weave, 2026-10-06).

For each board R x C (exhaustive, chessboard_weave_20261006_knight_hp.c and
..._tours.c), with D(u,v) = # Hamiltonian paths with ends u,v (undirected),
h(e) = # undirected Hamiltonian paths through edge e, c(e) = # closed tours
through e, H(x) = # Hamiltonian paths with an end at x:
  T1 (PROVED, Thomason lollipop, free first edge): for every x,
       sum_{z: deg z even} D(x,z) is even.
  T2 (PROVED, lollipop with fixed first edge x->y): for every directed edge,
       c(xy) == sum_z N(x,y,z) (deg z - 1)  (mod 2), N = # Ham paths x,y,...,z.
  plus empirical parity statistics (odd counts of D, h, c, totals).
Run: python3 chessboard_weave_20261006_knight_parity.py
"""
import os, subprocess, tempfile
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="knightpar_")
def build(name):
    exe = os.path.join(TMP, name)
    subprocess.run(["cc", "-O2", "-o", exe, os.path.join(HERE, f"chessboard_weave_20261006_knight_{name}.c")], check=True)
    return exe
HPX, TOURS = build("hp"), build("tours")

def graph(R, C):
    nb = defaultdict(set)
    for i in range(R):
        for j in range(C):
            for di, dj in [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]:
                a, b = i + di, j + dj
                if 0 <= a < R and 0 <= b < C:
                    nb[i * C + j].add(a * C + b)
    return nb

def gf2rank(rows):
    rows = [r for r in rows if r]
    rank = 0
    while rows:
        p = rows.pop()
        rank += 1
        lb = p & -p
        rows = [r ^ p if r & lb else r for r in rows]
        rows = [r for r in rows if r]
    return rank

def name(k, C):
    return "abcdefghijklmnop"[k % C] + str(k // C + 1)

BOARDS = [(3, 4), (3, 7), (3, 8), (3, 9), (3, 10), (3, 11), (3, 12), (4, 5), (4, 6), (4, 7), (4, 8),
          (5, 5), (5, 6), (5, 7), (6, 6)]
summary = []
for R, C in BOARDS:
    nb = graph(R, C)
    NV = R * C
    deg = {v: len(nb[v]) for v in range(NV)}
    f1, f2 = os.path.join(TMP, "a.txt"), os.path.join(TMP, "b.txt")
    out = subprocess.run([HPX, str(R), str(C), f1, f2], capture_output=True, text=True, check=True).stdout.strip()
    Nxyz = {}
    for l in open(f1):
        x, y, z, n_ = map(int, l.split()); Nxyz[(x, y, z)] = n_
    hdir = {}
    for l in open(f2):
        u, v, n_ = map(int, l.split()); hdir[(u, v)] = n_
    Tdir = sum(Nxyz.values())
    D = defaultdict(int)
    for (x, y, z), n_ in Nxyz.items():
        D[(x, z)] += n_
    assert all(D[(u, v)] == D[(v, u)] for (u, v) in list(D))
    assert all(n_ % 2 == 0 for n_ in hdir.values())
    h = {e: n_ // 2 for e, n_ in hdir.items()}
    H = {x: sum(D[(x, z)] for z in range(NV)) for x in range(NV)}
    # T1
    t1 = all(sum(D[(x, z)] for z in range(NV) if deg[z] % 2 == 0) % 2 == 0 for x in range(NV))
    # closed tours through each edge from the open-path table
    cxy = {}
    t2 = True
    for x in range(NV):
        for y in nb[x]:
            cyc = sum(Nxyz.get((x, y, z), 0) for z in nb[x])
            cxy[(x, y)] = cyc
            rhs = sum(n_ * (deg[z] - 1) for (a, b, z), n_ in Nxyz.items() if a == x and b == y)
            if (cyc - rhs) % 2:
                t2 = False
    ncl = None
    if R >= 3 and C >= 3:
        o = subprocess.run([TOURS, str(R), str(C)], capture_output=True, text=True, check=True).stdout.strip()
        ncl = int(o.split()[-1])
        assert sum(cxy[(x, y)] for x in range(NV) for y in nb[x] if x < y) == NV * ncl
    oddD = [(u, v) for (u, v), n_ in D.items() if u < v and n_ % 2]
    oddh = [e for e, n_ in h.items() if n_ % 2]
    oddc = [(x, y) for (x, y), n_ in cxy.items() if x < y and n_ % 2]
    nE = sum(len(nb[v]) for v in range(NV)) // 2
    rows = [sum(1 << z for z in range(NV) if D[(x, z)] % 2) for x in range(NV)]
    print(f"===== {R}x{C}: {out.split(': ')[1]} directed Ham paths = {Tdir//2} undirected "
          f"({'odd' if (Tdir//2) % 2 else 'even'}); closed tours {ncl}; degrees {dict(sorted(Counter(deg.values()).items()))}")
    print(f"  T1 (even-degree-end paths from each x even): {t1};  T2 (lollipop, fixed first edge): {t2}")
    print(f"  endpoint pairs with D(u,v)>0: {sum(1 for (u, v), n_ in D.items() if u < v and n_)};"
          f" odd D(u,v): {len(oddD)}; GF(2)-rank of [D mod 2]: {gf2rank(rows)}")
    print(f"  odd-D pairs by sorted (deg u, deg v): {dict(sorted(Counter(tuple(sorted((deg[u], deg[v]))) for u, v in oddD).items()))}")
    print(f"  H(x) (paths ending at x) odd at {sum(1 for x in range(NV) if H[x] % 2)} squares; "
          f"H mod 2 by degree: {dict(sorted(Counter((deg[x], H[x] % 2) for x in range(NV)).items()))}")
    print(f"  edges: {nE}; edges with odd h(e) (Ham paths through e): {len(oddh)}; "
          f"edges with h(e)=0: {sum(1 for n_ in h.values() if n_ == 0)}")
    if ncl:
        print(f"  edges with odd c(e) (closed tours through e): {len(oddc)}: "
              + " ".join(f"{name(x, C)}-{name(y, C)}:{cxy[(x, y)]}" for x, y in sorted(oddc))[:400])
    summary.append((R, C, Tdir // 2, ncl, t1, t2, len(oddD), len(oddh), nE, len(oddc) if ncl else None))

print("\nSUMMARY  board | undirected HPs | closed tours | T1 | T2 | #odd D(u,v) | #odd h(e)/edges | #odd c(e)")
for R, C, T, ncl, t1, t2, oD, oh, nE, oc in summary:
    print(f"  {R}x{C:<3d} {T:>10d} {str(ncl):>8s} {str(t1):>5s} {str(t2):>5s} {oD:>6d} {oh:>5d}/{nE:<4d} {str(oc):>5s}")

print("""
ANTI-REDEI (PROVED): for every R x C board with R, C >= 2 the number of undirected
Hamiltonian paths of the knight graph is EVEN.  Proof: the Klein group V4 = {id, left-right
mirror, up-down mirror, half-turn} acts faithfully on the squares.  If g != id maps an
undirected Hamiltonian path P to itself it must reverse P (preserving the direction would fix
every square).  Two distinct non-identity elements cannot both reverse P (their product would
preserve P, hence be the identity).  So every V4-orbit of paths has size 2 or 4.
(Consistent with every row of the table above.)""")

print("\n== closed tours through each edge, c(e) mod 2, on boards with closed tours ==")
for R, C in [(3, 10), (5, 6), (3, 12), (6, 6), (3, 14), (5, 8), (3, 16), (6, 7), (5, 10)]:
    o = subprocess.run([TOURS, str(R), str(C), "-e"], capture_output=True, text=True, check=True).stdout.splitlines()
    ncl = int(o[0].split()[-1])
    ce = [tuple(map(int, l.split())) for l in o[1:]]
    odd = [(u, v) for u, v, c_ in ce if c_ % 2]
    adj = defaultdict(list)
    for u, v in odd:
        adj[u].append(v); adj[v].append(u)
    seen, comps = set(), []
    for s0 in adj:
        if s0 in seen:
            continue
        comp, st = [s0], [s0]; seen.add(s0)
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); comp.append(y); st.append(y)
        comps.append(len(comp))
    print(f"  {R}x{C}: {ncl} closed tours, {len(ce)} edges, odd c(e) on {len(odd)} edges;"
          f" odd-c subgraph degrees {sorted(Counter(len(a) for a in adj.values()).items())},"
          f" component sizes {sorted(comps)}"
          + (f"; e.g. {name(odd[0][0], C)}-{name(odd[0][1], C)} c={dict(((u, v), c_) for u, v, c_ in ce)[odd[0]]}" if odd else ""))
print("  (the odd-c edge set is always an even subgraph: sum_{e at v} c(e) = 2 x #tours)")
