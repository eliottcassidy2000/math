#!/usr/bin/env python3
"""Tour-blocking number of the 6x6 TORUS knight graph vs the Hall (2-factor) bound
(chessboard-weave, 2026-10-06).

G = knight graph on Z6 x Z6 (8-regular, 144 edges, bipartite 18+18).
  beta_HC(G) = min #edges whose deletion leaves no Hamiltonian cycle.
  Hall/2-factor bound: min #edges whose deletion leaves no 2-factor = 7 (PROVED:
  a bipartite G has a 2-factor iff e(S, B\\T) >= 2|S| - 2|T| for all S in A, T in B;
  in an 8-regular graph e(S,B\\T) >= 8|S|-8|T|, so a violation costs
  >= 6(|S|-|T|)+1 >= 7 deletions, attained by S={a},T={} i.e. leaving a with degree 1).
Steps
 1. structure checks; edge-transitivity of the explicit group Z6^2 x| D4.
 2. C program (..._knight_torusblock.c), K=6: every 6-set containing e0 leaves a
    Hamiltonian cycle (=> beta_HC >= 7 by edge-transitivity); pool cycles re-verified here.
 3. K=7 (negative control + census): the 7-sets containing e0 that kill all Hamiltonian
    cycles; expected exactly the 14 vertex-isolating ones.
 4. Hamiltonian PATH blocking: each blocking 7-set still leaves a Hamiltonian path
    (CP-SAT), so beta_HP = 8 (delete all 8 edges at a vertex).
 5. independent CP-SAT spot checks on random and structured 6-sets.
Run: python3 chessboard_weave_20261006_knight_torus.py
"""
import os, random, subprocess, tempfile, time
from itertools import combinations
from ortools.sat.python import cp_model

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="knighttorus_")
EXE = os.path.join(TMP, "torusblock")
subprocess.run(["cc", "-O2", "-o", EXE, os.path.join(HERE, "chessboard_weave_20261006_knight_torusblock.c")], check=True)

# edge indexing identical to the C program
D8 = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
eid, EU = {}, []
for i in range(6):
    for j in range(6):
        u = 6 * i + j
        for di, dj in D8:
            v = 6 * ((i + di) % 6) + (j + dj) % 6
            if u < v and (u, v) not in eid:
                eid[(u, v)] = len(EU); EU.append((u, v))
E = EU
assert len(E) == 144
adj = {v: set() for v in range(36)}
for u, v in E:
    adj[u].add(v); adj[v].add(u)
print("1. 6x6 torus knight graph: 36 vertices, 144 edges, degrees", sorted({len(adj[v]) for v in adj}),
      "; bipartite by (i+j) mod 2:", all((u // 6 + u % 6 + v // 6 + v % 6) % 2 == 1 for u, v in E))
def act(f, e):
    a, b = e
    fa, fb = f(divmod(a, 6)), f(divmod(b, 6))
    x, y = 6 * (fa[0] % 6) + fa[1] % 6, 6 * (fb[0] % 6) + fb[1] % 6
    return (min(x, y), max(x, y))
gens = [lambda p: (p[0] + 1, p[1]), lambda p: (p[0], p[1] + 1), lambda p: (p[1], -p[0]), lambda p: (p[1], p[0])]
orb, st = {E[0]}, [E[0]]
while st:
    e = st.pop()
    for g in gens:
        f = act(g, e)
        assert f in eid
        if f not in orb:
            orb.add(f); st.append(f)
print("   edge orbit of e0 =", E[0], "under translations + D4 has size", len(orb), "-> edge-transitive:", len(orb) == 144)

def run(K, seed, pool=None):
    args = [EXE, str(K), "1000", str(seed)] + ([pool] if pool else [])
    t0 = time.time()
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    return out, time.time() - t0

print("\n2. K = 6: every 6-edge deletion containing e0")
pf = os.path.join(TMP, "pool6.txt")
out6, dt = run(6, 1, pf)
print("   " + "\n   ".join(l for l in out6.splitlines() if not l.startswith("  first free edge")), f"\n   ({dt:.0f}s)")
ok = True
masks = []
for l in open(pf):
    w = [int(x, 16) for x in l.split()]
    m = w[0] | (w[1] << 64) | (w[2] << 128)
    edges = [E[k] for k in range(144) if m >> k & 1]
    deg = {v: 0 for v in range(36)}
    a = {v: [] for v in range(36)}
    for u, v in edges:
        deg[u] += 1; deg[v] += 1; a[u].append(v); a[v].append(u)
    seen, stck = {0}, [0]
    while stck:
        x = stck.pop()
        for y in a[x]:
            if y not in seen:
                seen.add(y); stck.append(y)
    ok &= len(edges) == 36 and all(d == 2 for d in deg.values()) and len(seen) == 36
    masks.append(m)
print(f"   independent re-verification: all {len(masks)} pool masks are Hamiltonian cycles of G: {ok}")
# independently re-certify a random sample of 6-sets from the pool
rng = random.Random(7)
sample_ok = 0
for _ in range(20000):
    D = {0} | set(rng.sample(range(1, 144), 5))
    Dm = sum(1 << k for k in D)
    sample_ok += any(m & Dm == 0 for m in masks)
print(f"   pool alone certifies {sample_ok}/20000 random 6-sets containing e0 (rest handled by the C finder)")

print("\n3. K = 7 (negative control): 7-sets containing e0 that leave NO Hamiltonian cycle")
out7, dt = run(7, 2)
blk = [l for l in out7.splitlines() if l.startswith("BLOCKING") or l.startswith("UNRESOLVED")]
print("   " + "\n   ".join(l for l in out7.splitlines() if l.startswith("TOTAL") or l.startswith("finder")), f"\n   ({dt:.0f}s)")
for l in blk:
    print("   ", l)

def ham(D, path=False, tl=60):
    m = cp_model.CpModel()
    arcs = []
    for k, (u, v) in enumerate(E):
        if k in D:
            continue
        arcs.append((u, v, m.NewBoolVar(""))); arcs.append((v, u, m.NewBoolVar("")))
    if path:   # dummy vertex 36 joined to everything: Ham cycle in G+dummy <=> Ham path in G
        for v in range(36):
            arcs.append((36, v, m.NewBoolVar(""))); arcs.append((v, 36, m.NewBoolVar("")))
    m.AddCircuit(arcs)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl
    s.parameters.num_search_workers = 2
    return s.StatusName(s.Solve(m))

print("\n4. Hamiltonian PATHS after deleting each blocking 7-set (CP-SAT):")
for l in blk:
    pairs = [tuple(map(int, t.strip("{}").split(","))) for t in l.split("=")[1].split()]
    D = {eid[(min(p), max(p))] for p in pairs}
    print("   ", "HC:", ham(D), " HP:", ham(D, path=True), " D-edges", sorted(D))

print("\n5. independent CP-SAT checks of 6-sets (all must be FEASIBLE/OPTIMAL):")
tests = []
at0 = [k for k, e in enumerate(E) if 0 in e]
tests.append(("6 edges at one vertex", set(at0[:6])))
nb = sorted(adj[0])
at1 = [k for k, e in enumerate(E) if nb[0] in e and 0 not in e]
tests.append(("3+3 at two adjacent vertices", set(at0[:3]) | set(at1[:3])))
tests.append(("4+2 at two adjacent vertices", set(at0[:4]) | set(at1[:2])))
for t in range(30):
    tests.append((f"random #{t}", set(rng.sample(range(144), 6))))
res = {}
for nm, D in tests:
    r = ham(D)
    if nm.startswith("random"):
        res[r] = res.get(r, 0) + 1
    else:
        print(f"   {nm}: {r}")
print("   30 random 6-sets:", res)
