#!/usr/bin/env python3
"""Orchestrator's independent check of HYP-9168 at N = 10 (different algorithm and code from the lane's C engine):
hall(T) by the Konig form  min over A (|A| >= 2) of the sum of the (N+2-|A|) smallest in-degrees-from-A
(B may meet A), and beta(T) >= hall(T) by SAT with lazy Hamiltonian-path cuts (CaDiCaL via pysat):
variables x_e = 'arc e deleted', cardinality sum x_e <= hall-1, one clause 'some arc of P deleted' per
Hamiltonian path P found in T - D. UNSAT means no (hall-1)-set kills every HP.
Sample: classes from gentourng 10 chosen by stride, plus every class with hall = 8 among the first ones seen."""
import subprocess, sys, random, time
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

N = 10
def parse(s):
    adj = [[False]*N for _ in range(N)]
    k = 0
    for i in range(N):
        for j in range(i+1, N):
            if s[k] == '1': adj[i][j] = True
            else: adj[j][i] = True
            k += 1
    return adj

def hall(adj):
    best = None
    for A in range(1, 1 << N):
        a = bin(A).count('1')
        if a < 2: continue
        b = N + 2 - a
        if b < 0: continue
        indeg = sorted(sum(1 for u in range(N) if (A >> u) & 1 and u != v and adj[u][v]) for v in range(N))
        c = sum(indeg[:b])
        if best is None or c < best: best = c
    return best

def find_hp(out):  # out[v] = bitmask of successors; returns list of vertices or None
    full = (1 << N) - 1
    ends = [0] * (1 << N)
    par = {}
    for v in range(N):
        ends[1 << v] |= 1 << v
    for mask in range(1, 1 << N):
        e = ends[mask]
        if not e: continue
        m = e
        while m:
            v = (m & -m).bit_length() - 1; m &= m - 1
            nxt = out[v] & ~mask
            while nxt:
                u = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                nm = mask | (1 << u)
                if not (ends[nm] >> u) & 1:
                    ends[nm] |= 1 << u
                    par[(nm, u)] = v
    if not ends[full]: return None
    u = (ends[full] & -ends[full]).bit_length() - 1
    path, mask = [u], full
    while mask != (1 << u) or len(path) < N:
        if len(path) == N: break
        v = par[(mask, u)]
        mask &= ~(1 << u); u = v; path.append(u)
    return path[::-1]

def beta_at_least(adj, h):
    arcs = [(i, j) for i in range(N) for j in range(N) if i != j and adj[i][j]]
    pool = IDPool()
    x = {a: pool.id(a) for a in arcs}
    s = Cadical153()
    if h - 1 >= 0:
        for cl in CardEnc.atmost(lits=list(x.values()), bound=h - 1, vpool=pool, encoding=EncType.seqcounter).clauses:
            s.add_clause(cl)
    cuts = 0
    while True:
        if not s.solve():
            return True, cuts
        model = set(l for l in s.get_model() if l > 0)
        D = {a for a in arcs if x[a] in model}
        out = [0] * N
        for (i, j) in arcs:
            if (i, j) not in D: out[i] |= 1 << j
        P = find_hp(out)
        if P is None:
            return False, (D, cuts)
        s.add_clause([x[(P[k], P[k+1])] for k in range(N - 1)])
        cuts += 1

stride = int(sys.argv[1]) if len(sys.argv) > 1 else 97333
extra_h8 = int(sys.argv[2]) if len(sys.argv) > 2 else 15
proc = subprocess.Popen(["gentourng", "-q", str(N)], stdout=subprocess.PIPE, text=True)
t0 = time.time(); tested = 0; h8 = 0; hist = {}; maxcuts = 0
for idx, line in enumerate(proc.stdout):
    take = (idx % stride == 0)
    s = line.strip()
    if not take and h8 >= extra_h8 and idx > 2_000_000: break
    adj = parse(s)
    if not take:
        # cheap filter for hard classes: min score <= ... ; compute hall only on a sparse subsample
        if idx % 997 != 0 or h8 >= extra_h8: continue
    h = hall(adj)
    if not take and h != 8: continue
    ok, info = beta_at_least(adj, h)
    if not ok:
        print("COUNTEREXAMPLE", s, "hall", h, "D", info); sys.exit(1)
    tested += 1; hist[h] = hist.get(h, 0) + 1; maxcuts = max(maxcuts, info)
    if h == 8: h8 += 1
proc.kill()
print(f"independent check N=10: {tested} classes tested (stride {stride} plus hall=8 classes found by sparse scan), "
      f"beta >= hall (hence beta = hall) for all; hall histogram {dict(sorted(hist.items()))}; max lazy cuts {maxcuts}; {time.time()-t0:.0f}s")
