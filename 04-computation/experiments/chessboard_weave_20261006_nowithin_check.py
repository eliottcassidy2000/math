#!/usr/bin/env python3
"""Independent check (orchestrator, pysat + lazy subtour cuts; not CP-SAT):
the 8x8 knight graph minus its 24 within-ring edges has NO Hamiltonian cycle,
i.e. every closed knight's tour uses at least one move inside a ring.
Positive control: with the within-ring edges restored, a Hamiltonian cycle is found."""
from itertools import combinations
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

N = 8
def ring(i, j):
    c = (N - 1) / 2
    return int(max(abs(i - c), abs(j - c)) - 0.5)

V = [(i, j) for i in range(N) for j in range(N)]
E_all = set()
for (i, j) in V:
    for a, b in ((1, 2), (2, 1), (-1, 2), (-2, 1)):
        u, v = i + a, j + b
        if 0 <= u < N and 0 <= v < N:
            E_all.add(frozenset(((i, j), (u, v))))
E_all = sorted(tuple(sorted(e)) for e in E_all)
within = [e for e in E_all if ring(*e[0]) == ring(*e[1])]
assert len(E_all) == 168 and len(within) == 24

def hamiltonian(E, label):
    pool = IDPool()
    x = {e: pool.id(("e",) + e) for e in E}
    inc = {v: [] for v in V}
    for e in E:
        inc[e[0]].append(x[e]); inc[e[1]].append(x[e])
    s = Cadical153()
    for v in V:
        if len(inc[v]) < 2:
            print(f"  [{label}] vertex {v} has degree {len(inc[v])} < 2: no Hamiltonian cycle"); return False
        for cl in CardEnc.equals(lits=inc[v], bound=2, vpool=pool, encoding=EncType.seqcounter).clauses:
            s.add_clause(cl)
    cuts = 0
    while True:
        if not s.solve():
            print(f"  [{label}] UNSAT after {cuts} subtour cuts: no Hamiltonian cycle")
            return False
        m = set(l for l in s.get_model() if l > 0)
        chosen = [e for e in E if x[e] in m]
        # components of the 2-factor
        adj = {v: [] for v in V}
        for a, b in chosen:
            adj[a].append(b); adj[b].append(a)
        seen, comps = set(), []
        for v in V:
            if v in seen: continue
            stack, comp = [v], []
            seen.add(v)
            while stack:
                y = stack.pop(); comp.append(y)
                for z in adj[y]:
                    if z not in seen: seen.add(z); stack.append(z)
            comps.append(comp)
        if len(comps) == 1:
            print(f"  [{label}] Hamiltonian cycle found after {cuts} cuts")
            return True
        for comp in comps:
            cs = set(comp)
            cut = [x[e] for e in E if (e[0] in cs) != (e[1] in cs)]
            s.add_clause(cut)           # at least one edge leaves every proper subset (2 by parity)
            cuts += 1

print("8x8 knight graph:", len(E_all), "edges;", len(within), "within-ring edges")
ok_pos = hamiltonian(E_all, "all 168 edges (positive control)")
ok_neg = hamiltonian([e for e in E_all if e not in within], "minus the 24 within-ring edges")
print("RESULT: positive control found a tour:", ok_pos, "| tour avoiding every within-ring move exists:", ok_neg)
