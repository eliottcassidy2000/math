#!/usr/bin/env python3
"""Independent audit (written from scratch) of Section 5 of
05-knowledge/results/chessboard_weave_20261006.md: Theorem 5.1 (weave law), Hall form, 6x6 analogue, 10x10.
Solver: OR-tools CP-SAT AddCircuit, 2 workers. Every tour returned is re-validated by plain code.
"""
import sys
import itertools
from collections import Counter
from ortools.sat.python import cp_model

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


KN = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def ring(i, j, N):
    return (max(abs(2 * i - (N - 1)), abs(2 * j - (N - 1))) - 1) // 2


def knight_edges(N):
    E = set()
    for i in range(N):
        for j in range(N):
            for dx, dy in KN:
                a, b = i + dx, j + dy
                if 0 <= a < N and 0 <= b < N:
                    E.add(tuple(sorted([(i, j), (a, b)])))
    return sorted(E)


def validate_tour(N, tour):
    if len(tour) != N * N or len(set(tour)) != N * N:
        return False
    for t in range(N * N):
        u, v = tour[t], tour[(t + 1) % (N * N)]
        if sorted([abs(u[0] - v[0]), abs(u[1] - v[1])]) != [1, 2]:
            return False
    return True


def profile(N, tour):
    c = Counter()
    for t in range(N * N):
        u, v = tour[t], tour[(t + 1) % (N * N)]
        c[tuple(sorted((ring(*u, N), ring(*v, N))))] += 1
    return c


def build(N, extra=None, objective=None, time_limit=120.0, seed=0):
    V = [(i, j) for i in range(N) for j in range(N)]
    idx = {v: k for k, v in enumerate(V)}
    E = knight_edges(N)
    m = cp_model.CpModel()
    arcs, x = [], {}
    for (u, v) in E:
        for a, b in ((u, v), (v, u)):
            lit = m.NewBoolVar(f"x{idx[a]}_{idx[b]}")
            x[(a, b)] = lit
            arcs.append((idx[a], idx[b], lit))
    m.AddCircuit(arcs)
    pair_terms = {}
    for (u, v) in E:
        key = tuple(sorted((ring(*u, N), ring(*v, N))))
        pair_terms.setdefault(key, []).extend([x[(u, v)], x[(v, u)]])
    msum = {k: sum(vs) for k, vs in pair_terms.items()}
    if extra:
        extra(m, msum)
    if objective:
        sense, expr = objective(msum)
        (m.Minimize if sense == 'min' else m.Maximize)(expr)
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 2
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.random_seed = seed
    st = solver.Solve(m)
    tour = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        succ = {}
        for (a, b), lit in x.items():
            if solver.Value(lit):
                succ[a] = b
        tour = [V[0]]
        while len(tour) < N * N:
            tour.append(succ[tour[-1]])
    return solver.StatusName(st), tour, (solver.ObjectiveValue() if objective and tour else None), (solver.BestObjectiveBound() if objective else None)


# ---------------- structural facts used in the proof (8x8) ----------------
N = 8
E8 = knight_edges(8)
R = {(i, j): ring(i, j, 8) for i in range(8) for j in range(8)}
S = set(v for v in R if R[v] in (0, 3))
Sc = set(v for v in R if R[v] in (1, 2))
check(len(S) == 32 and len(Sc) == 32, "|S| = |S^c| = 32")
inside_S = [e for e in E8 if e[0] in S and e[1] in S]
print("knight moves with both ends in ring0 u ring3:", inside_S)
check(len(inside_S) == 8 and all(R[e[0]] == 3 and R[e[1]] == 3 for e in inside_S), "only the 8 ring-3 moves inside S")
corners = [(0, 0), (0, 7), (7, 0), (7, 7)]
check(all(min(max(abs(u[0] - c[0]), abs(u[1] - c[1])) for c in corners) <= 2 and min(max(abs(v[0] - c[0]), abs(v[1] - c[1])) for c in corners) <= 2 for (u, v) in inside_S), "ring-3 inside moves sit at corners")
check(len(set((u[0] + u[1]) % 2 for u in S)) == 2, "S has both colours")
# k = 0 colour lock: S-S^c knight edges split the board into the two classes (colour+side) and never join them
lock = Counter(((u[0] + u[1] + (u in S)) % 2, (v[0] + v[1] + (v in S)) % 2) for (u, v) in E8 if (u in S) != (v in S))
print("S-S^c knight edges by (colour+side) class of endpoints:", dict(lock))
check(all(a == b for (a, b) in lock), "colour+side invariant on S-S^c edges")

# ---------------- Hall form ----------------
within = set(e for e in E8 if R[e[0]] == R[e[1]])
check(len(within) == 24, "24 within-ring moves")
Gp = [e for e in E8 if e not in within]
for col in (0, 1):
    X = set(v for v in S if (v[0] + v[1]) % 2 == col)
    NX = set()
    for (u, v) in Gp:
        if u in X:
            NX.add(v)
        if v in X:
            NX.add(u)
    target = set(v for v in Sc if (v[0] + v[1]) % 2 != col)
    print(f"Hall form colour {col}: |X| = {len(X)}, |N(X)| = {len(NX)}, N(X) == other colour in rings 1-2: {NX == target}")
    check(len(X) == 16 and NX == target, f"Hall form colour {col}")

# ---------------- each k in 0..8 (CP-SAT) ----------------
print("=== 8x8: feasibility of m_33 = k ===")
found = {}
for k in range(0, 9):
    st, tour, _, _ = build(8, extra=lambda m, ms, k=k: m.Add(ms.get((3, 3), 0) == k), time_limit=300.0, seed=k)
    if tour:
        ok = validate_tour(8, tour)
        pr = profile(8, tour)
        ident = pr[(3, 3)] == pr[(1, 1)] + pr[(1, 2)] + pr[(2, 2)]
        check(ok and pr[(3, 3)] == k and ident, f"tour with k={k} valid & identity")
        found[k] = tour
        print(f"k={k}: {st}; valid closed tour {ok}; profile {dict(sorted(pr.items()))}; identity m33=m11+m12+m22: {ident}")
    else:
        print(f"k={k}: {st}")
        check(k == 0 and st == 'INFEASIBLE', f"k={k} status")
print("tours (as square sequences) for k=1..8:")
for k in sorted(found):
    print(f"  k={k}:", ' '.join(f"{i}{j}" for (i, j) in found[k]))

# ---------------- identity on many diverse tours ----------------
print("=== identity on tours with random objectives ===")
import random
rng = random.Random(20261006)
kseen = Counter()
for trial in range(40):
    w = {key: rng.randint(-5, 5) for key in [(0, 1), (0, 2), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]}
    st, tour, _, _ = build(8, objective=lambda ms, w=w: ('max', sum(w[k] * ms.get(k, 0) for k in w)), time_limit=8.0, seed=trial)
    if tour:
        pr = profile(8, tour)
        check(validate_tour(8, tour) and pr[(3, 3)] == pr[(1, 1)] + pr[(1, 2)] + pr[(2, 2)] and pr[(3, 3)] >= 1, f"identity random tour {trial}")
        kseen[pr[(3, 3)]] += 1
print("k values seen on 40 objective-randomised tours:", dict(sorted(kseen.items())))

# ---------------- 6x6 ----------------
print("=== 6x6 ===")
E6 = knight_edges(6)
R6 = {(i, j): ring(i, j, 6) for i in range(6) for j in range(6)}
print("6x6 ring sizes:", dict(Counter(R6.values())), " knight pair profile:", dict(sorted(Counter(tuple(sorted((R6[u], R6[v]))) for u, v in E6).items())))
res6 = {}
for sense in ('min', 'max'):
    st, tour, val, bnd = build(6, objective=lambda ms, s=sense: (s, ms.get((2, 2), 0)), time_limit=300.0)
    pr = profile(6, tour)
    inner = pr[(0, 1)] + pr[(1, 1)] + pr[(0, 0)]
    check(validate_tour(6, tour) and pr[(2, 2)] == 4 + inner, f"6x6 identity {sense}")
    res6[sense] = (st, val, bnd, dict(sorted(pr.items())))
    print(f"6x6 {sense} #moves inside ring 2: {st} value {val} bound {bnd}; profile {dict(sorted(pr.items()))}; =4+inner: {pr[(2, 2)] == 4 + inner}")
for trial in range(20):
    w = {key: rng.randint(-5, 5) for key in [(0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]}
    st, tour, _, _ = build(6, objective=lambda ms, w=w: ('max', sum(w[k] * ms.get(k, 0) for k in w)), time_limit=8.0, seed=trial)
    if tour:
        pr = profile(6, tour)
        check(validate_tour(6, tour) and pr[(2, 2)] == 4 + pr[(0, 1)] + pr[(1, 1)], f"6x6 identity random {trial}")

# ---------------- 10x10 and general balanced splits ----------------
print("=== balanced ring splits on 2m x 2m boards ===")
for m in range(2, 9):
    sizes = [4 * (2 * k + 1) for k in range(m)]
    half = 2 * m * m
    splits = [c for r in range(1, m) for c in itertools.combinations(range(m), r) if sum(sizes[k] for k in c) == half and 0 in c]
    print(f"{2*m}x{2*m}: ring sizes {sizes}, half {half}, balanced splits (containing ring 0): {splits}")
    if m == 5:
        check(splits == [], "10x10 no balanced split")
print()
print("TOTAL FAILURES:", len(fails))
