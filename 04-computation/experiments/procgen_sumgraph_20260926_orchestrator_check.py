#!/usr/bin/env python3
"""Orchestrator audit of lane `sumgraph` (sum graphs as unions of reflection matchings), written from
the note's statements; the lane's solver/theory code was not read.

G_S(n): vertices 1..n, x ~ y (x != y) iff x + y in S. M_t = {{x, t-x}}: a matching.
  1. Theorem A2 (two targets): brute force, n = 3..60, all pairs: a Hamiltonian path exactly for
     {n,n+1}, {n+1,n+2} (gap 1) and {n-1,n+1}, {n,n+2}, {n+1,n+3} with odd smaller target (gap 2); never a cycle.
  2. Theorem A3 (three targets): brute force, n = 3..24, all triples: Hamiltonian path iff (i) max degree <= 2,
     (ii) |M_a|+|M_b|+|M_c| = n-1, (iii) gcd(b-a, c-b) = 1, or 2 with a odd.
  3. Corollary A3''' (Pythagorean zigzags): every primitive triple s^2+t^2=u^2 (s<t, t <= 60): the three squares
     give a Hamiltonian path of [n] exactly for n in {t^2-1, t^2} (and 17 for (3,4,5)) among n in [t^2-3, t^2+3],
     with ends {s^2, t^2/2} (t even) or {s^2/2, s^2} (t odd) at n = t^2-1. Corollary A3'': three consecutive
     squares give a path only for j = 4, n in {15,16,17} (j <= 25, n <= 700).
  4. Lemma T: every target (2-powers or 3-powers) in (max(P,Q), 2n-1] is 2P or 3Q (n <= 5000).
  5. Propositions W1, W2, W3: the certificates (leaves, saturated neighbours, chokes, disjoint resolving sets),
     re-derived from neighbour lists, at every n of the ranges for a = 5, 7, 10 and at the ends plus 200
     random n of the ranges for a = 12; the rho-criteria for nonemptiness for a <= 200.
  6. Spot check of W_5 = [162, 243] by OR-tools CP-SAT (AddCircuit with a dummy node): PATH at 179, 180, 194, 224, 243 and NONE at 195, 200, 206, 207, 223.
"""
import math, random
from math import gcd


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


def matching(t, n):
    return [(x, t - x) for x in range(max(1, t - n), min(n, t - 1) + 1) if x < t - x]


def is_ham_path(n, edges):
    if len(edges) != n - 1:
        return False
    deg = [0] * (n + 1)
    par = list(range(n + 1))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (x, y) in edges:
        deg[x] += 1
        deg[y] += 1
        if deg[x] > 2 or deg[y] > 2:
            return False
        rx, ry = f(x), f(y)
        if rx == ry:
            return False
        par[rx] = ry
    return True


def ends_of(n, edges):
    deg = [0] * (n + 1)
    for (x, y) in edges:
        deg[x] += 1
        deg[y] += 1
    return sorted(v for v in range(1, n + 1) if deg[v] == 1)


# ---------------------------------------------------------------- 1. two targets
cnt = 0
for n in range(3, 61):
    for s in range(3, 2 * n):
        for t in range(s + 1, 2 * n):
            E = matching(s, n) + matching(t, n)
            # no cycle: union-find
            par = list(range(n + 1))

            def f(x):
                while par[x] != x:
                    par[x] = par[par[x]]
                    x = par[x]
                return x
            cyc = False
            for (x, y) in E:
                rx, ry = f(x), f(y)
                if rx == ry:
                    cyc = True
                par[rx] = ry
            assert not cyc, (n, s, t)
            ham = is_ham_path(n, E)
            d = t - s
            pred = (d == 1 and s in (n, n + 1)) or (d == 2 and s % 2 == 1 and s in (n - 1, n, n + 1))
            assert ham == pred, (n, s, t, ham, pred)
            cnt += ham
check(True, f"Theorem A2: two-target unions are acyclic and Hamiltonian exactly as predicted for 3 <= n <= 60 ({cnt} Hamiltonian pairs)")

# ---------------------------------------------------------------- 2. three targets
cnt = tot = 0
for n in range(3, 25):
    for a in range(3, 2 * n):
        for b in range(a + 1, 2 * n):
            for c in range(b + 1, 2 * n):
                E = matching(a, n) + matching(b, n) + matching(c, n)
                deg = [0] * (n + 1)
                for (x, y) in E:
                    deg[x] += 1
                    deg[y] += 1
                c1 = max(deg) <= 2
                c2 = len(E) == n - 1
                g = gcd(b - a, c - b)
                c3 = g == 1 or (g == 2 and a % 2 == 1)
                ham = is_ham_path(n, E)
                assert ham == (c1 and c2 and c3), (n, a, b, c)
                cnt += ham
                tot += 1
check(True, f"Theorem A3: over all {tot} triples with 3 <= n <= 24, the union is a Hamiltonian path iff (i) max degree <= 2, (ii) n-1 edges, (iii) gcd condition ({cnt} paths)")

# ---------------------------------------------------------------- 3. Pythagorean zigzags
trip = []
for m in range(2, 40):
    for k in range(1, m):
        if (m - k) % 2 == 1 and gcd(m, k) == 1:
            x, y, z = m * m - k * k, 2 * m * k, m * m + k * k
            s, t = min(x, y), max(x, y)
            if t <= 60:
                trip.append((s, t, z))
trip.sort()
for (s, t, u) in trip:
    A, Bq, C = s * s, t * t, u * u
    for n in range(t * t - 3, t * t + 4):
        E = matching(A, n) + matching(Bq, n) + matching(C, n)
        ham = is_ham_path(n, E)
        pred = n in (t * t - 1, t * t) or ((s, t, u) == (3, 4, 5) and n == 17)
        assert ham == pred, (s, t, u, n, ham)
    n = t * t - 1
    E = matching(A, n) + matching(Bq, n) + matching(C, n)
    en = ends_of(n, E)
    want = sorted([s * s, t * t // 2]) if t % 2 == 0 else sorted([s * s // 2, s * s])
    assert en == want, (s, t, u, en, want)
check(True, f"Corollary A3''': all {len(trip)} primitive triples with t <= 60: the squares s^2, t^2, u^2 give a Hamiltonian path of [n] exactly at n = t^2-1, t^2 (and 17 for (3,4,5)) within [t^2-3, t^2+3], with the stated ends; (3,4,5) gives Q_15 with ends 9 and 8")
bad = []
for j in range(2, 26):
    for n in range(3, 701):
        E = matching((j - 1) ** 2, n) + matching(j * j, n) + matching((j + 1) ** 2, n)
        if (j - 1) ** 2 >= 3 and is_ham_path(n, E):
            bad.append((j, n))
check(bad == [(4, 15), (4, 16), (4, 17)], f"Corollary A3'': three consecutive squares give a Hamiltonian union only at {bad} (j <= 25, n <= 700)")

# ---------------------------------------------------------------- 4. Lemma T
T23 = sorted(set([2 ** i for i in range(1, 40)] + [3 ** i for i in range(1, 26)]))
for n in range(3, 5001):
    P = 1 << (n.bit_length() - 1)
    Q = 1
    while Q * 3 <= n:
        Q *= 3
    M = max(P, Q)
    tops = [t for t in T23 if M < t <= 2 * n - 1]
    assert set(tops) <= {2 * P, 3 * Q}, (n, tops)
check(True, "Lemma T: for 3 <= n <= 5000 every target in (max(P,Q), 2n-1] is 2P or 3Q")


# ---------------------------------------------------------------- 5. W1, W2, W3
def nbrs(x, n):
    out = []
    for t in T23 + [2 ** i for i in range(40, 70)]:
        if x < t <= x + n and t != 2 * x and 1 <= t - x <= n:
            out.append(t - x)
    return sorted(set(out))


def W_ranges(a):
    B = 3 ** (a - 1)
    P = 1 << (B.bit_length())          # the power of two in (B, 2B)
    if not (B < P < 2 * B):
        P //= 2
    assert B < P < 2 * B
    F = lambda num, den: Fraction(num, den)
    from fractions import Fraction
    Pf, Bf = Fraction(P), Fraction(B)
    w1 = (max(5 * Pf / 4, 3 * Bf - 3 * Pf / 4), min(3 * Pf / 2, 3 * Bf - Pf / 2))
    w2 = (max(5 * Pf / 4, 3 * Bf - 3 * Pf / 4, 2 * Pf - Bf + Pf / 4, 2 * Bf + Pf / 4), min(7 * Pf / 4, 3 * Bf - Pf / 4))
    w3 = (max(7 * Pf / 4, 3 * Bf - 3 * Pf / 8, 2 * Bf + Pf / 4, 2 * Pf - Bf + Pf / 4), min(15 * Pf / 8, 3 * Bf - Pf / 4))
    return B, P, w1, w2, w3


def int_range(r):
    lo, hi = r
    return range(math.ceil(lo), math.ceil(hi))       # lo <= n < hi


def rigid_top(v, n, M, ends):
    return v > M and v not in ends and len(nbrs(v, n)) == 2


def certify_W1(n, B, P):
    M = max(P, B)
    if nbrs(P, n) != [3 * B - P]:
        return False
    if nbrs(P // 2, n) != [B - P // 2]:
        return False
    ends = {P, P // 2}
    if nbrs(P // 4, n) != sorted([3 * P // 4, B - P // 4]):
        return False
    # saturated: the two named partners of 3P/4 are rigid top vertices (other partners do not matter)
    sat = [5 * P // 4, 3 * B - 3 * P // 4]
    nb = nbrs(3 * P // 4, n)
    return all(v in nb and rigid_top(v, n, M, ends) for v in sat)


def certify_W2(n, B, P):
    M = max(P, B)
    if nbrs(P, n) != [3 * B - P]:
        return False
    ends = {P}
    if nbrs(P // 4, n) != sorted([3 * P // 4, B - P // 4]):
        return False
    s1 = [5 * P // 4, 3 * B - 3 * P // 4]
    s2 = [2 * P - B + P // 4, 2 * B + P // 4]
    nb1, nb2 = nbrs(3 * P // 4, n), nbrs(B - P // 4, n)
    ok1 = all(v in nb1 and rigid_top(v, n, M, ends) for v in s1)
    ok2 = all(v in nb2 and rigid_top(v, n, M, ends) for v in s2)
    return ok1 and ok2 and len(set(s1 + s2)) == 4


def certify_W3(n, B, P):
    """chokes at P/4 and P/8 with disjoint resolving sets (one forced end P)."""
    M = max(P, B)
    if nbrs(P, n) != [3 * B - P]:
        return False
    ends = {P}

    def resolving(d):
        """choke at d (need 2): rigid partners provide an edge; a partner saturated away from d (two rigid
        neighbours other than d) is unusable unless the free end is one of its saturators; other partners are
        usable. If at most one edge is available, the free end lies in R(d) = saturators (+ d itself if it has
        an available edge)."""
        R = set()
        usable = 0
        for v in nbrs(d, n):
            if rigid_top(v, n, M, ends):
                usable += 1
                continue
            sat = [w for w in nbrs(v, n) if w != d and rigid_top(w, n, M, ends)]
            if len(sat) >= (1 if v in ends else 2):
                R.update(sat)
            else:
                usable += 1
        if usable >= 2:
            return None
        if usable == 1:
            R.add(d)
        return R
    R4 = resolving(P // 4)
    R8 = resolving(P // 8)
    if R4 is None or R8 is None:
        return False
    return R4.isdisjoint(R8)


from fractions import Fraction
tested = {}
for a in (5, 7, 10, 12):
    B, P, w1, w2, w3 = W_ranges(a)
    for name, rng, cert in (("W1", w1, certify_W1), ("W2", w2, certify_W2), ("W3", w3, certify_W3)):
        ns = list(int_range(rng))
        if not ns:
            continue
        if a == 12:
            random.seed(a)
            ns = sorted(set([ns[0], ns[-1]] + random.sample(ns, min(200, len(ns)))))
        okc = all(cert(n, B, P) for n in ns)
        assert okc, (a, name)
        tested[(a, name)] = (int_range(rng).start, int_range(rng).stop - 1, len(ns))
check(True, "W1-W3 certificates re-derived from neighbour lists: " + "; ".join(f"a={a} {nm} [{lo}, {hi}] ({k} n)" for (a, nm), (lo, hi, k) in sorted(tested.items())))
# rho criteria for a <= 200
for a in range(3, 201):
    B = 3 ** (a - 1)
    P = 1 << (B.bit_length())
    if not (B < P < 2 * B):
        P //= 2
    rho = Fraction(P, B)
    _, _, w1, w2, w3 = W_ranges(a)
    assert (w1[0] < w1[1]) == (Fraction(4, 3) < rho < Fraction(12, 7)), a
    assert (w2[0] < w2[1]) == (Fraction(4, 3) < rho < Fraction(8, 5)), a
    assert (w3[0] < w3[1]) == (Fraction(4, 3) < rho < Fraction(3, 2)), a
check(True, "W1/W2/W3 ranges are nonempty exactly for rho in (4/3,12/7), (4/3,8/5), (4/3,3/2), all levels 3 <= a <= 200")


# ---------------------------------------------------------------- 6. spot check Hamiltonian n
def cpsat_ham_path(n, time_limit=120.0):
    """Hamiltonian path of C_n as a Hamiltonian circuit through an extra node 0 (OR-tools CP-SAT AddCircuit,
    2 workers). Returns ('PATH', path) with the path verified edge by edge, or ('NONE', None) if proved
    infeasible, or ('UNKNOWN', None)."""
    from ortools.sat.python import cp_model
    adj = {x: nbrs(x, n) for x in range(1, n + 1)}
    m = cp_model.CpModel()
    arcs, lit = [], {}
    for x in range(1, n + 1):
        for y in adj[x]:
            v = m.NewBoolVar("")
            lit[(x, y)] = v
            arcs.append((x, y, v))
        a0, b0 = m.NewBoolVar(""), m.NewBoolVar("")
        lit[(0, x)], lit[(x, 0)] = a0, b0
        arcs.append((0, x, a0))
        arcs.append((x, 0, b0))
    m.AddCircuit(arcs)
    sv = cp_model.CpSolver()
    sv.parameters.num_search_workers = 2
    sv.parameters.max_time_in_seconds = time_limit
    st = sv.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        nxt = {}
        for (x, y), v in lit.items():
            if sv.Value(v):
                nxt[x] = y
        path = []
        cur = nxt[0]
        while cur != 0:
            path.append(cur)
            cur = nxt[cur]
        assert len(path) == n and len(set(path)) == n and all(path[i + 1] in adj[path[i]] for i in range(n - 1))
        return "PATH", path
    if st == cp_model.INFEASIBLE:
        return "NONE", None
    return "UNKNOWN", None


import sys
sys.setrecursionlimit(10000)
res = {}
expect = {179: "PATH", 180: "PATH", 194: "PATH", 195: "NONE", 200: "NONE", 206: "NONE", 207: "NONE", 223: "NONE", 224: "PATH", 243: "PATH"}
for n in expect:
    res[n] = cpsat_ham_path(n)[0]
check(all(res[n] == expect[n] for n in expect),
      "CP-SAT (independent of both lanes' custom solvers; paths verified edge by edge) agrees with the Hamiltonian set of W_5 at n = " +
      ", ".join(f"{n}:{res[n]}" for n in expect))
