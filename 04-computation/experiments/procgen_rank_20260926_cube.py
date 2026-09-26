#!/usr/bin/env python3
"""
procgen_rank_20260926_cube.py -- the strategy cube side of the rank lane.

  (a) census, levels 2..5 (exhaustive): class (i) (no expanding cycle), class (i') (every SCC that contains an
      expanding cycle is that single cycle, i.e. finitely many expanding periodic points), other
  (b) the mechanism of Theorem R+ on a non-Collatz strategy: the level-5 strategy with the isolated expanding
      2-cycle {1/5, -1/5}: charges c > chi on the cycle plus a periodic h descend inside the shadow; the entry
      step from the uncharged preimage 2/5 fails with slope c
  (c) CP-SAT (ortools, 2 workers): an isolated expanding loop at -1 exists from level 6 on, but never together
      with "all other cycles contracting" for levels 4..10 (class (i') via the -1 loop is empty there)
"""
import math
import time
from collections import deque
from fractions import Fraction

from procgen_rank_20260926_lib import (KAPPA, LN2, LN3, check, say, v2, good_integer, classify_strategy, cube_graph,
                                       sccs)


def sig_from_bits(k, bits):
    odds = list(range(1, 1 << k, 2))
    return {r: (-1 if (bits >> i) & 1 else 1) for i, r in enumerate(odds)}


def part_census():
    say("== cube.a census of classes (i) and (i') at levels 2..5 (exhaustive) ==")
    expected_i = {2: 1, 3: 1, 4: 16, 5: 1052}
    example = None
    for k in range(2, 6):
        counts = {"i": 0, "i'": 0, "other": 0}
        n_iso = 0
        for bits in range(1 << (1 << (k - 1))):
            sig = sig_from_bits(k, bits)
            c, iso = classify_strategy(k, sig)
            counts[c] += 1
            n_iso += len(iso)
            if iso and example is None and k == 5:
                example = (bits, iso[0])
        say("   level %d: class (i) %d, class (i') %d, other %d; isolated expanding SCCs %d" %
            (k, counts["i"], counts["i'"], counts["other"], n_iso))
        check(counts["i"] == expected_i[k] and counts["i'"] == 0,
              "level %d: %d class-(i) strategies (as THM-4474) and no class-(i') strategy" % (k, counts["i"]))
        if k == 5:
            check(n_iso == 4, "level 5: exactly 4 isolated expanding SCCs occur (all inside strategies that also "
                  "have a non-isolated expanding cycle)")
    return example


def periodic_point_of_cycle(k, sig, nodes):
    A, B, C = 1, 1, 0
    for s in nodes:
        if s % 2:
            A, C = 3 * A, 3 * C + sig[s] * B
        B *= 2
    return Fraction(C, B - A)


def T_sig_frac(x, sig, k):
    x = Fraction(x)
    if x.numerator % 2 == 0:
        return x / 2
    r = (x.numerator * pow(x.denominator, -1, 1 << k)) % (1 << k)
    return (3 * x + sig[r]) / 2


def T_sig_int(n, sig, k):
    return (3 * n + sig[n % (1 << k)]) >> 1 if n & 1 else n >> 1


def part_Rplus(example):
    say("== cube.b Theorem R+ mechanism on a non-Collatz strategy (level 5, isolated expanding 2-cycle) ==")
    k = 5
    bits, comp = example
    sig = sig_from_bits(k, bits)
    adj = cube_graph(k, sig)
    # order the cycle nodes along the edges
    start = comp[0]
    order = [start]
    while True:
        nxt = [u for u in adj[order[-1]] if u in comp]
        assert len(nxt) == 1
        if nxt[0] == start:
            break
        order.append(nxt[0])
    x0 = periodic_point_of_cycle(k, sig, order)
    pts = [x0]
    for _ in range(len(order) - 1):
        pts.append(T_sig_frac(pts[-1], sig, k))
    check(T_sig_frac(pts[-1], sig, k) == x0 and len(order) == 2,
          "strategy #%d at level 5 has the isolated expanding cycle %s with periodic orbit %s (odd density 1, "
          "multiplier 9/4)" % (bits, order, [str(p) for p in pts]))
    word = [1 if (Fraction(p).numerator % 2) else 0 for p in pts]
    p = len(word)
    chi = (sum(word) * LN3 - p * LN2) / p
    c = chi + 0.02
    tau = [sum(v2(pts[j] - pts[i]) for i in range(p) if i != j) for j in range(p)]
    mubar = chi - c
    g = [0.0]
    for j in range(p - 1):
        w = KAPPA if word[j] else -LN2
        g.append(g[-1] + mubar - w + c)
    K = max(v2(pts[i] - pts[j]) for i in range(p) for j in range(p) if i != j) + 2
    K = max(K, k + 1)
    hmap = {}
    for j in range(p):
        r = (pts[j].numerator * pow(pts[j].denominator, -1, 1 << K)) % (1 << K)
        hmap[r] = -c * tau[j] + g[j]

    def Rk(n):
        s = math.log(n) + hmap.get(n % (1 << K), 0.0)
        for q in pts:
            s += c * v2(Fraction(n) - q)
        return s
    worst = -1e9
    steps = 0
    for D in (K + 40, K + 120):
        y = good_integer(pts[0], D)
        a = y
        for _ in range(D - K - 1):
            b = T_sig_int(a, sig, k)
            worst = max(worst, Rk(b) - Rk(a))
            a = b
            steps += 1
    check(worst <= mubar / 2, "inside the shadow of the cycle (%d steps): R(T n) - R(n) <= %.4f <= mubar/2 = %.4f "
          "with charge c = chi + 0.02 = %.4f on the two cycle points" % (steps, worst, mubar / 2, c))
    pre = 2 * pts[0]
    incs = []
    for D in (60, 120, 240):
        y = good_integer(pre, D)
        incs.append(Rk(T_sig_int(y, sig, k)) - Rk(y))
    check(incs[0] > 0 and abs((incs[1] - incs[0]) - c * 60) < 1e-6 and abs((incs[2] - incs[1]) - c * 120) < 1e-6,
          "the entry step from the uncharged even preimage %s increases R by %.2f, %.2f, %.2f at depths 60, 120, "
          "240 (slope c): charging the backward orbit is forced here too" % (pre, incs[0], incs[1], incs[2]))


def cpsat_loop_models(k, with_potential, tlimit=300):
    """sigma(-1) = +; the exit 2^(k-1)-1 of the loop cannot reach -1 (closed set); optionally every cycle
    avoiding -1 has density < c (integer potential at F = best lower approximation with denominator <= 2^k)"""
    from ortools.sat.python import cp_model
    M = 1 << k
    H = M >> 1
    C = math.log(2) / math.log(3)
    F = Fraction(0)
    for d in range(1, M + 1):
        a = math.floor(C * d)
        if 3 ** a < 2 ** d and Fraction(a, d) > F:
            F = Fraction(a, d)
    q, r = F.numerator, F.denominator
    m = cp_model.CpModel()
    sig = {s: m.NewBoolVar("p%d" % s) for s in range(1, M, 2)}
    m.Add(sig[M - 1] == 1)
    reach = [m.NewBoolVar("r%d" % v) for v in range(M)]
    m.Add(reach[H - 1] == 1)
    m.Add(reach[M - 1] == 0)
    psi = [m.NewIntVar(0, 4 * M * r, "psi%d" % v) for v in range(M)] if with_potential else None
    for s in range(M):
        if s % 2 == 0:
            opts = [(None, (s // 2) % H)]
        else:
            opts = [(sig[s], ((3 * s + 1) // 2) % H), (sig[s].Not(), ((3 * s - 1) // 2) % H)]
        wF = (r - q) if s % 2 else -q
        for lit, t in opts:
            for u in (t, t + H):
                if lit is None:
                    m.AddImplication(reach[s], reach[u])
                else:
                    m.AddBoolOr([reach[s].Not(), lit.Not(), reach[u]])
                if with_potential and s != M - 1 and u != M - 1:
                    if lit is None:
                        m.Add(psi[u] <= psi[s] - wF)
                    else:
                        m.Add(psi[u] <= psi[s] - wF).OnlyEnforceIf(lit)
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 2
    solver.parameters.max_time_in_seconds = tlimit
    st = solver.Solve(m)
    name = solver.StatusName(st)
    sol = None
    if name in ("OPTIMAL", "FEASIBLE"):
        sol = {s: (1 if solver.Value(sig[s]) else -1) for s in range(1, M, 2)}
    return name, sol


def part_cpsat():
    say("== cube.c the loop at -1: isolation is possible, isolation with a contracting rest is not (k <= 10) ==")
    t0 = time.time()
    name, sol = cpsat_loop_models(6, False)
    ok = False
    if sol is not None:
        adj = cube_graph(6, sol)
        seen = {31}
        dq = deque([31])
        while dq:
            v = dq.popleft()
            for u in adj[v]:
                if u not in seen:
                    seen.add(u)
                    dq.append(u)
        ok = 63 not in seen
        comps = [c for c in sccs(adj) if 63 in c]
        ok &= comps == [[63]]
    check(name in ("OPTIMAL", "FEASIBLE") and ok,
          "level 6: CP-SAT finds a strategy whose loop at -1 is an isolated SCC (re-verified by BFS: -1 is not "
          "reachable from its exit 31)")
    res = []
    for k in range(4, 11):
        nm, _ = cpsat_loop_models(k, True)
        res.append((k, nm))
    check(all(nm == "INFEASIBLE" for _, nm in res),
          "levels 4..10: CP-SAT proves INFEASIBLE 'loop at -1 isolated and every other cycle contracting' (%s; %.1f s)"
          % (", ".join("k=%d %s" % t for t in res), time.time() - t0))


def run():
    ex = part_census()
    part_Rplus(ex)
    part_cpsat()


if __name__ == "__main__":
    run()
