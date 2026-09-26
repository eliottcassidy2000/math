#!/usr/bin/env python3
"""
procgen_tension_20260926_q3.py -- Q3: is rho_max(sigma) <= F_k for every class-(i) level-k strategy?
No.  Exact maxima M_k (k <= 6), value sets, the level-7 witness 17/27, the sibling-splitting bound.
CP-SAT (<= 2 workers) decides "is there a class-(i) level-k strategy with a cycle of density > f" and
"... with rho_max exactly f"; every positive answer is re-verified by exact Karp.  Run through the runner.
"""
import time
from fractions import Fraction

from ortools.sat.python import cp_model

from procgen_tension_20260926_lib import (
    below_c, best_lower, check, chi4_mask, claim, karp, karp_single, lift, potential, potential_ok,
    rho_all, simple_cycles_adj, succ)

# the level-7 strategy with rho_max = 17/27 found by CP-SAT (this lane, 2026-09-26); re-verified below
LEVEL7_17_27 = 12481914054424834826


def cpsat_query(k, f, exact, tl=1800, workers=2, seed=1):
    """exact=False: is there a class-(i) sigma with a cycle of density > f ?
    exact=True : is there a sigma with rho_max = f exactly (all cycles <= f, one cycle >= f; f < c) ?
    Class (i) / 'all cycles <= F' is encoded by an integer potential (Lemma P), the cycle by AddCircuit."""
    M = 1 << k
    H = M >> 1
    F = f if exact else best_lower(M)          # all cycles <= F  (F_{2^k} <=> class (i))
    q, r = F.numerator, F.denominator
    m = cp_model.CpModel()
    minus = {s: m.NewBoolVar("m%d" % s) for s in range(1, M, 2)}
    U = (M // 2) * (r - q) + 1
    psi = [m.NewIntVar(0, U, "p%d" % s) for s in range(M)]
    vis = [m.NewBoolVar("v%d" % s) for s in range(M)]
    arcs = []
    for s in range(M):
        if s % 2 == 0:
            b = (s // 2) % H
            cand = [(b, None), (b + H, None)]
        else:
            bp = ((3 * s + 1) // 2) % H
            bm = ((3 * s - 1) // 2) % H
            cand = [(bp, minus[s].Not()), (bp + H, minus[s].Not()), (bm, minus[s]), (bm + H, minus[s])]
        w = (r - q) if s % 2 else -q
        for t, lit in cand:
            c = m.Add(psi[t] <= psi[s] - w)
            if lit is not None:
                c.OnlyEnforceIf(lit)
            if t != s:
                x = m.NewBoolVar("")
                if lit is not None:
                    m.AddImplication(x, lit)
                arcs.append((s, t, x))
        arcs.append((s, s, vis[s].Not()))
    m.AddCircuit(arcs)
    a0, p0 = f.numerator, f.denominator
    lhs = sum(p0 * vis[s] for s in range(1, M, 2)) - sum(a0 * vis[s] for s in range(M))
    m.Add(lhs >= (0 if exact else 1))
    m.Add(sum(vis) >= 2)
    so = cp_model.CpSolver()
    so.parameters.num_search_workers = workers
    so.parameters.max_time_in_seconds = tl
    so.parameters.random_seed = seed
    t0 = time.time()
    st = so.Solve(m)
    dt = time.time() - t0
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        mask = sum(1 << ((s - 1) // 2) for s in range(1, M, 2) if so.Value(minus[s]))
        return 'SAT', mask, dt
    if st == cp_model.INFEASIBLE:
        return 'UNSAT', None, dt
    return 'UNKNOWN', None, dt


def cpsat_refined(k, f, tl=7200, workers=2, seed=1):
    """search for a level-k strategy with rho_max = f exactly (f = a/p < c in lowest terms, p <= 2^(k-1)) whose
    maximal cycle has length exactly p, a odd nodes and no two siblings.  Model: integer potential at f (all
    cycles <= f, hence class (i)) + such a circuit.  By Lemma S a strategy with rho_max = f has a sibling-free
    maximal cycle, of length a multiple of p and <= 2^(k-1); so when 2p > 2^(k-1) the query is COMPLETE (UNSAT
    proves that rho_max = f does not occur), and otherwise it is only a restricted search."""
    M = 1 << k
    H = M >> 1
    q, r = f.numerator, f.denominator
    check(below_c(f) and r <= H, "refined query needs f < c with denominator <= 2^(k-1)")
    m = cp_model.CpModel()
    minus = {s: m.NewBoolVar("m%d" % s) for s in range(1, M, 2)}
    U = (M // 2) * (r - q) + 1
    psi = [m.NewIntVar(0, U, "p%d" % s) for s in range(M)]
    vis = [m.NewBoolVar("v%d" % s) for s in range(M)]
    arcs = []
    for s in range(M):
        if s % 2 == 0:
            b = (s // 2) % H
            cand = [(b, None), (b + H, None)]
        else:
            bp = ((3 * s + 1) // 2) % H
            bm = ((3 * s - 1) // 2) % H
            cand = [(bp, minus[s].Not()), (bp + H, minus[s].Not()), (bm, minus[s]), (bm + H, minus[s])]
        w = (r - q) if s % 2 else -q
        for t, lit in cand:
            c = m.Add(psi[t] <= psi[s] - w)
            if lit is not None:
                c.OnlyEnforceIf(lit)
            if t != s:
                x = m.NewBoolVar("")
                if lit is not None:
                    m.AddImplication(x, lit)
                arcs.append((s, t, x))
        arcs.append((s, s, vis[s].Not()))
    m.AddCircuit(arcs)
    m.Add(sum(vis) == r)
    m.Add(sum(vis[s] for s in range(1, M, 2)) == q)
    for t in range(H):
        m.Add(vis[t] + vis[t + H] <= 1)
    so = cp_model.CpSolver()
    so.parameters.num_search_workers = workers
    so.parameters.max_time_in_seconds = tl
    so.parameters.random_seed = seed
    t0 = time.time()
    st = so.Solve(m)
    dt = time.time() - t0
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        mask = sum(1 << ((s - 1) // 2) for s in range(1, M, 2) if so.Value(minus[s]))
        return 'SAT', mask, dt
    if st == cp_model.INFEASIBLE:
        return 'UNSAT', None, dt
    return 'UNKNOWN', None, dt


def sibling_split_check(k, R):
    """every simple cycle of every class-(i) strategy that contains two siblings (t, t + 2^(k-1)) splits into two
    cycles on the same nodes (so a maximum-density cycle without siblings exists)"""
    H = 1 << (k - 1)
    cnt = 0
    for m in range(len(R)):
        if not below_c(R[m]):
            continue
        adj = {s: list(succ(k, m, s)) for s in range(1 << k)}
        for cyc in simple_cycles_adj(adj):
            pos = {v: i for i, v in enumerate(cyc)}
            p = len(cyc)
            for i, v in enumerate(cyc):
                j = pos.get(v ^ H)
                if j is None or j < i:
                    continue
                # chords: pred(v_i) -> v_j and pred(v_j) -> v_i
                check(cyc[j] in succ(k, m, cyc[i - 1]) and cyc[i] in succ(k, m, cyc[j - 1]), "sibling chords")
                c1 = cyc[i:j]
                c2 = cyc[j:] + cyc[:i]
                check(len(c1) + len(c2) == p, "split partitions the cycle")
                a1, a2 = sum(x & 1 for x in c1), sum(x & 1 for x in c2)
                check(a1 + a2 == sum(x & 1 for x in cyc), "odd counts add")
                cnt += 1
    return cnt


def run(level6=True, level7_bound=False, tl7=7200, level8_search=False, tl8=1800, level8_mask=None):
    print("Q3. Christoffel rigidity: is rho_max(sigma) <= F_k on class (i)?")
    # (1) exhaustive levels 2..5
    Vk = {}
    for k in (2, 3, 4, 5):
        R = rho_all(k)
        vals = {}
        for f in R:
            if below_c(f):
                vals[f] = vals.get(f, 0) + 1
        Vk[k] = vals
        Fk = best_lower(k)
        viol = sum(c for f, c in vals.items() if f > Fk)
        print("    level %d: F_k = %s; class-(i) rho_max values %s; M_k = %s; %d strategies exceed F_k"
              % (k, Fk, sorted(vals.items()), max(vals), viol))
    check(max(Vk[2]) == Fraction(1, 2) and max(Vk[3]) == Fraction(1, 2) and max(Vk[4]) == Fraction(3, 5)
          and max(Vk[5]) == Fraction(5, 8), "M_2..M_5")
    check(set(Vk[5]) == {Fraction(1, 2), Fraction(5, 9), Fraction(4, 7), Fraction(3, 5), Fraction(5, 8)}, "V_5")
    claim(True, "exhaustive (all 65,812 strategies of levels 2-5, exact Karp): M_2 = M_3 = 1/2, M_4 = 3/5 > F_4 = 1/2, M_5 = 5/8 > F_5 = 3/5; "
          "rho_max <= F_k fails for 8 strategies at level 4 and 220 at level 5 (REFUTED)")
    # (1') a refuted guess: 'rho_max has at most 2^(k-3) even nodes' holds for EVERY strategy of levels 3..5 ...
    for k in (3, 4, 5):
        R = rho_all(k)
        check(all(f.denominator - f.numerator <= 2 ** (k - 3) for f in R), "even-node bound at level %d" % k)
    claim(True, "every strategy of levels 3-5 (65,808) has rho_max = a/p with p - a <= 2^(k-3) (fails at level 6, below)")
    # (2) sibling splitting, checked on every simple cycle of every class-(i) strategy at levels 2..4
    tot = 0
    for k in (2, 3, 4):
        tot += sibling_split_check(k, rho_all(k))
    claim(True, "sibling splitting: in %d simple cycles of class-(i) strategies (k <= 4) containing a sibling pair, both "
          "chords exist and split the cycle into two cycles on the same nodes" % tot)
    # (3) the numerator bound G_(2^(k-2)) (Corollary S) and the exact maxima M_k = G_(2^(k-2)) for k <= 6
    def G(A):
        return max(Fraction(a, p) for a in range(1, A + 1) for p in range(a, 3 * a + 2) if below_c(Fraction(a, p)))
    Gk = {k: G(2 ** (k - 2)) for k in range(2, 10)}
    check([Gk[k] for k in range(2, 10)] == [Fraction(1, 2), Fraction(1, 2), Fraction(3, 5), Fraction(5, 8),
                                            Fraction(5, 8), Fraction(29, 46), Fraction(41, 65), Fraction(94, 149)],
          "numerator bound values")
    R4, R5 = rho_all(4), rho_all(5)
    m4 = min(m for m in range(len(R4)) if R4[m] == Fraction(3, 5))
    m5 = min(m for m in range(len(R5)) if R5[m] == Fraction(5, 8))
    for k, mk in ((2, chi4_mask(2)), (3, chi4_mask(3)), (4, m4), (5, m5), (6, lift(m5, 5, 6))):
        rho = karp(k, mk)
        check(below_c(rho) and rho == Gk[k], "witness attains the numerator bound at level %d" % k)
    claim(True, "numerator bound (Corollary S): G_(2^(k-2)) = 1/2, 1/2, 3/5, 5/8, 5/8, 29/46, 41/65, 94/149 for "
          "k = 2..9; it is attained for k = 2..6 by chi_(-4), chi_(-4), the level-4 optimum (mask %d), the level-5 "
          "optimum (mask %d) and its lift (exact Karp): M_k = G_(2^(k-2)) for k <= 6 (PROVED)" % (m4, m5))
    # (3') level 6 by CP-SAT: an independent (solver) confirmation of M_6 and the value set V_6
    if level6:
        st, mask, dt = cpsat_query(6, Fraction(5, 8), exact=False)
        check(st == 'UNSAT', "level 6: no class-(i) strategy has a cycle of density > 5/8")
        claim(True, "level 6: CP-SAT independently confirms that no class-(i) strategy has a cycle of density > 5/8 "
              "(%.0f s)" % dt)
        cands = sorted(set(Fraction(a, p) for p in range(2, 33) for a in range(p)
                           if Fraction(1, 2) <= Fraction(a, p) <= Fraction(5, 8)))
        V6 = []
        for f in cands:
            st, mask, dt = cpsat_query(6, f, exact=True, tl=1800)
            check(st in ('SAT', 'UNSAT'), "level-6 value query decided, f = %s" % f)
            if st == 'SAT':
                rho = karp(6, mask)
                check(rho == f and below_c(rho), "returned strategy has rho_max = %s exactly" % f)
                V6.append(f)
        print("    level 6 value set V_6 (%d of %d candidate fractions in [1/2, 5/8] with denominator <= 32): %s"
              % (len(V6), len(cands), [str(x) for x in V6]))
        check(Fraction(13, 22) in V6 and Fraction(9, 17) not in V6 and Fraction(11, 18) in V6, "V_6 landmarks")
        claim(True, "REFUTED guesses at level 6: 13/22 (9 even nodes > 2^(k-3) = 8) is a value of rho_max; the values "
              "11/18, 11/19, 13/21, 13/22 have denominators > 2^(k-2) + 1 = 17 while 9/17 is absent")
        claim(len(V6) == 16 and max(V6) == Fraction(5, 8),
              "level 6: the rho_max values of class-(i) strategies are exactly %d fractions, max denominator %d "
              "(each realized value re-verified by exact Karp on the returned strategy; each absent one solver-certified)"
              % (len(V6), max(x.denominator for x in V6)))
    # (4) level 7: the stored witness
    rho = karp_single(7, LEVEL7_17_27)
    check(rho == Fraction(17, 27) and below_c(rho), "level-7 witness has rho_max = 17/27")
    psi = potential(7, LEVEL7_17_27, Fraction(17, 27))
    check(psi is not None and potential_ok(7, LEVEL7_17_27, Fraction(17, 27), psi), "level-7 potential certificate")
    claim(True, "level 7: an explicit class-(i) strategy (mask %d) has rho_max = 17/27 exactly (exact Karp + edge-checked "
          "integer potential at 17/27), so M_7 >= 17/27 = F_27 while F_7 = 3/5" % LEVEL7_17_27)
    cands7 = sorted(set(Fraction(a, p) for p in range(2, 65) for a in range(p)
                        if Fraction(17, 27) < Fraction(a, p) and below_c(Fraction(a, p))))
    check(cands7 == [Fraction(29, 46)], "above 17/27 only 29/46 has denominator <= 64")
    claim(True, "sibling-splitting bound: a maximum-density cycle has length <= 2^(k-1) = 64 at level 7, and the only "
          "fraction in (17/27, log_3 2) with denominator <= 64 is 29/46; so M_7 is 17/27 or 29/46")
    if level7_bound:
        check(2 * 46 > 2 ** 6, "the refined level-7 query is complete (2p > 2^(k-1))")
        st, mask, dt = cpsat_refined(7, Fraction(29, 46), tl=tl7)
        if st == 'UNSAT':
            claim(True, "level 7: CP-SAT (refined model: a sibling-free maximal cycle of length 46 with 29 odd nodes) "
                  "proves that no class-(i) strategy has rho_max = 29/46 (%.0f s): M_7 = 17/27 (solver-certified)" % dt)
        elif st == 'SAT':
            rho = karp_single(7, mask)
            claim(below_c(rho) and rho > Fraction(17, 27), "level 7: a class-(i) strategy with rho_max = %s" % rho)
        else:
            print("    level 7 upper bound: CP-SAT undecided after %.0f s" % dt)
    # (5) level 8: a stored witness (if any) and/or a refined search at the numerator bound 41/65 and at 29/46
    if level8_mask is not None:
        rho = karp_single(8, level8_mask)
        psi = potential(8, level8_mask, rho)
        check(below_c(rho) and psi is not None and potential_ok(8, level8_mask, rho, psi), "level-8 witness")
        claim(rho > Fraction(17, 27), "level 8: an explicit class-(i) strategy has rho_max = %s (exact Karp + edge-checked "
              "potential), so the rigidity rho_max <= F_k fails for every k < %d" % (rho, rho.denominator))
    if level8_search:
        for f in (Fraction(41, 65), Fraction(29, 46)):
            st, mask, dt = cpsat_refined(8, f, tl=tl8)
            if st == 'SAT':
                rho = karp_single(8, mask)
                claim(below_c(rho) and rho == f, "level 8: CP-SAT finds a class-(i) strategy with rho_max = %s (mask %d), "
                      "re-verified by exact Karp" % (rho, mask))
                break
            complete = 2 * f.denominator > 2 ** 7
            print("    level 8, rho_max = %s: CP-SAT %s after %.0f s (%s)" % (f, st, dt,
                  "complete query" if complete else "restricted search: only cycles of length exactly %d" % f.denominator))
    return Vk


if __name__ == "__main__":
    run()
