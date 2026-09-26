#!/usr/bin/env python3
"""
procgen_drift_20260926_minmax.py -- the min-max cycle density rho*(q,k) = min over level-k sign strategies of the
maximum odd density of a cycle of G_sigma (exact, with certificates); drift lane 2026-09-26.

Class (i) is nonempty at level k iff rho*(q,k) < log_q 2 (THM-4474 Theorem A).  rho*(q,k) is non-increasing in k
(lifting keeps the map).

Method (exact): descent on the value.
  * Start from a strategy, compute its exact rho_max F0 (Dinkelbach iteration with a verified cycle of density F0
    and a verified integer potential at threshold F0).
  * Ask for a strategy all of whose cycles have density < F0, i.e. density <= F' = the largest fraction < F0 with
    denominator <= 2^k (every simple cycle has length <= 2^k): SAT with lazily generated no-goods.  A no-good is an
    expanding-at-threshold cycle (density >= F0) with the signs on its odd nodes; any strategy agreeing with those
    signs keeps the cycle.  No-goods for threshold F0 stay valid for every smaller threshold.
  * UNSAT  =>  rho*(q,k) = F0 (the no-goods are the certificate);  SAT  =>  new strategy with smaller rho_max.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_drift_20260926_lib as L    # noqa: E402


def farey_pred(x, D):
    """largest fraction < x with denominator <= D (x = a/b in lowest terms, b <= D, x > 0)"""
    a, b = x.numerator, x.denominator
    L.check(b <= D and a > 0, "farey_pred domain")
    if b == 1:
        # x = a/1: predecessor is (a*D - 1)/D
        return Fraction(a * D - 1, D)
    p1 = pow(a, -1, b)                   # a p1 = 1 mod b
    p1 += ((D - p1) // b) * b             # largest p1' = p1 mod b with p1' <= D
    a1 = (a * p1 - 1) // b
    L.check(a * p1 - b * a1 == 1 and 0 < p1 <= D and p1 + b > D, "farey_pred certificate")
    return Fraction(a1, p1)


def greedy_flips(k, q):
    """max-halving strategy: sigma(r) = sign maximizing v_2(q r + sigma) (ties impossible: exactly one of q r +- 1 is
    divisible by 4); returned as the flip set relative to all-plus"""
    return set(r for r in range(1, 1 << k, 2) if (q * r + 1) % 4 != 0)


def harvest(k, q, flipset, Fp, M=60, P=None, cap=4000):
    """no-goods (cycle, minus-odd-set) for cycles of G_sigma with density > Fp (engine at threshold Fp)"""
    fl = L.flip_array(k, flipset)
    E = L.Engine(k, q, fl, Fp)
    out = []
    for _ in range(M):
        if E.solve():
            break
        cyc = E.cycle()
        a, p, _ = L.verify_expanding_walk(k, q, fl, cyc)
        L.check(Fraction(a, p) > Fp, "harvested cycle below threshold")
        out.append(cyc)
        E.kill(cyc)
    if P:
        E2 = L.Engine(k, q, fl, Fp)
        for cyc in E2.short_cycles(P, cap, needflip=False):
            a, p, _ = L.verify_expanding_walk(k, q, fl, cyc)
            L.check(Fraction(a, p) > Fp, "short cycle below threshold")
            out.append(cyc)
    return [(c, frozenset(v for v in c if (v & 1) and fl[v])) for c in out]


def cycle_ok(k, q, cyc, mo, F0):
    """re-check a stored no-good: closed walk under the recorded signs with density >= F0"""
    N = 1 << k
    H = N >> 1
    p = len(cyc)
    for i in range(p):
        s, t = cyc[i], cyc[(i + 1) % p]
        u = (q * s + (-1 if s in mo else 1)) // 2 if s & 1 else s // 2
        if t % H != u % H:
            return False
    a = sum(1 for v in cyc if v & 1)
    return Fraction(a, p) >= F0


def minmax_density(k, q, time_limit=3600.0, log=None, P=None):
    from pysat.solvers import Glucose4
    t0 = time.time()
    D = 1 << k
    odd = list(range(1, D, 2))
    var = {v: i + 1 for i, v in enumerate(odd)}
    cur = greedy_flips(k, q)
    F0, cyc0, _ = L.rho_max(k, q, L.flip_array(k, cur))
    best = (F0, sorted(cur))
    nogoods = []                      # (cycle, minus-odd set, threshold F0 at which it was generated)
    solver = Glucose4()
    seen = set()
    P = P if P is not None else k + 3
    rounds = 0
    while True:
        if time.time() - t0 > time_limit:
            solver.delete()
            return {'status': 'timeout', 'best': best, 'nogoods': nogoods}
        if F0 == 0:
            solver.delete()
            return {'status': 'optimal', 'rho': F0, 'best': best, 'nogoods': nogoods}
        Fp = farey_pred(F0, D)
        # SAT loop at threshold F0 (forbid cycles with density >= F0)
        while True:
            rounds += 1
            ng = harvest(k, q, cur, Fp, P=P)
            added = 0
            for cyc, mo in ng:
                lits = tuple(sorted((var[v], 0 if v in mo else 1) for v in set(u for u in cyc if u & 1)))
                if lits in seen:
                    continue
                seen.add(lits)
                nogoods.append((cyc, mo, F0))
                solver.add_clause([vv if val == 1 else -vv for vv, val in lits])
                added += 1
            if not ng:
                break                      # cur has all cycles < F0
            if not solver.solve():
                solver.delete()
                # UNSAT: every strategy keeps a cycle of density >= F0
                return {'status': 'optimal', 'rho': F0, 'best': best, 'nogoods': nogoods,
                        'seconds': time.time() - t0, 'rounds': rounds}
            m = solver.get_model()
            val = {abs(l): l > 0 for l in m}          # variables in no clause are absent: treat as False
            cur = set(v for v in odd if val.get(var[v], False))
        F1, cyc1, _ = L.rho_max(k, q, L.flip_array(k, cur))
        L.check(F1 < F0, "descent did not decrease rho_max")
        F0 = F1
        best = (F0, sorted(cur))
        if log:
            log(f"    q={q} k={k}: strategy with rho_max {F0} = {float(F0):.5f} found ({len(nogoods)} no-goods, "
                f"{time.time() - t0:.0f}s)")


def verify_minmax(k, q, res):
    """re-check: best strategy has exact rho_max = rho (Dinkelbach, verified certificates) and the stored no-goods
    (all genuine cycles with density >= rho under the recorded signs) are jointly unsatisfiable (Glucose4)."""
    from pysat.solvers import Glucose4
    rho = res['rho']
    F, _, _ = L.rho_max(k, q, L.flip_array(k, set(res['best'][1])))
    L.check(F == rho, "best strategy rho_max mismatch")
    odd = list(range(1, 1 << k, 2))
    var = {v: i + 1 for i, v in enumerate(odd)}
    with Glucose4() as s:
        for cyc, mo, F0 in res['nogoods']:
            L.check(F0 >= rho, "no-good generated below the final value")
            L.check(cycle_ok(k, q, cyc, mo, rho), "stored no-good is not a cycle of density >= rho")
            s.add_clause([var[v] if v not in mo else -var[v] for v in set(u for u in cyc if u & 1)])
        L.check(not s.solve(), "no-goods are satisfiable: rho* not certified")
    return True


if __name__ == '__main__':
    q = int(sys.argv[1])
    for k in [int(a) for a in sys.argv[2:]]:
        res = minmax_density(k, q, log=lambda *a: print(*a, flush=True))
        if res['status'] == 'optimal':
            verify_minmax(k, q, res)
            c = L.log_c(q)
            print(f"q={q} k={k}: rho* = {res['rho']} = {float(res['rho']):.5f} (log_q 2 = {c:.5f}; class (i) "
                  f"{'NONEMPTY' if q ** res['rho'].numerator < 2 ** res['rho'].denominator else 'EMPTY'}); "
                  f"{len(res['nogoods'])} no-goods, {res.get('seconds', 0):.0f}s", flush=True)
        else:
            print(f"q={q} k={k}: timeout; best {res['best'][0]}", flush=True)
