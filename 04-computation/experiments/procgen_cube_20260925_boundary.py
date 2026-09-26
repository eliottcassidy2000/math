#!/usr/bin/env python3
"""
procgen_cube_20260925_boundary.py -- how far is Collatz from the provable strategies?

delta_k := least number of odd residues mod 2^k at which sigma must differ from Collatz (all +)
           for T_sigma to be in class (i) (every cycle of G_sigma contracting).
Haar distance = delta_k / 2^(k-1).

Method: lazy-constraint CP-SAT.  Variables x_r (flip residue r), objective min sum x_r.  Each
iteration solves to optimality, builds G_sigma, and asks the C engine ('expcyc' mode) for expanding
simple cycles; every returned cycle is re-verified exactly in Python (closed walk of G_sigma, 3^a > 2^p)
and turned into the valid no-good clause  OR_{odd r on the cycle} (x_r != current x_r)
(the cycle exists whenever sigma agrees with the current one on its odd nodes).  The clauses are
necessary conditions, so the CP-SAT optimum is a lower bound; the first optimum without expanding
cycles is therefore optimal (solver-certified), and it is re-verified with the exact descent
certificate + finite check of procgen_cube_20260925_core.py.

Called from procgen_cube_20260925_run.py (section H3) or standalone:  python3 ... 3 4 5 6 7 8
"""
import os
import sys
import time
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_cube_20260925_core as core  # noqa: E402

try:
    from ortools.sat.python import cp_model
except ImportError:            # pragma: no cover
    cp_model = None


def engine_path():
    root = os.path.abspath(os.path.join(HERE, '..', '..'))
    return os.path.join(root, 'scratch', 'procgen_cube', 'strategy_cube', 'engine')


def expanding_cycles(k, mask, M=40):
    out = subprocess.run([engine_path(), 'expcyc', str(k), str(M)], input=format(mask, 'x') + '\n',
                         capture_output=True, text=True, check=True).stdout
    cycles = []
    for line in out.splitlines():
        if line.startswith('C '):
            cycles.append([int(v) for v in line[2:].split(',')])
    # exact verification
    succ = core.graph(k, mask)
    for cyc in cycles:
        p = len(cyc)
        for i in range(p):
            if cyc[(i + 1) % p] not in succ[cyc[i]]:
                raise SystemExit("CHECK FAILED: engine cycle is not a cycle of G_sigma")
        a = sum(1 for v in cyc if v & 1)
        if not 3 ** a > 2 ** p:
            raise SystemExit("CHECK FAILED: engine cycle not expanding")
    return cycles


def min_flips_to_I(k, timeout=3000, workers=2, verbose=False):
    res = list(range(1, 1 << k, 2))
    model = cp_model.CpModel()
    x = {r: model.NewBoolVar(f"x{r}") for r in res}
    model.Minimize(sum(x.values()))
    t0 = time.time()
    it = ncl = 0
    while True:
        it += 1
        solver = cp_model.CpSolver()
        solver.parameters.num_workers = workers
        solver.parameters.max_time_in_seconds = max(10.0, timeout - (time.time() - t0))
        st = solver.Solve(model)
        if st != cp_model.OPTIMAL:
            return {'status': 'timeout', 'k': k, 'iterations': it, 'clauses': ncl, 'lb': int(solver.BestObjectiveBound())}
        flips = [r for r in res if solver.Value(x[r])]
        mask = 0
        for r in flips:
            mask |= 1 << ((r - 1) // 2)
        cycs = expanding_cycles(k, mask)
        if not cycs:
            return {'status': 'optimal', 'k': k, 'delta': len(flips), 'flips': flips, 'mask': mask,
                    'iterations': it, 'clauses': ncl, 'seconds': time.time() - t0}
        for cyc in cycs:
            lits = [x[v].Not() if solver.Value(x[v]) else x[v] for v in cyc if v & 1]
            model.AddBoolOr(lits)
            ncl += 1
        if verbose:
            print(f"   k={k} it {it}: lower bound {len(flips)}, clauses {ncl}, {time.time() - t0:.0f} s", flush=True)
        if time.time() - t0 > timeout:
            return {'status': 'timeout', 'k': k, 'iterations': it, 'clauses': ncl, 'lb': len(flips)}


def min_flips_to_I_rc2(k, M=40, timeout=3600, verbose=False):
    """Same lazy scheme with the exact MaxSAT solver RC2 (pysat): soft unit clauses (not x_r), weight 1."""
    from pysat.formula import WCNF
    from pysat.examples.rc2 import RC2
    res = list(range(1, 1 << k, 2))
    var = {r: i + 1 for i, r in enumerate(res)}
    hard = []
    t0 = time.time()
    it = 0
    while True:
        it += 1
        w = WCNF()
        for c in hard:
            w.append(c)
        for r in res:
            w.append([-var[r]], weight=1)
        with RC2(w) as rc2:
            model = rc2.compute()
        flips = [r for r in res if model[var[r] - 1] > 0]
        mask = 0
        for r in flips:
            mask |= 1 << ((r - 1) // 2)
        cycs = expanding_cycles(k, mask, M)
        if verbose:
            print(f"   k={k} it {it}: lower bound {len(flips)}, clauses {len(hard)}, {time.time() - t0:.0f} s", flush=True)
        if not cycs:
            return {'status': 'optimal', 'k': k, 'delta': len(flips), 'flips': flips, 'mask': mask,
                    'iterations': it, 'clauses': len(hard), 'seconds': time.time() - t0}
        fl = set(flips)
        for cyc in cycs:
            hard.append([(-var[v] if v in fl else var[v]) for v in cyc if v & 1])
        if time.time() - t0 > timeout:
            return {'status': 'timeout', 'k': k, 'iterations': it, 'clauses': len(hard), 'lb': len(flips)}


def cycles_in_subset(k, mask, Y, sign, M=40):
    text = format(mask, 'x') + '\n' + ''.join('1' if v in Y else '0' for v in range(1 << k)) + '\n'
    out = subprocess.run([engine_path(), 'cycsub', str(k), str(M), str(sign)], input=text,
                         capture_output=True, text=True, check=True).stdout
    cycles = [[int(v) for v in line[2:].split(',')] for line in out.splitlines() if line.startswith('C ')]
    succ = core.graph(k, mask)
    for cyc in cycles:
        p = len(cyc)
        for i in range(p):
            if cyc[(i + 1) % p] not in succ[cyc[i]] or cyc[i] not in Y:
                raise SystemExit("CHECK FAILED: engine subset cycle invalid")
        a = sum(1 for v in cyc if v & 1)
        if (3 ** a > 2 ** p) != (sign > 0):
            raise SystemExit("CHECK FAILED: engine subset cycle has the wrong sign")
    return cycles


def min_flips_to_II_rc2(k, M=40, timeout=3600, verbose=False):
    """delta^(ii)_k: least number of flips from Collatz for a closed node set Y all of whose cycles expand.
    Variables x_r (flip), y_s (s in Y); hard: Y closed under the sigma-successors, Y nonempty;
    lazy no-goods: a contracting cycle inside Y under the current flips."""
    from pysat.formula import WCNF
    from pysat.examples.rc2 import RC2
    K = 1 << k
    H = K >> 1
    res = list(range(1, K, 2))
    xv = {r: i + 1 for i, r in enumerate(res)}
    yv = {s: len(res) + 1 + s for s in range(K)}
    base = []
    for s in range(K):
        if s % 2 == 0:
            for u in (s // 2 % H, s // 2 % H + H):
                base.append([-yv[s], yv[u]])
        else:
            tp = ((3 * s + 1) // 2) % H
            tm = ((3 * s - 1) // 2) % H
            for u in (tp, tp + H):
                base.append([-yv[s], xv[s], yv[u]])
            for u in (tm, tm + H):
                base.append([-yv[s], -xv[s], yv[u]])
    base.append([yv[s] for s in range(K)])
    lazy = []
    t0 = time.time()
    it = 0
    while True:
        it += 1
        w = WCNF()
        for c in base + lazy:
            w.append(c)
        for r in res:
            w.append([-xv[r]], weight=1)
        with RC2(w) as rc2:
            model = rc2.compute()
        if model is None:
            return {'status': 'infeasible', 'k': k}
        pos = set(v for v in model if v > 0)
        flips = [r for r in res if xv[r] in pos]
        Y = set(s for s in range(K) if yv[s] in pos)
        mask = 0
        for r in flips:
            mask |= 1 << ((r - 1) // 2)
        cycs = cycles_in_subset(k, mask, Y, -1, M)
        if verbose:
            print(f"   (ii) k={k} it {it}: lower bound {len(flips)}, |Y| {len(Y)}, clauses {len(lazy)}, {time.time() - t0:.0f} s", flush=True)
        if not cycs:
            return {'status': 'optimal', 'k': k, 'delta': len(flips), 'flips': flips, 'mask': mask, 'Y': sorted(Y),
                    'iterations': it, 'clauses': len(lazy), 'seconds': time.time() - t0}
        fl = set(flips)
        for cyc in cycs:
            cl = [-yv[v] for v in cyc]
            cl += [(-xv[v] if v in fl else xv[v]) for v in cyc if v & 1]
            lazy.append(cl)
        if time.time() - t0 > timeout:
            return {'status': 'timeout', 'k': k, 'iterations': it, 'lb': len(flips)}


def verify_I(k, mask):
    succ = core.graph(k, mask)
    rmax = core.karp_density(succ, range(1 << k), True)
    ok = core.density_vs_critical(rmax) < 0
    cert = core.descent_certificate(k, mask, 400) if ok else None
    cycles = sorted(core.cycles_below(k, mask, cert[1])) if cert else None
    return ok, rmax, cert, cycles


if __name__ == '__main__':
    which = sys.argv[1]
    for k in map(int, sys.argv[2:]):
        if which == 'i':
            r = min_flips_to_I_rc2(k, verbose=True)
            print(r, flush=True)
            if r['status'] == 'optimal':
                print('  verify:', verify_I(k, r['mask'])[:3], flush=True)
        else:
            r = min_flips_to_II_rc2(k, verbose=True)
            print({kk: v for kk, v in r.items() if kk != 'Y'}, flush=True)
