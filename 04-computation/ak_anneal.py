#!/usr/bin/env python3
"""
ak_anneal.py — klein-S691: focused annealer for strict-mode AK certificates.

Starts from the run-2 record (13/7 on [2,2,2]) and cold starts; anneals
labels/wildcards with greedy-pruned seed completion; exact min-seed
branch-and-bound squeeze on improvements.  Targets the rung ladder
11/6, 9/5, 7/4, 5/3.  All records re-verified from scratch (exact
rational closure) before printing.
"""
import random, time, sys
from fractions import Fraction
from itertools import product, combinations
from ak_forcing_engine import AKInstance, force
from ak_strict_search import greedy_seeds, XFULL, POOL3, POOL5, X0

RUNG = [Fraction(13, 7), Fraction(11, 6), Fraction(9, 5),
        Fraction(7, 4), Fraction(5, 3)]


def try_force(inst):
    ok, T, _ = force(inst, "strict")
    return ok, T


def minseed_exact(dims, fs, T0, pool, ub, time_cap=25.0):
    """Branch & bound: minimal #seeds (< ub) making (dims, fs, T0) force.
    Branch on the first unfired vertex' seed options plus 'wildcardize'
    nothing (T0 fixed).  Returns best R or None."""
    t0 = time.time()
    best = [None]

    def rec(R):
        if time.time() - t0 > time_cap:
            return
        if best[0] is not None and len(R) >= len(best[0]):
            return
        inst = AKInstance(XFULL, dims, fs, T0, R)
        ok, T = try_force(inst)
        if ok:
            best[0] = list(R)
            return
        if len(R) + 1 >= ub:
            return
        unf = [v for v in inst.verts if inst.vid[v] not in T]
        # branch on seeds at the lattice-most-constrained vertex: heuristic =
        # first unfired; try each pool label
        v = unf[0]
        for lab in pool:
            rec(R + [(v, lab)])
        # also allow branching at a second unfired vertex to escape local traps
        if len(unf) > 1:
            v2 = unf[len(unf) // 2]
            if v2 != v:
                for lab in pool:
                    rec(R + [(v2, lab)])

    rec([])
    return best[0]


def anneal(dims, steps=1500, seed=0, start=None, T0_max=2, log=print):
    rng = random.Random(seed)
    k = len(dims)
    verts = list(product(*[range(1, d + 1) for d in dims]))

    def rand_fs():
        fs = []
        for i0 in range(k):
            dom = list(product(*([range(1, d + 1) for d in dims[:i0]]
                                 + [range(1, dims[i0])])))
            f = {}
            for pre in dom:
                if rng.random() < 0.6:
                    f[pre] = rng.choice(POOL5)
            fs.append(f)
        return fs

    if start:
        fs, T0 = [dict(f) for f in start[0]], list(start[1])
    else:
        fs, T0 = rand_fs(), []
    cur = greedy_seeds(dims, fs, T0, POOL3)
    cur_s = cur.score() if cur else Fraction(10)
    best, best_s = cur, cur_s
    temp0 = 0.25
    for step in range(steps):
        temp = temp0 * (1 - step / steps) + 0.01
        nfs = [dict(f) for f in fs]
        nT0 = list(T0)
        move = rng.random()
        if move < 0.55:
            i0 = rng.randrange(k)
            dom = list(product(*([range(1, d + 1) for d in dims[:i0]]
                                 + [range(1, dims[i0])])))
            if not dom:
                continue
            pre = rng.choice(dom)
            opts = [X0] + POOL5
            cv = nfs[i0].get(pre, X0)
            nv = rng.choice([o for o in opts if o != cv])
            if nv == X0:
                nfs[i0].pop(pre, None)
            else:
                nfs[i0][pre] = nv
        elif move < 0.8 and len(nT0) < T0_max:
            v = rng.choice(verts)
            if v not in nT0:
                nT0.append(v)
        elif nT0:
            nT0.pop(rng.randrange(len(nT0)))
        cand = greedy_seeds(dims, nfs, nT0, POOL3)
        if cand is None:
            continue
        s = cand.score()
        if s is None:
            continue
        import math
        if s <= cur_s or rng.random() < math.exp(float(cur_s - s) / temp):
            fs, T0, cur, cur_s = nfs, nT0, cand, s
            if s < best_s:
                # exact seed squeeze
                mR = minseed_exact(dims, fs, T0, POOL3, ub=len(cand.R))
                if mR is not None and len(mR) < len(cand.R):
                    cand = AKInstance(XFULL, dims, fs, T0, mR)
                    s = cand.score()
                # fresh re-verification
                chk = AKInstance(XFULL, dims, fs, T0, cand.R)
                ok, _ = try_force(chk)
                assert ok
                best, best_s = cand, s
                log(f"  [{'-'.join(map(str,dims))} seed{seed} step{step}] "
                    f"NEW BEST {s} = {float(s):.4f}  "
                    f"m={cand.m()} r={cand.r()} n={cand.n} t={cand.t()}")
                for rung in RUNG:
                    if s <= rung:
                        log(f"    RUNG {rung} REACHED")
                        log(f"    fs={cand.fs} T0={cand.T0} R={cand.R}")
                        break
    return best, best_s


if __name__ == "__main__":
    t0 = time.time()
    RECORD = ([{(1,): (1, 2)}, {(1, 1): (1, 0), (2, 1): (1, 0)},
               {(2, 1, 1): (2, 1), (2, 2, 1): (0, 1)}], [(1, 1, 2)])
    overall = (Fraction(10), None)
    plans = [([2, 2, 2], 1200, RECORD), ([2, 2, 2], 1200, None),
             ([2, 2, 2, 2], 700, None), ([3, 3], 900, None),
             ([2, 2, 2, 2], 700, None), ([2, 3, 2], 800, None)]
    for i, (dims, steps, start) in enumerate(plans):
        print(f"== anneal {dims} steps={steps} start={'record' if start else 'cold'}")
        b, s = anneal(dims, steps=steps, seed=100 + i, start=start)
        if b is not None and s < overall[0]:
            overall = (s, b)
        print(f"   -> best {s} ({float(s):.4f})  [{time.time()-t0:.0f}s]")
    s, b = overall
    print("\n==== OVERALL BEST ====")
    if b is not None:
        print(f"score={s} = {float(s):.4f} m={b.m()} r={b.r()} n={b.n} t={b.t()}")
        print("dims=", b.dims)
        print("fs  =", b.fs)
        print("T0  =", b.T0)
        print("R   =", b.R)
