#!/usr/bin/env python3
"""audit H: THM-4611 statement 5 (census on Z_5): independent class sizes, property (i) for every sampled map, and an
independent check of every certificate printed in block_census.out (adaptive: greedy exact refinement with the printed
form, levels 3..5; fixed k = 2: every 2-step block exactly balanced by the printed form).
Usage: python3 i_census_check.py"""
import itertools, math, re, sys, time
import numpy as np
from hcore import ratio_group, primes_of, expvec, lattice_rank, factor
from hdp import Tables, balance_exact

def class_sizes():
    cands = [x for x in range(1, 40) if x % 5]
    ev = {}
    pr = sorted({p for x in cands for p in factor(x)})
    for x in cands:
        f = factor(x); ev[x] = tuple(f.get(p, 0) for p in pr)
    sizes = {}
    for ms in itertools.product(cands, repeat=4):
        if math.prod(ms) >= 3125: continue
        m = (1,) + ms
        G = ratio_group(list(m), 5)
        if len(G) == 1: continue
        rk = int(np.linalg.matrix_rank(np.array([ev[x] for x in ms], dtype=float)))
        if rk < 3: continue
        sizes[(rk, len(G))] = sizes.get((rk, len(G)), 0) + 1
    return dict(sorted(sizes.items()))

def parse(path):
    out = []
    for line in open(path):
        mm = re.search(r"m = \[([^\]]*)\]", line)
        if not mm or 'CERTIFIED' not in line: continue
        m = [int(x) for x in mm.group(1).split(',')]
        Q = eval(re.search(r"Q = (\[\[.*?\]\])", line).group(1))
        kind = 'fixed2' if 'fixed k = 2' in line else 'adaptive'
        out.append((m, Q, kind))
    return out

def prop_i(m, d=5):
    G = sorted(ratio_group(m, d))
    return not any(all(m[(a * j + b) % d] == m[j] for j in range(d)) for a in G for b in range(d) if (a, b) != (1, 0))

def fixed_check(m, Q, k=2):
    tb = Tables(5, m)
    T = tb.table(k).reshape(-1, tb.rho, tb.rho)
    nz = np.any(T.reshape(len(T), -1) != 0, axis=1)
    U = np.unique(T[nz].reshape(len(T[nz]), -1), axis=0).reshape(-1, tb.rho, tb.rho)
    return bool(balance_exact(Q, U).all()), len(U), int((~nz).sum())

if __name__ == '__main__':
    t0 = time.time()
    print("independent class sizes (rank, |G|):", class_sizes(), f"[{time.time() - t0:.1f}s]", flush=True)
    import h_greedy
    for m, Q, kind in parse('../block_census.out'):
        pi = prop_i(m)
        indep = None
        if kind == 'fixed2':
            ok, nU, nz0 = fixed_check(m, Q)
            print(f"m = {m}: property (i) {pi}; fixed k = 2 form {Q}: all {nU} nonzero blocks exactly balanced: {ok} (zero states {nz0})", flush=True)
        else:
            h_greedy.MAPS['tmp'] = (5, m, Q, 3, 5)
            print(f"m = {m}: property (i) {pi}; adaptive 3..5, form {Q}:", flush=True)
            h_greedy.run('tmp')
