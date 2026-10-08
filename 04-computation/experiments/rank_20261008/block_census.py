#!/usr/bin/env python3
"""Census of block-balance certificates on Z_5: contracting maps T(x) = (m_i x + r_i)/5 with m_0 = 1, multipliers in
[1, 39] prime to 5, standard constants r_i = -m_i i mod 5, debt rank >= 3, and an even-order coupling group
(<m_i/m_j mod 5> = {1, 4} or all units), i.e. exactly the maps out of reach of one-step forms (THM-4609 (5)).
A seeded random sample from each (rank, |G|) class.  Version 2: for each map try the fixed block length k = 2 (all
2-step block matrices balanced by one exactly verified integer form); otherwise an adaptive certificate (block_adaptive:
state-dependent block lengths K0..KMAX).  (Version 1, fixed k <= 4 only, is block_census_fixed_k4.out.)
Usage: python3 block_census.py PER_CLASS K0 KMAX"""
import itertools, math, random, sys, time
import numpy as np
from block_balance import int_coords, ratio_group, std_r
from block_balance2 import integer_forms
from block_lmi2 import last_level_unique, active_set_opt, exact_check_vec
from block_adaptive import adaptive_certify

def certify(d, m, r, kmax):
    rho, v = int_coords(m)
    for k in range(2, kmax + 1):
        U, nstates, frozen = last_level_unique(d, m, r, v, rho, k)
        rk = np.linalg.matrix_rank(U.astype(float))
        if rk.min() < 3: continue
        Uf = U.astype(float); Uf /= np.trace(Uf, axis1=1, axis2=2)[:, None, None]
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            P, mg = active_set_opt(Uf, rho)
        if mg <= 0:
            last = (k, mg); continue
        Q = np.linalg.inv(P)
        for Qi in integer_forms(Q, rho):
            ok, _ = exact_check_vec(Qi, U)
            if ok: return k, mg, Qi, len(U), frozen
        last = (k, mg)
    return None, last[1] if 'last' in locals() else None, None, None, None

if __name__ == '__main__':
    per, K0, KMAX = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    d = 5
    cands = [x for x in range(1, 40) if x % 5]
    classes = {}
    for ms in itertools.product(cands, repeat=4):
        m = (1,) + ms
        if math.prod(m) >= 5 ** 5: continue
        G = ratio_group(list(m), 5)
        if len(G) == 1: continue
        rho, _ = int_coords(list(m))
        if rho < 3: continue
        classes.setdefault((rho, len(G)), []).append(m)
    rnd = random.Random(20261008)
    print("class sizes:", {k: len(v) for k, v in sorted(classes.items())}, flush=True)
    tally = {}
    for key in sorted(classes):
        sample = rnd.sample(classes[key], min(per, len(classes[key])))
        for m in sample:
            m = list(m); r = std_r(m)
            Lam = sum(math.log(x / d) for x in m) / d
            t0 = time.time()
            k, mg, Qi, nU, frozen = certify(d, m, r, 2)
            if k:
                status = f"CERTIFIED fixed k = 2, margin {mg:+.4f}, Q = {Qi}, {nU} blocks"; tag = 'k=2'
            else:
                res = adaptive_certify(d, m, r, K0, KMAX, 0.02, log=lambda x: None)
                if res:
                    Qi, nU, counts, mgl = res
                    status = f"CERTIFIED adaptive lengths {K0}..{KMAX}, nodes refined per level {counts}, leaf margin {mgl:+.4f}, Q = {Qi}, {nU} leaf blocks"; tag = f'adaptive {K0}..{KMAX}'
                else:
                    status = f"not certified (fixed k = 2; adaptive {K0}..{KMAX})"; tag = 'none'
            print(f"rank {key[0]} |G| = {key[1]}  m = {m} r = {r} Lambda = {Lam:+.3f}: {status}  [{time.time() - t0:.0f}s]", flush=True)
            tally.setdefault(key, []).append(tag)
    print("summary (class: outcomes):", {k: v for k, v in sorted(tally.items())}, flush=True)
