#!/usr/bin/env python3
"""audit H, task C: independent re-verification of THM-4611 statement 3 (fixed-length certificates).

For each map: branch condition, Lambda, rank, independence, coupling group, property (i), minimal one-step coupling rank,
then ALL k-step block matrices A_k(h), h in H_k x Z/d^k (own DP, level k in chunks of M), with
  - zero blocks listed, distinct nonzero blocks counted (exact dedup), minimal exact rank,
  - exact balance by the stated integer form (charpoly-sign PD test of tr(A adj Q) Q - 2 det(Q) A, Python ints),
  - the float margin min (tr - 2 lmax)/tr and alpha_max = min tr/lmax - 2 of the stated form.
Usage: python3 c_certificates.py NAME [NAME ...]   (names below)"""
import math, sys, time, itertools
import numpy as np
from hdp import Tables, balance_exact, margins_float, exact_rank_vec
from hcore import ratio_group, det_int

CERTS = {
    'Z5_123711': (5, [1, 2, 3, 7, 11], 2, [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]),
    'Z5_12371': (5, [1, 2, 3, 7, 1], 4, [[20, -4, -8], [-4, 15, -3], [-8, -3, 20]]),
    'Z5_11237': (5, [1, 1, 2, 3, 7], 5, [[8, -2, -3], [-2, 6, -1], [-3, -1, 8]]),
    'Z7_1111235': (7, [1, 1, 1, 1, 2, 3, 5], 4, [[6, -1, -1], [-1, 5, -1], [-1, -1, 6]]),
    'Z7_11123511': (7, [1, 1, 1, 2, 3, 5, 11], 2, [[5, -1, -2, -1], [-1, 5, 0, -1], [-2, 0, 6, -1], [-1, -1, -1, 5]]),
}

def rank_float(A):
    return int(np.linalg.matrix_rank(np.array(A, dtype=float))) if np.any(A) else 0

def one_step_info(tb):
    d = tb.d
    G = sorted(ratio_group(tb.m, d))
    ranks = {}
    zero_nonid = []
    for a in G:
        for b in range(d):
            if (a, b) == (1, 0):
                continue
            Dm = tb.Dab[a, b]
            rk = rank_float(Dm)
            ranks[(a, b)] = rk
            if rk == 0:
                zero_nonid.append((a, b))
    # property (i) directly from multipliers as well
    prop_i = not any(all(tb.m[(a * j + b) % d] == tb.m[j] for j in range(d)) for a in G for b in range(d) if (a, b) != (1, 0))
    return G, min(ranks.values()), [k for k, v in ranks.items() if v == min(ranks.values())], zero_nonid, prop_i

def independence(m):
    # independent configuration: values equal on a set Z (the most common value), all others multiplicatively independent
    from hcore import primes_of, expvec, lattice_rank
    d = len(m)
    vals = {}
    for x in m:
        vals[x] = vals.get(x, 0) + 1
    c = max(vals, key=lambda x: vals[x])
    others = [x for x in m if x != c]
    pr = primes_of(m)
    vecs = [tuple(a - b for a, b in zip(expvec(x, pr), expvec(c, pr))) for x in others]
    return lattice_rank(vecs) == len(others), c, len(others)

def run(name, chunk=40, maxchunks=None):
    d, m, k, Q = CERTS[name]
    t0 = time.time()
    tb = Tables(d, m)
    Lam = sum(math.log(x / d) for x in m) / d
    indep, c, nonunit = independence(m)
    G, min1, argmin1, zero_nonid, prop_i = one_step_info(tb)
    print(f"{name}: d = {d}, m = {m}, r = {tb.r}, residues {[x % d for x in m]}, Lambda = {Lam:+.4f}, rank {tb.rho}, "
          f"independent configuration {indep} (repeated value {c}, {nonunit} others), coords den {tb.den}", flush=True)
    print(f"   coupling group {G}; property (i) {prop_i} (zero-root non-identity couplings {zero_nonid}); "
          f"min one-step coupling rank {min1} at {argmin1[:6]}{'...' if len(argmin1) > 6 else ''}", flush=True)
    detQ = det_int(Q)
    print(f"   form Q = {Q}, det {detQ}, leading minors {[det_int([row[:t] for row in Q[:t]]) for t in range(1, len(Q) + 1)]}", flush=True)
    for s in range(1, k):
        tb.table(s)
    H = tb.Hs(k)
    E = np.arange(d ** k, dtype=np.int64)
    nstates = 0; zero_states = []; keys = []; allbal = True; nbad = 0
    minmg = np.inf; minal = np.inf; minrank = 99
    rho = tb.rho
    iu = np.triu_indices(rho)
    for c0 in range(0, len(H) if maxchunks is None else min(len(H), maxchunks * chunk), chunk):
        Hc = H[c0:c0 + chunk]
        Mv = np.repeat(Hc, d ** k); Ev = np.tile(E, len(Hc))
        A = tb.blocks_at(k, Mv, Ev)
        nstates += len(A)
        flat = A.reshape(len(A), -1)
        nz = np.any(flat != 0, axis=1)
        for ix in np.nonzero(~nz)[0]:
            zero_states.append((int(Mv[ix]), int(Ev[ix])))
        An = A[nz]
        # dedup within chunk on the upper triangle
        U = np.unique(An[:, iu[0], iu[1]], axis=0)
        assert np.abs(U).max() < 2 ** 15
        keys.append(U.astype(np.int16))
        # rebuild symmetric blocks
        B = np.zeros((len(U), rho, rho), dtype=np.int64)
        B[:, iu[0], iu[1]] = U
        B[:, iu[1], iu[0]] = U
        ok = balance_exact(Q, B)
        if not ok.all():
            allbal = False; nbad += int((~ok).sum())
        mg, al = margins_float(Q, B)
        minmg = min(minmg, float(mg.min())); minal = min(minal, float(al.min()))
        if rho == 3:
            minrank = min(minrank, int(exact_rank_vec(B).min()))
        else:
            minrank = min(minrank, min(rank_float(x) for x in B))
    allU = np.unique(np.concatenate(keys), axis=0)
    print(f"   k = {k}: {nstates} states (|H_k| = {len(H)}), zero blocks at {zero_states[:5]}{'...' if len(zero_states) > 5 else ''} "
          f"({len(zero_states)} zero), {len(allU)} distinct nonzero blocks, min exact rank {minrank}", flush=True)
    print(f"   EXACT balance by Q of every nonzero block: {allbal} ({nbad} failures); float margin of this Q {minmg:+.5f}, "
          f"alpha_max of this Q {minal:.4f}  [{time.time() - t0:.1f}s]", flush=True)

if __name__ == '__main__':
    args = sys.argv[1:]
    mc = None
    if args and args[0].startswith('--maxchunks='):
        mc = int(args[0].split('=')[1]); args = args[1:]
    for nm in args:
        run(nm, maxchunks=mc)
