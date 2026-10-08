#!/usr/bin/env python3
"""Large block-balance certificates: chunked DP at the last level, active-set ellipsoid optimum, rigorous vectorized
exact check.  Same criterion as block_balance.py / block_lmi.py.

Exact check of S' = tr(A adj Q) Q - 2 det(Q) A > 0 (Sylvester) for every distinct block matrix A (rho = 3 or 4):
entries of S' are computed exactly in int64 (bounds asserted); leading minors are evaluated in float64 together with an
explicit rounding-error bound (sum of absolute values of the expansion terms times 8 rho! eps); a minor whose float value
does not clear its bound is recomputed in exact Python integers.  So every verdict is exact.
Usage: python3 block_lmi2.py MAPNAME k"""
import itertools, math, sys, time
import numpy as np
from block_balance import int_coords, ratio_group, coupling_D, std_r
from block_balance2 import block_arrays, integer_forms, adj3, det_int
from block_lmi import ellipsoid_opt, psd_sqrt

MAPS = {
    'Z5_12371': (5, [1, 2, 3, 7, 1], [0, 3, 4, 4, 1]),
    'Z5_11237': (5, [1, 1, 2, 3, 7], std_r([1, 1, 2, 3, 7])),
    'Z5_123711': (5, [1, 2, 3, 7, 11], std_r([1, 2, 3, 7, 11])),
    'Z7_1111235': (7, [1, 1, 1, 1, 2, 3, 5], std_r([1, 1, 1, 1, 2, 3, 5])),
    'Z7_11123511': (7, [1, 1, 1, 2, 3, 5, 11], std_r([1, 1, 1, 2, 3, 5, 11])),
}

def last_level_unique(d, m, r, v, rho, k, chunk=64):
    """Distinct nonzero k-step block matrices, computing level k in chunks of M-residues."""
    H1, A1 = block_arrays(d, m, r, v, rho, k - 1)
    D = np.zeros((d, d, rho, rho), dtype=np.int64)
    for a in range(1, d):
        if math.gcd(a, d) != 1: continue
        for b in range(d): D[a, b] = np.array(coupling_D(v, a, b, d, rho), dtype=np.int64)
    mod = d ** k; pmod = d ** (k - 1)
    H = np.array(ratio_group(m, mod), dtype=np.int64); E = np.arange(mod, dtype=np.int64)
    pos = -np.ones(pmod, dtype=np.int64); pos[H1] = np.arange(len(H1))
    mm = np.array(m, dtype=np.int64); rr = np.array(r, dtype=np.int64)
    inv = np.array([pow(x, -1, mod) for x in m], dtype=np.int64)
    uniq = []; frozen = 0; nstates = 0
    for c0 in range(0, len(H), chunk):
        Hc = H[c0:c0 + chunk]
        Mg = Hc[:, None] * np.ones((1, mod), dtype=np.int64); Eg = np.ones((len(Hc), 1), dtype=np.int64) * E[None, :]
        a = Mg % d; b = Eg % d
        A = D[a, b] * pmod
        for j in range(d):
            i = (a * j + b) % d
            t = (mm[i] * inv[j]) % mod
            t = (t * Mg) % mod
            t = (t * rr[j]) % mod
            N = (mm[i] * Eg + rr[i] - t) % mod
            assert np.all(N % d == 0)
            e2 = (N // d) % pmod
            M2 = (((mm[i] * inv[j]) % mod) * Mg) % pmod
            idx = pos[M2]; assert np.all(idx >= 0)
            A += A1[idx, e2]
        flat = A.reshape(-1, rho * rho); nstates += len(flat)
        nz = np.any(flat != 0, axis=1); frozen += int((~nz).sum())
        uniq.append(np.unique(flat[nz], axis=0))
    U = np.unique(np.concatenate(uniq), axis=0).reshape(-1, rho, rho)
    return U, nstates, frozen

def margins_P(P, Ah, As):
    X = Ah @ P @ Ah
    w = np.linalg.eigvalsh(X)
    tr = np.einsum('nij,ij->n', As, P)
    return 1 - 2 * w[:, -1] / tr

def active_set_opt(Uf, rho, rounds=8, act=40000, iters=700):
    Ah_all = psd_sqrt(Uf)
    sel = np.random.default_rng(1).choice(len(Uf), size=min(act, len(Uf)), replace=False)
    best = None
    for rd in range(rounds):
        P, mg = ellipsoid_opt(Uf[sel], iters=iters)
        full = np.concatenate([margins_P(P, Ah_all[s:s + 500000], Uf[s:s + 500000]) for s in range(0, len(Uf), 500000)])
        fm = full.min()
        print(f"   round {rd}: active {len(sel)}, active-set optimum {mg:+.5f}, full minimum {fm:+.5f}", flush=True)
        if best is None or fm > best[1]: best = (P, fm)
        if fm >= mg - 1e-6: break
        worst = np.argsort(full)[:act // 2]
        sel = np.unique(np.concatenate([sel, worst]))
    return best

def det3_float(S):
    a, b, c = S[:, 0, 0], S[:, 0, 1], S[:, 0, 2]; dd, e, f = S[:, 1, 0], S[:, 1, 1], S[:, 1, 2]; g, h, i = S[:, 2, 0], S[:, 2, 1], S[:, 2, 2]
    terms = [a * e * i, -a * f * h, -b * dd * i, b * f * g, c * dd * h, -c * e * g]
    val = sum(terms); bound = sum(np.abs(t) for t in terms) * 8 * 6 * np.finfo(float).eps
    return val, bound

def exact_check_vec(Qi, U):
    rho = len(Qi)
    if not all(det_int([row[:t] for row in Qi[:t]]) > 0 for t in range(1, rho + 1)): return False, 0   # form not positive definite
    adj = np.array(adj3(Qi), dtype=np.int64); dq = det_int(Qi); Qa = np.array(Qi, dtype=np.int64)
    assert np.abs(U).max() < 2 ** 20 and np.abs(adj).max() < 2 ** 20
    t = np.einsum('nxy,yx->n', U, adj)
    S = t[:, None, None] * Qa[None] - 2 * dq * U
    assert np.abs(S).max() < 2 ** 50
    Sf = S.astype(float)
    ok = (S[:, 0, 0] > 0)
    m2 = S[:, 0, 0].astype(object) * S[:, 1, 1].astype(object) - S[:, 0, 1].astype(object) * S[:, 1, 0].astype(object)
    ok &= np.array([x > 0 for x in m2], dtype=bool)
    if rho == 3:
        val, bound = det3_float(Sf)
        sure_pos = val > bound; sure_neg = val < -bound
        unsure = ~(sure_pos | sure_neg)
        res = sure_pos.copy()
        for n in np.nonzero(unsure)[0]:
            res[n] = det_int(S[n].tolist()) > 0
        ok &= res
        return bool(ok.all()), int(unsure.sum())
    # rho = 4: exact integer minors for the 3x3 and 4x4 leading blocks (object arithmetic)
    for n in range(len(S)):
        if not ok[n]: continue
        Sn = S[n].tolist()
        if not (det_int([row[:3] for row in Sn[:3]]) > 0 and det_int(Sn) > 0): ok[n] = False
    return bool(ok.all()), 0

if __name__ == '__main__':
    name = sys.argv[1]; k = int(sys.argv[2])
    d, m, r = MAPS[name]
    for i in range(d): assert (m[i] * i + r[i]) % d == 0
    rho, v = int_coords(m)
    Lam = sum(math.log(x / d) for x in m) / d
    t0 = time.time()
    U, nstates, frozen = last_level_unique(d, m, r, v, rho, k)
    print(f"{name} m = {m} r = {r} rank {rho} Lambda = {Lam:+.4f}; k = {k}: {nstates} states, {frozen} frozen, "
          f"{len(U)} distinct nonzero blocks  [{time.time() - t0:.1f}s]", flush=True)
    rk = np.linalg.matrix_rank(U.astype(float))
    if rk.min() < 3:
        print(f"   impossible: a block of rank {rk.min()}", flush=True); sys.exit(0)
    Uf = U.astype(float); Uf /= np.trace(Uf, axis1=1, axis2=2)[:, None, None]
    P, mg = active_set_opt(Uf, rho)
    amax = 2 / (1 - mg) - 2 if mg > 0 else float('nan')
    Q = np.linalg.inv(P)
    print(f"   optimal margin {mg:+.5f} (alpha_max {amax:.4f}); Q ~ {np.round(Q / np.abs(Q).max(), 4).tolist()}  [{time.time() - t0:.1f}s]", flush=True)
    if mg > 0:
        for Qi in integer_forms(Q, rho):
            ok, unsure = exact_check_vec(Qi, U)
            if ok:
                print(f"   EXACT: integer form Q = {Qi} balances all {len(U)} distinct block matrices "
                      f"({unsure} minors needed exact integers)  [{time.time() - t0:.1f}s]", flush=True)
                break
        else:
            print("   exact rationalization failed at the tried scales", flush=True)
