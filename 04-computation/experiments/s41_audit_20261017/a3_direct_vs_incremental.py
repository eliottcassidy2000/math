#!/usr/bin/env python3
"""S41 audit A3.
 (i)   An own incremental implementation (pure Python dictionaries, own P and seed maps) for every N < 2^14: exact
       energies E_d(N) for all d, octave traces T_d(X) for X = 2^1..2^13, pointwise maxima, D_X, shell sums;
       compared with the committed root_transport_incremental_trace_22_20261017.json.
 (ii)  An own direct transport (vectorised, exact int64; a different algorithm) for all N <= 512 and ~330 sampled N
       (octaves 2^12, 2^13, 2^16, 2^17 and the note's named N), compared depth by depth with (i) and with the session's
       incremental numba kernels (build_tree/apply, driven here to N < 2^18), and with the committed _17 .npz arrays.
 (iii) ROOT shares and participation ratios at N = 7910, 10870, 66463, 116599 (note, section 5).
No session code is used in (i), (ii)-direct, (iii). Numba caches are redirected outside the repository.
Usage: python a3_direct_vs_incremental.py
"""
import os, sys, time, json, random
sys.dont_write_bytecode = True
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', '/tmp'), 's41_audit_numba_cache'))
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
RES = os.path.join(EXP, '..', '..', '05-knowledge', 'results')


# ---------------- own primitives ----------------
def seed(m):
    while m % 8 == 5:
        m = (m - 1) // 4
    return m


NEXT = {}


def P(s):
    r = NEXT.get(s)
    if r is None:
        y = 3 * s + 1
        a = 0
        while y % 2 == 0:
            y //= 2
            a += 1
        assert a in (1, 2)
        r = (seed(y), 1 if a == 2 else -1)
        NEXT[s] = r
    return r


def w(m):
    return 2 if m % 3 == 2 else -1


def thr(m):            # -3 chi(m) - w(m)
    return (1, -2, 1)[m % 3]


def my_incremental(Nmax, DMAX, targets):
    c = [dict() for _ in range(DMAX)]
    E = [0] * DMAX
    Vr = [0] * DMAX
    per = {}
    full = {}
    for N in range(1, Nmax):
        for m, amt in ((4 * N - 3, w(4 * N - 3)), (4 * N - 1, w(4 * N - 1)), (2 * N - 1, thr(2 * N - 1))):
            u = seed(m)
            d = 0
            sign = 1
            while u != 1:
                delta = sign * amt
                cd = c[d]
                old = cd.get(u, 0)
                E[d] += 2 * old * delta + delta * delta
                if old + delta:
                    cd[u] = old + delta
                else:
                    del cd[u]
                u, sg = P(u)
                sign *= sg
                d += 1
                assert d < DMAX
            delta = sign * amt
            for dd in range(d, DMAX):
                old = Vr[dd]
                E[dd] += 2 * old * delta + delta * delta
                Vr[dd] = old + delta
        best = max(range(DMAX), key=lambda i: (E[i], -i))
        per[N] = (E[best] / N, best, Vr[-1], max(v * v for v in Vr) / N)
        if N in targets:
            full[N] = (list(E), list(Vr))
        yield N, E, per[N]
    my_incremental.per = per
    my_incremental.full = full


# ---------------- own direct transport ----------------
def seed_vec(u):
    u = u.copy()
    m = (u % 8) == 5
    while m.any():
        u[m] = (u[m] - 1) // 4
        m = (u % 8) == 5
    return u


def aggregate(keys, vals):
    order = np.argsort(keys, kind='stable')
    k = keys[order]
    v = vals[order]
    if len(k) == 0:
        return k, v
    starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    s = np.add.reduceat(v, starts)
    uk = k[starts]
    keep = s != 0
    return uk[keep], s[keep]


def initial_vec(N):
    m = np.arange(1, 4 * N + 1, 2, dtype=np.int64)
    r = m % 3
    chi = np.where(r == 0, 0, np.where(r == 1, 1, -1))
    wv = np.where(r == 2, 2, -1)
    k = np.where(m < 2 * N, -3 * chi, wv).astype(np.int64)
    return aggregate(seed_vec(m), k)


def step_vec(S, C):
    y = 3 * S + 1
    low = y & (-y)
    assert np.all((low == 2) | (low == 4))
    a = np.where(low == 4, 2, 1)
    t = seed_vec(y >> a)
    return aggregate(t, np.where(a == 2, C, -C))


def profile_direct(N):
    S, C = initial_vec(N)
    E, V = [], []
    while True:
        E.append(int(np.dot(C, C)))
        V.append(int(C[S == 1].sum()) if len(S) else 0)
        if len(S) == 0 or (len(S) == 1 and S[0] == 1):
            break
        S, C = step_vec(S, C)
    return E, V


def main():
    t0 = time.time()
    committed = json.load(open(os.path.join(RES, 'root_transport_incremental_trace_22_20261017.json')))
    crow = {r['k']: r for r in committed['octaves']}
    random.seed(20261017)
    targets = sorted(set(list(range(1, 513)) + random.sample(range(2 ** 12, 2 ** 13), 100)
                         + random.sample(range(2 ** 13, 2 ** 14), 100)
                         + [7910, 10870, 8191, 8192, 16383]))
    big = sorted(set(random.sample(range(2 ** 16, 2 ** 17), 40) + random.sample(range(2 ** 17, 2 ** 18), 10)
                     + [66463, 116599, 169309]))
    report = dict(octaves=[], sampled=[])
    # ---------- (i) own incremental, every N < 2^14 ----------
    KM = 13
    DMAX = 300
    T = {k: [0] * DMAX for k in range(0, KM + 1)}
    DX = {k: 0 for k in range(0, KM + 1)}
    for N, E, (sup, dsup, R, shellmax) in my_incremental(2 ** (KM + 1), DMAX, set(targets)):
        k = N.bit_length() - 1
        Tk = T[k]
        for d in range(DMAX):
            Tk[d] += E[d]
        DX[k] += R * R
    per = my_incremental.per
    print('own incremental done for every N < 2^%d [%.0f s]' % (KM + 1, time.time() - t0), flush=True)
    for k in range(1, KM + 1):
        X = 2 ** k
        prof = [t / X ** 2 for t in T[k]]
        Ns = range(X, 2 * X)
        sup = max(per[N][0] for N in Ns)
        Nsup = max(Ns, key=lambda N: per[N][0])
        frac0 = sum(per[N][1] == 0 for N in Ns) / X
        shell = sum(per[N][3] * N for N in Ns) / X ** 2
        e6 = None
        c = crow[k]
        dprof = max(abs(a - b) for a, b in zip(prof, c['profile']))
        row = dict(k=k, T0=prof[0], Tmax=max(prof), d_Tmax=int(np.argmax(prof)), Tfinal=prof[-1], DX=DX[k] / X ** 2,
                   final_equals_DX=(T[k][-1] == DX[k]), sup_pointwise=sup, N_at_sup=Nsup, frac_sup_at_0=frac0,
                   shell_max_sum=shell, max_abs_diff_profile_vs_committed=dprof,
                   committed=dict((key, c[key]) for key in ('T0', 'Tmax', 'd_Tmax', 'DX', 'sup_pointwise', 'N_at_sup',
                                                            'frac_sup_at_0', 'shell_max_sum')))
        if k == 1:
            row['T5'] = prof[5]
        report['octaves'].append(row)
        print('X = 2^%2d: T0 %.4f, max %.4f at d = %d, final %.4f = D_X %.4f (%s); sup_N %.4f at N = %d (committed %.4f at %d); '
              'peak at d=0: %.4f (committed %.4f); shell %.4f (committed %.4f); max |profile diff| %.1e'
              % (k, prof[0], max(prof), row['d_Tmax'], prof[-1], DX[k] / X ** 2, row['final_equals_DX'], sup, Nsup,
                 c['sup_pointwise'], c['N_at_sup'], frac0, c['frac_sup_at_0'], shell, c['shell_max_sum'], dprof), flush=True)
    allN = [N for N in per if N >= 1]
    best = max(allN, key=lambda N: per[N][0])
    top = sorted(allN, key=lambda N: -per[N][0])[:6]
    print('pointwise: max over 1 <= N < 2^14 of sup_d E_d(N)/N = %.6f at N = %d; top: %s'
          % (per[best][0], best, [(N, round(per[N][0], 4)) for N in top]), flush=True)
    r2n = max((per[N][2] ** 2 / N, N) for N in allN if N >= 2 ** 10)
    print('largest R^2/N for 2^10 <= N < 2^14: %.4f at N = %d' % r2n)
    report['pointwise_max_below_2^14'] = dict(N=best, value=per[best][0], top=[(N, per[N][0]) for N in top])

    # ---------- (ii) direct transport vs own incremental and session kernels ----------
    sys.path.insert(0, EXP)
    from root_transport_incremental_trace_20261017 import build_tree, apply, wfun, thr as thr_s
    KMAX = 17
    Nmax = 2 ** (KMAX + 1)
    level0, succ, sgn, sizes = build_tree(4 * Nmax)
    DM2 = len(sizes) + 1
    V = np.zeros(len(succ), dtype=np.int64)
    Vroot = np.zeros(DM2, dtype=np.int64)
    Einc = np.zeros(DM2, dtype=np.int64)
    tset = set(targets) | set(big)
    rec = {}
    for N in range(1, Nmax):
        for m, amt in ((4 * N - 3, wfun(4 * N - 3)), (4 * N - 1, wfun(4 * N - 1)), (2 * N - 1, thr_s(2 * N - 1))):
            s = m
            while s % 8 == 5:
                s = (s - 1) // 4
            apply(s, amt, level0, succ, sgn, V, Vroot, Einc, DM2)
        if N in tset:
            rec[N] = (Einc.copy(), Vroot.copy())
    print('session incremental kernels driven to N < 2^%d (DMAX %d) [%.0f s]' % (KMAX + 1, DM2, time.time() - t0), flush=True)
    npz = np.load(os.path.join(RES, 'root_transport_incremental_trace_17_20261017.npz'))
    nbad = 0
    for N in sorted(tset):
        E, Vd = profile_direct(N)
        Ei, Vi = rec[N]
        L = len(Ei)
        Ed = np.array(E + [E[-1]] * (L - len(E)), dtype=np.int64)
        Vdd = np.array(Vd + [Vd[-1]] * (L - len(Vd)), dtype=np.int64)
        same_session = bool(np.array_equal(Ed, Ei) and np.array_equal(Vdd, Vi))
        same_own = True
        if N in my_incremental.full:
            Eo, Vo = my_incremental.full[N]
            same_own = (Eo[:len(E)] == E and all(e == E[-1] for e in Eo[len(E):]) and Vo[:len(Vd)] == Vd)
        sup = max(E) / N
        dsup = int(np.argmax(E))
        npz_ok = (abs(float(npz['supE'][N]) - sup) <= 1e-6 * max(1.0, sup) and int(npz['dsup'][N]) == dsup
                  and int(npz['R'][N]) == Vd[-1])
        ok = same_session and same_own and npz_ok
        nbad += not ok
        if not ok:
            print('  MISMATCH N = %d: session %s own %s npz %s' % (N, same_session, same_own, npz_ok))
        report['sampled'].append(dict(N=N, supE=sup, d_sup=dsup, R=Vd[-1], depth=len(E) - 1, ok=ok))
    print('direct transport compared at %d values of N, all depths: mismatches %d [%.0f s]'
          % (len(tset), nbad, time.time() - t0), flush=True)

    # ---------- (iii) ROOT shares ----------
    shares = {}
    for N in (7910, 10870, 66463, 116599):
        S, C = initial_vec(N)
        d = 0
        out = {}
        while True:
            if d in (20, 40, 60):
                Cf = C.astype(np.float64)
                Ef = float(Cf @ Cf)
                root = float(C[S == 1].sum()) ** 2
                pr = Ef ** 2 / float(np.sum(Cf ** 4))
                rest = Ef - root
                Cr = Cf[S != 1]
                prest = rest ** 2 / float(np.sum(Cr ** 4)) if rest > 0 else 0.0
                out[d] = dict(root_share=root / Ef, participation=pr, E_over_N=Ef / N, nonroot_E_over_N=rest / N,
                              nonroot_participation=prest)
            if len(S) == 0 or (len(S) == 1 and S[0] == 1):
                break
            S, C = step_vec(S, C)
            d += 1
        R = int(C[0]) if len(C) else 0
        shares[N] = dict(R2_over_N=R * R / N, depths=out)
        print('N = %d (R^2/N = %.3f): %s' % (N, R * R / N,
                                              {dd: {kk: round(vv, 3) for kk, vv in v.items()} for dd, v in out.items()}))
    report['root_shares'] = shares
    json.dump(report, open(os.path.join(HERE, 'a3_direct_vs_incremental.json'), 'w'), indent=0, default=str)
    print('[%.0f s]' % (time.time() - t0))


if __name__ == '__main__':
    main()
