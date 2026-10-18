#!/usr/bin/env python3
"""S42 audit a4: clock_duality's certified closed frontiers V_100, V_1000 (frontier_renewal.json) and the closed-frontier scan.

1. Structure: 1 in V, every element a seed, P(V) subset V (P(1) = 1), recorded rank strictly descends along P, and V equals
   the forward P-closure of {seeds <= Y}.
2. Brute force from the definition (fold n in (N, 4N] to seeds; transport until the FIRST ENTRY into V) for every N < 2^12
   and at the reported arg-max N (4120, 4810, 48192687), compared with the session's closed_frontier_map + run octave
   statistics (numba disk cache disabled) and with the .out values.
3. Independent incremental scan (own code) of every N < 2^26 for both sets: octave max/argmax/mean/min of E_V(N)/N.
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
sys.dont_write_bytecode = True
import numpy as np
from numba import njit, prange


def P(s):
    v = 3 * s + 1
    while v % 2 == 0:
        v //= 2
    while v % 8 == 5:
        v = (v - 1) // 4
    return v


@njit(cache=False)
def fpar(n):
    c = 0
    while n != 1:
        if n & 1:
            n = (3 * n + 1) >> 1
        else:
            n >>= 1
        c += 1
    return 1 - 2 * (c & 1)


@njit(cache=False)
def brute_V(N, inV, maxV):
    c = np.zeros(maxV + 1, dtype=np.int64)
    rdir = 0
    for n in range(N + 1, 4 * N + 1):
        wn = 2 if n % 3 == 2 else -1
        m = n
        k = 0
        while m % 2 == 0:
            m //= 2
            k += 1
        sg = 1 - 2 * (k & 1)
        x = m
        while x % 8 == 5:
            x = (x - 1) // 4
        while not (x <= maxV and inV[x]):
            v = 3 * x + 1
            a = 0
            while v % 2 == 0:
                v //= 2
                a += 1
            if a & 1:
                sg = -sg
            while v % 8 == 5:
                v = (v - 1) // 4
            x = v
        c[x] += wn * sg
        rdir += wn * fpar(n)
    E = 0
    R = 0
    for u in range(1, maxV + 1):
        if c[u] != 0:
            E += c[u] * c[u]
            R += c[u] * fpar(u)
    return E, R, rdir


@njit(cache=False)
def parity_table(M):
    par = np.zeros(M + 1, dtype=np.int8)
    for n in range(2, M + 1):
        if n % 2 == 0:
            par[n] = par[n // 2] ^ 1
        else:
            x = n
            st = 0
            while x >= n:
                x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
                st += 1
            par[n] = (st + par[x]) & 1
    return par


@njit(cache=False)
def v_dest(x, inV, maxV):
    while x % 8 == 5:
        x = (x - 1) // 4
    s = 1
    while True:
        if x <= maxV:
            if inV[x]:
                break
        v = 3 * x + 1
        a = 0
        while v % 2 == 0:
            v //= 2
            a += 1
        if a & 1:
            s = -s
        while v % 8 == 5:
            v = (v - 1) // 4
        x = v
    return x, s


@njit(parallel=True, cache=False)
def vmap(M, inV, maxV):
    n_odd = (M + 1) // 2
    dest = np.zeros(n_odd, dtype=np.int32)
    sgn = np.zeros(n_odd, dtype=np.int8)
    for i in prange(n_odd):
        x, s = v_dest(2 * i + 1, inV, maxV)
        dest[i] = x
        sgn[i] = s
    return dest, sgn


@njit(cache=False)
def scan(Nmax, maxV, dest, sgn, par):
    c = np.zeros(maxV + 1, dtype=np.int64)
    E = np.int64(0); R = np.int64(0); H4 = np.int64(0); H1 = np.int64(0)
    bad = 0
    K = 40
    emax = np.zeros(K); earg = np.zeros(K, dtype=np.int64); esum = np.zeros(K); emin = np.full(K, 1e300)
    cnt = np.zeros(K, dtype=np.int64)
    for N in range(1, Nmax):
        for t in range(3):
            m = 4 * N - 3 if t == 0 else (4 * N - 1 if t == 1 else 2 * N - 1)
            r = m % 3
            w = 2 if r == 2 else -1
            if t == 2:
                amt = -3 * (0 if r == 0 else (1 if r == 1 else -1)) - w
            else:
                amt = w
            i = (m - 1) // 2
            u = dest[i]
            d = sgn[i] * amt
            E += d * (2 * c[u] + d)
            c[u] += d
            R += (1 - 2 * par[u]) * d
        for kk in range(4 * N - 3, 4 * N + 1):
            H4 += (2 if kk % 3 == 2 else -1) * (1 - 2 * par[kk])
        H1 += (2 if N % 3 == 2 else -1) * (1 - 2 * par[N])
        if R != H4 - H1:
            bad += 1
        k = 0
        x = N
        while x > 1:
            x >>= 1
            k += 1
        e = E / N
        if e > emax[k]:
            emax[k] = e
            earg[k] = N
        if e < emin[k]:
            emin[k] = e
        esum[k] += e
        cnt[k] += 1
    return bad, emax, earg, esum, emin, cnt


def main():
    rec = json.load(open(os.path.join(EXP, 'energy_extraction_20261017i', 'frontier_renewal.json')))
    report = {}
    sets = {}
    for cf in rec['closed_frontiers']:
        Y = cf['Y']
        V = [v[0] for v in cf['vertices']]; rank = {v[0]: v[1] for v in cf['vertices']}
        Vs = set(V)
        seeds_ok = all(v % 2 == 1 and v % 8 != 5 for v in V)
        closed = all((P(v) if v != 1 else 1) in Vs for v in V)
        rank_desc = all(rank[P(v)] < rank[v] for v in V if v != 1)
        # forward closure of {seeds <= Y}
        clo = set(s for s in range(1, Y + 1, 2) if s % 8 != 5)
        stack = list(clo)
        while stack:
            x = stack.pop()
            if x == 1:
                continue
            y = P(x)
            if y not in clo:
                clo.add(y); stack.append(y)
        report[Y] = dict(size=len(V), max=max(V), contains_1=1 in Vs, all_seeds=seeds_ok, forward_closed=closed,
                         rank_strictly_descends=rank_desc, equals_forward_closure_of_seeds_le_Y=(clo == Vs),
                         recorded_seeds=cf['seeds'], recorded_max=cf['max_seed'])
        print('Y = %d: |V| = %d (recorded %d), max %d (recorded %d); 1 in V %s; all seeds %s; P(V) in V %s; rank descends %s; '
              'V = forward closure of seeds <= Y: %s' % (Y, len(V), cf['seeds'], max(V), cf['max_seed'], 1 in Vs, seeds_ok,
                                                          closed, rank_desc, clo == Vs), flush=True)
        if Y >= 100:
            sets[Y] = V
    sys.path.insert(0, EXP); sys.path.insert(0, HERE)
    import root_transport_closed_frontier_energy_20261017 as SC
    import root_transport_frontier_energy_big_20261017 as S2
    from _nocache import nocache
    nocache(SC); nocache(S2)
    t0 = time.time()
    for Y, V in sets.items():
        maxV = max(V)
        inV = np.zeros(maxV + 1, dtype=np.bool_); inV[V] = True
        # brute force all N < 2^12, octave stats
        st = {}
        mism_R = 0
        for N in range(1, 2 ** 12):
            E, R, rdir = brute_V(N, inV, maxV)
            if R != rdir:
                mism_R += 1
            k = N.bit_length() - 1
            s = st.setdefault(k, [0.0, 0, 0.0, 1e300, 0])
            e = E / N
            if e > s[0]:
                s[0] = e; s[1] = N
            s[2] += e; s[3] = min(s[3], e); s[4] += 1
        par = S2.clock_parity8(max(4 * 2 ** 12, maxV + 1))  # run() reads par[u] for u in V (up to maxV)
        fu, fs = SC.closed_frontier_map(4 * 2 ** 12, inV, maxV)
        ok, emax, eargmax, esum, emin, cnt = SC.run(2 ** 12, maxV, fu, fs, par)
        agree = all(abs(st[k][0] - emax[k]) < 1e-12 and st[k][1] == eargmax[k] and abs(st[k][2] / st[k][4] - esum[k] / cnt[k]) < 1e-9
                    and abs(st[k][3] - emin[k]) < 1e-12 for k in range(8, 12))
        print('V_%d: brute force (first entry into V) for every N < 4096: R = direct sum at every N: %s; octave stats 2^8..2^11 '
              'agree with the session run(): %s [%.0f s]' % (Y, mism_R == 0, agree, time.time() - t0), flush=True)
        spot = {}
        for N in (4120, 4810, 48192687, 47930133):
            E, R, rdir = brute_V(N, inV, maxV)
            spot[N] = dict(E_over_N=E / N, R=int(R), R_direct=int(rdir))
            print('   brute force N = %d: E_V/N = %.6f, R = %d (direct %d)' % (N, E / N, R, rdir), flush=True)
        report[Y].update(brute_all_N_lt_4096_R_ok=(mism_R == 0), octave_agree_2_8_to_2_11=bool(agree),
                         spot={str(k): v for k, v in spot.items()})
    # independent incremental scan of every N < 2^26
    KMAX = 25
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    par = parity_table(M)
    for Y, V in sets.items():
        maxV = max(V)
        inV = np.zeros(maxV + 1, dtype=np.bool_); inV[V] = True
        dest, sgn = vmap(M, inV, maxV)
        bad, emax, earg, esum, emin, cnt = scan(Nmax, maxV, dest, sgn, par)
        rows = [dict(k=k, max=float(emax[k]), argmax=int(earg[k]), mean=float(esum[k] / cnt[k]), min=float(emin[k]))
                for k in range(8, KMAX + 1)]
        report[Y]['scan_invariant_failures'] = int(bad)
        report[Y]['scan'] = rows
        print('V_%d independent scan N < 2^26: invariant failures %d [%.0f s]' % (Y, bad, time.time() - t0), flush=True)
        for r in rows:
            print('   2^%2d: max %.6f (N=%d) mean %.6f min %.6f' % (r['k'], r['max'], r['argmax'], r['mean'], r['min']), flush=True)
    json.dump({str(k): v for k, v in report.items()}, open(os.path.join(HERE, 'a4_closed_frontier_check.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
