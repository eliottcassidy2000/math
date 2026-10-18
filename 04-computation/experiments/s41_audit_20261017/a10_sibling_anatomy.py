#!/usr/bin/env python3
"""S41 audit A10 (independent referee code): the sibling anatomy of the note (section 3) and THM-4667 statement 7.

Own incremental transport for every N < 2^13 that also maintains INTRA_d(N) = -2 sum_(s = 1 mod 8) c_d(s) c_d(2s+1)
(ROOT 1 is a head with tail 3): a change delta of c_d(u) changes INTRA_d by -2 delta c_d(partner(u)). Octave sums for
2^10, 2^11, 2^12: E_d, dE_d, INTRA_d, INTER_d, the trace minimum, INTRA/INTER before and after it, and the INTRA share
of dE at d = 0..5. Then the 12 sampled N near 2^16 of the session output, by an own direct transport.
Usage: python a10_sibling_anatomy.py
"""
import os, sys, json, time
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


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
        r = (seed(y), 1 if a == 2 else -1)
        NEXT[s] = r
    return r


def w(m):
    return 2 if m % 3 == 2 else -1


def thr(m):
    return (1, -2, 1)[m % 3]


def partner(u):
    if u % 8 == 1:
        return 2 * u + 1
    if u % 16 == 3:
        return (u - 1) // 2
    return None


def run(Nmax, DMAX, octaves):
    c = [dict() for _ in range(DMAX)]
    E = [0] * DMAX
    IN = [0] * DMAX
    Vr = [0] * DMAX
    acc = {k: dict(E=np.zeros(DMAX), IN=np.zeros(DMAX)) for k in octaves}

    def get(d, u):
        return Vr[d] if u == 1 else c[d].get(u, 0)

    for N in range(1, Nmax):
        for m, amt in ((4 * N - 3, w(4 * N - 3)), (4 * N - 1, w(4 * N - 1)), (2 * N - 1, thr(2 * N - 1))):
            u = seed(m)
            d = 0
            sign = 1
            while u != 1:
                delta = sign * amt
                p = partner(u)
                if p is not None:
                    IN[d] += -2 * delta * get(d, p)
                old = c[d].get(u, 0)
                E[d] += 2 * old * delta + delta * delta
                if old + delta:
                    c[d][u] = old + delta
                else:
                    del c[d][u]
                u, sg = P(u)
                sign *= sg
                d += 1
            delta = sign * amt
            for dd in range(d, DMAX):
                IN[dd] += -2 * delta * c[dd].get(3, 0)
                old = Vr[dd]
                E[dd] += 2 * old * delta + delta * delta
                Vr[dd] = old + delta
        k = N.bit_length() - 1
        if k in acc:
            acc[k]['E'] += np.array(E, dtype=np.float64)
            acc[k]['IN'] += np.array(IN, dtype=np.float64)
    return acc


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
    starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    s = np.add.reduceat(v, starts)
    keep = s != 0
    return k[starts][keep], s[keep]


def direct_with_intra(N):
    m = np.arange(1, 4 * N + 1, 2, dtype=np.int64)
    r = m % 3
    k = np.where(m < 2 * N, -3 * np.where(r == 0, 0, np.where(r == 1, 1, -1)), np.where(r == 2, 2, -1)).astype(np.int64)
    S, C = aggregate(seed_vec(m), k)
    E, IN = [], []
    while True:
        E.append(float(np.dot(C, C)))
        heads = S % 8 == 1
        t = 2 * S[heads] + 1
        pos = np.searchsorted(S, t)
        pos[pos >= len(S)] = 0
        hit = S[pos] == t
        IN.append(-2.0 * float(np.dot(C[heads][hit], C[pos[hit]])))
        if len(S) == 1 and S[0] == 1:
            break
        y = 3 * S + 1
        a = np.where((y & (-y)) == 4, 2, 1)
        S, C = aggregate(seed_vec(y >> a), np.where(a == 2, C, -C))
    return np.array(E), np.array(IN[:-1])


def main():
    t0 = time.time()
    DMAX = 200
    acc = run(2 ** 13, DMAX, (10, 11, 12))
    out = {}
    for k in (10, 11, 12):
        X2 = float(2 ** k) ** 2
        E = acc[k]['E']
        IN = acc[k]['IN']
        dE = np.diff(E)
        INd = IN[:-1]
        last = int(np.flatnonzero(dE)[-1]) + 1
        dmin = int(np.argmin(E[:last + 1]))
        tot_in = INd.sum() / X2
        tot_inter = (dE - INd).sum() / X2
        bef_in = INd[:dmin].sum() / X2
        bef_inter = (dE[:dmin] - INd[:dmin]).sum() / X2
        shares = [float(INd[d] / dE[d]) for d in range(6)]
        out[k] = dict(total_INTRA=tot_in, total_INTER=tot_inter, d_min=dmin, E_min=float(E[dmin] / X2),
                      E_final=float(E[-1] / X2), INTRA_after=tot_in - bef_in, INTER_after=tot_inter - bef_inter,
                      INTRA_share_d0_5=shares)
        print('2^%d: total INTRA %.4f INTER %.4f; trace minimum %.4f at d = %d (final %.4f); after the minimum INTRA %+.4f '
              'INTER %+.4f; INTRA share of dE at d = 0..5: %s'
              % (k, tot_in, tot_inter, E[dmin] / X2, dmin, E[-1] / X2, tot_in - bef_in, tot_inter - bef_inter,
                 ' '.join('%.3f' % x for x in shares)), flush=True)
    samp = [109060, 128067, 119852, 115863, 111937, 73791, 82194, 119733, 87418, 121180, 77100, 88459]
    rows = []
    for N in samp:
        E, IN = direct_with_intra(N)
        dE = np.diff(E)
        dmin = int(np.argmin(E))
        rows.append(dict(N=N, INTRA_before=IN[:dmin].sum() / N, INTER_before=(dE[:dmin] - IN[:dmin]).sum() / N,
                         INTRA_after=IN[dmin:].sum() / N, INTER_after=(dE[dmin:] - IN[dmin:]).sum() / N))
    print('sampled N near 2^16: INTRA before the minimum in [%.3f, %.3f]; INTER before in [%.3f, %.3f] [%.0f s]'
          % (min(r['INTRA_before'] for r in rows), max(r['INTRA_before'] for r in rows),
             min(r['INTER_before'] for r in rows), max(r['INTER_before'] for r in rows), time.time() - t0))
    out['sampled'] = rows
    json.dump(out, open(os.path.join(HERE, 'a10_sibling_anatomy.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
