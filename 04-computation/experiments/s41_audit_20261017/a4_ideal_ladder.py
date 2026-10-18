#!/usr/bin/env python3
"""S41 audit A4 (independent referee code): the ideal energy ladder of THM-4667 statement 6.

(a) K_0(r, b) = h_r(3/2^b) is rebuilt from the prefix/rim definition (not copied), and K_d is computed by a FORWARD
    push in the integer gauge, (s, b') -> (F_o(s), b' + o) with sign (-1)^(o+1), F_o(s) = (3s + c_A) 2^(-A) mod 3^(d+1),
    A = o + 1, c_A = (3 + 2^a - 2^A)/3, on the full band range [0, P_6) with no periodicity assumed. Checks: R8.1
    (period P_d = 2 3^(d+1)) and R8.2 (last-digit balance) where the range allows; exact rational gamma_0..gamma_3 and
    E_0..E_3; gamma_d and E_d for d <= 6 from one period; the diagonal identity diag_d = E_(d-1).
(b) Float push in the positive gauge, residue-major layout, band cap B, for d <= DMAXF: E_d, gamma_d, cross terms.
(c) The band-tail bound of root_transport_energy_ladder_banded_20261017.py, replicated, and a corrected version whose
    initial band marginal includes the bands >= BMAX and whose convolution horizon is extended.
Usage: python a4_ideal_ladder.py [DMAXF] [B]
"""
import os, sys, time, json
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', '/tmp'), 's41_audit_numba_cache'))
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))


def w(n):
    return 2 if n % 3 == 2 else -1


def chi(n):
    return (0, 1, -1)[n % 3]


def h(r, x):
    """Ideal initial coefficient: prefix weight -3 chi(r+k) while 4^k x <= 2, rim weight w(r+k) when 2 < 4^k x <= 4."""
    tot = 0
    k = 0
    while True:
        y = x * Fr(4) ** k
        if y <= 2:
            tot += -3 * chi(r + k)
        elif y <= 4:
            tot += w(r + k)
        else:
            break
        k += 1
    return tot


def K0_table(B):
    return np.array([[h(r, Fr(3, 2 ** b)) for r in range(3)] for b in range(B)], dtype=np.int64)   # [b][r]


def maps(Q, A):
    """F_A(s) = (3 s + c_A) 2^(-A) mod Q as (mult, add)."""
    a = 1 if A % 2 == 1 else 2
    num = 3 + 2 ** a - 2 ** A
    assert num % 3 == 0
    cA = num // 3
    inv = pow(2, -A, Q)
    return (3 * inv) % Q, (cA % Q) * inv % Q


@njit(cache=True)
def push_int(prev, Q, mult, add, out, diag_acc, wts):
    """prev: [q, B] int64 (residue-major), out: [Q, B]. For each A = o+1 <= B: out[F(s), b'+o] += (-1)^(o+1) prev[s, b']."""
    q, B = prev.shape
    for o in range(B):
        sg = 1 if (o + 1) % 2 == 0 else -1
        m = mult[o]
        ad = add[o]
        for s in range(q):
            t = (m * s + ad) % Q
            for bp in range(B - o):
                v = prev[s, bp]
                if v != 0:
                    out[t, bp + o] += sg * v
                    diag_acc[0] += wts[bp + o] * v * v
    return out


@njit(cache=True)
def push_pos(prev, Q, mult, add, out):
    q, B = prev.shape
    p = 0.5
    for o in range(B):
        m = mult[o]
        ad = add[o]
        for s in range(q):
            t = (m * s + ad) % Q
            for bp in range(B - o):
                out[t, bp + o] += p * prev[s, bp]
        p *= 0.5
    return out


def energies_from_period(K, d):
    """K: [3^(d+1), >= P_d] exact. gamma_d and E_d by the periodic formula (float) and the truncated full sum."""
    q = 3 ** (d + 1)
    Pd = 2 * q
    wts = 0.5 ** np.arange(Pd)
    absb = np.abs(K[:, :Pd]).sum(axis=0).astype(np.float64)
    sqb = (K[:, :Pd].astype(np.float64) ** 2).sum(axis=0)
    norm = 4.0 * 2.0 ** d * (1.0 - 2.0 ** (-Pd))
    return float((wts * absb).sum() / norm), float((wts * sqb).sum() / norm)


def exact_rationals(K, d):
    q = 3 ** (d + 1)
    Pd = 2 * q
    g = Fr(0)
    e = Fr(0)
    for b in range(Pd):
        col = K[:, b]
        g += Fr(int(np.abs(col).sum()), 2 ** b)
        e += Fr(int((col.astype(object) ** 2).sum()), 2 ** b)
    norm = Fr(4 * 2 ** d) * (1 - Fr(1, 2 ** Pd))
    return g / norm, e / norm


def main():
    DMAXF = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    t0 = time.time()
    out = dict()
    # ---------------- (a) exact integer push, D = 6 ----------------
    D = 6
    BR = 2 * 3 ** (D + 1)                   # P_6 = 4374
    K0 = K0_table(BR)                       # [b][r]
    assert np.array_equal(K0[:6], np.array([[-1, -1, 2], [0, -3, 3], [-1, -1, 2], [-3, 0, 3], [-1, -1, 2], [0, 0, 0]]))
    assert all(np.array_equal(K0[b], K0[b % 6]) for b in range(BR)), 'K_0 has period 6'
    K = K0.T.copy()                         # [r][b]
    rows = []
    prevE = None
    wts = 0.5 ** np.arange(BR)
    for d in range(0, D + 1):
        if d > 0:
            Q = 3 ** (d + 1)
            mult = np.zeros(BR, dtype=np.int64)
            add = np.zeros(BR, dtype=np.int64)
            for o in range(BR):
                mult[o], add[o] = maps(Q, o + 1)
            new = np.zeros((Q, BR), dtype=np.int64)
            acc = np.zeros(1)
            new = push_int(K, Q, mult, add, new, acc, wts)
            diag = acc[0] / (4.0 * 2.0 ** d)
            K = new
        q = 3 ** (d + 1)
        Pd = 2 * q
        period_ok = None
        if 2 * Pd <= BR:
            period_ok = bool(np.array_equal(K[:, :BR - Pd], K[:, Pd:]))
        third = q // 3
        balance_ok = bool(np.all(K[:third] + K[third:2 * third] + K[2 * third:] == 0)) if d >= 1 else bool(np.all(K.sum(axis=0) == 0))
        g, e = energies_from_period(K, d)
        # full truncated sum (no periodicity used) for the energy
        e_trunc = float((wts * (K.astype(np.float64) ** 2).sum(axis=0)).sum() / (4.0 * 2.0 ** d))
        row = dict(d=d, gamma=g, E=e, E_truncated_sum=e_trunc, period_ok=period_ok, balance_ok=balance_ok,
                   max_abs=int(np.abs(K).max()))
        if d > 0:
            row['diag'] = diag
            row['diag_minus_Eprev'] = diag - prevE
            row['cross'] = e_trunc - prevE
        if d <= 3:
            ge, ee = exact_rationals(K, d)
            row['gamma_exact'] = str(ge)
            row['E_exact'] = str(ee)
        rows.append(row)
        prevE = e_trunc
        print('(a) d = %d: gamma %.10f  E %.9f (truncated sum %.9f)  period %s balance %s max|K| %d%s [%.0f s]'
              % (d, g, e, e_trunc, period_ok, balance_ok, row['max_abs'],
                 ('  diag - E_(d-1) = %.2e, cross %.6f' % (row['diag_minus_Eprev'], row['cross'])) if d else '',
                 time.time() - t0), flush=True)
        if d <= 3:
            print('      exact: gamma = %s, E = %s' % (row['gamma_exact'], row['E_exact']), flush=True)
    out['integer'] = rows
    del K
    # ---------------- (b) positive gauge, float, band cap B ----------------
    nu = np.zeros((3, B))
    for b in range(B):
        for r in range(3):
            nu[r, b] = (-1) ** b * K0[b % 6, r] / (4.0 * 2.0 ** b)
    nu0_band = np.abs(nu).sum(axis=0)
    rowsf = []
    prevE = None
    bw = 2.0 ** np.arange(B)
    for d in range(0, DMAXF + 1):
        if d > 0:
            Q = 3 ** (d + 1)
            mult = np.zeros(B, dtype=np.int64)
            add = np.zeros(B, dtype=np.int64)
            for o in range(B):
                mult[o], add[o] = maps(Q, o + 1)
            new = np.zeros((Q, B))
            nu = push_pos(nu, Q, mult, add, new)
        E = 4.0 * 2.0 ** d * float((nu ** 2 * bw).sum())
        g = float(np.abs(nu).sum())
        cross = None if prevE is None else E - prevE
        rowsf.append(dict(d=d, E=E, gamma=g, cross=cross))
        print('(b) d = %2d: E %.9f gamma %.12f cross %s [%.0f s]' % (d, E, g, '%+.6f' % cross if cross is not None else '-',
                                                                    time.time() - t0), flush=True)
        prevE = E
    out['positive_gauge'] = rowsf
    ratios = [rowsf[i + 1]['cross'] / rowsf[i]['cross'] for i in range(1, len(rowsf) - 1)]
    print('    cross-term ratios:', ' '.join('%.4f' % x for x in ratios))
    out['ratios'] = ratios
    # ---------------- (c) band-tail bounds ----------------
    def tail(mu0, d, BMAX, LONG):
        mu = np.zeros(LONG)
        mu[:len(mu0)] = mu0[:LONG]
        g = 0.5 ** (np.arange(LONG) + 1)
        for _ in range(d):
            mu = np.convolve(mu, g)[:LONG]
        bb = np.arange(LONG)
        return 4 * 2.0 ** d * float(np.sum((2.0 ** bb[BMAX:]) * mu[BMAX:] ** 2)), float(np.sum(mu[BMAX:]))
    BMAX = 200
    mu_session = np.zeros(BMAX)
    mu_full = np.zeros(4000)
    for b in range(4000):
        v = np.abs(K0[b % 6]).sum() / (4.0 * 2.0 ** b)
        mu_full[b] = v
        if b < BMAX:
            mu_session[b] = v
    tails = []
    for d in (1, 6, 10, 12):
        s_e, s_g = tail(mu_session, d, BMAX, 1000)
        f_e, f_g = tail(mu_full, d, BMAX, 4000)
        tails.append(dict(d=d, session_E_tail=s_e, corrected_E_tail=f_e, session_gamma_tail=s_g, corrected_gamma_tail=f_g))
        print('(c) d = %2d: session tail bound E %.3e gamma %.3e | corrected (initial bands >= 200 included, horizon 4000) '
              'E %.3e gamma %.3e' % (d, s_e, s_g, f_e, f_g))
    out['tails'] = tails
    json.dump(out, open(os.path.join(HERE, 'a4_ideal_ladder.json'), 'w'), indent=1)
    print('[%.0f s]' % (time.time() - t0))


if __name__ == '__main__':
    main()
