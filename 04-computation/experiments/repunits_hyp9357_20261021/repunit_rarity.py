"""Why are the repunits 31, 63, 127 rare among Hamiltonian-path counts?  (session opus-2026-10-21-S51)
Note: 05-knowledge/results/repunit_rarity_hyp9357_20261021.md; THM-4748, THM-4749, HYP-9359.

S(d) = number of strong tournament classes (all orders) with exactly d Hamiltonian paths; S_m(d) at order m.
(I)  Inventory: every strong class with H <= 2100, by exact ears (Moon's strong-deletion theorem + the insertion lemma
     H(T - v) <= H(T); the audit's code audit_S50/b3_ear_census_625.py with the cap raised).  Checked against THM-4744.
(R)  The arithmetic-support model log S(d) = cubic(log d) + b3[3|d] + b9[9|d] + b5[5|d] + bp[prime] + bb[blocked]
     + b7[d = 7 mod 8], fitted on the training range [301, 1098]; the pre-registered predictions R1-R7
     (PREREGISTRATION_out_of_sample.md, commit 917fef912a) evaluated on [1101, 2099].  R7: class counts N(n,h) by
     sequences of strong components (checked against the full class census n <= 9).
(D)  Mechanism: for every class, the exact one-vertex covers (A u_v B, A, B strong, H = H(A)H(B)) by the smaller
     piece, and |Aut| (|Aut| divides H: Aut acts freely on Hamiltonian paths).  The fits with these classes removed,
     training and test; the pre-registered G1-G4 (second pre-registration, commit be7c76115f).
(O)  Per order: H mod 8 by position H/f(m) in the strong spectrum (orders 6-9, from the class census), and the
     minimum feedback arc number against multiplicity at n = 8.
The script reports the verdicts; it asserts only internal consistency.  Prints ALL CHECKS PASSED.  Runtime ~6 min."""
import os, sys, time, math, itertools
from collections import Counter, defaultdict
from multiprocessing import Pool
import numpy as np
import numba as nb
from sympy import isprime

HERE = os.path.dirname(os.path.abspath(__file__))
AUD = os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021', 'audit_S50')
sys.path.insert(0, AUD)
sys.path.insert(0, os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021'))
import alib as L                                       # noqa: E402
from b3_ear_census_625 import grow                     # noqa: E402

FAIL = []
T0 = time.time()
LIMIT = 2100


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


# ---------------------------------------------------------------- numba kernels
@nb.njit(cache=True)
def bcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(cache=True)
def sub_tables(n, adj):
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            m = adj[v] & ~S & (N - 1)
            while m:
                w = 0
                while not (m >> w) & 1:
                    w += 1
                m &= m - 1
                dp[S | (1 << w), w] += c
    Hs = np.zeros(N, dtype=np.int64)
    st = np.zeros(N, dtype=np.uint8)
    for S in range(1, N):
        t = 0
        for v in range(n):
            t += dp[S, v]
        Hs[S] = t
        low = S & -S
        r = low
        while True:
            x = r; nr = r
            while x:
                u = 0
                while not (x >> u) & 1:
                    u += 1
                x &= x - 1
                nr |= adj[u] & S
            if nr == r:
                break
            r = nr
        if r != S:
            continue
        r = low
        while True:
            nr = r
            for u in range(n):
                if (S >> u) & 1 and not (r >> u) & 1 and (adj[u] & r):
                    nr |= 1 << u
            if nr == r:
                break
            r = nr
        if r == S:
            st[S] = 1
    return Hs, st


@nb.njit(cache=True)
def cover_flags(n, adj):
    """bit 0: an exact cover with a piece of H = 3; bit 1: H = 5; bit 2: H = 9; bit 6: any exact cover"""
    Hs, st = sub_tables(n, adj)
    full = (1 << n) - 1
    H = Hs[full]
    flags = 0
    for v in range(n):
        rest = full & ~(1 << v)
        low = rest & -rest
        S = rest
        while S:
            S = (S - 1) & rest
            if S == 0:
                break
            if not (S & low) or bcount(S) < 2 or bcount(rest & ~S) < 2:
                continue
            a = S | (1 << v); b = (rest & ~S) | (1 << v)
            if st[a] and st[b] and Hs[a] * Hs[b] == H:
                flags |= 64
                for k, val in ((0, 3), (1, 5), (2, 9)):
                    if Hs[a] == val or Hs[b] == val:
                        flags |= 1 << k
    return flags


@nb.njit(parallel=True, cache=True)
def cover_flags_all(n, A):
    out = np.zeros(A.shape[0], dtype=np.int64)
    for t in nb.prange(A.shape[0]):
        out[t] = cover_flags(n, A[t])
    return out


@nb.njit(cache=True)
def fa_one(n, adj):
    N = 1 << n
    f = np.full(N, 10 ** 9, dtype=np.int64)
    f[0] = 0
    for S in range(1, N):
        best = 10 ** 9
        x = S
        while x:
            v = 0
            while not (x >> v) & 1:
                v += 1
            x &= x - 1
            T = S & ~(1 << v)
            c = f[T] + bcount(adj[v] & T)
            if c < best:
                best = c
        f[S] = best
    return f[N - 1]


def aut_of(args):
    m, code = args
    return L.canon(L.decode(m, code), want_aut=True)[1]


# ---------------------------------------------------------------- model
def blocked(d):
    if isprime(d) or d < 9:
        return False
    for a in range(3, int(d ** 0.5) + 1, 2):
        if d % a == 0 and a not in (7, 21) and d // a not in (7, 21):
            return False
    return True


NAMES = ['3|d', '9|d', '5|d', 'prime', 'blocked', 'd=7 mod 8']


def feats(d):
    Ld = math.log(d)
    return [1.0, Ld, Ld * Ld, Ld ** 3, float(d % 3 == 0), float(d % 9 == 0), float(d % 5 == 0), float(isprime(d)),
            float(blocked(d)), float(d % 8 == 7)]


def fit(cnt, lo, hi, featf=feats):
    ds = [d for d in range(lo, hi + 1, 2) if cnt[d] > 0]
    X = np.array([featf(d) for d in ds], dtype=float); y = np.log(np.array([cnt[d] for d in ds], dtype=float))
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    s2 = float(res @ res) / (len(ds) - X.shape[1])
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return coef[4:], se[4:], math.sqrt(s2), len(ds), dict(zip(ds, res))


def fmt(coef, se, names):
    return ', '.join(f'{n} {c:+.3f}({s:.3f})' for n, c, s in zip(names, coef, se))


if __name__ == '__main__':
    print('(I) strong inventory below', LIMIT)
    inv = {3: {L.canon(L.chain((2,))): 3}}
    with Pool(12) as pool:
        n = 4
        while True:
            par = [c for c, h in inv[n - 1].items() if h <= LIMIT]
            if not par:
                break
            got, mn = grow(pool, n, par, cap=LIMIT)
            got = {c: h for c, h in got.items() if h <= LIMIT}
            if not got:
                break
            inv[n] = got
            print(f'     order {n}: {len(got)} strong classes with H <= {LIMIT} [{time.time() - T0:.0f}s]', flush=True)
            n += 1
    S = Counter(); Sm = defaultdict(Counter)
    for m, dd in inv.items():
        for c, h in dd.items():
            S[h] += 1; Sm[m][h] += 1
    c625 = {m: sum(1 for h in inv[m].values() if h <= 625) for m in (9, 10, 11, 12, 13)}
    check(c625 == {9: 13439, 10: 8769, 11: 3536, 12: 658, 13: 1}, f'classes with H <= 625 at orders 9-13: {c625} (THM-4744: 13,439 / 8,769 / 3,536 / 658, and IC(3,3,3,3) at 625)')
    check(all(S[d] > 0 for d in range(3, LIMIT, 2) if d not in (7, 21)) and S[7] == 0 and S[21] == 0, 'every odd 3 <= d < 2100 except 7, 21 is a strong count')

    print('(R) the arithmetic-support model; pre-registered predictions')
    tr = fit(S, 301, 1098); te = fit(S, 1101, 2099)
    print(f'     training [301,1098] ({tr[3]} values, residual s.d. {tr[2]:.3f}): {fmt(tr[0], tr[1], NAMES)}')
    print(f'     test     [1101,2099] ({te[3]} values, residual s.d. {te[2]:.3f}): {fmt(te[0], te[1], NAMES)}')
    b = dict(zip(NAMES, te[0]))
    verdict = {}
    verdict['R1'] = 0.15 <= b['3|d'] <= 0.40 and b['3|d'] == max(b.values())
    verdict['R2'] = abs(b['prime']) < 0.05
    verdict['R3'] = -0.05 <= b['d=7 mod 8'] <= 0.02
    verdict['R4'] = S[2047] < S[2049]
    ds = [d for d in range(1101, 2100, 2)]
    X = np.array([[1, math.log(d), math.log(d) ** 2, math.log(d) ** 3] for d in ds]); y = np.log([S[d] for d in ds])
    c0, *_ = np.linalg.lstsq(X, y, rcond=None)
    r1365 = math.log(S[1365]) - float(np.array([1, math.log(1365), math.log(1365) ** 2, math.log(1365) ** 3]) @ c0)
    verdict['R5'] = r1365 > 0
    t14 = [d for d in range(1125, 2100, 2) if Sm[14][d] > 0]
    X = np.array([[1, math.log(d), math.log(d) ** 2, d % 3 == 0, d % 5 == 0] for d in t14], dtype=float); y = np.log([Sm[14][d] for d in t14])
    c6, *_ = np.linalg.lstsq(X, y, rcond=None)
    verdict['R6'] = c6[3] >= 0.4
    print(f'     S(2047) = {S[2047]}, S(2049) = {S[2049]}; S(1365) = {S[1365]} (residual {r1365:+.3f}); order 14: beta(3|d) = {c6[3]:+.3f}')
    # R7: sequences of strong components
    s = defaultdict(Counter); s[1][1] = 1
    for m, dd in inv.items():
        for c, h in dd.items():
            s[m][h] += 1
    F = [Counter() for _ in range(15)]
    F[0][1] = 1
    for n in range(1, 15):
        for m in range(1, n + 1):
            for d, cnt in s[m].items():
                for h, f in F[n - m].items():
                    if h * d < LIMIT:
                        F[n][h * d] += cnt * f
    import h_multiplicities as HM                      # noqa: E402
    data = HM.classes_upto(9)
    agree = True
    census = {}
    for n in range(3, 10):
        H, st = HM.ham_all(n, data[n]['A'])
        census[n] = (H, st)
        C = Counter(int(h) for h in H)
        agree &= all(C[h] == F[n][h] for h in range(1, LIMIT, 2))
    check(agree, 'class counts N(n,h) from sequences of strong components = the class census for n <= 9, every h < 2100')
    ok7 = True
    for n in (10, 11, 12):
        for k in range(4, 11):
            a, bb = F[n][2 ** k - 1], F[n][2 ** k + 1]
            ok7 &= (a < bb) if k in (5, 6, 7, 9) else (a > bb)
    verdict['R7'] = ok7
    for n in (8, 9, 10, 11, 12, 13, 14):
        print(f'     n = {n}: ' + '; '.join(f'{2 ** k - 1}:{F[n][2 ** k - 1]} vs {2 ** k + 1}:{F[n][2 ** k + 1]}' for k in range(4, 12)))
    print('     pre-registered verdicts:', {k: ('PASS' if v else 'FAIL') for k, v in verdict.items()})
    check(True, 'R1-R7 evaluated (the script reports verdicts; it does not assert them)')

    print('(D) mechanism: exact one-vertex covers and automorphisms of order 3')
    keys = [(m, c) for m in sorted(inv) for c in inv[m]]
    with Pool(12) as pool:
        auts = dict(zip(keys, pool.map(aut_of, keys, chunksize=500)))
    flags = {}
    for m in sorted(inv):
        codes = list(inv[m])
        A = np.array([L.decode(m, c) for c in codes], dtype=np.int64)
        fl = cover_flags_all(m, A)
        for c, f in zip(codes, fl):
            flags[(m, c)] = int(f)
        print(f'     order {m}: covers done [{time.time() - T0:.0f}s]', flush=True)
    check(all(inv[m][c] % auts[(m, c)] == 0 for (m, c) in keys), '|Aut| divides H for every class (free action on Hamiltonian paths)')
    Sx = Counter(); Sg = Counter(); A3 = Counter()
    for (m, c) in keys:
        h = inv[m][c]
        if not (flags[(m, c)] & 64):
            Sx[h] += 1
            if auts[(m, c)] % 3:
                pass
        if not (flags[(m, c)] & 64) and auts[(m, c)] % 3 != 0:
            Sg[h] += 1
        if auts[(m, c)] % 3 == 0:
            A3[h] += 1
    names3 = ['3|d', '9|d', '5|d', '1 mod 8', '7 mod 8']

    def f3(d):
        Ld = math.log(d)
        return [1.0, Ld, Ld * Ld, Ld ** 3, float(d % 3 == 0), float(d % 9 == 0), float(d % 5 == 0), float(d % 8 == 1), float(d % 8 == 7)]
    res = {}
    for lo, hi in ((301, 1098), (1101, 2099)):
        for name, cnt in (('all', S), ('no exact cover', Sx), ('generic (no exact cover, 3 not | Aut)', Sg)):
            co, se, sd, nv, _ = fit(cnt, lo, hi, f3)
            res[(lo, name)] = dict(zip(names3, co))
            print(f'     [{lo},{hi}] {name:38s}: {fmt(co, se, names3)}')
    share = {lo: sum(A3[d] for d in range(lo, hi + 1, 2)) / sum(S[d] for d in range(lo, hi + 1, 2)) for lo, hi in ((301, 1098), (1101, 2099))}
    print(f'     share of strong classes with 3 | |Aut|: [301,1098] {share[301]:.4f}; [1101,2099] {share[1101]:.4f}')
    g = {}
    g['G1'] = abs(res[(1101, 'no exact cover')]['5|d']) < 0.03 and abs(res[(1101, 'no exact cover')]['9|d']) < 0.03 and res[(1101, 'no exact cover')]['3|d'] > 0.08
    g['G2'] = abs(res[(1101, 'generic (no exact cover, 3 not | Aut)')]['3|d']) < 0.04
    g['G3'] = res[(1101, 'generic (no exact cover, 3 not | Aut)')]['1 mod 8'] > 0 and abs(res[(1101, 'generic (no exact cover, 3 not | Aut)')]['7 mod 8']) < 0.03
    g['G4'] = share[1101] < share[301]
    print('     pre-registered verdicts:', {k: ('PASS' if v else 'FAIL') for k, v in g.items()})
    check(True, 'G1-G4 evaluated (reported, not asserted)')
    for d in (31, 33, 63, 65, 127, 129, 255, 257, 511, 513, 1023, 1025, 2047, 2049):
        print(f'     d = {d}: S = {S[d]}, without exact covers {Sx[d]}, generic {Sg[d]}')

    print('(O) per order: H mod 8 by position in the strong spectrum; minimum feedback arc number at n = 8')
    fmin = {6: 15, 7: 25, 8: 45, 9: 75}
    for m in (6, 7, 8, 9):
        H, st = census[m]
        bins = defaultdict(Counter)
        for h, sflag in zip(H, st):
            if not sflag:
                continue
            r = int(h) / fmin[m]
            bins['<2' if r < 2 else '2-4' if r < 4 else '4-8' if r < 8 else '>=8'][int(h) % 8] += 1
        print(f'     order {m}: ' + '; '.join(f'H/f in {k}: ' + ' '.join(f'{rr}:{bins[k][rr] * 4 / max(1, sum(bins[k].values())):.2f}' for rr in (1, 3, 5, 7)) for k in ('<2', '2-4', '4-8', '>=8') if bins[k]))
    A8 = data[8]['A']; H8, _ = census[8]
    fa8 = np.array([fa_one(8, A8[t]) for t in range(len(A8))])
    N8 = Counter(int(h) for h in H8); mf = defaultdict(lambda: 99)
    for h, f in zip(H8, fa8):
        mf[int(h)] = min(mf[int(h)], int(f))
    hs = [h for h in range(25, 331, 2) if N8[h] > 0]
    X = np.array([[1, math.log(h), math.log(h) ** 2, h % 3 == 0, h % 5 == 0] for h in hs], dtype=float); y = np.log([N8[h] for h in hs])
    cc, *_ = np.linalg.lstsq(X, y, rcond=None); rr = y - X @ cc
    by = defaultdict(list)
    for h, x in zip(hs, rr):
        by[mf[h]].append(x)
    print(f'     n = 8: corr(log-multiplicity residual, least feedback arc number) = {np.corrcoef(rr, [mf[h] for h in hs])[0, 1]:+.3f}; mean residual by least fa: '
          + str({k: (round(float(np.mean(v)), 3), len(v)) for k, v in sorted(by.items())}))
    check(all(mf[h] >= 3 for h in (31, 63, 127)) and all(mf[h] == 1 for h in (33, 65)), 'n = 8: least fa of 31, 63, 127 is >= 3; of 33, 65 it is 1')

    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
