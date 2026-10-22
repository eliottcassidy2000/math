"""S51 audit c5: HYP-9359 and the eleven pre-registered predictions, recomputed from the auditor's own inventory
(c4_inventory.py, cached in the scratchpad) with the auditor's own features, fits, exact one-vertex covers and |Aut|.
  (R) the arithmetic-support model on [301, 1098] and [1101, 2099]; R1-R7 as worded in PREREGISTRATION_out_of_sample.md,
      with robustness variants where the wording leaves a choice (R5: which smooth fit; R6: which order-14 model);
  (D) exact one-vertex covers A u_v B (A, B strong, |A|, |B| >= 3, H = H(A)H(B)) and 3 | |Aut|; the fits without them;
      G1-G4; the per-value table printed by the session;
  (X) the HYP-9359 residuals quoted for 63 and 127 (-1.68, -0.55), searched over natural fit ranges.
Prints ALL CHECKS PASSED or the failures."""
import os, sys, math, time, itertools
from collections import Counter, defaultdict
import numpy as np
import numba as nb
from sympy import isprime, divisors
import alib51 as L

FAIL = []
T0 = time.time()
SCR = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\audit_S51'


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def close(a, b, tol=0.0005):
    return abs(a - b) <= tol


# ---------------------------------------------------------------- inventory
Z = np.load(os.path.join(SCR, 'inventory_2100.npz'))
orders = sorted(int(k[1:]) for k in Z.files if k.startswith('H'))
inv = {n: dict(H=Z[f'H{n}'], aut=Z[f'aut{n}'], A=Z[f'A{n}']) for n in orders}
S = Counter(); Sm = defaultdict(Counter)
for n in orders:
    for h in inv[n]['H']:
        S[int(h)] += 1; Sm[n][int(h)] += 1
print(f'     inventory: {sum(len(inv[n]["H"]) for n in orders)} strong classes, orders {orders} [{time.time() - T0:.0f}s]')


# ---------------------------------------------------------------- features (auditor's own)
def w_len(h, memo={}):
    """longest factorization of h into odd factors >= 3 other than 7, 21 (-1 if none; 0 for h = 1)"""
    if h == 1:
        return 0
    if h in memo:
        return memo[h]
    best = -1
    for d in divisors(h):
        if d >= 3 and d % 2 == 1 and d not in (7, 21):
            r = w_len(h // d)
            if r >= 0:
                best = max(best, r + 1)
    memo[h] = best
    return best


def blocked(d):
    return (not isprime(d)) and d >= 9 and w_len(d) == 1


def F6(d):
    Ld = math.log(d)
    return [1.0, Ld, Ld ** 2, Ld ** 3, d % 3 == 0, d % 9 == 0, d % 5 == 0, bool(isprime(d)), blocked(d), d % 8 == 7]


def F5(d):
    Ld = math.log(d)
    return [1.0, Ld, Ld ** 2, Ld ** 3, d % 3 == 0, d % 9 == 0, d % 5 == 0, d % 8 == 1, d % 8 == 7]


def ols(cnt, ds, feat):
    ds = [d for d in ds if cnt[d] > 0]
    X = np.array([feat(d) for d in ds], dtype=float)
    y = np.log(np.array([cnt[d] for d in ds], dtype=float))
    # solve the normal equations by QR (independent of lstsq's SVD path)
    Q, R = np.linalg.qr(X)
    beta = np.linalg.solve(R, Q.T @ y)
    res = y - X @ beta
    dof = len(ds) - X.shape[1]
    s2 = float(res @ res) / dof
    Rinv = np.linalg.inv(R)
    se = np.sqrt(np.diag(s2 * (Rinv @ Rinv.T)))
    return beta, se, math.sqrt(s2), len(ds), dict(zip(ds, res))


TR = list(range(301, 1099, 2)); TE = list(range(1101, 2100, 2))
N6 = ['3|d', '9|d', '5|d', 'prime', 'blocked', '7 mod 8']
print('(R) the arithmetic-support model')
btr, setr, sdtr, ntr, _ = ols(S, TR, F6)
bte, sete, sdte, nte, _ = ols(S, TE, F6)
fm = lambda b, s: ', '.join(f'{n} {x:+.3f}({e:.3f})' for n, x, e in zip(N6, b[4:], s[4:]))
print(f'     training [301,1098]: {ntr} values, s.d. {sdtr:.3f}: {fm(btr, setr)}')
print(f'     test [1101,2099]:   {nte} values, s.d. {sdte:.3f}: {fm(bte, sete)}')
sess_tr = [0.239, 0.088, 0.065, -0.002, 0.028, -0.027]; sess_tr_se = [0.014, 0.018, 0.014, 0.014, 0.023, 0.012]
sess_te = [0.219, 0.060, 0.054, 0.005, 0.008, -0.015]; sess_te_se = [0.006, 0.008, 0.006, 0.006, 0.010, 0.005]
check(ntr == 399 and nte == 500 and close(sdtr, 0.100) and close(sdte, 0.047), 'value counts 399 / 500 and residual s.d. 0.100 / 0.047')
check(all(close(a, b) for a, b in zip(btr[4:], sess_tr)) and all(close(a, b) for a, b in zip(setr[4:], sess_tr_se)), 'training coefficients and standard errors = HYP-9359 table (to 0.0005)')
check(all(close(a, b) for a, b in zip(bte[4:], sess_te)) and all(close(a, b) for a, b in zip(sete[4:], sess_te_se)), 'test coefficients and standard errors = HYP-9359 table (to 0.0005)')
b = dict(zip(N6, bte[4:]))
R = {}
R['R1'] = 0.15 <= b['3|d'] <= 0.40 and b['3|d'] == max(b.values())
R['R2'] = abs(b['prime']) < 0.05
R['R3'] = -0.05 <= b['7 mod 8'] <= 0.02
R['R4'] = S[2047] < S[2049]
cub = lambda d: [1.0, math.log(d), math.log(d) ** 2, math.log(d) ** 3]
r5 = {}
for name, rng in (('cubic on [1101,2099] (session)', TE), ('cubic on [1101,2099] without 1365', [d for d in TE if d != 1365]),
                  ('cubic on [301,2099]', list(range(301, 2100, 2))), ('cubic on training [301,1098], extrapolated', TR)):
    bb, *_ = ols(S, rng, cub)
    r5[name] = math.log(S[1365]) - float(np.array(cub(1365)) @ bb)
R['R5'] = r5['cubic on [1101,2099] (session)'] > 0
print('     R5 residual of S(1365) = %d: ' % S[1365] + '; '.join(f'{k}: {v:+.3f}' for k, v in r5.items()))
t14 = [d for d in range(1125, 2100, 2) if Sm[14][d] > 0]
r6 = {}
for name, feat in (('session model [1, L, L^2, 3|d, 5|d]', lambda d: [1.0, math.log(d), math.log(d) ** 2, d % 3 == 0, d % 5 == 0]),
                   ('[1, L, L^2, 3|d]', lambda d: [1.0, math.log(d), math.log(d) ** 2, d % 3 == 0]),
                   ('full F6 (cubic + six indicators)', F6),
                   ('F5 (cubic + 3,9,5 | d, 1 and 7 mod 8)', F5)):
    bb, se, *_ = ols(Sm[14], t14, feat)
    r6[name] = (float(bb[4] if feat in (F6, F5) else bb[3]), float(se[4] if feat in (F6, F5) else se[3]))
R['R6'] = r6['session model [1, L, L^2, 3|d, 5|d]'][0] >= 0.4
print(f'     R6: order-14 classes at {len(t14)} values d in [1125, 2099]; beta(3|d): ' + '; '.join(f'{k}: {v[0]:+.3f} ({v[1]:.3f})' for k, v in r6.items()))
# R7 by sequences of strong components (auditor's convolution)
s = defaultdict(Counter); s[1][1] = 1
for n in orders:
    for h in inv[n]['H']:
        s[n][int(h)] += 1
NM = 14
Fn = [Counter() for _ in range(NM + 1)]
Fn[0][1] = 1
for n in range(1, NM + 1):
    for m in range(1, n + 1):
        for d, c in s[m].items():
            for h, f in Fn[n - m].items():
                if h * d < 2100:
                    Fn[n][h * d] += c * f
ok7 = all((Fn[n][2 ** k - 1] < Fn[n][2 ** k + 1]) if k in (5, 6, 7, 9) else (Fn[n][2 ** k - 1] > Fn[n][2 ** k + 1]) for n in (10, 11, 12) for k in (4, 5, 6, 7, 8, 9, 10))
R['R7'] = ok7
check(Fn[14][2047] == 7426 and Fn[14][2049] == 10168 and Fn[10][1023] == 602 and Fn[10][1025] == 543 and Fn[12][255] == 538,
      'N(n, h) by the auditor\'s convolution reproduces the session table (spot values n = 10, 12, 14)')
print('     pre-registered R1-R7 recomputed:', {k: 'PASS' if v else 'FAIL' for k, v in R.items()})
check(all(R.values()), 'R1-R7 all PASS on the auditor\'s data with the session\'s operationalisation')
check(all(v > 0 for v in r5.values()), 'R5 is robust: S(1365) lies above every natural cubic smooth fit tried')
check(all(v[0] >= 0.4 for v in r6.values()), 'R6 is robust: the order-14 beta(3|d) >= 0.4 under every model tried')

print('(D) mechanism: exact one-vertex covers and 3 | |Aut|')


@nb.njit(cache=True)
def popc(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(cache=True)
def cover(n, adj):
    Hs = L.ham_subsets(n, adj)
    full = (1 << n) - 1
    H = Hs[full]
    flags = 0
    for v in range(n):
        rest = full ^ (1 << v)
        low = rest & -rest
        S_ = rest
        while True:
            S_ = (S_ - 1) & rest
            if S_ == 0:
                break
            if not (S_ & low):
                continue
            if popc(S_) < 2 or popc(rest ^ S_) < 2:
                continue
            a = S_ | (1 << v); b2 = (rest ^ S_) | (1 << v)
            if Hs[a] * Hs[b2] == H:
                if L.is_strong_mask(n, adj, a) and L.is_strong_mask(n, adj, b2):
                    flags |= 64
                    for k, val in ((0, 3), (1, 5), (2, 9)):
                        if Hs[a] == val or Hs[b2] == val:
                            flags |= 1 << k
    return flags


@nb.njit(parallel=True, cache=True)
def cover_all(n, A):
    out = np.zeros(A.shape[0], dtype=np.int64)
    for t in nb.prange(A.shape[0]):
        out[t] = cover(n, A[t])
    return out


Sx = Counter(); Sg = Counter(); A3 = Counter(); ncov = 0
for n in orders:
    fl = cover_all(n, inv[n]['A'])
    inv[n]['flags'] = fl
    for h, f, a in zip(inv[n]['H'], fl, inv[n]['aut']):
        h = int(h)
        if not (f & 64):
            Sx[h] += 1
            if a % 3:
                Sg[h] += 1
        else:
            ncov += 1
        if a % 3 == 0:
            A3[h] += 1
    print(f'     order {n}: covers done [{time.time() - T0:.0f}s]', flush=True)
N5 = ['3|d', '9|d', '5|d', '1 mod 8', '7 mod 8']
res = {}
for lo, rng in ((301, TR), (1101, TE)):
    for name, cnt in (('all', S), ('no exact cover', Sx), ('generic', Sg)):
        bb, se, *_ = ols(cnt, rng, F5)
        res[(lo, name)] = (dict(zip(N5, bb[4:])), dict(zip(N5, se[4:])))
        print(f'     [{lo}..] {name:15s}: ' + ', '.join(f'{k} {res[(lo, name)][0][k]:+.3f}({res[(lo, name)][1][k]:.3f})' for k in N5))
share = {lo: sum(A3[d] for d in rng) / sum(S[d] for d in rng) for lo, rng in ((301, TR), (1101, TE))}
print(f'     share of classes with 3 | |Aut|: training {share[301]:.4f}, test {share[1101]:.4f}; classes with an exact cover: {ncov}')
sess = {(301, 'all'): [0.237, 0.088, 0.064, 0.038, -0.015], (301, 'no exact cover'): [0.139, 0.005, -0.003, 0.045, -0.011],
        (301, 'generic'): [-0.015, 0.024, 0.002, 0.046, -0.011], (1101, 'all'): [0.217, 0.060, 0.052, 0.013, -0.011],
        (1101, 'no exact cover'): [0.135, -0.009, -0.005, 0.014, -0.010], (1101, 'generic'): [0.010, -0.007, -0.005, 0.017, -0.012]}
check(all(all(close(res[key][0][k], v) for k, v in zip(N5, vals)) for key, vals in sess.items()),
      'every mechanism-fit coefficient printed in repunit_rarity.out (D) is reproduced (to 0.0005)')
check(close(share[301], 0.0514, 0.00005) and close(share[1101], 0.0447, 0.00005), '3 | |Aut| shares 0.0514 / 0.0447')
G = {}
nx_ = res[(1101, 'no exact cover')][0]; gg = res[(1101, 'generic')][0]
G['G1'] = abs(nx_['5|d']) < 0.03 and abs(nx_['9|d']) < 0.03 and nx_['3|d'] > 0.08
G['G2'] = abs(gg['3|d']) < 0.04
G['G3'] = gg['1 mod 8'] > 0 and abs(gg['7 mod 8']) < 0.03
G['G4'] = share[1101] < share[301]
print('     pre-registered G1-G4 recomputed:', {k: 'PASS' if v else 'FAIL' for k, v in G.items()})
check(all(G.values()), 'G1-G4 all PASS on the auditor\'s data')
gse = res[(301, 'generic')][1]['1 mod 8']
print(f'     prereg baseline quoted "+0.049 +- 0.016" for the training generic 1 mod 8; recomputed {res[(301, "generic")][0]["1 mod 8"]:+.3f} +- {gse:.3f} (F5 model)')
tab = {d: (S[d], Sx[d], Sg[d]) for d in (31, 33, 63, 65, 127, 129, 255, 257, 511, 513, 1023, 1025, 2047, 2049)}
want = {31: (4, 4, 4), 33: (8, 8, 6), 63: (2, 2, 2), 65: (14, 14, 14), 127: (12, 12, 12), 129: (40, 38, 37), 255: (107, 75, 67),
        257: (73, 73, 73), 511: (164, 164, 164), 513: (265, 230, 192), 1023: (663, 593, 500), 1025: (537, 510, 510),
        2047: (1599, 1587, 1587), 2049: (1968, 1828, 1622)}
check(tab == want, 'per-value table S / without exact cover / generic reproduced for the 14 repunit-side values')
# sanity of the mechanism facts themselves
check(all(int(h) % int(a) == 0 for n in orders for h, a in zip(inv[n]['H'], inv[n]['aut'])), '|Aut| divides H for every class (free action)')
import random
rng_ = random.Random(51)


def sub_adj(adj, vs):
    return np.array([sum(1 << vs.index(w) for w in vs if (int(adj[u]) >> w) & 1) for u in vs], dtype=np.int64)


def cover_slow(n, adj):
    H = int(L.ham(n, adj))
    for v in range(n):
        rest = [u for u in range(n) if u != v]
        for r in range(2, n - 2):
            for part in itertools.combinations(rest, r):
                if rest[0] not in part:
                    continue
                a = sorted(list(part) + [v]); b_ = sorted([u for u in rest if u not in part] + [v])
                A_, B_ = sub_adj(adj, a), sub_adj(adj, b_)
                if int(L.ham(len(a), A_)) * int(L.ham(len(b_), B_)) == H and L.is_strong_mask(len(a), A_, (1 << len(a)) - 1) and L.is_strong_mask(len(b_), B_, (1 << len(b_)) - 1):
                    return True
    return False


okm = True; ntest = 0; npos = 0
for n in (5, 6, 7, 8, 9):
    idx = rng_.sample(range(len(inv[n]['H'])), min(60, len(inv[n]['H'])))
    pos = list(np.nonzero(inv[n]['flags'] & 64)[0][:20])
    for t in idx + pos:
        f = cover_slow(n, inv[n]['A'][t])
        okm &= f == bool(inv[n]['flags'][t] & 64)
        ntest += 1; npos += int(f)
check(okm, f'exact-cover flags agree with a slow independent re-implementation on {ntest} classes of orders 5-9 ({npos} with a cover)')

print('(X) the quoted post-model residuals of 63 (-1.68) and 127 (-0.55)')
found = []
for lo in (25, 27, 29, 31, 45, 47, 49, 51, 61, 63, 75, 101, 125):
    for hi in (299, 301, 623, 625, 1097, 1098, 1099, 2099):
        for featn, feat in (('F6', F6), ('F5', F5)):
            rng = list(range(lo if lo % 2 else lo + 1, hi + 1, 2))
            bb, se, sd, nv, rr = ols(S, rng, feat)
            if 63 in rr and 127 in rr and close(rr[63], -1.68, 0.006) and close(rr[127], -0.55, 0.006):
                found.append((lo, hi, featn, round(rr[63], 3), round(rr[127], 3)))
print(f'     fit ranges reproducing (-1.68, -0.55): {found[:6]}')
check(len(found) > 0, 'the quoted residuals of 63 and 127 are reproducible by some natural fit range (the committed scripts do not print them)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
