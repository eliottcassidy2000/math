"""j_inline_checks.py - the auditor's small numpy analyses (run from audit/), one function per output file.
Q is always built from its defining recurrence; tail values are scanned over every m up to the complete horizon
(h + (R - d) m <= H).  usage: python3 j_inline_checks.py NAME   with NAME in
  D_identities | exceptions_summary | shared_halves | pairwise_excess | halves_product | intensity_weighted
"""
import sys, math, itertools, json, collections
import numpy as np
from math import gcd, log

LO, HI = 1 << 22, 1 << 23
NODD = (HI - LO) // 2

def qarr(mmax):
    Q = np.zeros(mmax + 1, dtype=np.int64); lo = 1
    while lo <= mmax:
        hi = min(2 * lo, mmax + 1); idx = np.arange(lo, hi); k = idx >> 1; Q[idx] = Q[k] + (idx & 1) * (k + (k & 1)); lo = hi
    return Q

def load(fn):
    return [tuple(map(int, l.split())) for l in open(fn) if l.strip()]

def image_bits(Q, s, par=None):
    h, R, e, d = s
    m = np.arange(0 if par is None else par, (HI - 1 - h) // (R - d) + 1, 1 if par is None else 2)
    H = h + R * m - d * Q[m] - e * (m & 1); H = H[(H >= LO) & (H < HI)]
    b = np.zeros(NODD, dtype=bool); b[(H - LO) >> 1] = True
    return b

def D_identities():
    N = 1 << 20; Q = qarr(N - 1); m = np.arange(N); Ht = 1 + 2 * m - Q
    bank = load('../bankA.txt'); bad = cnt = 0
    for (h, R, e, d) in bank:
        if R != 2 * d: continue
        cnt += 1
        if not np.array_equal(h + R * m - d * Q - e * (m & 1), (h - d) + d * Ht - e * (m & 1)): bad += 1
    print('(D)3 dilated-root identity on', cnt, 'R=2d bank-A seeds, m<2^20: failures', bad)
    po = [(a, b) for a, b in itertools.combinations(bank, 2) if (a[0] - a[2], a[1], a[3]) == (b[0] - b[2], b[1], b[3])]
    pe = [(a, b) for a, b in itertools.combinations(bank, 2) if (a[0], a[1], a[3]) == (b[0], b[1], b[3])]
    ok = True
    for a, b in po:
        ok &= bool(np.array_equal((a[0] + a[1] * m - a[3] * Q - a[2] * (m & 1))[1::2], (b[0] + b[1] * m - b[3] * Q - b[2] * (m & 1))[1::2]))
    for a, b in pe:
        ok &= bool(np.array_equal((a[0] + a[1] * m - a[3] * Q - a[2] * (m & 1))[0::2], (b[0] + b[1] * m - b[3] * Q - b[2] * (m & 1))[0::2]))
    print('(D)2: bank-A pairs sharing (D-E,U,d):', len(po), ' sharing (D,U,d):', len(pe), ' identical halves:', ok)
    print('   odd-sharing pairs:', po)
    Z = lambda D, E, U, d: (D - E + U, d, 2 * U - d, d); O = lambda D, E, U, d: (D, 0, 2 * U, d)
    print('   same Z-child:', all(Z(a[0], a[2], a[1], a[3]) == Z(b[0], b[2], b[1], b[3]) for a, b in po),
          ' same O-child (even pairs):', all(O(a[0], a[2], a[1], a[3]) == O(b[0], b[2], b[1], b[3]) for a, b in pe))
    h, R, e, d = 81, 99, 27, 45; c = log(R / (R - d)) / d
    for P in (3, 9):
        print('(81,99,27,45) one class mod', P, 'prediction %.5f' % ((1 / P) * (1 - math.exp(-P * 2 * c))))

def exceptions_summary():
    rows = {}
    for l in open('out_exceptions_w22_M63.txt'):
        t = l.split(); s = tuple(map(int, t[1:5])); c, n, p, v = map(int, t[6:10]); rows.setdefault(s, []).append((n, p, v))
    print('window [2^22,2^23), classes mod 63 with p_c>0; kappa=1 Poisson law miss_c=exp(-lambda_c)')
    for s, cl in rows.items():
        live = [(n, p, v) for n, p, v in cl if p > 0]
        dev = [(1 - v / n) - math.exp(-p / n) for n, p, v in live]; dinj = [(1 - v / n) - (1 - p / n) for n, p, v in live]
        lam = sum(p for n, p, v in live) / sum(n for n, p, v in live)
        h, R, e, d = s; kind = 'd=1' if d == 1 else ('R=2d' if R == 2 * d else ('balanced' if h == R and e == d else 'generic'))
        print(f'{str(s):18s} {kind:9s} lambda_live {lam:.3f}  mean signed dev {sum(dev)/len(dev):+.4f}  max|dev| {max(abs(x) for x in dev):.4f}  (no-collision model mean dev {sum(dinj)/len(dinj):+.4f})')

def shared_halves():
    s1, s2 = (25, 30, 10, 15), (27, 30, 12, 15)
    Q = qarr((HI - 25) // 15 + 2)
    A, B = image_bits(Q, s1), image_bits(Q, s2)
    Ao, Bo, Ae, Be = image_bits(Q, s1, 1), image_bits(Q, s2, 1), image_bits(Q, s1, 0), image_bits(Q, s2, 0)
    odd = np.arange(LO + 1, HI, 2)
    print('|A|', int(A.sum()), '|B|', int(B.sum()), '|A&B|', int((A & B).sum()))
    print('odd-m halves identical:', bool(np.array_equal(Ao, Bo)), ' size', int(Ao.sum()), '; even halves intersect:', int((Ae & Be).sum()),
          '; even-half residues mod 15:', set((odd[Ae] % 15).tolist()), set((odd[Be] % 15).tolist()), 'odd-half residues mod 15:', set((odd[Ao] % 15).tolist()))
    for M in (3, 15):
        cls = odd % M; nc = np.bincount(cls, minlength=M).astype(float)
        ca, cb = np.bincount(cls[A], minlength=M), np.bincount(cls[B], minlength=M)
        print('independent baseline resolved mod', M, ':', round(float(np.sum(ca * cb / nc)), 1))

def pairwise_excess():
    seeds = load('bankA_nobal.txt'); Q = qarr(max((HI - 1 - h) // (R - d) for h, R, e, d in seeds) + 1)
    odd = np.arange(LO + 1, HI, 2); M = 315; cls = odd % M; nc = np.bincount(cls, minlength=M).astype(float)
    bm = [image_bits(Q, s) for s in seeds]; cc = [np.bincount(cls[b], minlength=M).astype(float) for b in bm]
    rows = []
    for i, j in itertools.combinations(range(len(seeds)), 2):
        inter = int(np.count_nonzero(bm[i] & bm[j])); exp = float(np.sum(cc[i] * cc[j] / np.where(nc > 0, nc, 1)))
        a, b = seeds[i], seeds[j]; kind = []
        if (a[0] - a[2], a[1], a[3]) == (b[0] - b[2], b[1], b[3]): kind.append('shared odd half (D)2')
        if (a[0], a[1], a[3]) == (b[0], b[1], b[3]): kind.append('shared even half (D)2')
        if a[1] == 2 * a[3] and b[1] == 2 * b[3]: kind.append('both R=2d (D)3')
        rows.append((inter - exp, inter, exp, a, b, ', '.join(kind) or '-'))
    rows.sort(key=lambda r: -r[0])
    print('window [2^22,2^23), 50 non-balanced bank-A seeds, 1225 pairs; independent baseline resolved mod 315')
    for r in rows[:15]: print('excess %9.0f  inter %7d  indep %9.0f  %s %s  %s' % r)
    top = rows[:11]
    print('of the top 11 pairs, (D)2 shared-half pairs:', sum('(D)2' in r[5] for r in top), '; (D)3 R=2d pairs:', sum('(D)3' in r[5] for r in top))
    print('all 11 (D)2 pairs ranks:', [k + 1 for k, r in enumerate(rows) if '(D)2' in r[5]])

def halves_product():
    """product of per-class miss rates over DISTINCT half-images (key (h or h-e, R, d, parity)) versus the union."""
    allb, gen, r2d = load('bankA_nobal.txt'), load('bankA_generic34.txt'), load('bankA_r2d16.txt')
    Q = qarr(HI); odd = np.arange(LO + 1, HI, 2)
    def halves(seeds):
        hv = {}
        for s in seeds:
            h, R, e, d = s
            for par in (0, 1):
                key = ((h if par == 0 else h - e), R, d, par)
                if key not in hv: hv[key] = image_bits(Q, s, par)
        return list(hv.values())
    def union(bs):
        u = np.zeros(NODD, dtype=bool)
        for b in bs: u |= b
        return u
    def prod(members, M):
        cls = odd % M; nc = np.bincount(cls, minlength=M).astype(float); ok = nc > 0; lp = np.zeros(M)
        for b in members:
            v = np.bincount(cls[b], minlength=M).astype(float); lp[ok] += np.log(np.clip(1 - v[ok] / nc[ok], 1e-300, None))
        return float(np.sum(nc[ok] * np.exp(lp[ok])) / nc.sum())
    print('method: union miss of the seeds; product over seeds; product over distinct half-images; window [2^22,2^23)')
    for name, seeds in (('bank A without balanced (50)', allb), ('34 generic', gen), ('16 R=2d (root incl.)', r2d)):
        imgs = [image_bits(Q, s) for s in seeds]; hs = halves(seeds); U = union(imgs)
        for M in (3, 945):
            print('%-30s mod %4d: union miss %.5f | product over seeds %.5f | product over %d distinct half-images %.5f'
                  % (name, M, 1 - U.mean(), prod(imgs, M), len(hs), prod(hs, M)))
    g, r = union(halves(gen)), union(halves(r2d)); M = 945
    cls = odd % M; nc = np.bincount(cls, minlength=M).astype(float); ok = nc > 0
    mg = 1 - np.bincount(cls[g], minlength=M) / np.where(ok, nc, 1); mr = 1 - np.bincount(cls[r], minlength=M) / np.where(ok, nc, 1)
    print('exact union(generic) x exact union(R=2d) as two independent blocks, mod %d: %.5f  vs union of all 50: %.5f'
          % (M, float(np.sum(nc[ok] * mg[ok] * mr[ok]) / nc.sum()), 1 - (g | r).mean()))

def intensity_weighted():
    bank = load('../bankA.txt'); bal = {tuple(s) for s in json.load(open('out_states_U4096.json'))['bankA_lists']['balanced']}
    c = lambda h, R, e, d: log(R / (R - d)) / d
    tot = sum(2 * c(*s) for s in bank)
    print('bank A nominal intensity 2 sum c_s = %.3f' % tot)
    print('with the root counted once per involution pair (subtract ln 2): %.3f' % (tot - log(2)))
    nob = sum(2 * c(*s) for s in bank if s not in bal)
    print('without the 35 balanced seeds (root halved): %.3f ; (root not halved): %.3f' % (nob - log(2), nob))
    print('-ln(0.0219) = %.3f' % (-log(0.02187)))
    c0 = lambda h, R, e, d: log(R / (R - d / 2)) / d; c1 = lambda h, R, e, d: log((R - d / 2) / (R - d)) / d
    for label, skip in (('', set()), ('same without the 35 nested balanced seeds: ', bal)):
        worst = (1e9, None)
        for p in range(3, 106, 2):
            for a in range(p):
                W = log(2) / (2 * p)
                for s in bank[1:]:
                    if s in skip: continue
                    h, R, e, d = s; g = gcd(gcd(R, d), p)
                    W += (g / p) * (c0(*s) * ((a - h) % g == 0) + c1(*s) * ((a - (h - e)) % g == 0))
                if W * 2 * p < worst[0]: worst = (W * 2 * p, (p, a))
        print(label + 'min over odd p<=105 and classes a of (weighted capacity)/(class size) = %.3f at (p,a)=%s' % worst)
    for lab, lam, frac in [('41.411 (S82 smallest class mean, no loss)', 41.411, 1 / 6), ('41.411/2 ("loss of half")', 41.411 / 2, 1 / 6),
                           ('41.411*3.82/8.29 (bank A nominal->observed)', 41.411 * 3.82 / 8.29, 1 / 6),
                           ('41.411*3.82/5.15 (bank A after nesting)', 41.411 * 3.82 / 5.15, 1 / 6),
                           ('48.2*3.82/5.15 (overall mean, all odd N)', 48.2 * 3.82 / 5.15, 1 / 2)]:
        print('lambda_eff = %-45s %.2f -> expected first miss near N ~ %.1e ; expected misses below 2^35: %.3g'
              % (lab, lam, 1 / (frac * math.exp(-lam)), (2 ** 35) * frac * math.exp(-lam)))

if __name__ == '__main__':
    {'D_identities': D_identities, 'exceptions_summary': exceptions_summary, 'shared_halves': shared_halves,
     'pairwise_excess': pairwise_excess, 'halves_product': halves_product, 'intensity_weighted': intensity_weighted}[sys.argv[1]]()
