#!/usr/bin/env python3
"""
collatz_procgen_20260922_order_laws.py

Procedural search for SIGN-SPECIFIC ORDER LAWS of odd-only Syracuse orbits
    U_b(x) = oddpart(3x + b),  b = +1 (plus sheet) or b = -1 (minus sheet),
on the positive odd integers.  Orbits are STOPPED at the first element of
the cycle set C_+ = {1}, C_- = {1, 5, 7, 17, 25, 37, 41, 55, 61, 91}.

Sections of the output
  S0  C census engine (collatz_procgen_20260922_order_laws.c) at N = 20001
      (check), 10^6 + 1 (search) and 10^7 + 1 (extension), both sheets
  S1  independent pure-Python recomputation of EVERY C output at N = 20001
  S2  the grammar (seven families, 6,6xx candidates) and its evaluation
  S3  automatic classification of the survivors, with witnesses
  S4  the deviation spectrum (gate crossings): exact big-integer gates,
      clocks, computer-assisted completeness, anatomy (dip and return)
  S5  the coefficient stopping time law sigma = tau on both sheets
  S6  extension of every survivor to 10^7
  S7  the three best candidates

Residues in the main grammar are SIGNED (b*x mod M), so a word function is
the same function on both sheets (the word bijection of the minus-sheet
control lane sends the plus class r to the minus class -r).  Raw residues
and the raw map n -> 4n+1 appear only as a labelled control group.

Reproduce:
  python3 04-computation/experiments/collatz_procgen_20260922_order_laws.py \
      > 05-knowledge/results/collatz_procgen_20260922_order_laws.out
(about 2.5 minutes on one core, the C runs dominate; peak RSS about 1.05 GB; timing
lines go to stderr; set ORDER_LAWS_RERUN=1 to force the C runs even if the
scratch outputs exist.)
"""
import decimal
import itertools
import math
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
CSRC = os.path.join(HERE, 'collatz_procgen_20260922_order_laws.c')
SCR = os.path.join(ROOT, 'scratch', 'procgen_order')
BIN = os.path.join(SCR, 'order_laws')
N_CHECK = 20001
N_SEARCH = 1000001
N_EXT = 10000001
LAG_SEARCH = 64
LAG_EXT = 128
LOG23 = math.log2(3)
CYC = {1: {1}, -1: {1, 5, 7, 17, 25, 37, 41, 55, 61, 91}}
SHEET = {1: 'plus', -1: 'minus'}
T0 = time.time()


def log(msg):
    sys.stderr.write('[%7.1fs] %s\n' % (time.time() - T0, msg))
    sys.stderr.flush()


def check(cond, msg):
    if not cond:
        raise RuntimeError('CHECK FAILED: ' + msg)


# ----------------------------------------------------------------------------
# dynamics, words, patterns
# ----------------------------------------------------------------------------
def U(x, b):
    y = 3 * x + b
    k = (y & -y).bit_length() - 1
    return y >> k, k


def stopped_orbit(n, b, extra=0):
    xs = [n]
    ks = []
    x = n
    while x not in CYC[b]:
        x, k = U(x, b)
        ks.append(k)
        xs.append(x)
    T = len(xs) - 1
    for _ in range(extra):
        x, k = U(x, b)
        ks.append(k)
        xs.append(x)
    return xs, ks, T


def sres(x, M, b):
    return (b * x) % M


def lehmer(vals):
    w = len(vals)
    idx = 0
    for i in range(w):
        c = 0
        vi = vals[i]
        for j in range(i + 1, w):
            if vals[j] < vi:
                c += 1
        idx += c * math.factorial(w - 1 - i)
    return idx


def perm_of(idx, w):
    """rank tuple (0 = smallest) of a Lehmer index."""
    digits = []
    for i in range(w):
        f = math.factorial(w - 1 - i)
        digits.append(idx // f)
        idx %= f
    avail = list(range(w))
    return tuple(avail.pop(d) for d in digits)


def ptxt(idx, w):
    return ''.join(str(r) for r in perm_of(idx, w))


def growth(l, K):
    """True iff 3^l > 2^K (exact)."""
    return 3 ** l > 2 ** K


def generic_pattern(word):
    K = 0
    vals = [0.0]
    for i, k in enumerate(word):
        K += k
        vals.append((i + 1) * LOG23 - K)
    return lehmer(vals)


def generic_sets():
    """w -> {generic pattern index -> words (k_i <= 9) realizing it}.  k_i <= 9 suffices: a step with
    k >= ceil(5 log2 3) + 1 = 9 already makes every window containing it a decay window."""
    out = {}
    for w in range(3, 7):
        d = defaultdict(list)
        for word in itertools.product(range(1, 10), repeat=w - 1):
            d[generic_pattern(word)].append(word)
        for v in d.values():
            v.sort(key=lambda wd: (sum(wd), wd))
        out[w] = d
    return out


def word_class(word, b):
    """(r, 2^(K+1)): the odd class of n whose first exponents are `word` on sheet b (backward 2-adic solve)."""
    k = word[-1]
    mod = 2 ** (k + 1)
    a = ((2 ** k - b) * pow(3, -1, mod)) % mod
    for kk in reversed(word[:-1]):
        mod2 = mod * 2 ** kk
        a = ((2 ** kk * a - b) * pow(3, -1, mod2)) % mod2
        mod = mod2
    return a, mod


def carry(word, b):
    """B with 2^K x_L = 3^L x_0 + B."""
    Bv = 0
    K = 0
    L = len(word)
    for i, k in enumerate(word):
        Bv += 3 ** (L - 1 - i) * 2 ** K
        K += k
    return b * Bv


def gate_bound(l, b, kextra=80):
    """G_b(l) = max over words of length l with a POSITIVE gate on sheet b of |B|/|2^K - 3^l|, using
    B_max(l,K) = 3^(l-1) + 2^(K-l+1)(3^(l-1) - 2^(l-1)) (all later exponents 1, first exponent K-l+1).
    plus: decay words (2^K > 3^l, K unbounded: the tail K > Kc + kextra is below the K -> infinity limit
    (3^(l-1) - 2^(l-1))/2^(l-1), which is included); minus: growth words (l <= K < l log2 3)."""
    best = Fraction(0)
    Kc = math.floor(l * LOG23)
    Ks = range(Kc + 1, Kc + 1 + kextra) if b > 0 else range(l, Kc + 1)
    for K in Ks:
        D = 2 ** K - 3 ** l
        if (b > 0 and D <= 0) or (b < 0 and D >= 0):
            continue
        Bm = 3 ** (l - 1) + 2 ** (K - l + 1) * (3 ** (l - 1) - 2 ** (l - 1))
        g = Fraction(Bm, abs(D))
        if g > best:
            best = g
    if b > 0:
        # tail K >= Ks = Kc + kextra + 1: B_max/D = A 2^K/(2^(l-1)(2^K - 3^l)) + 3^(l-1)/(2^K - 3^l) with
        # A = 3^(l-1) - 2^(l-1); both terms decrease in K, so the value at K = Ks bounds the whole tail
        Ks_ = Kc + kextra + 1
        A = 3 ** (l - 1) - 2 ** (l - 1)
        D = 2 ** Ks_ - 3 ** l
        tail = Fraction(A * 2 ** Ks_, 2 ** (l - 1) * D) + Fraction(3 ** (l - 1), D)
        best = max(best, tail)
    return best


def cf_log23(nterms=16):
    decimal.getcontext().prec = 150
    x = decimal.Decimal(3).ln() / decimal.Decimal(2).ln()
    a = []
    for _ in range(nterms):
        q = int(x)
        a.append(q)
        x = 1 / (x - q)
    return a


def convergents(a):
    out = []
    p0, q0, p1, q1 = 1, 0, a[0], 1
    out.append((p1, q1))
    for ai in a[1:]:
        p0, q0, p1, q1 = p1, q1, ai * p1 + p0, ai * q1 + q0
        out.append((p1, q1))
    return out


def sturm_excluded(g):
    """True iff a GENERIC first passage above the start at lag g >= 2 is impossible:
    floor(g log2 3) - floor((g-1) log2 3) = 1 (exact via integer powers)."""
    fl = lambda m: max(K for K in range(0, 2 * m + 2) if 2 ** K < 3 ** m) if m > 0 else 0
    return g >= 2 and fl(g) - fl(g - 1) == 1


def sturm_dp(g, kcap=14):
    """independent check: does a word exist with 2^K_m > 3^m for 1 <= m < g and 2^K_g < 3^g ?"""
    S = {0}
    for m in range(1, g):
        S = {K + k for K in S for k in range(1, kcap) if 2 ** (K + k) > 3 ** m}
        S = {K for K in S if K < m * LOG23 + kcap}
    return any(2 ** (K + 1) < 3 ** g for K in S)


def clock_type(K, L, conv):
    for i, (p, q) in enumerate(conv):
        if (p, q) == (K, L):
            return 'convergent p%d/q%d' % (i, i)
    for i in range(1, len(conv) - 1):
        (pm, qm), (p, q) = conv[i - 1], conv[i]
        for j in range(1, 60):
            if (pm + j * p, qm + j * q) == (K, L):
                return 'intermediate (p%d+%dp%d)/(q%d+%dq%d)' % (i - 1, j, i, i - 1, j, i)
    return 'not a best approximation'


# ----------------------------------------------------------------------------
# S0: C engine
# ----------------------------------------------------------------------------
def build():
    os.makedirs(SCR, exist_ok=True)
    if (not os.path.exists(BIN)) or os.path.getmtime(BIN) < os.path.getmtime(CSRC):
        subprocess.check_call(['cc', '-O2', '-o', BIN, CSRC, '-lm'])


def run_c(b, N, lagmax, tag):
    pre = os.path.join(SCR, '%s_b%d' % (tag, b))
    if (not os.path.exists(pre + '_pern.bin')) or os.environ.get('ORDER_LAWS_RERUN') == '1' \
            or os.path.getmtime(pre + '_pern.bin') < os.path.getmtime(BIN):
        t = time.time()
        subprocess.check_call([BIN, str(b), str(N), pre, str(lagmax)], stderr=subprocess.DEVNULL)
        log('C run b=%+d N=%d in %.1fs' % (b, N, time.time() - t))
    return pre


PERN_DT = np.dtype([('sigma', '<i2'), ('tau', '<i2'), ('T', '<i2'), ('argmax', '<i2'), ('rho', '<i2'),
                    ('rmax', '<i2'), ('rmin', '<i2'), ('Ksig', '<i2'), ('fmr', '<i2'), ('nD', '<i2'),
                    ('nU', '<i2'), ('KT', '<i2'), ('entry', 'i1'), ('cst', 'i1'), ('pad', 'V6'),
                    ('M', '<u8')])
PERN_FIELDS = ['sigma', 'tau', 'T', 'argmax', 'rho', 'rmax', 'rmin', 'Ksig', 'fmr', 'nD', 'nU', 'KT', 'entry',
               'cst', 'M']


def load_c(pre, b, N):
    D = {'b': b, 'N': N}
    pats, wins = {}, {}
    with open(pre + '_patterns.txt') as f:
        D['head'] = f.readline().strip()
        for line in f:
            p = line.split()
            if p[0] == 'W':
                wins[(int(p[1]), int(p[2]))] = tuple(int(q.split('=')[1]) for q in p[3:6])
            elif p[0] == 'P':
                pats[(int(p[1]), int(p[2]), int(p[3]))] = tuple(int(q) for q in p[4:8])
    D['pats'], D['wins'] = pats, wins
    lt = {}
    with open(pre + '_lags.txt') as f:
        for line in f:
            p = line.split()
            if p[0] == 'L':
                lt[tuple(int(q) for q in p[1:5])] = (int(p[5]), int(p[6]))
    D['lags'] = lt
    D['lagmax'] = max(k[0] for k in lt)
    devs, mar, cst, rec = {}, {}, {}, {}
    with open(pre + '_devs.txt') as f:
        for line in f:
            p = line.split()
            if p[0] == 'D':
                devs[(int(p[1]), int(p[2]))] = (int(p[3]), int(p[4]), int(p[5]), int(p[6]))
    with open(pre + '_margins.txt') as f:
        for line in f:
            p = line.split()
            if p[0] == 'M':
                mar[(int(p[1]), int(p[2]))] = (int(p[3]), int(p[4]), float(p[5]), int(p[6]), int(p[7]))
    with open(pre + '_cst.txt') as f:
        for line in f:
            p = line.split()
            if p[0] == 'C':
                cst[(int(p[1]), int(p[2]))] = tuple(int(v) for v in p[3:])
    with open(pre + '_records.txt') as f:
        for line in f:
            p = line.split()
            if p[0] == 'R':
                rec[(int(p[1]), int(p[2]), p[3], int(p[4]))] = (int(p[5]), int(p[6]))
    D['devs'], D['margins'], D['cst'], D['records'] = devs, mar, cst, rec
    D['pern'] = np.fromfile(pre + '_pern.bin', dtype=PERN_DT)
    D['nn'] = np.arange(1, N + 1, 2, dtype=np.int64)
    check(len(D['pern']) == len(D['nn']), 'pern size')
    D['notC'] = ~np.isin(D['nn'], sorted(CYC[b]))
    return D


# ----------------------------------------------------------------------------
# S1: independent pure-Python recomputation (exact integer comparisons, Fractions)
# ----------------------------------------------------------------------------
def python_census(b, N, lagmax):
    pats = defaultdict(lambda: [0, 0, 0, 0])
    wins = defaultdict(lambda: [0, 0, 0])
    lags = defaultdict(lambda: [0, 0])
    devs, mar, cst = {}, {}, {}
    rec = defaultdict(lambda: [0, 0])
    pern = []
    pow3 = [3 ** i for i in range(2000)]
    for n in range(1, N + 1, 2):
        xs, ks, T = stopped_orbit(n, b, extra=16)
        Kp = [0]
        for k in ks:
            Kp.append(Kp[-1] + k)
        sig = next((j for j in range(1, T + 17) if xs[j] < n), -1)
        tau = next((j for j in range(1, T + 17) if 2 ** Kp[j] > pow3[j]), -1)
        orb = xs[:T + 1]
        M = max(orb)
        am = orb.index(M)
        rho = sum(1 for v in orb[1:] if v > n)
        cmax = cmin = n
        lastmax = lastmin = 0
        rmax = rmin = 0
        fmr = -1
        for j in range(1, T + 1):
            for typ in (0, 1):
                if not (orb[j] > cmax if typ == 0 else orb[j] < cmin):
                    continue
                last = lastmax if typ == 0 else lastmin
                fv = {'val_s16': sres(orb[j], 16, b), 'val_s27': sres(orb[j], 27, b),
                      'pred_s16': sres(orb[j - 1], 16, b), 'pred_s27': sres(orb[j - 1], 27, b),
                      'val_r16': orb[j] % 16, 'val_r27': orb[j] % 27,
                      'pred_r16': orb[j - 1] % 16, 'pred_r27': orb[j - 1] % 27,
                      'kprev': min(ks[j - 1], 63), 'gap': min(j - last, 255), 'idx': min(j, 1023)}
                for s in (0, 1):
                    if s == 1 and j == T:
                        continue
                    for f, v in fv.items():
                        e = rec[(typ, s, f, v)]
                        if e[0] == 0:
                            e[1] = n
                        e[0] += 1
                if typ == 0:
                    cmax, lastmax = orb[j], j
                    rmax += 1
                    if fmr < 0:
                        fmr = j
                else:
                    cmin, lastmin = orb[j], j
                    rmin += 1
        for w in range(3, 7):
            for t in range(0, T - w + 2):
                win = orb[t:t + w]
                ia = lehmer(win)
                Kmax = Kp[t + w - 1] - Kp[t]
                gv = [pow3[i] * 2 ** (Kmax - (Kp[t + i] - Kp[t])) for i in range(w)]   # exact generic values
                ig = lehmer(gv)
                mn = min(win)
                scopes = (True, t == 0, t + w - 1 < T, mn >= 1000)
                for s in range(4):
                    if not scopes[s]:
                        continue
                    e = pats[(w, s, ia)]
                    if e[0] == 0:
                        e[1], e[2] = n, t
                    e[0] += 1
                    e[3] = n
                    ww = wins[(w, s)]
                    ww[0] += 1
                    if ia != ig:
                        if ww[1] == 0:
                            ww[2] = n
                        ww[1] += 1
        nD = nU = 0
        for i in range(T):
            xi = orb[i]
            sig_i = tau_i = -1
            for j in range(i + 1, T + 1):
                l = j - i
                K = Kp[j] - Kp[i]
                p2 = 1 << K
                g = pow3[l] > p2
                xj = orb[j]
                up = xj > xi
                if sig_i < 0 and xj < xi:
                    sig_i = j
                if tau_i < 0 and not g:
                    tau_i = j
                if l <= lagmax:
                    e = lags[(l, 0, int(g), int(up))]
                    if e[0] == 0:
                        e[1] = n
                    e[0] += 1
                    if i == 0:
                        e = lags[(l, 1, int(g), int(up))]
                        if e[0] == 0:
                            e[1] = n
                        e[0] += 1
                if g != up:
                    if g:
                        nD += 1
                    else:
                        nU += 1
                    key = (xi, l)
                    if key not in devs:
                        devs[key] = [xj, K, 1 if j == T else 0, n]
                    if i == 0:
                        devs[key][2] |= 2
                if ((not g) if b > 0 else g):
                    den = p2 * xj - pow3[l] * xi          # = B' (sign b)
                    num = (xi - xj) * p2
                    below = (num < den) if den > 0 else (num > den)   # mu = num/den < 1
                    if below and (xi, l) not in mar:
                        mar[(xi, l)] = (xj, K, Fraction(num, den), 1 if j == T else 0, n)
            if sig_i >= 0 and tau_i >= 0 and sig_i != tau_i:
                key = (xi, sig_i - i)
                if key not in cst:
                    cst[key] = (orb[sig_i], tau_i - i, 2 if i == 0 else 0, n)
        ent = 0 if orb[T] == 1 else (1 if orb[T] in (5, 7) else 2)
        pern.append((sig, tau, T, am, rho, rmax, rmin, Kp[sig] if sig > 0 else -1, fmr,
                     min(nD, 32767), min(nU, 32767), Kp[T], ent, 1 if sig == tau else 0, M))
    return dict(pats=pats, wins=wins, lags=lags, devs=devs, margins=mar, cst=cst, records=rec, pern=pern)


def compare_census(b, C, P):
    bad = []
    for key, v in C['pats'].items():
        pv = P['pats'].get(key, [0, 0, 0, 0])
        if v[0] != pv[0] or (v[0] and tuple(v[1:]) != tuple(pv[1:])):
            bad.append(('pat', key))
    for key, v in C['wins'].items():
        pv = P['wins'].get(key, [0, 0, 0])
        if v[0] != pv[0] or v[1] != pv[1] or (v[1] and v[2] != pv[2]):
            bad.append(('win', key))
    for key, v in C['lags'].items():
        if key[0] > LAG_SEARCH:
            continue
        pv = P['lags'].get(key, [0, 0])
        if v[0] != pv[0] or (v[0] and v[1] != pv[1]):
            bad.append(('lag', key))
    if set(C['devs']) != set(P['devs']):
        bad.append(('devset',))
    else:
        for k, v in C['devs'].items():
            if tuple(v) != tuple(P['devs'][k]):
                bad.append(('dev', k))
    if set(C['margins']) != set(P['margins']):
        bad.append(('marset', len(C['margins']), len(P['margins'])))
    else:
        for k, v in C['margins'].items():
            pv = P['margins'][k]
            if abs(v[2] - float(pv[2])) > 1e-6 or v[0] != pv[0] or v[1] != pv[1] or v[4] != pv[4]:
                bad.append(('mar', k))
    if set(C['cst']) != set(P['cst']):
        bad.append(('cst',))
    nz = {k: v for k, v in P['records'].items() if v[0]}
    if set(C['records']) != set(nz):
        bad.append(('recset',))
    else:
        for k, v in C['records'].items():
            if tuple(v) != tuple(nz[k]):
                bad.append(('rec', k))
    pn = C['pern']
    if len(pn) != len(P['pern']):
        bad.append(('pernlen',))
    else:
        cols = {f: pn[f].astype(np.int64) if f != 'M' else pn[f] for f in PERN_FIELDS}
        for i, row in enumerate(P['pern']):
            for f, v in zip(PERN_FIELDS, row):
                if int(cols[f][i]) != int(v):
                    bad.append(('pern', i, f))
                    break
    return bad


# ----------------------------------------------------------------------------
# S2: grammar
# ----------------------------------------------------------------------------
SCOPES = ['all windows', 'initial windows (t=0)', 'interior windows (x_T excluded)',
          'large windows (all values >= 1000)']
RRANGE = {'val_s16': 16, 'val_s27': 27, 'pred_s16': 16, 'pred_s27': 27, 'val_r16': 16, 'val_r27': 27,
          'pred_r16': 16, 'pred_r27': 27, 'kprev': 41, 'gap': 65, 'idx': 129}
SF = ['sigma', 'tau', 'T', 'argmax', 'rho', 'rmax', 'rmin', 'fmr']
XF = ['sigma', 'tau', 'T', 'argmax', 'rho', 'rmax', 'rmin', 'fmr', 'M']
XP = {'SP': 'm = n + 2^(K_sigma(n)+1) (same word through the stopping time of n)',
      'SPT': 'm = n + 2^(K_T(n)+1) (same stopped word as n)',
      'R': 'm = R_b(n) = 4n + b (U_b(m) = U_b(n))',
      'RAW': 'm = 4n + 1 on both sheets (NON-equivariant control)',
      'ADJ': 'm = n + 2'}
CST_STATEMENTS = [
    ('sigma=tau', 'sigma(n) = tau(n): the first descent below n happens exactly at the first decay prefix'),
    ('sigma>=tau', 'sigma(n) >= tau(n): no descent below n before the first decay prefix'),
    ('sigma<=tau', 'sigma(n) <= tau(n): the first decay prefix always descends below n'),
    ('sigma<inf', 'sigma(n) < infinity for n not in C_b'),
    ('tau<inf', 'tau(n) < infinity for n not in C_b'),
    ('orbit-points', 'sigma = tau at every point x_i of every stopped orbit where both occur inside it')]


def grammar(gen):
    G = []
    for w in range(3, 7):
        for s in range(4):
            for i in range(math.factorial(w)):
                G.append(dict(fam='P', id='P w%d s%d %s' % (w, s, ptxt(i, w)), w=w, s=s, idx=i,
                              generic=(i in gen[w]),
                              text='no %s of %d consecutive odd iterates has ordinal pattern %s'
                                   % (SCOPES[s], w, ptxt(i, w))))
    for l in range(1, LAG_SEARCH + 1):
        for fs in (0, 1):
            where = '(x_l vs the start n)' if fs else '(x_(i+l) vs x_i, any i)'
            for (g, u, txt) in ((1, 0, 'every GROWTH window (3^l > 2^K) of lag %d goes up' % l),
                                (0, 1, 'every DECAY window (2^K > 3^l) of lag %d goes down' % l),
                                (None, 1, 'no lag-%d comparison goes up' % l),
                                (None, 0, 'no lag-%d comparison goes down' % l)):
                G.append(dict(fam='L', id='L l%d %s %s' % (l, 'start' if fs else 'any',
                                                           {(1, 0): 'growth-up', (0, 1): 'decay-down',
                                                            (None, 1): 'never-up', (None, 0): 'never-down'}[(g, u)]),
                              l=l, fs=fs, g=g, u=u, text=txt + ' ' + where))
    for typ in ('D', 'U'):
        for fs in (0, 1):
            G.append(dict(fam='L', id='L all %s %s' % ('start' if fs else 'any',
                                                       'growth-up' if typ == 'D' else 'decay-down'),
                          l=0, fs=fs, g=(1 if typ == 'D' else 0), u=(0 if typ == 'D' else 1),
                          text=('every growth window (any lag) goes up' if typ == 'D'
                                else 'every decay window (any lag) goes down') +
                               (' (from the start)' if fs else ' (anywhere)')))
    for c in (0.0, 0.01, 0.03, 0.1, 0.3, 0.5, 1.0):
        G.append(dict(fam='G', id='G mu>=%g' % c, c=c,
                      text='every window with a POSITIVE gate g = B/(2^K-3^l) has x_i %s (1+%g) g'
                           % ('>' if c == 0 else '>=', c)))
    for typ in (0, 1):
        for s in (0, 1):
            for f, rng in RRANGE.items():
                for v in range(rng):
                    if f in ('kprev', 'gap', 'idx') and v == 0:
                        continue
                    if f.endswith('16') and v % 2 == 0:
                        continue
                    if f.endswith('27') and v % 3 == 0:
                        continue
                    G.append(dict(fam='R', id='R %s %s %s=%d' % (('max', 'min')[typ], ('all', 'int')[s], f, v),
                                  typ=typ, s=s, f=f, v=v,
                                  text='no running-%s record x_j (%s) has %s = %d'
                                       % (('max', 'min')[typ], ('j >= 1', '1 <= j < T')[s], f, v)))
    for f in SF:
        for v in range(-1, 65):
            G.append(dict(fam='S', id='S %s=%d' % (f, v), f=f, v=v,
                          text='no odd n outside C_b has %s(n) = %s' % (f, 'infinity' if v == -1 and f in
                                                                         ('sigma', 'tau') else v)))
    for p in XP:
        for f in XF:
            for rel in ('le', 'ge', 'eq'):
                G.append(dict(fam='X', id='X %s %s %s' % (p, f, rel), ptype=p, f=f, rel=rel,
                              text='%s: %s(m) %s %s(n)' % (XP[p], f, {'le': '<=', 'ge': '>=', 'eq': '='}[rel], f)))
    for cid, txt in CST_STATEMENTS:
        G.append(dict(fam='C', id='C ' + cid, cid=cid, text=txt))
    return G


class Evaluator:
    """evaluates grammar entries on one C dataset (one sheet, one N)."""

    def __init__(self, D):
        self.D = D
        self.b = D['b']
        P = D['pern']
        nn, notC = D['nn'], D['notC']
        # S family: first occurrence of every value
        self.sfirst = {}
        for f in SF:
            vals = P[f][notC]
            u, first = np.unique(vals, return_index=True)
            nsel = nn[notC]
            self.sfirst[f] = {int(a): int(nsel[i]) for a, i in zip(u, first)}
        # X family
        self.x = {}
        N = D['N']
        for p in XP:
            if p == 'SP':
                ok = (P['sigma'] > 0) & (P['Ksig'] + 1 <= 40)
                sh = np.clip(P['Ksig'].astype(np.int64) + 1, 0, 40)
                m = nn + np.where(ok, np.left_shift(np.int64(1), sh), 0)
            elif p == 'SPT':
                ok = (P['KT'] + 1 <= 40)
                sh = np.clip(P['KT'].astype(np.int64) + 1, 0, 40)
                m = nn + np.where(ok, np.left_shift(np.int64(1), sh), 0)
            elif p == 'R':
                ok = np.ones(len(nn), bool)
                m = 4 * nn + self.b
            elif p == 'RAW':
                ok = np.ones(len(nn), bool)
                m = 4 * nn + 1
            else:
                ok = np.ones(len(nn), bool)
                m = nn + 2
            ok = ok & (m <= N) & notC
            mi = np.clip((m - 1) // 2, 0, len(nn) - 1)
            ok = ok & notC[mi]
            for f in XF:
                okf = ok
                if f in ('sigma', 'tau', 'fmr'):
                    okf = ok & (P[f] > 0) & (P[f][mi] > 0)
                a = P[f].astype(np.int64) if f != 'M' else P[f].astype(np.float64)
                am = a[mi]
                for rel in ('le', 'ge', 'eq'):
                    good = (am <= a) if rel == 'le' else (am >= a) if rel == 'ge' else (am == a)
                    viol = okf & ~good
                    if viol.any():
                        j = int(np.argmax(viol))
                        self.x[(p, f, rel)] = (False, (int(nn[j]), int(m[j])), int(okf.sum()), int(viol.sum()))
                    else:
                        self.x[(p, f, rel)] = (True, None, int(okf.sum()), 0)
        # C family
        self.c = {}
        sg, ta = P['sigma'], P['tau']
        tests = {
            'sigma=tau': sg != ta,
            'sigma>=tau': (sg > 0) & (((ta > 0) & (sg < ta)) | (ta < 0)),
            'sigma<=tau': (ta > 0) & ((sg < 0) | (sg > ta)),
            'sigma<inf': sg < 0,
            'tau<inf': ta < 0}
        for cid, viol in tests.items():
            viol = viol & notC
            self.c[cid] = (not viol.any(), int(nn[np.argmax(viol)]) if viol.any() else None, int(viol.sum()))
        cs = D['cst']
        self.c['orbit-points'] = (len(cs) == 0, min(cs)[0] if cs else None, len(cs))

    def __call__(self, g):
        D, b = self.D, self.b
        fam = g['fam']
        if fam == 'P':
            c = D['pats'][(g['w'], g['s'], g['idx'])]
            return (c[0] == 0, None if c[0] == 0 else (c[1], c[2]), c[0])
        if fam == 'L':
            if g['l'] == 0:
                typ = 'D' if g['g'] == 1 else 'U'
                sel = sorted(k for k, v in D['devs'].items()
                             if growth(k[1], v[1]) == (typ == 'D') and (not g['fs'] or k[0] <= D['N']))
                return (not sel, sel[0] if sel else None, len(sel))
            lt = D['lags']
            if g['g'] is None:
                cs = [lt[(g['l'], g['fs'], gg, g['u'])] for gg in (0, 1)]
                cnt = cs[0][0] + cs[1][0]
                fw = min([c[1] for c in cs if c[0]] or [0])
            else:
                cnt, fw = lt[(g['l'], g['fs'], g['g'], g['u'])]
            return (cnt == 0, fw if cnt else None, cnt)
        if fam == 'G':
            c = g['c']
            bad = sorted((v[2], k) for k, v in D['margins'].items() if (v[2] <= 0 if c == 0 else v[2] < c))
            return (not bad, bad[0][1] if bad else None, len(bad))
        if fam == 'R':
            e = D['records'].get((g['typ'], g['s'], g['f'], g['v']))
            return (e is None, None if e is None else e[1], 0 if e is None else e[0])
        if fam == 'S':
            fn = self.sfirst[g['f']].get(g['v'])
            return (fn is None, fn, 0 if fn is None else 1)
        if fam == 'X':
            h, w, npairs, nv = self.x[(g['ptype'], g['f'], g['rel'])]
            return (h, w, nv)
        if fam == 'C':
            return self.c[g['cid']]
        raise ValueError(fam)


# ----------------------------------------------------------------------------
# witness construction for generic patterns
# ----------------------------------------------------------------------------
def window_ok(n, b, w, scope, target):
    xs, ks, T = stopped_orbit(n, b)
    if T < w - 1 or (scope == 2 and T == w - 1):
        return False
    win = xs[:w]
    if scope == 3 and min(win) < 1000:
        return False
    return lehmer(win) == target


def construct_witness(w, scope, idx, b, gen, maxwords=2000):
    best = None
    for word in gen[w].get(idx, [])[:maxwords]:
        r, mod = word_class(list(word), b)
        for j in range(0, 400):
            n = r + j * mod
            if best is not None and n >= best[0]:
                break
            if window_ok(n, b, w, scope, idx):
                best = (n, word)
                break
    return best


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    print('=' * 110)
    print('collatz_procgen_20260922_order_laws: procedural search for sign-specific ORDER laws of Syracuse orbits')
    print('U_b(x) = oddpart(3x+b), b = +1 (plus sheet) / -1 (minus sheet), positive odd x; stopped orbits x_0..x_T')
    print('(T = first index with x_T in C_b; C_+ = {1}; C_- = {1,5,7,17,25,37,41,55,61,91}); word k_i = v_2(3x_i+b),')
    print('K_j = k_0+...+k_(j-1); signed residues b*x mod M.  A window (i, i+l) is GROWTH if 3^l > 2^K, DECAY if 2^K > 3^l')
    print('(K = K_(i+l) - K_i); its GATE is g = B/(2^K - 3^l) where 2^K x_(i+l) = 3^l x_i + B, sign(B) = b.')
    print('=' * 110)
    build()
    R = {}
    for tag, N, lag in (('chk', N_CHECK, LAG_SEARCH), ('srch', N_SEARCH, LAG_SEARCH), ('ext', N_EXT, LAG_EXT)):
        for b in (1, -1):
            R[(tag, b)] = load_c(run_c(b, N, lag, tag), b, N)
    log('C outputs loaded')

    print('\n## S0. C census runs')
    for tag in ('chk', 'srch', 'ext'):
        for b in (1, -1):
            print('  %-4s %-5s %s' % (tag, SHEET[b], R[(tag, b)]['head']))

    # ---------------- S1 ----------------
    print('\n## S1. Independent pure-Python recomputation at N = %d (exact integers / Fractions)' % N_CHECK)
    for b in (1, -1):
        t = time.time()
        P = python_census(b, N_CHECK, LAG_SEARCH)
        bad = compare_census(b, R[('chk', b)], P)
        check(not bad, 'census mismatch b=%d: %s' % (b, bad[:10]))
        print('  %s: patterns (w=3..6 x 4 scopes, counts + first/last witness), lag tables (l <= %d), %d gate '
              'crossings, %d near-gate windows (mu < 1), %d sigma/tau disagreements, %d record-histogram cells, '
              '%d per-n rows: IDENTICAL' % (SHEET[b], LAG_SEARCH, len(P['devs']), len(P['margins']), len(P['cst']),
                                             len([1 for v in P['records'].values() if v[0]]), len(P['pern'])))
        log('python census b=%d %.1fs' % (b, time.time() - t))
    # random-sample recomputation of the 10^7 per-n rows (second code path for the extension data)
    import random
    rng = random.Random(20260922)
    for b in (1, -1):
        D = R[('ext', b)]
        sample = sorted(set(2 * rng.randrange(N_EXT // 2 + 1) + 1 for _ in range(3000)))
        bad = 0
        for n in sample:
            xs, ks, T = stopped_orbit(n, b, extra=16)
            Kp = [0]
            for k in ks:
                Kp.append(Kp[-1] + k)
            sig = next((j for j in range(1, T + 17) if xs[j] < n), -1)
            tau = next((j for j in range(1, T + 17) if 2 ** Kp[j] > 3 ** j), -1)
            orb = xs[:T + 1]
            row = D['pern'][(n - 1) // 2]
            if (int(row['sigma']), int(row['tau']), int(row['T']), int(row['M']), int(row['argmax']),
                    int(row['KT'])) != (sig, tau, T, max(orb), orb.index(max(orb)), Kp[T]):
                bad += 1
        check(bad == 0, 'random sample of 10^7 rows b=%d' % b)
        print('  %s: %d random odd n <= 10^7: sigma, tau, T, max, argmax, K_T recomputed in Python: identical'
              % (SHEET[b], len(sample)))
    # the crude gate bound dominates the true gate of every word (brute force, l <= 7, k_i <= 10)
    for b in (1, -1):
        for l in range(1, 7):
            G_ = gate_bound(l, b)
            worst = Fraction(0)
            for word in itertools.product(range(1, 12), repeat=l):
                K = sum(word)
                D_ = 2 ** K - 3 ** l
                if (b > 0 and D_ > 0) or (b < 0 and D_ < 0):
                    worst = max(worst, Fraction(abs(carry(word, b)), abs(D_)))
            check(worst <= G_, 'gate bound l=%d b=%d' % (l, b))
    print('  gate bound G_b(l) >= max |B|/|2^K - 3^l| over all words with k_i <= 11, l <= 6, both sheets: checked')
    # cross-orbit identities used in S3 (independent spot checks, n <= N_CHECK)
    for b in (1, -1):
        nsp = nr = 0
        for n in range(3, N_CHECK + 1, 2):
            if n in CYC[b]:
                continue
            m = 4 * n + b
            check(U(m, b)[0] == U(n, b)[0], 'R_b shares the image')
            nr += 1
            xs, ks, T = stopped_orbit(n, b, extra=16)
            sig = next((j for j in range(1, len(xs)) if xs[j] < n), None)
            if sig is None:
                continue
            K = sum(ks[:sig])
            m = n + 2 ** (K + 1)
            ym, km, Tm = stopped_orbit(m, b, extra=16)
            check(km[:sig] == ks[:sig], 'SP pair shares the word through sigma(n)')
            Kj = 0
            for j in range(1, sig + 1):
                Kj += ks[j - 1]
                check(ym[j] - xs[j] == 3 ** j * (m - n) // 2 ** Kj and ym[j] > xs[j], 'SP order preserved')
            nsp += 1
        print('  %s: U_b(4n+b) = U_b(n) for %d n; SP pairs (n, n + 2^(K_sigma+1)): same word through sigma(n) and '
              'x_j(m) - x_j(n) = 3^j (m-n)/2^K_j > 0 for all j <= sigma(n), %d pairs' % (SHEET[b], nr, nsp))

    gen = generic_sets()
    cfa = cf_log23(16)
    conv = convergents(cfa)
    print('  log2 3 = [%s; %s]; convergents %s' % (cfa[0], ', '.join(map(str, cfa[1:])),
                                                  ', '.join('%d/%d' % pq for pq in conv[:11])))
    ngen = {w: len(gen[w]) for w in range(3, 7)}
    print('  generic (word-realizable, k_i <= 9 suffices) ordinal patterns: %s of %s'
          % (ngen, {w: math.factorial(w) for w in range(3, 7)}))

    # ---------------- S2 ----------------
    G = grammar(gen)
    EV = {(tag, b): Evaluator(R[(tag, b)]) for tag in ('srch', 'ext') for b in (1, -1)}
    for g in G:
        for tag in ('srch', 'ext'):
            for b in (1, -1):
                g[(tag, b)] = EV[(tag, b)](g)
    log('grammar evaluated: %d candidates' % len(G))

    def side(g, tag):
        hp, hm = g[(tag, 1)][0], g[(tag, -1)][0]
        return 'both' if hp and hm else 'plus-only' if hp else 'minus-only' if hm else 'neither'

    FAMS = ['P', 'L', 'G', 'R', 'S', 'X', 'C']
    FAMNAME = {'P': 'ordinal patterns of windows, w = 3..6, 4 scopes',
               'L': 'lag comparisons x_(i+l) vs x_i (l <= 64; from start / anywhere; word-conditioned or not) + all-lag',
               'G': 'gate margins mu = x_i/g - 1 on positive-gate windows (equivariant)',
               'R': 'running max/min records: signed/raw residues mod 16, 27; k_(j-1); gap; index',
               'S': 'per-n features NEVER(v): sigma, tau, T, argmax, rho, #max/#min records, first record',
               'X': 'cross-orbit pairs (same prefix SP/SPT, R_b(n) = 4n+b, raw 4n+1, n+2) x 9 features x <=,>=,=',
               'C': 'coefficient stopping time sigma vs tau; finiteness'}
    print('\n## S2. Grammar and evaluation at N = %d (all odd n <= N, both sheets)' % N_SEARCH)
    print('  family | description | candidates | plus-only | minus-only | both | neither')
    tot = Counter()
    for fam in FAMS:
        cnt = Counter(side(g, 'srch') for g in G if g['fam'] == fam)
        n = sum(cnt.values())
        tot.update(cnt)
        tot['all'] += n
        print('  %s | %s | %d | %d | %d | %d | %d' % (fam, FAMNAME[fam], n, cnt['plus-only'], cnt['minus-only'],
                                                      cnt['both'], cnt['neither']))
    print('  TOTAL | | %d | %d | %d | %d | %d' % (tot['all'], tot['plus-only'], tot['minus-only'], tot['both'],
                                                  tot['neither']))
    print('  grammar parameters: P: w in {3,4,5,6}, all w! patterns, scopes {all, t=0, x_T excluded, values>=1000};')
    print('  L: l = 1..64, scope {any i, i = 0}, conditions {growth->up, decay->down, never up, never down}, plus the')
    print('  all-lag growth/decay statements; G: c in {0, .01, .03, .1, .3, .5, 1}; R: record type {max, min} x scope')
    print('  {j>=1, j<T} x feature {signed / raw value and predecessor mod 16 (odd v) and mod 27 (3 !| v), k_(j-1) =')
    print('  1..40, gap = 1..64, index = 1..128}; S: 8 features x v = -1..64 (n outside C_b); X: 5 pair types x 9')
    print('  features x {<=, >=, =}; C: 6 statements.')
    for s in range(4):
        for w in range(3, 7):
            for b in (1, -1):
                realized = set(i for i in range(math.factorial(w)) if R[('srch', b)]['pats'][(w, s, i)][0] > 0)
                check(realized == set(gen[w]), 'pattern family: realized set = generic set (w=%d s=%d b=%d)'
                      % (w, s, b))
    print('  P family: for every w and scope, the set of realized patterns on EACH sheet equals the set of generic')
    print('  (word-realizable) patterns exactly (checked): %s generic patterns for w = 3..6.'
          % [len(gen[w]) for w in range(3, 7)])

    # ---------------- S3: classification ----------------
    print('\n## S3. Classification of the survivors (search range N = %d; extension N = %d)' % (N_SEARCH, N_EXT))
    for g in G:
        g['side'] = side(g, 'srch')
        g['side7'] = side(g, 'ext')
    labels = {}
    for g in G:
        sd = g['side']
        if sd == 'neither':
            continue
        fam = g['fam']
        lab = None
        if fam == 'P':
            if sd == 'both':
                lab = ('WORD-FORBIDDEN' if not g['generic'] else 'SMALL-N-BOTH')
            else:
                zero = 1 if sd == 'plus-only' else -1
                if g['generic']:
                    wtn = construct_witness(g['w'], g['s'], g['idx'], zero, gen)
                    g['constructed'] = wtn
                    lab = 'SMALL-N'
                else:
                    lab = 'UNEXPLAINED'
        elif fam == 'L':
            if sd == 'both':
                lab = 'NO-CROSSING-AT-THIS-LAG' if g['g'] is not None else 'TRIVIAL'
            elif sd == 'plus-only' and g['g'] == 1:
                lab = 'SIGN-LAW'
            elif sd == 'minus-only' and g['g'] == 0:
                lab = 'SIGN-LAW'
            else:
                lab = 'UNEXPLAINED'
        elif fam == 'G':
            lab = 'UNEXPLAINED' if sd != 'both' else 'BOTH'
        elif fam == 'R':
            if sd == 'both':
                lab = 'BOTH'
            elif g['f'].startswith(('val_r', 'pred_r')):
                lab = 'NORMALIZATION'
            elif g['f'] == 'gap' and g['typ'] == 0 and sturm_excluded(g['v']):
                lab = 'SIGN-LAW-STURMIAN'
            else:
                lab = 'SMALL-N' if g['side7'] != sd else 'UNRESOLVED'
        elif fam == 'S':
            if sd == 'both':
                lab = 'BOTH'
            elif g['f'] == 'fmr' and sturm_excluded(g['v']):
                lab = 'SIGN-LAW-STURMIAN'
            elif g['f'] == 'rmin' and g['v'] == 0:
                lab = 'TERMINAL'
            else:
                lab = 'SMALL-N' if g['side7'] != sd else 'UNRESOLVED'
        elif fam == 'X':
            if sd == 'both':
                lab = 'BOTH'
            elif g['ptype'] == 'RAW':
                lab = 'NORMALIZATION'
            elif g['ptype'] == 'SPT' and g['f'] in ('sigma', 'tau'):
                other = -1 if sd == 'plus-only' else 1
                n0 = g[('srch', other)][1][0]
                P0 = R[('srch', other)]['pern'][(n0 - 1) // 2]
                lab = 'TERMINAL' if int(P0['sigma']) > int(P0['T']) else 'UNRESOLVED'
            else:
                lab = 'SMALL-N' if g['side7'] != sd else 'UNRESOLVED'
        elif fam == 'C':
            lab = 'BOTH' if sd == 'both' else 'UNRESOLVED'
        g['label'] = lab
        labels.setdefault((sd, fam, lab), []).append(g)
    print('  side | family | label | count')
    for k in sorted(labels):
        print('  %s | %s | %s | %d' % (k[0], k[1], k[2], len(labels[k])))

    print('\n  label key: SIGN-LAW = holds on one sheet by sign(B) = b applied to a window (PROVED there), fails on')
    print('  the other at gate crossings (S4); SIGN-LAW-STURMIAN = a record gap / first-record index g that no WORD')
    print('  can produce (floor(g log2 3) - floor((g-1) log2 3) = 1, S4b), realized only through a crossing;')
    print('  TERMINAL = fails only because a minus orbit enters C_- above its start (names the cycles); SMALL-N = the')
    print('  other sheet realizes it by 10^7 (explicit n); NORMALIZATION = raw residue / raw 4n+1, whose signed')
    print('  (equivariant) form is sheet-blind; WORD-FORBIDDEN = pattern produced by no word, and windows this short')
    print('  are generic (S4); UNRESOLVED = still one-sided at 10^7 with no explanation; BOTH = holds on both sheets.')
    unres = [g for g in G if g.get('label') == 'UNRESOLVED' or g.get('label') == 'UNEXPLAINED']
    print('  UNRESOLVED / UNEXPLAINED survivors: %d %s' % (len(unres), [g['id'] for g in unres]))
    # verification of the STURMIAN and TERMINAL labels on their witnesses
    devx = {b: R[('ext', b)]['devs'] for b in (1, -1)}
    for g in G:
        if g.get('label') == 'SIGN-LAW-STURMIAN':
            other = -1 if g['side'] == 'plus-only' else 1
            n0 = g[('srch', other)][1]
            xs, ks, T = stopped_orbit(n0, other)
            cm, recs = xs[0], [0]
            for j in range(1, T + 1):
                if xs[j] > cm:
                    cm = xs[j]
                    recs.append(j)
            gv = g['v']
            if g['fam'] == 'R':
                if g['s'] == 1:
                    pairs = [(i, j) for i, j in zip(recs, recs[1:]) if j - i == gv and j < T]
                else:
                    pairs = [(i, j) for i, j in zip(recs, recs[1:]) if j - i == gv]
            else:
                pairs = [(0, recs[1])] if len(recs) > 1 and recs[1] == gv else []
            check(pairs, 'sturmian witness has the gap')
            i, j = pairs[0]
            check((xs[i], gv) in devx[other], 'sturmian witness window is a gate crossing')
            g['sturm_window'] = (n0, xs[i], gv, xs[j])
        if g.get('label') == 'TERMINAL':
            other = -1 if g['side'] == 'plus-only' else 1
            check(other == -1, 'terminal failures are minus-sheet')
    st = [g for g in G if g.get('label') == 'SIGN-LAW-STURMIAN']
    print('  SIGN-LAW-STURMIAN witnesses verified (witness n, record x_i, gap g, next record): %s'
          % sorted(set(g['sturm_window'] for g in st)))
    P = R[('ext', -1)]['pern']
    nn = R[('ext', -1)]['nn']
    tset = nn[(P['rmin'] == 0) & R[('ext', -1)]['notC']]
    print('  TERMINAL: the minus orbits (n <= 10^7) with no new minimum before entering C_- are n = %s; each enters '
          'C_- at an element larger than n (entries %s)'
          % (list(map(int, tset)), [stopped_orbit(int(v), -1)[0][-1] for v in tset]))
    for v in tset:
        xs, ks, T = stopped_orbit(int(v), -1)
        check(xs[-1] > v and min(xs) == v, 'terminal explanation')

    for sd in ('plus-only', 'minus-only'):
        other = -1 if sd == 'plus-only' else 1
        zero = -other
        print('\n  --- %s survivors (true on the %s sheet for all odd n <= %d, false on the %s sheet) ---'
              % (sd, SHEET[zero], N_SEARCH, SHEET[other]))
        for fam in FAMS:
            gs = [g for g in G if g['side'] == sd and g['fam'] == fam]
            if not gs:
                continue
            print('  [%s] %d' % (fam, len(gs)))
            if fam == 'P':
                for g in gs:
                    wt = g.get('constructed')
                    e7 = g[('ext', zero)]
                    print('    %-24s %s witness n=%s(t=%s); %s sheet: 10^7 first n=%s, constructed n=%s word=%s'
                          % (g['id'], SHEET[other], g[('srch', other)][1][0], g[('srch', other)][1][1], SHEET[zero],
                             e7[1][0] if e7[1] else 'none', wt[0] if wt else None, wt[1] if wt else None))
            elif fam in ('L', 'G', 'C'):
                for g in gs:
                    print('    %-26s %-10s first %s witness %s (%d failures); 10^7: %s'
                          % (g['id'], g['label'], SHEET[other], g[('srch', other)][1], g[('srch', other)][2],
                             g['side7']))
            elif fam == 'X':
                for g in gs:
                    print('    %-22s %-14s %s witness (n,m)=%s; %d violations; 10^7: %s'
                          % (g['id'], g['label'], SHEET[other], g[('srch', other)][1], g[('srch', other)][2],
                             g['side7']))
            else:
                lab = Counter(g['label'] for g in gs)
                print('    labels: %s' % dict(lab))
                for g in [g for g in gs if g['label'] not in ('NORMALIZATION',)][:40]:
                    print('    %-26s %-12s %s witness n=%s; 10^7 on the %s sheet: %s'
                          % (g['id'], g['label'], SHEET[other], g[('srch', other)][1], SHEET[zero],
                             g[('ext', zero)][1]))
    print('\n  --- both-sheet survivors that are order statements, not word functions by construction ---')
    print('    (X: pairs tested at 10^6 / 10^7 per sheet; SP/SPT pairs exist only when m <= N, i.e. short shared words)')
    for g in G:
        if g['side'] == 'both' and g['fam'] in ('C', 'X', 'G'):
            extra = ''
            if g['fam'] == 'X':
                extra = ' [pairs %s]' % ', '.join('%s %d/%d' % (SHEET[b], EV[('srch', b)].x[(g['ptype'], g['f'],
                                                                                           g['rel'])][2],
                                                               EV[('ext', b)].x[(g['ptype'], g['f'], g['rel'])][2])
                                                  for b in (1, -1))
            print('    %-22s 10^7: %-9s %s%s' % (g['id'], g['side7'], g['text'], extra))
    bl = [g for g in G if g['side'] == 'both' and g['fam'] == 'L' and g['g'] is not None and g['l'] > 0]
    print('    [L] %d word-conditioned lag statements hold on both sheets (all lags except the crossing lags)' % len(bl))
    bp = Counter(g['label'] for g in G if g['side'] == 'both' and g['fam'] == 'P')
    print('    [P] %s' % dict(bp))

    # ---------------- S4 ----------------
    print('\n## S4. Deviation spectrum (gate crossings on transients) with exact gates')
    print('  Direction lemma (PROVED): x_(i+l) - x_i = (B - (2^K - 3^l) x_i)/2^K and sign(B) = b.  Hence on the plus')
    print('  sheet a GROWTH window always goes up and only a DECAY window can fail (go up); on the minus sheet a DECAY')
    print('  window always goes down and only a GROWTH window can fail (go down); a failure needs 0 < x_i < gate.')
    for b in (1, -1):
        for tag in ('srch', 'ext'):
            ds = R[(tag, b)]['devs']
            kinds = Counter(('D' if growth(l, v[1]) else 'U', l, v[1]) for (xi, l), v in ds.items())
            print('  %s N=%d: %d distinct crossing windows; (type, l, K): count = %s'
                  % (SHEET[b], R[(tag, b)]['N'], len(ds), dict(sorted(kinds.items()))))
        ds = R[('ext', b)]['devs']
        check(all((not growth(l, v[1])) if b > 0 else growth(l, v[1]) for (xi, l), v in ds.items()),
              'direction lemma violated')
        nex = 0
        for (xi, l), v in sorted(ds.items()):
            xs, ks, T = stopped_orbit(xi, b)
            check(T >= l, 'window outside the stopped orbit of x_i')
            word = ks[:l]
            K = sum(word)
            check(K == v[1] and xs[l] == v[0], 'window data')
            Bw = carry(word, b)
            check(2 ** K * xs[l] == 3 ** l * xi + Bw, 'carry identity')
            gt = Fraction(Bw, 2 ** K - 3 ** l)
            check(gt > xi > 0, 'not a gate crossing')
            nex += 1
        print('  %s: all %d crossings re-verified in exact arithmetic (word, 2^K x_l = 3^l x_0 + B, 0 < x_i < gate)'
              % (SHEET[b], nex))
        for (K, l) in sorted(set((v[1], l) for (xi, l), v in ds.items()), key=lambda t: t[1]):
            ws = sorted(xi for (xi, ll), v in ds.items() if ll == l and v[1] == K)
            gb = gate_bound(l, b)
            print('    clock %d/%d (%s): 2^K/3^l = %.6f; %d windows, x_i in [%d, %d]; G_b(l) = %.4g -> %s'
                  % (K, l, clock_type(K, l, conv), 2 ** K / 3 ** l, len(ws), ws[0], ws[-1], float(gb),
                     'list COMPLETE (PROVED)' if gb <= N_EXT else 'completeness not certified'))
            print('      x_i -> x_(i+l): %s' % ', '.join('%d->%d' % (xi, ds[(xi, l)][0]) for xi in ws))
        Lc = 0
        while gate_bound(Lc + 1, b) <= N_EXT:
            Lc += 1
        lagset = sorted(set(l for (xi, l) in ds))
        first = lagset[0]
        print('  %s: G_b(l) <= 10^7 for every l <= %d (G_b(%d) = %.4g), so the 10^7 census PROVES that the gate '
              'crossings of lag <= %d are exactly the lags %s listed above'
              % (SHEET[b], Lc, Lc + 1, float(gate_bound(Lc + 1, b)), Lc, [l for l in lagset if l <= Lc]))
        print('  %s: short-lag genericity (PROVED): for l < %d every comparison x_(i+l) vs x_i on a transient equals '
              'its word prediction (max G_b(l), l < %d, is %.4g)'
              % (SHEET[b], first, first, float(max(gate_bound(l, b) for l in range(1, first)))))
    print('\n  near-gate windows (0 < mu < 1: not crossings) at 10^7: plus %d, minus %d'
          % (sum(1 for v in R[('ext', 1)]['margins'].values() if v[2] > 0),
             sum(1 for v in R[('ext', -1)]['margins'].values() if v[2] > 0)))
    for b in (1, -1):
        mar = R[('ext', b)]['margins']
        top = sorted((v[2], k, v) for k, v in mar.items() if v[2] > 0)[:4]
        print('    %s smallest positive margins: %s' % (SHEET[b], ', '.join('(%d, l=%d, K=%d, mu=%.4f)'
                                                                          % (k[0], k[1], v[1], m) for m, k, v in top)))
    print('\n  anatomy of the smallest crossings (dip = minimum inside the window):')
    for b in (1, -1):
        ds = R[('ext', b)]['devs']
        for (xi, l) in sorted(ds)[:3]:
            xs, ks, T = stopped_orbit(xi, b)
            win = xs[:l + 1]
            print('    %s: %d -(%d steps, K=%d)-> %d  dip %d at step %d, peak %d; word %s'
                  % (SHEET[b], xi, l, sum(ks[:l]), xs[l], min(win), win.index(min(win)), max(win), ks[:l]))
        dips = Counter()
        for (xi, l) in ds:
            xs, ks, T = stopped_orbit(xi, b)
            win = xs[:l + 1]
            dips[min(win)] += 1
            check(min(win) < xi and min(win) < xs[l], 'crossing without an interior dip')
        print('    %s: every crossing window dips strictly below both endpoints; dip values: %s'
              % (SHEET[b], dict(sorted(dips.items()))))
    xs27 = set(stopped_orbit(27, 1)[0])
    ds = R[('ext', 1)]['devs']
    on27 = sum(1 for (xi, l) in ds if set(stopped_orbit(xi, 1)[0][:l + 1]) & xs27)
    print('  plus: %d of %d crossing windows share a value with the orbit of 27' % (on27, len(ds)))
    for b in (1, -1):
        P = R[('ext', b)]['pern']
        print('  %s: %d of %d odd n <= 10^7 have a gate crossing inside their stopped orbit'
              % (SHEET[b], int(((P['nD'] > 0) | (P['nU'] > 0)).sum()), len(P)))
    # nested crossings: sub-windows of a crossing window that are themselves crossings
    for b in (1, -1):
        ds = R[('ext', b)]['devs']
        nest = Counter()
        for (xi, l) in ds:
            xs, ks, T = stopped_orbit(xi, b)
            subs = tuple(sorted(set((ll) for a in range(0, l) for ll in range(1, l - a + 1)
                                    if (a, ll) != (0, l) and (xs[a], ll) in ds)))
            nest[(l, subs)] += 1
        print('  %s: (lag, lags of crossing sub-windows) -> count: %s' % (SHEET[b], dict(sorted(nest.items()))))

    # ---------------- S4b: Sturmian ladder epochs ----------------
    print('\n## S4b. Ladder epochs: which record gaps can a WORD produce?')
    excl = [g for g in range(2, 129) if sturm_excluded(g)]
    check(all(sturm_excluded(g) == (not sturm_dp(g)) for g in range(2, 65)), 'sturmian DP cross-check')
    print('  Lemma (PROVED): if x_j > x_i and x_m < x_i for i < m < j (a first passage above x_i at lag g = j - i),')
    print('  and no pair in the window is a crossing, then 2^(K_m) > 3^m for 0 < m < g and 2^(K_g) < 3^g with the last')
    print('  exponent 1, which forces floor(g log2 3) - floor((g-1) log2 3) = 2 (g >= 2).  Excluded lags g <= 128')
    print('  (density 2 - log2 3 = 0.415): %s' % excl)
    print('  (cross-checked for g <= 64 by a dynamic program over exponent words: identical)')
    for b in (1, -1):
        gaps = {k[3]: v for k, v in R[('ext', b)]['records'].items() if k[:3] == (0, 0, 'gap')}
        fm = R[('ext', b)]['pern']['fmr']
        nc = R[('ext', b)]['notC']
        fvals = set(int(v) for v in np.unique(fm[nc]))
        print('  %s (10^7): max-record gaps realized: %d values up to %d; EXCLUDED gaps realized: %s; first-record '
              'indices realized at excluded lags: %s'
              % (SHEET[b], len(gaps), max(gaps), {g: gaps[g] for g in sorted(gaps) if g in excl},
                 sorted(v for v in fvals if v in excl)))
    print('  minus sheet, PROVED (computer-assisted): an excluded gap g needs a crossing (x_i, m), m < g, whose start')
    print('  x_i is a running maximum; on the minus sheet the crossings with m <= 35 are exactly (165|309|549, 12)')
    print('  (S4), and a running maximum x_i forces the start n <= x_i <= 549, inside the census.  Hence NO minus')
    print('  orbit (any n) has a max-record gap or first-record index g <= 36 that is excluded; to 10^7 no excluded')
    print('  gap of any size occurs on the minus sheet.  Plus sheet: excluded gaps occur exactly at the crossing')
    print('  lags 17, 29, 41, 46 (witnesses verified in S3).')

    # ---------------- S5 ----------------
    print('\n## S5. Coefficient stopping time: sigma(n) = tau(n) on both sheets')
    for b in (1, -1):
        for tag in ('srch', 'ext'):
            D = R[(tag, b)]
            P, nn, nc = D['pern'], D['nn'], D['notC']
            eq = int(((P['sigma'] == P['tau']) & nc).sum())
            print('  %s N=%d: sigma = tau for %d of %d odd n outside C_b; sigma = infinity for %d of them; '
                  'orbit-point disagreements %d; max sigma %d at n = %d'
                  % (SHEET[b], D['N'], eq, int(nc.sum()), int(((P['sigma'] < 0) & nc).sum()), len(D['cst']),
                     int(P['sigma'].max()), int(nn[np.argmax(P['sigma'])])))
        info = []
        for n0 in sorted(CYC[b]):
            xs, ks, T = stopped_orbit(n0, b, extra=70)
            K = 0
            tau = None
            for j, k in enumerate(ks, 1):
                K += k
                if 2 ** K > 3 ** j:
                    tau = j
                    break
            sig = next((j for j in range(1, len(xs)) if xs[j] < n0), None)
            info.append('%d:(%s,%s)' % (n0, sig, tau))
        print('  %s cycle elements n:(sigma, tau) over 70 steps: %s' % (SHEET[b], ' '.join(info)))
    Lc = {}
    for b in (1, -1):
        L = 0
        while gate_bound(L + 1, b) <= N_EXT:
            L += 1
        Lc[b] = L
    print('  PROVED (computer-assisted): a failure of sigma = tau at n is a gate crossing from the start at lag')
    print('  l = tau(n) (plus: the first decay prefix does not descend) or l = sigma(n) (minus: a growth prefix')
    print('  descends), hence n < G_b(l); the 10^7 census has none, so sigma = tau holds for every odd n >= 3')
    print('  with tau(n) <= %d on the plus sheet and for every odd n outside C_- with sigma(n) <= %d on the minus'
          % (Lc[1], Lc[-1]))
    print('  sheet.  The minus cycles are compatible with the law because their words never become decay words')
    print('  (tau = sigma = infinity at every cycle element); a plus cycle would violate it at its minimum')
    print('  (tau finite because 2^K > 3^L, sigma infinite).')
    # CST margins: how close does a start come to the gate of its stopping window?
    for b in (1, -1):
        D = R[('ext', b)]
        P = D['pern']
        rows = []
        for (x, l), v in D['margins'].items():
            if x > D['N']:
                continue
            r = P[(x - 1) // 2]
            if b > 0 and l == int(r['tau']):
                rows.append((v[2], x, l, v[1], v[0]))
            if b < 0 and (int(r['sigma']) < 0 or l < int(r['sigma'])):
                rows.append((v[2], x, l, v[1], v[0]))
        rows.sort()
        L2 = 0
        while 2 * gate_bound(L2 + 1, b) <= N_EXT:
            L2 += 1
        what = 'first-decay window (0, tau(n))' if b > 0 else 'growth prefix (0, l), l < sigma(n)'
        print('  %s CST margin: windows %s with mu = n/gate - 1 < 1, n <= 10^7: %s; complete for l <= %d '
              '(2 G_b(l) <= 10^7)' % (SHEET[b], what,
                                     ['n=%d l=%d K=%d x_l=%d mu=%.4f' % (x, l, K, xl, m) for m, x, l, K, xl in rows]
                                     or 'none', L2))

    # ---------------- S6 ----------------
    print('\n## S6. Extension of every one-sided survivor to N = %d' % N_EXT)
    agg = Counter()
    for g in G:
        if g['side'] in ('plus-only', 'minus-only'):
            agg[(g['side'], g['fam'], g['label'], 'unchanged' if g['side7'] == g['side'] else 'now ' + g['side7'])] += 1
    for k in sorted(agg):
        print('  %s | %s | %s | %s | %d' % (k + (agg[k],)))
    newly = Counter((g['fam'], g['side7']) for g in G if g['side'] == 'neither' and g['side7'] != 'neither')
    print('  statements that were two-sided failures at 10^6 and one-sided at 10^7: %s (impossible: failures persist)'
          % dict(newly))
    check(not newly, 'monotonicity of failures')
    both7 = Counter((g['fam'], g['side7']) for g in G if g['side'] == 'both' and g['side7'] != 'both')
    print('  both-sheet survivors at 10^6 that fail somewhere at 10^7: %s' % dict(both7))
    for g in G:
        if g['side'] == 'both' and g['side7'] in ('plus-only', 'minus-only'):
            other = -1 if g['side7'] == 'plus-only' else 1
            lab = 'SIGN-LAW-STURMIAN' if ((g['fam'] == 'R' and g['f'] == 'gap' and g['typ'] == 0) or
                                          (g['fam'] == 'S' and g['f'] == 'fmr')) and sturm_excluded(g['v']) \
                else 'TAIL (first realized between 10^6 and 10^7, on one sheet so far)'
            print('    %-24s -> %-10s at 10^7 (first %s witness n = %s)  %s'
                  % (g['id'], g['side7'], SHEET[other], g[('ext', other)][1], lab))
    for b in (1, -1):
        P7 = R[('ext', b)]['pern']
        print('  %s 10^7 maxima: rmax %d, rmin %d, fmr %d, argmax %d, rho %d, T %d'
              % (SHEET[b], P7['rmax'].max(), P7['rmin'].max(), P7['fmr'].max(), P7['argmax'].max(),
                 P7['rho'].max(), P7['T'].max()))

    # ---------------- S7 ----------------
    print('\n## S7. The three best candidates')
    print('  (1) WINDOW SIGN LAW (plus-only).  "Every growth window grows": if 3^l > 2^(K_(i+l)-K_i) then')
    print('      x_(i+l) > x_i.  PLUS: PROVED for all n (2^K x_(i+l) - 3^l x_i = B > 0).  MINUS: false exactly at the')
    print('      gate crossings; smallest witness 165 -(12 odd steps, clock 19/12)-> 163; complete list at lag 12:')
    print('      165, 309, 549 (PROVED complete); lag 53 (clock 84/53): 12 windows, 2981 ... 44433 (to 10^7); no other')
    print('      lag <= 35 (PROVED).  Plus sheet 10^7: 0 failures.  STATUS: PROVED on plus, a restatement of the sign')
    print('      law sign(B) = b on windows; the new content is the anatomy (transient near-cycles at the growth-side')
    print('      convergents 19/12, 84/53; the minus cycles occupy 1/1, 3/2, 11/7).')
    gw = {g['v']: g[('ext', 1)][1] for g in G if g['fam'] == 'R' and g['f'] == 'gap' and g['typ'] == 0
          and g['s'] == 0 and g.get('side7') == 'minus-only'}
    fw = {g['v']: g[('ext', 1)][1] for g in G if g['fam'] == 'S' and g['f'] == 'fmr'
          and g.get('side7') == 'minus-only' and sturm_excluded(g['v'])}
    print('  (2) STURMIAN RECORD LAW (minus-only).  "Consecutive running-max records are never g apart, and the first')
    print('      record is never at index g, for g with floor(g log2 3) - floor((g-1) log2 3) = 1 (g = 3, 5, 8, 10,')
    print('      13, 15, 17, ...)".  MINUS: PROVED for g <= 36 (all n), FINITE-EXACT for every g to 10^7.  PLUS: false;')
    print('      first plus witnesses (gap g: n) %s and (first-record index g: n) %s; every witness window is a plus'
          % (dict(sorted(gw.items())), dict(sorted(fw.items()))))
    print('      gate crossing.  STATUS: the word part is PROVED; the sheet asymmetry is the sign law again (plus')
    print('      crossings are decay windows that go up, and only they can create an excluded first passage).')
    print('  (3) COEFFICIENT STOPPING TIME sigma(n) = tau(n) (BOTH sheets; order-based, not a word function).')
    print('      FINITE-EXACT: all 5000000 plus and 4999991 minus non-cycle odd n <= 10^7, and every orbit point.')
    print('      PROVED (computer-assisted) for tau(n) <= %d (plus) and sigma(n) <= %d (minus).  OPEN in general'
          % (Lc[1], Lc[-1]))
    print('      (Terras\'s coefficient stopping time conjecture on the plus sheet).  Sheet-blind as a statement;')
    print('      combined with the sign law it excludes every nontrivial positive 3n+1 cycle, while the minus cycles')
    print('      satisfy it (their words never turn decay).  It is the order law that the minus cycles do not')
    print('      violate and a plus cycle would.')
    return G, R, gen, conv


if __name__ == '__main__':
    main()
    log('done')
