"""Confirmatory run for THM-4883 / HYP-9476 (kind-pasteur-2026-11-05-S12), predictions P1-P6 of PREREGISTRATION.md.

    python3 verify.py --out VALIDATION_NORMAL.json
    python3 -O verify.py --out VALIDATION_OPTIMIZED.json

The two receipts must agree byte for byte.  Failures of a *prediction* are recorded in the receipt (they are outcomes);
failures of an *identity that the proof asserts* raise RuntimeError (explicit gates, not removable asserts).
Requires numpy.  P6 reads witness_length_3_5.txt and witness_length_8_9.txt produced by witness_length.c.
"""
import sys, json, hashlib, math
from fractions import Fraction
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'dependencies'))
sys.setrecursionlimit(10000)
HEADMAX = 1 << 20


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


# ---------------------------------------------------------------- g, two independent evaluations

def g_array(N):
    """g[1..N] by the defining recurrence, level by level (numpy)."""
    g = np.zeros(N + 1, dtype=np.int64)
    g[1] = 1
    lo = 2
    while lo <= N:
        hi = min(2 * lo - 1, N)
        idx = np.arange(lo, hi + 1, dtype=np.int64)
        ev = idx[idx % 2 == 0]
        g[ev] = g[ev // 2]
        od = idx[idx % 2 == 1]
        m = (od + 1) // 2
        g[od] = g[m] + m - (m % 2)
        lo *= 2
    return g


def g_closed_array(N):
    """Closed bit formula g(n) = f(n-1), f(x) = 1 + sum_{i < top(x), bit_i(x) = 0} (floor(x/2^(i+1)) + bit_(i+1)(x))."""
    x = np.arange(0, N, dtype=np.int64)          # x = n - 1 for n = 1..N
    f = np.ones(N, dtype=np.int64)
    for i in range(int(N).bit_length() + 1):
        q = x >> (i + 1)
        take = (q > 0) & (((x >> i) & 1) == 0)
        f += np.where(take, q + (q & 1), 0)
    out = np.zeros(N + 1, dtype=np.int64)
    out[1:] = f
    return out


def g_word(n):
    """Big-integer g(n) by running the physical O/Z transfer along the word of n (n >= 2)."""
    D, E, U = 1, 0, 2
    for ch in bin(n - 1)[3:]:
        if ch == '1':
            U, E = 2 * U, 0
        else:
            D, U, E = D - E + U, 2 * U - 1, 1
    need(U == n, 'word reconstruction')
    return D


def g_closed(n):
    """Big-integer g(n) by the closed bit formula (independent of the transfer)."""
    x = n - 1
    if x == 0:
        return 1
    bits = bin(x)[2:][::-1]
    top = len(bits) - 1
    total, q = 1, x >> 1
    for i in range(top):
        if bits[i] == '0':
            total += q + (q & 1)
        q >>= 1
    return total


# ---------------------------------------------------------------- slopes, digits, automaton

def v2(x):
    s = 0
    while x % 2 == 0:
        x //= 2
        s += 1
    return s


def digits(a, b, count):
    out, num = [], a
    for _ in range(count):
        num *= 2
        if num >= b:
            out.append(1)
            num -= b
        else:
            out.append(0)
    return out


HALF, ONE = Fraction(1, 2), Fraction(1)


def gen_pieces(al, be, c, tau, star, E, out, extra_j=2, window=None):
    """Append a description of c * S^(*)(al, be) + tau (al, be > 0).  Pieces:
    ('FULL', c, tau, al, be, star, E)  odd denominator, interior slope: a full parity class scaled by c, shifted by tau;
    ('ONE',  c, tau, al, be, star, E)  one-signed slope (<= 1/2 or >= 1);
    ('PT',   c, tau, al, be, True, E)  the single n = 1 value c (al - be) + tau of an interior S* piece.
    E: every original n in the piece satisfies n < 2^E * (argument).  Branch lemma (PROOF.md, Lemma 3):
    S(al,be) = (S(2al-be, be) + be - al)  u  U_{j>=1} (S*(2^j (2al-be), be) - al).
    The infinite tail of j is cut after the pieces are positive and j >= v2(be) + extra_j (they then share one
    residue class and one sign, PROOF.md section 4); with window=X the cut is instead where the piece minimum exceeds X.
    """
    gg = gcd(al, be)
    al, be, c = al // gg, be // gg, c * gg
    x = Fraction(al, be)
    if HALF < x < ONE:
        if be % 2 == 1:
            out.append(('FULL', c, tau, al, be, star, E))
            if star:
                out.append(('PT', c, tau, al, be, True, E))
            return
        if star:
            out.append(('PT', c, tau, al, be, True, E))
        al1 = 2 * al - be
        gen_pieces(al1, be, c, tau + c * (be - al), False, E + 1, out, extra_j, window)
        s_loc = v2(be)
        j = 1
        while True:
            slope = Fraction(2 ** j * al1, be)
            if window is None:
                if slope > ONE and j >= s_loc + extra_j:
                    break
            else:
                if slope > ONE and j >= s_loc + extra_j and c * (2 ** j * al1 - be) + tau - c * al > window:
                    break
            gen_pieces(2 ** j * al1, be, c, tau - c * al, True, E + j + 1, out, extra_j, window)
            j += 1
        return
    out.append(('ONE', c, tau, al, be, star, E))


def piece_residue(p):
    kind, c, tau, al, be, star, E = p
    return (tau + c * ((al - be) % 2)) % (2 * c), 2 * c


def piece_bound(p):
    """('POS', lower) or ('NEG', upper) for a ONE piece; exact extremum of the piece (attained at n = 1 or n = 3)."""
    kind, c, tau, al, be, star, E = p
    x = Fraction(al, be)
    ext = (al - be) if star else 3 * (al - be)
    if x >= ONE:
        return 'POS', tau + c * ext
    need(x <= HALF, 'one-signed piece slope')
    return 'NEG', tau + c * ext


def classify(a, b):
    """The automaton of PROOF.md section 4: classes mod 2^(s+1) of the parity of a - b, each FULL or signed."""
    s = v2(b)
    M = 2 ** (s + 1)
    out = []
    gen_pieces(a, b, 1, 0, False, 0, out)
    classes = {}
    for cl in range(M):
        if (cl - (a - b)) % 2 == 0:
            classes[cl] = {'full': False, 'pos': [], 'neg': [], 'pts': []}
    for p in out:
        if p[0] == 'PT':
            val = p[2] + p[1] * (p[3] - p[4])
            classes[val % M]['pts'].append(val)
            continue
        res, mod = piece_residue(p)
        if p[0] == 'FULL':
            need(mod == M, 'full piece modulus')
            classes[res]['full'] = True
            continue
        sign, bd = piece_bound(p)
        targets = [res % M] if mod > M else [cl for cl in classes if (cl - res) % mod == 0]
        for cl in targets:
            classes[cl]['pos' if sign == 'POS' else 'neg'].append(bd)
    result = {}
    for cl, info in sorted(classes.items()):
        if info['full']:
            result[cl] = ('FULL', None)
        elif info['pos'] and not info['neg']:
            result[cl] = ('POS', min(info['pos'] + info['pts']))
        elif info['neg'] and not info['pos']:
            result[cl] = ('NEG', max(info['neg'] + info['pts']))
        elif not info['pos'] and not info['neg']:
            result[cl] = ('FINITE', (min(info['pts']), max(info['pts'])) if info['pts'] else None)
        else:
            result[cl] = ('MIXED', None)
    return s, M, result


def apply_letters(m, word):
    for ch in word:
        m = 2 * m if ch == 'O' else 2 * m - 1
    return m


def tail_data(a, b):
    """T2: L, the chronological tail word d_(L-1)..d_1, the head slope frac(2^(L-1) a/b) (odd denominator)."""
    s = v2(b)
    ds = digits(a, b, s + 4 * b + 8)
    L = next(i + 1 for i in range(s, len(ds)) if ds[i] == 1)
    tail = ''.join('Z' if ds[i] else 'O' for i in range(L - 2, -1, -1))
    xh = Fraction((a * 2 ** (L - 1)) % b, b)
    return L, tail, xh, ds


def t3_data(a, b, ds):
    s = v2(b)
    plus = s >= 1 and all(d == 0 for d in ds[1:s + 1])
    minus = s >= 1 and all(d == 1 for d in ds[0:s + 1])
    F = lambda n: a * n - b * g_word(n)
    cplus = None
    if plus:
        cplus = 5 * (a - b) if s == 1 else max(5 * (a - b), F(2 ** s - 1))
    cminus = (2 ** s + 1) * (a - b) if minus else None
    return plus, minus, cplus, cminus


def clock(a, b):
    """Session 198 clock for an interior slope with odd denominator: word W and drop K on strong inputs."""
    k = 1
    while (2 ** k - 1) % b:
        k += 1
    A = a * (2 ** k - 1) // b
    W = ''.join('Z' if (A >> i) & 1 else 'O' for i in range(k))
    drops = set()
    for m in (3, 5, 7, 9, 11, 13, 15, 17, 19, 21):
        mm = apply_letters(m, W)
        drops.add((a * m - b * g_word(m)) - (a * mm - b * g_word(mm)))
    need(len(drops) == 1, 'clock drop independent of the strong input')
    K = drops.pop()
    need(K > 0 and K % 2 == 0, 'clock drop positive and even')
    return W, K


# ---------------------------------------------------------------- literal tournaments (subset DP, independent)

def one_apex_graph(n):
    """Literal one-apex source for cover count n >= 2: root 0 -> chain vertex 1, then O = universal sink,
    Z = vertex beating only the root.  Returns out-neighbour bitmasks."""
    T = [2, 0]
    for ch in bin(n - 1)[3:]:
        v = len(T)
        T.append(0)
        for i in range(v):
            if ch == '0' and i == 0:
                T[v] |= 1
            else:
                T[i] |= 1 << v
    return T


_LAYERS = {}


def _layers(n):
    if n not in _LAYERS:
        N = 1 << n
        pop = np.zeros(N, dtype=np.int64)
        for i in range(n):
            pop += (np.arange(N, dtype=np.int64) >> i) & 1
        _LAYERS[n] = [np.nonzero(pop == k)[0] for k in range(n + 1)]
    return _LAYERS[n]


def path_tables(T, start=None):
    """end[M, v] = number of directed paths covering exactly the vertex set M and ending at v (subset DP);
    with start given, only paths beginning at that vertex."""
    n = len(T)
    N = 1 << n
    adj = np.array([[(T[u] >> v) & 1 for v in range(n)] for u in range(n)], dtype=np.int64)
    dp = np.zeros((N, n), dtype=np.int64)
    if start is None:
        for v in range(n):
            dp[1 << v, v] = 1
    else:
        dp[1 << start, start] = 1
    layers = _layers(n)
    for k in range(2, n + 1):
        masks = layers[k]
        if start is not None:
            masks = masks[(masks >> start) & 1 == 1]
        for w in range(n):
            if start is not None and w == start:
                continue
            sel = masks[(masks >> w) & 1 == 1]
            if len(sel):
                dp[sel, w] = dp[sel ^ (1 << w), :] @ adj[:, w]
    return dp


def literal_profile(T, root):
    """(D, E, U, d, s) of a rooted tournament by literal subset counting (independent of the g recurrence)."""
    n = len(T)
    N = 1 << n
    full = N - 1
    end = path_tables(T)
    conv = [full ^ (1 << v) ^ T[v] for v in range(n)]
    endc = path_tables(conv)                      # ending at v in the converse = starting at v in T
    H = end.sum(axis=1)
    H[0] = 1
    D = int(H[full])
    E = int(end[full, root])
    s = int(endc[full, root])
    Ms = np.arange(N, dtype=np.int64)
    Ms = Ms[(Ms >> root) & 1 == 1]
    U = int((H[full ^ Ms] * endc[Ms, root]).sum())
    d = int(H[full ^ (1 << root)])
    return D, E, U, d, s


def ham_counts_np(T):
    """Hamiltonian paths and Hamiltonian cycles of a tournament by layered numpy subset DP (independent counter)."""
    n = len(T)
    N = 1 << n
    dp = path_tables(T)
    H = int(dp[N - 1].sum())
    del dp
    dp0 = path_tables(T, start=0)
    C = int(sum(int(dp0[N - 1, u]) for u in range(1, n) if (T[u] >> 0) & 1))
    return H, C


# ---------------------------------------------------------------- main

def slopes_upto(bmax):
    out = []
    for b in range(3, bmax + 1):
        for a in range(1, b):
            if gcd(a, b) == 1 and 2 * a > b:
                out.append((a, b))
    return out


def main():
    outpath = HERE / 'VALIDATION.json'
    if '--out' in sys.argv:
        outpath = HERE / sys.argv[sys.argv.index('--out') + 1]
    receipt = {}
    global HEADMAX
    smoke = '--smoke' in sys.argv
    N22 = 1 << (16 if smoke else 22)
    N24 = 1 << (18 if smoke else 24)
    BMAX1, BMAX2, Rw, MQ = (12, 12, 200, 1 << 8) if smoke else (64, 48, 2000, 1 << 12)
    HEADMAX = 1 << (14 if smoke else 20)
    receipt['parameters'] = {'n_census': N22, 'n_direct': N24, 'bmax_P1': BMAX1, 'bmax_P2': BMAX2,
                             'target_window_P2': Rw, 'head_search': HEADMAX, 'q_cover_bound_P4': MQ}

    # ---------------- g arrays, cross-checked
    g24 = g_array(N24)
    gc = g_closed_array(N24)
    need(np.array_equal(g24[1:], gc[1:]), 'recurrence and closed bit formula agree for n < 2^24')
    for n in list(range(2, 2000)) + [2 ** 30 + 12345, 3 ** 40 + 7, 2 ** 200 - 1, 2 ** 200 + 1]:
        need(g_word(n) == g_closed(n), 'big-integer g evaluations agree')
        if n < N24:
            need(g_word(n) == int(g24[n]), 'big-integer g agrees with array')
    receipt['g_crosscheck'] = 'recurrence == closed bit formula for all n < 2^24; transfer == closed formula on 2002 big/small n'
    nodd22 = np.arange(3, N22, 2, dtype=np.int64)
    godd22 = g24[3:N22:2]
    nodd24 = np.arange(3, N24, 2, dtype=np.int64)
    godd24 = g24[3:N24:2]

    # ---------------- Lemma checks: g odd, g <= n, g >= (n+3)/2 on odd n >= 3; branch identity
    need(bool(np.all(g24[1:] % 2 == 1)), 'g odd')
    need(bool(np.all(g24[1:] <= np.arange(1, N24 + 1))), 'g <= n')
    need(bool(np.all(2 * godd24 >= nodd24 + 3)), 'g >= (n+3)/2 on odd n >= 3')
    receipt['lemma1'] = 'g odd, g <= n, 2 g(n) >= n + 3 for odd n >= 3: all n < 2^24'

    # ---------------- P1 census and P3
    P1 = {'slopes': 0, 'full_count_ok': 0, 'violations': 0, 'classes_nonfull': 0,
          'attained_classes': 0, 'not_attained_classes': [], 'group_attained_ok': 0, 'group_checks': 0,
          'mixed_or_unlabelled': 0, 'cstar_digit_agree': 0, 'suffix_class_ok': 0}
    P3 = {'slopes': 0, 'plus_cases': 0, 'minus_cases': 0, 'neither_cases': 0, 'ok': 0, 'fail': []}
    t2_table = []
    for (a, b) in slopes_upto(BMAX1):
        P1['slopes'] += 1
        s, M, cls = classify(a, b)
        bp = b >> s
        nfull = sum(1 for v in cls.values() if v[0] == 'FULL')
        if nfull == (1 if bp > 1 else 0):
            P1['full_count_ok'] += 1
        F = a * nodd22 - b * godd22
        res = F % M
        mins = {cl: None for cl in cls}
        maxs = {cl: None for cl in cls}
        order = np.argsort(res, kind='stable')
        rs = res[order]
        Fs = F[order]
        bounds = np.searchsorted(rs, np.arange(M + 1))
        for cl in cls:
            seg = Fs[bounds[cl]:bounds[cl + 1]]
            if len(seg):
                mins[cl], maxs[cl] = int(seg.min()), int(seg.max())
        for cl, (kind, bd) in cls.items():
            if kind == 'FULL':
                continue
            P1['classes_nonfull'] += 1
            if kind == 'POS':
                if mins[cl] is not None and mins[cl] < bd:
                    P1['violations'] += 1
                if mins[cl] == bd:
                    P1['attained_classes'] += 1
                else:
                    P1['not_attained_classes'].append([a, b, cl, 'POS', bd, mins[cl]])
            elif kind == 'NEG':
                if maxs[cl] is not None and maxs[cl] > bd:
                    P1['violations'] += 1
                if maxs[cl] == bd:
                    P1['attained_classes'] += 1
                else:
                    P1['not_attained_classes'].append([a, b, cl, 'NEG', bd, maxs[cl]])
            else:
                P1['mixed_or_unlabelled'] += 1
        # supplementary (not pre-registered): attainment per group of classes sharing the same bound label
        groups = {}
        for cl, (kind, bd) in cls.items():
            if kind in ('POS', 'NEG'):
                groups.setdefault((kind, bd), []).append(cl)
        for (kind, bd), members in groups.items():
            P1['group_checks'] += 1
            vals = [mins[c] if kind == 'POS' else maxs[c] for c in members if mins[c] is not None]
            ext = (min(vals) if kind == 'POS' else max(vals)) if vals else None
            if ext == bd:
                P1['group_attained_ok'] += 1
        # T2: explicit C*
        if bp > 1:
            L, tail, xh, ds = tail_data(a, b)
            cstar = [cl for cl, v in cls.items() if v[0] == 'FULL'][0]
            # every odd n >= 3 whose word ends with Z + tail: n - 1 = (m - 1) 2^(L-1) + tailbits with m odd >= 3
            low = sum((1 - ds[i]) << i for i in range(L))            # bits of n - 1 below position L
            mask = (nodd22 - 1) % (1 << L) == low
            mask &= (nodd22 - 1) >= (1 << L)
            vals = np.unique(F[mask] % M)
            if len(vals) == 1 and int(vals[0]) == cstar:
                P1['suffix_class_ok'] += 1
            mstar = apply_letters(3, tail)
            if (a * mstar - b * g_word(mstar)) % M == cstar:
                P1['cstar_digit_agree'] += 1
            if b <= 24:
                t2_table.append([a, b, s, L, 'Z' + tail, cstar, M])
            # P3
            if s >= 1:
                P3['slopes'] += 1
                plus, minus, cplus, cminus = t3_data(a, b, ds)
                outside = F[res != cstar]
                ok = True
                if plus:
                    P3['plus_cases'] += 1
                    ok = ok and int(outside.max()) == cplus
                    ok = ok and all(v[0] in ('NEG', 'FULL') for v in cls.values())
                if minus:
                    P3['minus_cases'] += 1
                    ok = ok and int(outside.min()) == cminus
                    ok = ok and all(v[0] in ('POS', 'FULL') for v in cls.values())
                if not plus and not minus:
                    P3['neither_cases'] += 1
                    ok = ok and int(outside.max()) > 200 and int(outside.min()) < -200
                if ok:
                    P3['ok'] += 1
                else:
                    P3['fail'].append([a, b])
        else:
            # dyadic: no full class; T2 not applicable
            if nfull == 0:
                P1['suffix_class_ok'] += 1
                P1['cstar_digit_agree'] += 1
    P1['not_attained_count'] = len(P1['not_attained_classes'])
    P1['not_attained_sample'] = P1['not_attained_classes'][:12]
    del P1['not_attained_classes']
    P1['prediction_structure'] = (P1['full_count_ok'] == P1['slopes'] and P1['violations'] == 0
                                  and P1['mixed_or_unlabelled'] == 0)
    P1['prediction_attained_every_class'] = (P1['not_attained_count'] == 0)
    receipt['P1'] = P1
    receipt['P3'] = P3
    receipt['P3_prediction'] = (len(P3['fail']) == 0)
    receipt['T2_table_b_le_24'] = t2_table

    # ---------------- P2 fullness in |r| <= 2000 for b <= 48, b' > 1
    P2 = {'slopes': 0, 'targets': 0, 'direct': 0, 'recipe': 0, 'failures': [], 'max_recipe_bits': 0,
          'max_clock_reps': 0, 'literal_n_checked': 0, 'tail_identity_checks': 0}
    literal_ns = set()
    for (a, b) in slopes_upto(BMAX2):
        s = v2(b)
        bp = b >> s
        if bp == 1:
            continue
        P2['slopes'] += 1
        _, M, cls = classify(a, b)
        cstar = [cl for cl, v in cls.items() if v[0] == 'FULL'][0]
        F = a * nodd24 - b * godd24
        sel = np.abs(F) <= Rw
        Fsel, nsel = F[sel], nodd24[sel]
        order = np.lexsort((nsel, Fsel))
        Fsel, nsel = Fsel[order], nsel[order]
        first = np.ones(len(Fsel), dtype=bool)
        first[1:] = Fsel[1:] != Fsel[:-1]
        direct = {int(f): int(n) for f, n in zip(Fsel[first], nsel[first])}
        L, tail, xh, ds = tail_data(a, b)
        a2, b2 = xh.numerator, xh.denominator
        need(b2 == bp and HALF < xh < ONE, 'head slope odd denominator, interior')
        # tail lemma: F(tail(m)) = 2^s F''(m) + kappa for strong heads m
        kap = set()
        for m in range(3, 400, 2):
            nn = apply_letters(m, tail)
            kap.add((a * nn - b * g_word(nn)) - (2 ** s) * (a2 * m - b2 * g_word(m)))
        need(len(kap) == 1, 'tail lemma constant')
        kappa = kap.pop()
        P2['tail_identity_checks'] += 1
        W2, K2 = clock(a2, b2)
        mh = np.arange(3, HEADMAX, 2, dtype=np.int64)
        Fh = a2 * mh - b2 * g24[3:HEADMAX:2]
        # per residue class mod K2: values sorted ascending, ties by smaller head
        resh = Fh % K2
        ordh = np.lexsort((mh, Fh, resh))
        Fh_s, mh_s, rs_s = Fh[ordh], mh[ordh], resh[ordh]
        cut = np.searchsorted(rs_s, np.arange(K2 + 1))
        byclass = {cl: (Fh_s[cut[cl]:cut[cl + 1]], mh_s[cut[cl]:cut[cl + 1]]) for cl in range(K2)}
        for r in range(-Rw, Rw + 1):
            if r % M != cstar:
                continue
            P2['targets'] += 1
            if r in direct:
                n = direct[r]
                P2['direct'] += 1
            else:
                t2 = r - kappa
                need(t2 % (2 ** s) == 0, 'head target integral')
                t2 //= 2 ** s
                vals, heads = byclass[t2 % K2]
                pos = int(np.searchsorted(vals, t2))
                need(pos < len(vals), 'a head value above the target in its class (modular seed)')
                v0, m0 = int(vals[pos]), int(heads[pos])
                h = (v0 - t2) // K2
                need((v0 - t2) % K2 == 0, 'clock congruence')
                m = apply_letters(m0, W2 * h)
                n = apply_letters(m, tail)
                P2['recipe'] += 1
                P2['max_clock_reps'] = max(P2['max_clock_reps'], h)
                P2['max_recipe_bits'] = max(P2['max_recipe_bits'], n.bit_length())
            ok = n % 2 == 1 and n >= 3 and a * n - b * g_word(n) == r and a * n - b * g_closed(n) == r
            if not ok:
                P2['failures'].append([a, b, r])
            if n < 2 ** 12:
                literal_ns.add(n)
    for n in sorted(literal_ns):
        T = one_apex_graph(n)
        D, E, U, d, st = literal_profile(T, 0)
        need((D, E, U, d, st) == (int(g24[n]), 1, n, 1, 1), 'literal one-apex profile (g(n), 1, n, 1), start one')
        P2['literal_n_checked'] += 1
    P2['prediction'] = (len(P2['failures']) == 0)
    receipt['P2'] = P2

    # ---------------- P4: C2 offsets for every left source
    import bilinear
    bank = bilinear.source_bank(300)
    srcs = sorted(v for v in bank if v != (1, 1, 1, 1))
    qs = [(int(g24[m]), 1, m, 1) for m in range(3, 83, 2)]
    P4 = {'bank_sources': len(srcs), 'pairs_checked': 0, 'formula_mismatch': 0, 'slope_claims_ok': 0,
          'literal_pairs': []}

    def offset_formula(v, q):
        Dv, Ev, Uv, dv = v
        al = Dv + 4 * Uv - 2 * dv
        be = al + Ev - dv
        c = -Dv + Uv + dv + 2 * Ev
        return al * q[0] - be * q[2] + c, al, be, c

    for v in srcs:
        Dv, Ev, Uv, dv = v
        for q in qs:
            D, E, U, d = q
            Oq = (D, 0, 2 * U, d)
            Jq = (D - E + U, D - E, D + U, D)
            m_ = bilinear.value('M', v, Jq)
            n_ = bilinear.value('N', v, Oq)
            P4['pairs_checked'] += 1
            if m_ - 2 * n_ != offset_formula(v, q)[0]:
                P4['formula_mismatch'] += 1
        _, al, be, c = offset_formula(v, qs[0])
        ok = al % 2 == 1
        if Ev < dv:
            x = Fraction(be, al)
            ok = ok and HALF < x < ONE and x.denominator % 2 == 1
        else:
            ok = ok and be == al
        if ok:
            P4['slope_claims_ok'] += 1

    def word_graph(word):
        T = [2, 0]
        root = 0
        for ch in word:
            n = len(T)
            T.append(0)
            for i in range(n):
                if ch in ('Z', 'J') and i == root:
                    T[n] |= 1 << i
                else:
                    T[i] |= 1 << n
            if ch == 'J':
                root = n
        return T, root

    # Left sources: all O/Z/J words of length <= 5 ending in O (so E = 0 < d), profile counted literally, with
    # gcd(alpha, d - E) = 1.  For each, the least strong unit cover m < MQ with offset -1, resp. +1.  Keep the first
    # occurrence of each (profile, offset) whose child graph has order <= 19 (subset-DP feasible), in the order
    # (child order, word, offset).
    import itertools
    cands = []
    for Lw in range(0, 6):
        for wt in itertools.product('OZJ', repeat=Lw):
            vw = ''.join(wt) + 'O'
            V, c0 = word_graph(vw)
            v = literal_profile(V, c0)[:4]
            if not v[1] < v[3]:
                continue
            _, al, be, c = offset_formula(v, qs[0])
            if gcd(al, v[3] - v[1]) != 1:
                continue
            got = set()
            for m in range(3, MQ, 2):
                dlt = al * int(g24[m]) - be * m + c
                if dlt in (-1, 1) and dlt not in got:
                    got.add(dlt)
                    cands.append((len(V) + (m - 1).bit_length() + 4, vw, v, m, dlt))
                if len(got) == 2:
                    break
    cands.sort()
    chosen, used = [], set()
    for (order_child, vw, v, m, dlt) in cands:
        if order_child <= 19 and (v, dlt) not in used:
            used.add((v, dlt))
            chosen.append((vw, v, m, dlt))
    P4['candidate_pairs'] = len(cands)
    for (vw, v, m, dlt) in chosen:
        V, c0 = word_graph(vw)
        if True:
            qword = ''.join('O' if ch == '1' else 'Z' for ch in bin(m - 1)[3:])
            Vq_O, cO = word_graph(qword + 'O')
            Vq_J, cJ = word_graph(qword + 'J')
            parent = bilinear.construct('N', V, c0, Vq_O, cO)[0]
            child = bilinear.construct('M', V, c0, Vq_J, cJ)[0]
            Hp, Cp = ham_counts_np(parent)
            Hc, Cc = ham_counts_np(child)
            q = (int(g24[m]), 1, m, 1)
            ok = (Cp == 2 and Cc == 2 and Hc == 2 * Hp + dlt
                  and Hp == bilinear.value('N', v, (q[0], 0, 2 * q[2], 1))
                  and Hc == bilinear.value('M', v, (q[0] - 1 + q[2], q[0] - 1, q[0] + q[2], q[0])))
            P4['literal_pairs'].append({'left_word': vw, 'left_profile': list(v), 'q_cover': m, 'offset': dlt,
                                        'parent_order': len(parent), 'child_order': len(child),
                                        'parent_H': Hp, 'child_H': Hc, 'cycles': [Cp, Cc], 'ok': ok})
    lp = P4['literal_pairs']
    P4['literal_distinct_left_profiles'] = len({tuple(p['left_profile']) for p in lp})
    P4['literal_offsets'] = sorted({p['offset'] for p in lp})
    P4['prediction'] = (P4['formula_mismatch'] == 0 and P4['slope_claims_ok'] == len(srcs)
                        and len(lp) >= 6 and all(p['ok'] for p in lp)
                        and P4['literal_distinct_left_profiles'] >= 3 and P4['literal_offsets'] == [-1, 1])
    receipt['P4'] = P4

    # ---------------- P5: densities of non-full classes on their unbounded side
    P5 = {'decided': [], 'lower_bounds_only': [], 'exceptions': []}
    X = 20000
    for (a, b) in [(3, 4), (5, 6), (7, 10), (9, 10), (7, 12), (11, 12), (9, 14), (13, 14), (5, 8), (7, 8),
                   (13, 24), (17, 24), (19, 24), (23, 24)]:
        s, M, cls = classify(a, b)
        pieces = []
        gen_pieces(a, b, 1, 0, False, 0, pieces, extra_j=2, window=X)
        F = a * nodd22 - b * godd22
        res = F % M
        for cl, (kind, bd) in cls.items():
            if kind not in ('POS', 'NEG'):
                continue
            decided = True
            for p in pieces:
                if p[0] != 'ONE':
                    continue
                pr, pm = piece_residue(p)
                if not ((pm <= M and (cl - pr) % pm == 0) or (pm > M and pr % M == cl)):
                    continue
                pk, pc, ptau, pal, pbe, pstar, pE = p
                x = Fraction(pal, pbe)
                if x == ONE or x == HALF:
                    decided = False
                    break
                if x > ONE:
                    mmax = Fraction(X - ptau, pc * (pal - pbe))
                else:
                    mmax = Fraction(X + ptau, pc * (Fraction(pbe, 2) - pal))
                if mmax > 0 and 2 ** pE * (mmax + 1) > N22:
                    decided = False
                    break
            sel = F[res == cl]
            if kind == 'POS':
                window = sel[(sel >= 500) & (sel <= X)]
                tot = len([r for r in range(500, X + 1) if r % M == cl])
            else:
                window = sel[(sel <= -500) & (sel >= -X)]
                tot = len([r for r in range(-X, -499) if r % M == cl])
            hit = len(np.unique(window))
            if hit < 50:
                continue
            frac = round(hit / tot, 6)
            row = [a, b, cl, kind, hit, tot, frac]
            if decided:
                P5['decided'].append(row)
                if not (0.02 < frac < 0.98):
                    P5['exceptions'].append(row)
            else:
                P5['lower_bounds_only'].append(row)
    P5['prediction'] = (len(P5['exceptions']) == 0)
    receipt['P5'] = P5

    # ---------------- T5: clock uniqueness on all words of length <= 12, b <= 31
    T5 = {'slopes': 0, 'words': 0, 'ok': 0}
    for (a, b) in slopes_upto(31 if not smoke else 9):
        T5['slopes'] += 1
        if b % 2:
            Wc, _ = clock(a, b)
        else:
            Wc = None
        for k in range(1, 13 if not smoke else 9):
            for A in range(2 ** k):
                w = ''.join('Z' if (A >> i) & 1 else 'O' for i in range(k))
                lam = a * (2 ** k - 1) - b * A
                is_power = Wc is not None and k % len(Wc) == 0 and w == Wc * (k // len(Wc))
                ch = {(a * apply_letters(m, w) - b * int(g24[apply_letters(m, w)])) - (a * m - b * int(g24[m]))
                      for m in (3, 5, 7, 9)}
                T5['words'] += 1
                if (lam == 0) == is_power and (len(ch) == 1) == (lam == 0):
                    T5['ok'] += 1
    T5['all_ok'] = (T5['ok'] == T5['words'])
    need(T5['all_ok'], 'clock uniqueness lemma')
    receipt['T5'] = T5

    # ---------------- P6: witness lengths (outputs of witness_length.c)
    P6 = {}
    for fname, bound_j in (('witness_length_3_5.txt', (4, 9)), ('witness_length_8_9.txt', (4, 9))):
        text = (HERE / fname).read_text().strip().splitlines()
        rows = {}
        for line in text:
            if line.startswith('window'):
                parts = line.split()
                w, tot, miss, ml = int(parts[1]), int(parts[3]), int(parts[5]), int(parts[7])
                rows[w] = (tot, miss, ml)
            if line.startswith('trailing-clock'):
                hist = {int(t.split('=')[0][1:]): int(t.split('=')[1]) for t in line.split(':')[1].split()}
        tot, miss, _ = rows[1024]
        a_ok = (tot - miss) / tot >= 0.99
        b_rows = {j: rows[2 ** j] for j in range(bound_j[0], bound_j[1] + 1)}
        b_ok = all(miss == 0 and ml <= 3 * j + 8 for j, (tot, miss, ml) in b_rows.items())
        b_fail = [j for j, (tot, miss, ml) in b_rows.items() if not (miss == 0 and ml <= 3 * j + 8)]
        t1 = sum(v for t, v in hist.items() if t >= 1)
        t2 = sum(v for t, v in hist.items() if t >= 2)
        c_ok = 2 * t2 <= t1
        P6[fname] = {'found_fraction_1024': round((tot - miss) / tot, 6), 'a': a_ok, 'b': b_ok,
                     'b_failing_j': b_fail, 'maxlen_by_j': {j: v[2] for j, v in sorted(b_rows.items())},
                     't_ge_1': t1, 't_ge_2': t2, 'c': c_ok,
                     'sha256': hashlib.sha256((HERE / fname).read_bytes()).hexdigest()}
    P6['prediction'] = all(v['a'] and v['b'] and v['c'] for v in P6.values() if isinstance(v, dict))
    receipt['P6'] = P6

    receipt['dependency_hashes'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in sorted((HERE / 'dependencies').glob('*.py'))}
    receipt['scope'] = ('THM-4883 confirmatory checks; the theorem is proved in PROOF.md, these finite checks '
                        'corroborate it and decide the pre-registered predictions P1-P6.')
    text = json.dumps(receipt, sort_keys=True, indent=1) + '\n'
    outpath.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
