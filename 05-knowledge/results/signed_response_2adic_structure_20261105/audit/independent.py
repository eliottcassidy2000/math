"""independent.py -- independent adversarial audit of THM-4883 / HYP-9476 (own code only).

Nothing is imported from the packet (no verify.py, no dependencies/bilinear.py or core.py).  The C2 graphs of T4 are
built here directly from the text of BILINEAR_C2.md (Section 2 graft, Section 3 pair forms, Section 4 converse), the
one-apex sources directly from the O/Z/J birth rules, and every Hamilton count by the subset DP in audit_census.c.

    python independent.py --exe <path/to/audit_census binary> --out INDEPENDENT.json

If the binary does not exist it is compiled from audit_census.c with gcc -O2 (default location: the system temp
directory, never the repository).  Output is deterministic (fixed seeds, sorted keys); python -O gives the same bytes.
About 5 minutes (the depth-35 witness search for 8/9 is about 2.5 of them).
"""
import sys, json, subprocess, math, random, itertools, time
from fractions import Fraction
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
HALF, ONE = Fraction(1, 2), Fraction(1)
FAILS = []          # hard failures of a claim (anything here makes the run report a problem)


def check(ok, msg):
    if not ok:
        FAILS.append(msg)
    return ok


# ====================================================================== g: three independent big-integer evaluations

def g_topdown(n):
    """g(n) by unwinding the recurrence from n down to 1."""
    if n < 1:
        raise ValueError
    total = 1
    while n > 1:
        if n & 1:
            m = (n + 1) >> 1
            total += m - (m & 1)
            n = m
        else:
            n >>= 1
    return total


def g_bits(n):
    """g(n) = f(n - 1), f(A) = 1 + sum_{i: e_i = 1} Z_i + #{i >= 1: e_(i-1) = 0, e_i = 1},
    Z_i = sum_{h < i, e_h = 0} 2^(i-1-h) (DIGIT_BLOCK's bit expansion, summed over set bits, own implementation)."""
    A = n - 1
    total, Z, prev, i = 1, 0, 1, 0
    while A >> i:
        e = (A >> i) & 1
        if e:
            total += Z
            if i >= 1 and prev == 0:
                total += 1
        Z = 2 * Z + (1 - e)
        prev = e
        i += 1
    return total


def g_transfer(n):
    """g(n) as the D coordinate of the seed (1,0,2,1) after the O/Z word of n (O: bit 1, Z: bit 0 of n-1)."""
    D, E, U = 1, 0, 2
    for ch in bin(n - 1)[3:]:
        if ch == '1':
            U, E = 2 * U, 0
        else:
            D, E, U = D - E + U, 1, 2 * U - 1
    if U != n:
        raise RuntimeError('transfer cover count')
    return D


def g_list(N):
    g = [0] * (N + 1)
    g[1] = 1
    for n in range(2, N + 1):
        if n & 1:
            m = (n + 1) >> 1
            g[n] = g[m] + m - (m & 1)
        else:
            g[n] = g[n >> 1]
    return g


def F3(a, b, n):
    """a n - b g(n) by all three evaluations; returns None if they disagree."""
    v = {a * n - b * g_topdown(n), a * n - b * g_bits(n), a * n - b * g_transfer(n)}
    return v.pop() if len(v) == 1 else None


def word_of(n):
    return ''.join('O' if ch == '1' else 'Z' for ch in bin(n - 1)[3:])


def apply_word(m, w):
    for ch in w:
        m = 2 * m if ch == 'O' else 2 * m - 1
    return m


def v2(x):
    s = 0
    while x % 2 == 0:
        x //= 2
        s += 1
    return s


def digits(a, b, cnt):
    out, r = [], a % b
    for _ in range(cnt):
        r *= 2
        out.append(1 if r >= b else 0)
        if r >= b:
            r -= b
    return out


# ====================================================================== C helper

class Helper:
    def __init__(self, exe):
        self.exe = str(exe)

    def census(self, log2n, rows):
        inp = ''.join(' '.join(str(x) for x in r) + '\n' for r in rows)
        out = subprocess.run([self.exe, 'census', str(log2n)], input=inp, capture_output=True, text=True, check=True).stdout
        res, cur = {}, None
        head = None
        for line in out.splitlines():
            p = line.split()
            if p[0] == 'LEMMA1_BAD':
                head = (int(p[1]), int(p[3]))
            elif p[0] == 'SLOPE':
                cur = {'classes': {}, 'M': int(p[4])}
                res[(int(p[1]), int(p[2]))] = cur
            elif p[0] == 'C':
                c = int(p[1])
                if int(p[2]) == 0:
                    cur['classes'][c] = None
                else:
                    cur['classes'][c] = dict(zip(['count', 'min', 'argmin', 'max', 'argmax', 'hit_pos', 'hit_neg', 'hit_win'],
                                                 [int(x) for x in p[2:10]]))
            elif p[0] == 'SUFFIX':
                cur['suffix'] = (int(p[1]), int(p[2]))
            elif p[0] == 'MISS':
                cur['miss'] = [int(x) for x in p[1:]]
        return head, res

    def graphs(self, items):
        """items: list of (out_masks, root).  Returns list of (D, E, U, d, s, C)."""
        inp = []
        for out, root in items:
            inp.append('%d %d %s\n' % (len(out), root, ' '.join(str(x) for x in out)))
        res = subprocess.run([self.exe, 'graph'], input=''.join(inp), capture_output=True, text=True, check=True).stdout
        return [tuple(int(x) for x in line.split()) for line in res.splitlines()]

    def wl(self, a, b, depth, R):
        res = subprocess.run([self.exe, 'wl', str(a), str(b), str(depth), str(R)], capture_output=True, text=True,
                             check=True).stdout
        best = {}
        for line in res.splitlines():
            p = line.split()
            if p[0] == 'W':
                best[int(p[1])] = (int(p[2]), p[3] if len(p) > 3 else '')
        return best


# ====================================================================== literal tournaments (out-neighbour bitmasks)

def tour_from_word(word):
    """Rooted tournament of the O/Z/J word applied to the two-vertex seed (root 0 -> vertex 1).
    O: newborn is a universal sink (root kept); Z: newborn beats only the current root (root kept);
    J: newborn beats only the current root and becomes the root (DIGIT_BLOCK / BILINEAR_C2 Section 3)."""
    out = [0b10, 0]
    root = 0
    for ch in word:
        v = len(out)
        out.append(0)
        for i in range(v):
            if ch in 'ZJ' and i == root:
                out[v] |= 1 << i
            else:
                out[i] |= 1 << v
        if ch == 'J':
            root = v
    return out, root


def converse(out):
    n = len(out)
    full = (1 << n) - 1
    return [full ^ (1 << v) ^ out[v] for v in range(n)]


def is_tournament(out):
    n = len(out)
    for u in range(n):
        if (out[u] >> u) & 1:
            return False
        for v in range(u + 1, n):
            if ((out[u] >> v) & 1) == ((out[v] >> u) & 1):
                return False
    return True


def strongly_connected(out):
    n = len(out)
    seen, st = 1, [0]
    while st:
        u = st.pop()
        nb = out[u] & ~seen
        while nb:
            w = (nb & -nb).bit_length() - 1
            nb &= nb - 1
            seen |= 1 << w
            st.append(w)
    if seen != (1 << n) - 1:
        return False
    rin = [0] * n
    for u in range(n):
        for v in range(n):
            if (out[u] >> v) & 1:
                rin[v] |= 1 << u
    seen, st = 1, [0]
    while st:
        u = st.pop()
        nb = rin[u] & ~seen
        while nb:
            w = (nb & -nb).bit_length() - 1
            nb &= nb - 1
            seen |= 1 << w
            st.append(w)
    return seen == (1 << n) - 1


def j_birth(out, root):
    """J_root(V): newborn whose only old out-neighbour is root; every other old vertex beats it."""
    v = len(out)
    out = list(out) + [1 << root]
    for i in range(v):
        if i != root:
            out[i] |= 1 << v
    return out, v


def graft(T, a, b, Gout, r):
    """BILINEAR_C2 Section 2: every cross edge T -> G except r -> a and r -> b."""
    nT, nG = len(T), len(Gout)
    out = [0] * (nT + nG)
    for i in range(nT):
        out[i] = T[i]
    for j in range(nG):
        m = 0
        for k in range(nG):
            if (Gout[j] >> k) & 1:
                m |= 1 << (nT + k)
        out[nT + j] = m
    for i in range(nT):
        for j in range(nG):
            if j == r and i in (a, b):
                out[nT + j] |= 1 << i
            else:
                out[i] |= 1 << (nT + j)
    return out


def pair_N(V, c):
    T, a = j_birth(V, c)
    return T, a, c


def pair_M(V, c):
    S, bb = j_birth(V, c)
    T, u = j_birth(S, bb)
    return T, u, c


# ====================================================================== transfers on profiles

def O_(p):
    D, E, U, d = p
    return (D, 0, 2 * U, d)


def Z_(p):
    D, E, U, d = p
    return (D - E + U, d, 2 * U - d, d)


def J_(p):
    D, E, U, d = p
    return (D - E + U, D - E, D + U, D)


def run_word(word, p=(1, 0, 2, 1)):
    for ch in word:
        p = {'O': O_, 'Z': Z_, 'J': J_}[ch](p)
    return p


def N_form(v, w):
    D, E, U, dv = v
    g, s, L, d = w
    return (D - E + U) * g - (D - E + dv) * s + (D + 3 * U - dv) * L - (U + E) * d


def M_form(v, w):
    D, E, U, dv = v
    g, s, L, d = w
    return (D + 2 * U) * g - (U + dv) * s + (2 * D - E + 6 * U - dv) * L - (U + E) * d


def t4_coeffs(v):
    D, E, U, d = v
    al = D + 4 * U - 2 * d
    be = al + E - d
    c = -D + U + d + 2 * E
    return al, be, c


# ====================================================================== slopes, clocks, recursion of the proof

def interior_slopes(bmin, bmax):
    return [(a, b) for b in range(bmin, bmax + 1) for a in range(1, b) if gcd(a, b) == 1 and 2 * a > b]


def clock(al, be):
    """Lemma 6 clock for 1/2 < al/be < 1, be odd: word W (chronological), K."""
    k = 1
    while (pow(2, k, be) - 1) % be:
        k += 1
    A = al * (2 ** k - 1) // be
    W = ''.join('Z' if (A >> i) & 1 else 'O' for i in range(k))
    return W, k, A


def clock_drop_exact(al, be):
    """K = be * Delta, Delta = sum_i e_i (u_i - e_(i-1)) cyclic, u_0 = R/(B - 1), u_(i+1) = 2 u_i - e_i (exact)."""
    W, k, A = clock(al, be)
    eps = [(A >> i) & 1 for i in range(k)]
    R = sum(eps[i] << (k - 1 - i) for i in range(k))
    u = Fraction(R, 2 ** k - 1)
    us = []
    for i in range(k):
        us.append(u)
        u = 2 * u - eps[i]
    periodic = (u == us[0]) and all(0 < x < 1 and eps[i] == math.floor(2 * x) for i, x in enumerate(us))
    Delta = sum(eps[i] * (us[i] - eps[i - 1]) for i in range(k))
    # run decomposition: each maximal cyclic run of j ones followed by a zero with orbit value v there
    runs = Fraction(0)
    start = next(i for i in range(k) if eps[i] == 0)       # rotate to a zero
    i = 0
    while i < k:
        idx = (start + i) % k
        if eps[idx] == 1:
            j = 0
            while eps[(start + i) % k] == 1:
                j += 1
                i += 1
            vz = us[(start + i) % k]
            runs += vz + (1 - vz) / Fraction(2 ** j)
        else:
            i += 1
    return W, k, A, be * Delta, periodic, runs == Delta


def structure(al, be, scale, shift, star, out, extra=3):
    """Atoms of scale * S^(*)(al, be) + shift (be >= 1, al > 0) by Lemma 3 / Lemma 2 (own implementation).
    ('FULL', mod, res) | ('ONE', mod, res, sign, ext) | ('PT', value)."""
    gg = gcd(al, be)
    al, be, scale = al // gg, be // gg, scale * gg
    x = Fraction(al, be)
    mod = 2 * scale
    res = (shift + scale * ((al - be) % 2)) % mod
    if HALF < x < ONE:
        if star:
            out.append(('PT', shift + scale * (al - be)))
        if be % 2 == 1:
            out.append(('FULL', mod, res))
            return
        al1 = 2 * al - be
        structure(al1, be, scale, shift + scale * (be - al), False, out, extra)
        sl = v2(be)
        j = 1
        tail = []
        while True:
            sub = []
            structure((2 ** j) * al1, be, scale, shift - scale * al, True, sub, extra)
            if j >= sl and Fraction((2 ** j) * al1, be) > 1:
                tail.append(sub)
                if len(tail) > extra:
                    break
            else:
                out.extend(sub)
            j += 1
        # the cut family: all pieces j >= first cut index are one-signed POS in one class; keep the first, and check
        first = tail[0]
        out.extend(first)
        f_atoms = [t for t in first if t[0] == 'ONE']
        for later in tail[1:]:
            l_atoms = [t for t in later if t[0] == 'ONE']
            ok = (len(l_atoms) == 1 and len(f_atoms) == 1 and l_atoms[0][3] == 'POS' and f_atoms[0][3] == 'POS'
                  and l_atoms[0][4] >= f_atoms[0][4] and l_atoms[0][1] == f_atoms[0][1] and l_atoms[0][2] == f_atoms[0][2])
            check(ok, 'family cut: later pieces in the same class with larger minima (%s/%s)' % (al, be))
            for t in later:
                if t[0] == 'PT':
                    out.append(t)
        return
    ext = (al - be) if star else 3 * (al - be)
    if x >= ONE:
        out.append(('ONE', mod, res, 'POS', shift + scale * ext, str(x)))
    else:
        out.append(('ONE', mod, res, 'NEG', shift + scale * ext, str(x)))


def classify(a, b):
    """Per class of P_s (odd residues mod 2^(s+1) of the parity of a - b): FULL / POS bound / NEG bound / MIXED."""
    s = v2(b)
    M = 2 ** (s + 1)
    atoms = []
    structure(a, b, 1, 0, False, atoms)
    cls = {}
    for c in range(M):
        if (c - (a - b)) % 2:
            continue
        full, pos, neg, keys = False, [], [], []
        for t in atoms:
            if t[0] == 'FULL':
                check(t[1] == M, 'FULL atom has modulus exactly 2^(s+1) (%d/%d)' % (a, b))
                if t[1] <= M and (c - t[2]) % t[1] == 0:
                    full = True
            elif t[0] == 'ONE':
                hit = (c - t[2]) % t[1] == 0 if t[1] <= M else t[2] % M == c
                if hit:
                    (pos if t[3] == 'POS' else neg).append(t[4])
                    keys.append(t)
        pts = [t[1] for t in atoms if t[0] == 'PT' and t[1] % M == c]
        if full:
            cls[c] = ('FULL', None, None)
        elif pos and not neg:
            cls[c] = ('POS', min(pos + pts), tuple(sorted(map(str, keys))))
        elif neg and not pos:
            cls[c] = ('NEG', max(neg + pts), tuple(sorted(map(str, keys))))
        elif not pos and not neg:
            cls[c] = ('FINITE', tuple(pts), None)
        else:
            cls[c] = ('MIXED', None, None)
    return s, M, cls


def t2_data(a, b):
    s = v2(b)
    ds = digits(a, b, s + 4 * b + 10)
    L = next(i + 1 for i in range(s, len(ds)) if ds[i] == 1)
    tail = ''.join('Z' if ds[i] else 'O' for i in range(L - 2, -1, -1))   # chronological d_(L-1) .. d_1
    low = sum((1 - ds[i]) << i for i in range(L))                         # n - 1 mod 2^L
    num = (a * 2 ** (L - 1)) % b
    a2, b2 = num >> s, b >> s
    return s, ds, L, tail, low, a2, b2


def cstar_from_digits(a, b):
    s, ds, L, tail, low, a2, b2 = t2_data(a, b)
    M = 2 ** (s + 1)
    m = apply_word(3, tail)                     # head 3 = seed . Z, then the tail
    return (a * m - b * g_topdown(m)) % M


def t3_pred(a, b):
    s = v2(b)
    ds = digits(a, b, s + 2)
    plus = all(d == 0 for d in ds[1:s + 1])
    minus = all(d == 1 for d in ds[0:s + 1])
    F = lambda n: a * n - b * g_topdown(n)
    cplus = (5 * (a - b) if s == 1 else max(5 * (a - b), F(2 ** s - 1))) if plus else None
    cminus = (2 ** s + 1) * (a - b) if minus else None
    return plus, minus, cplus, cminus


# ====================================================================== JOINT_MODULAR_GADGET Theorem 1 (own implementation)

def jmg_L(M, u, z, extra_loops=1):
    """Return L with L = u, g(L) = z (mod M), M even, z odd, following Session 197 Theorem 1 literally."""
    k, q = v2(M), M >> v2(M)
    a = 1
    while pow(2, a, q) != 1 % q or a < 2:
        a += 1
    Delta = (2 ** a - 1) * (2 ** k - 1) + 1
    r0 = (u - 1) % M
    width = max(1, (M - 1).bit_length())
    # f(r0) with f(x) = g(x + 1)
    fr = g_topdown(r0 + 1)
    jj = None
    for j in range(M // 2):
        if (fr + j * Delta - z) % M == 0:
            jj = j
            break
    if jj is None:
        return None
    jj += extra_loops * (M // 2)
    word = ('1' * a + '0' * k) * jj + format(r0, '0%db' % width)
    L = int(word, 2) + 1
    return L


def egcd(x, y):
    if y == 0:
        return (x, 1, 0)
    g0, p, q = egcd(y, x % y)
    return (g0, q, p - (x // y) * q)


GARR = None      # numpy g[0 .. 2^20], set in main


def a0_witness(al, be, r, use_jmg=False):
    """Theorem A0's route: residues (u, z) of (L, g(L)) mod K solving the congruence, head n0 = 2^(j+1) L - 1 (Lemma 3
    branch with m = 2^j L), then (v - r)/K clock periods.  L is the least odd L < 2^20 with those residues and
    value >= r (its existence for some L is the JOINT_MODULAR_GADGET theorem), or, with use_jmg, the gadget's own L."""
    W, k, A = clock(al, be)
    K = (al * 3 - be * 3) - (al * apply_word(3, W) - be * g_topdown(apply_word(3, W)))
    j = 1
    while (2 ** j) * (2 * al - be) <= be:
        j += 1
    sigma = (2 ** j) * (2 * al - be)
    check(gcd(sigma, be) == 1, 'gcd(sigma, beta) = 1 for %d/%d' % (al, be))
    need2 = r - (sigma - be - al)                  # value = sigma - be - al + 2 (sigma y - be x) (mod K)
    if need2 % 2:
        return None
    half = K // 2
    c = (need2 // 2) % half if half > 1 else 0
    g0, p, q = egcd(sigma, be)                     # sigma p + be q = 1
    y, x = c * p, -c * q                           # sigma y - be x = c
    u, z = (1 + 2 * y) % K, (1 + 2 * x) % K
    L = None
    if use_jmg:
        for extra in range(1, 200):
            Lc = jmg_L(K, u, z, extra)
            if Lc is not None and sigma * Lc - be * g_topdown(Lc) - al >= r:
                L = Lc
                break
    else:
        Lo = np.arange(1, 1 << 20, 2, dtype=np.int64)
        gL = GARR[1:1 << 20:2]
        val = sigma * Lo - be * gL - al
        ok = np.nonzero((Lo % K == u) & (gL % K == z) & (val >= r))[0]
        if len(ok):
            best = ok[np.argmin(val[ok])]          # least branch value >= r: fewest clock periods
            L = int(Lo[best])
    if L is None:
        return 'no L below 2^20 with the residues'
    check(L % K == u and g_topdown(L) % K == z, 'A0 head residues')
    n0 = 2 ** (j + 1) * L - 1
    v = al * n0 - be * g_topdown(n0)
    check(v == sigma * L - be * g_topdown(L) - al, 'Lemma 3 branch value at n0')
    if v < r or (v - r) % K:
        return 'branch value congruence failed'
    h = (v - r) // K
    if h * len(W) > LETTER_CAP:
        return 'over letter cap'
    return apply_word(n0, W * h), K, L, h


LETTER_CAP = 40000
_HEADS = {}


def direct_witness(al, be, r):
    """Practical recipe (same clock): least head value F(m0) >= r with F(m0) = r (mod K) over all odd m0 < 2^20, then
    (F(m0) - r)/K clock periods.  Returns (n, K, m0, h) or a reason string."""
    if (al, be) not in _HEADS:
        W, k, A = clock(al, be)
        K = (al * 3 - be * 3) - (al * apply_word(3, W) - be * g_topdown(apply_word(3, W)))
        mh = np.arange(3, 1 << 20, 2, dtype=np.int64)
        Fh = al * mh - be * GARR[3:1 << 20:2]
        resid = Fh % K
        order = np.lexsort((mh, Fh, resid))
        Fs, ms, rs = Fh[order], mh[order], resid[order]
        cut = np.searchsorted(rs, np.arange(K + 1))
        if len(_HEADS) > 64:
            _HEADS.clear()
        _HEADS[(al, be)] = (W, K, Fs, ms, cut)
    W, K, Fs, ms, cut = _HEADS[(al, be)]
    cl = r % K
    seg_F, seg_m = Fs[cut[cl]:cut[cl + 1]], ms[cut[cl]:cut[cl + 1]]
    pos = int(np.searchsorted(seg_F, r))
    if pos >= len(seg_F):
        return 'no head'
    v0, m0 = int(seg_F[pos]), int(seg_m[pos])
    h = (v0 - r) // K
    if h * len(W) > LETTER_CAP:
        return 'over letter cap'
    return apply_word(m0, W * h), K, m0, h


# ====================================================================== main

def main():
    args = sys.argv[1:]
    import tempfile                        # default binary location is outside the repository
    exe = (Path(args[args.index('--exe') + 1]) if '--exe' in args
           else Path(tempfile.gettempdir()) / ('audit_census_thm4883' + ('.exe' if sys.platform == 'win32' else '')))
    outp = Path(args[args.index('--out') + 1]) if '--out' in args else HERE / 'INDEPENDENT.json'
    quick = '--quick' in args
    if not exe.exists():
        subprocess.run(['gcc', '-O2', '-o', str(exe), str(HERE / 'audit_census.c')], check=True)
    H = Helper(exe)
    R = {}
    t_start = time.time()
    random.seed(20261105)

    def log(msg):
        print('[%7.1f s] %s' % (time.time() - t_start, msg), file=sys.stderr, flush=True)

    # ------------------------------------------------------------------ 0. g cross-checks
    log('section 0: g cross-checks')
    global GARR
    G20 = g_list(1 << 20)
    GARR = np.array(G20, dtype=np.int64)
    bad = 0
    for n in range(1, 1 << 14):
        vals = {g_topdown(n), g_bits(n), G20[n]}
        if n >= 2:
            vals.add(g_transfer(n))
        bad += len(vals) != 1
    bigs = [2 ** 61 - 1, 2 ** 200 + 1, 3 ** 300 + 2, 7 ** 500 + 3] + [random.getrandbits(1500) | 1 for _ in range(20)]
    for n in bigs:
        bad += len({g_topdown(n), g_bits(n), g_transfer(n)}) != 1
    R['g_three_evaluations_disagreements'] = bad
    check(bad == 0, 'g evaluations disagree')

    log('dictionary')
    # ------------------------------------------------------------------ 1. physical dictionary: literal one-apex profiles
    items, ns = [], []
    for n in list(range(2, 300)) + sorted(random.sample(range(300, 1 << 13), 60)):
        out, root = tour_from_word(word_of(n))
        items.append((out, root))
        ns.append(n)
    prof = H.graphs(items)
    bad = [n for n, p in zip(ns, prof) if p != (G20[n], n % 2, n, 1, 1, n % 2)]
    R['dictionary_literal_profiles'] = {'checked': len(ns), 'max_order': max(len(o) for o, _ in items),
                                        'claim': '(D,E,U,d,s,C) = (g(n), n mod 2, n, 1, 1, n mod 2)', 'bad': bad}
    check(not bad, 'literal one-apex profile mismatch')
    # converse: every rooted tournament (orders 2..6, root 0) with s = d = 1 is one-apex; s = E = d = 1 <=> odd n >= 3
    conv_bad, conv_count, conv_strong = 0, 0, 0
    for nv in range(2, 7):
        pairs = [(i, j) for i in range(nv) for j in range(i + 1, nv)]
        batch = []
        for bits in range(1 << len(pairs)):
            out = [0] * nv
            for t, (i, j) in enumerate(pairs):
                if (bits >> t) & 1:
                    out[i] |= 1 << j
                else:
                    out[j] |= 1 << i
            batch.append((out, 0))
        res = H.graphs(batch)
        for (out, _), (D, E, U, d, s1, C) in zip(batch, res):
            if s1 == 1 and d == 1:
                conv_count += 1
                ok = (D == G20[U] and E == U % 2 and U >= 2)
                if E == 1:
                    conv_strong += 1
                    ok = ok and U % 2 == 1 and U >= 3 and strongly_connected(out)
                conv_bad += not ok
    R['dictionary_converse_all_rooted_tournaments_orders_2_6'] = {'s_eq_d_eq_1': conv_count, 'with_E_1': conv_strong,
                                                                   'bad': conv_bad}
    check(conv_bad == 0, 'converse dictionary')

    log('lemma 3')
    # ------------------------------------------------------------------ 2. Lemma 3 (branch decomposition) numerically
    l3_bad = 0
    for (al, be) in [(3, 4), (5, 6), (-7, 3), (2, 9), (11, 4), (0, 5), (13, 24), (1, 1)]:
        lhs = sorted(al * n - be * G20[n] for n in range(3, 1 << 14, 2))
        rhs = []
        for m in range(3, 1 << 13, 2):          # n = 2m - 1, m odd >= 3
            rhs.append((2 * al - be) * m - be * G20[m] + be - al)
        for jx in range(1, 14):
            for mp in range(1, 1 << 13, 2):
                n = 2 ** (jx + 1) * mp - 1
                if n < (1 << 14):
                    rhs.append((2 ** jx) * (2 * al - be) * mp - be * G20[mp] - al)
        l3_bad += lhs != sorted(rhs)
    R['lemma3_multiset_identity_bad'] = l3_bad
    check(l3_bad == 0, 'Lemma 3')

    log('lemma 6')
    # ------------------------------------------------------------------ 3. Lemma 6 clocks (odd denominators <= 101)
    lem6 = {'slopes': 0, 'bad': [], 'K_min': None, 'K_max': None}
    for (al, be) in interior_slopes(3, 101):
        if be % 2 == 0:
            continue
        W, k, A, Kex, periodic, runs_ok = clock_drop_exact(al, be)
        drops = set()
        for m in [3, 5, 7, 9, 11, 13, 15, 101, 1001, 2 ** 40 + 1, 2 ** 70 - 1, 3 ** 50]:
            mm = apply_word(m, W)
            drops.add((al * m - be * g_topdown(m)) - (al * mm - be * g_topdown(mm)))
        lem6['slopes'] += 1
        ok = (len(drops) == 1 and Kex.denominator == 1 and int(Kex) == next(iter(drops)) and Kex > 0
              and int(Kex) % 2 == 0 and W[-1] == 'Z' and 'O' in W and periodic and runs_ok
              and W[0] == ('O' if al % 2 == 0 else 'Z'))
        if not ok:
            lem6['bad'].append([al, be])
        K = int(Kex)
        lem6['K_min'] = K if lem6['K_min'] is None else min(lem6['K_min'], K)
        lem6['K_max'] = K if lem6['K_max'] is None else max(lem6['K_max'], K)
    R['lemma6_clock'] = lem6
    check(not lem6['bad'], 'Lemma 6')

    log('theorem A0')
    # ------------------------------------------------------------------ 4. Theorem A0 by its own route (modular gadget head)
    a0 = {'slopes': 0, 'targets': 0, 'bad': [], 'max_bits': 0, 'max_reps': 0, 'jmg_residue_bad': 0}
    for M in range(2, 41, 2):                      # JOINT_MODULAR_GADGET Theorem 1 on its own
        for u in range(M):
            for z in range(1, M, 2):
                L = jmg_L(M, u, z)
                a0['jmg_residue_bad'] += not (L is not None and L % M == u and g_topdown(L) % M == z)
    check(a0['jmg_residue_bad'] == 0, 'JMG Theorem 1 residues')
    big = [-2001, -999, 1001, 4096, 12345]
    small = [-101, -11, -7, -3, -1, 0, 1, 2, 3, 5, 6, 7, 8, 15, 98, 99, 100]
    a0['jmg_head_witnesses'] = 0
    a0['branch_head_witnesses'] = 0
    a0['branch_route_over_letter_cap_then_direct_head'] = []
    a0['letter_cap'] = LETTER_CAP
    for (al, be) in interior_slopes(3, 31 if not quick else 11):
        if be % 2 == 0:
            continue
        a0['slopes'] += 1
        tset = small + (big if be <= 15 else [])
        K = int(clock_drop_exact(al, be)[3])
        for r in tset:
            if (r - (al - be)) % 2:
                continue
            a0['targets'] += 1
            use_jmg = (K == 2 and abs(r) <= 101)          # the gadget's own L is small enough when K = 2
            got = a0_witness(al, be, r, use_jmg=use_jmg)
            if isinstance(got, str):
                a0['branch_route_over_letter_cap_then_direct_head'].append([al, be, r, K, got])
                got = direct_witness(al, be, r)
                if isinstance(got, str):
                    a0['bad'].append([al, be, r, got])
                    continue
            else:
                a0['jmg_head_witnesses'] += use_jmg
                a0['branch_head_witnesses'] += 1
            n, K2, L, h = got
            ok = n % 2 == 1 and n >= 3 and F3(al, be, n) == r
            if not ok:
                a0['bad'].append([al, be, r])
            a0['max_bits'] = max(a0['max_bits'], n.bit_length())
            a0['max_reps'] = max(a0['max_reps'], h)
    R['theorem_A0_literal_route'] = a0
    check(not a0['bad'], 'Theorem A0 route')

    log('census')
    # ------------------------------------------------------------------ 5. census of S(a, b) (C, odd n < 2^LOGN)
    LOGN = 22 if quick else 24
    BMAX = 32 if quick else 64
    slopes = interior_slopes(3, BMAX)
    rows, meta = [], {}
    for (a, b) in slopes:
        s = v2(b)
        M = 2 ** (s + 1)
        if (b >> s) > 1 and s >= 1:
            _, ds, L, tail, low, a2, b2 = t2_data(a, b)
            cs = cstar_from_digits(a, b)
            rows.append((a, b, M, cs, L, low, 2000, 20000))
        else:
            cs = ((a - b) % 2) if (b >> s) > 1 else -1
            rows.append((a, b, M, cs, 0, 0, 2000, 20000))
        meta[(a, b)] = cs
    t0 = time.time()
    lemma_head, cen = H.census(LOGN, rows)
    R['census_parameters'] = {'odd_n_below': 2 ** LOGN, 'interior_slopes_b_le': BMAX, 'slopes': len(slopes),
                              'two_sided_threshold': 1 << 18}
    R['lemma1_violations_n_below_2^LOGN'] = lemma_head[0]
    check(lemma_head[0] == 0, 'Lemma 1')
    X = 1 << 18
    summary = {'slopes': len(slopes), 'b_prime_gt_1': 0, 'dyadic': 0, 'two_sided_class_count_ok': 0,
               'two_sided_equals_T2_class': 0, 'mixed_or_violations': [], 'nonfull_classes': 0,
               'class_attains_bound': 0, 'class_misses_bound': 0, 'first_miss': None,
               'atomgroups': 0, 'atomgroups_attained': 0, 'suffix_n': 0, 'suffix_violations': 0,
               'odd_b_slopes': 0, 'even_b_nondyadic': 0}
    t3 = {'slopes': 0, 'plus': 0, 'minus': 0, 'neither': 0, 'bad': []}
    for (a, b) in slopes:
        s = v2(b)
        bp = b >> s
        M = 2 ** (s + 1)
        dat = cen[(a, b)]['classes']
        two = [c for c, v in dat.items() if v and v['min'] < -X and v['max'] > X]
        if bp > 1:
            summary['b_prime_gt_1'] += 1
            if s == 0:
                summary['odd_b_slopes'] += 1
            else:
                summary['even_b_nondyadic'] += 1
        else:
            summary['dyadic'] += 1
        if len(two) == (1 if bp > 1 else 0):
            summary['two_sided_class_count_ok'] += 1
        else:
            check(False, 'two-sided class count %d/%d: %s' % (a, b, two))
        if bp > 1 and two == [meta[(a, b)]]:
            summary['two_sided_equals_T2_class'] += 1
        if s == 0:
            continue
        _, _, cls = classify(a, b)
        nfull = sum(1 for v in cls.values() if v[0] == 'FULL')
        check(nfull == (1 if bp > 1 else 0), 'automaton full count %d/%d' % (a, b))
        if bp > 1:
            check(cls[meta[(a, b)]][0] == 'FULL', 'automaton FULL class equals T2 class %d/%d' % (a, b))
        groups = {}
        lgroups = {}
        for c, (kind, bd, key) in cls.items():
            if kind == 'FULL':
                continue
            summary['nonfull_classes'] += 1
            v = dat[c]
            if kind not in ('POS', 'NEG'):
                summary['mixed_or_violations'].append([a, b, c, kind])
                continue
            if v is None:
                summary['mixed_or_violations'].append([a, b, c, 'empty'])
                continue
            if (kind == 'POS' and v['min'] < bd) or (kind == 'NEG' and v['max'] > bd):
                summary['mixed_or_violations'].append([a, b, c, kind, bd, v['min'], v['max']])
            ext = v['min'] if kind == 'POS' else v['max']
            if ext == bd:
                summary['class_attains_bound'] += 1
            else:
                summary['class_misses_bound'] += 1
                if summary['first_miss'] is None:
                    summary['first_miss'] = [a, b, c, kind, bd, ext]
            groups.setdefault((kind, bd, key), []).append(ext)
            lgroups.setdefault((kind, bd), []).append(ext)
        for (kind, bd, key), exts in groups.items():
            summary['atomgroups'] += 1
            gext = min(exts) if kind == 'POS' else max(exts)
            summary['atomgroups_attained'] += (gext == bd)
        for (kind, bd), exts in lgroups.items():
            summary['labelgroups'] = summary.get('labelgroups', 0) + 1
            gext = min(exts) if kind == 'POS' else max(exts)
            summary['labelgroups_attained'] = summary.get('labelgroups_attained', 0) + (gext == bd)
        if s >= 2 and bp > 1 and all(d == 0 for d in digits(a, b, s + 2)[1:s + 1]):
            Fm = a * (2 ** s - 1) - b * g_topdown(2 ** s - 1)
            summary['t3_cplus_max_redundant'] = summary.get('t3_cplus_max_redundant', 0) + (max(5 * (a - b), Fm) == Fm)
            summary['t3_cplus_s_ge_2_cases'] = summary.get('t3_cplus_s_ge_2_cases', 0) + 1
        if bp > 1:
            sc, sv = cen[(a, b)]['suffix']
            summary['suffix_n'] += sc
            summary['suffix_violations'] += sv
            # T3
            t3['slopes'] += 1
            plus, minus, cplus, cminus = t3_pred(a, b)
            cs = meta[(a, b)]
            outside = [v for c, v in dat.items() if c != cs and v]
            omax = max(v['max'] for v in outside)
            omin = min(v['min'] for v in outside)
            ok = True
            if plus:
                t3['plus'] += 1
                ok = ok and omax == cplus
            else:
                ok = ok and omax > X
            if minus:
                t3['minus'] += 1
                ok = ok and omin == cminus
            else:
                ok = ok and omin < -X
            if not plus and not minus:
                t3['neither'] += 1
            if plus and minus:
                ok = False
            if not ok:
                t3['bad'].append([a, b, plus, minus, cplus, cminus, omax, omin])
    R['census_T1_T2'] = summary
    R['census_T3'] = t3
    check(not summary['mixed_or_violations'], 'class bound violations')
    check(summary['suffix_violations'] == 0, 'T2 suffix violations')
    check(not t3['bad'], 'T3')

    # P2-style window counts (b <= 48, b' > 1, |r| <= 2000, direct n < 2^LOGN)
    p2 = {'slopes': 0, 'targets': 0, 'direct_hits': 0}
    for (a, b) in slopes:
        if b > 48:
            continue
        s = v2(b)
        if (b >> s) == 1:
            continue
        p2['slopes'] += 1
        M = 2 ** (s + 1)
        cs = meta[(a, b)]
        tgt = sum(1 for r in range(-2000, 2001) if r % M == cs)
        p2['targets'] += tgt
        p2['direct_hits'] += tgt - len(cen[(a, b)]['miss'])
    R['P2_window_counts'] = p2

    log('witnesses')
    # ------------------------------------------------------------------ 6. witnesses for unhit C* targets (T2 recipe)
    wit = {'slopes': 0, 'attempted': 0, 'ok': 0, 'bad': [], 'max_bits': 0, 'max_reps': 0, 'via_A0_head': 0,
           'A0_branch_over_cap_fell_back_to_direct_head': 0}
    rng = random.Random(7)
    wit['odd_b_slopes'] = 0
    for (a, b) in slopes:
        if b > 48:
            continue
        s = v2(b)
        if (b >> s) == 1:
            continue
        miss = cen[(a, b)]['miss']
        if not miss:
            continue
        wit['slopes'] += 1
        wit['odd_b_slopes'] += (s == 0)
        nsm = 10 if s >= 1 else 4
        sample = sorted(set(sorted(miss, key=lambda r: (abs(r), r))[:nsm] + rng.sample(miss, min(nsm, len(miss)))))
        _, ds, L, tail, low, a2, b2 = t2_data(a, b)
        check(gcd(a2, b2) == 1 and b2 % 2 == 1 and HALF < Fraction(a2, b2) < ONE, 'T2 head slope %d/%d' % (a, b))
        kap = {(a * apply_word(m, tail) - b * g_topdown(apply_word(m, tail))) - (2 ** s) * (a2 * m - b2 * g_topdown(m))
               for m in range(3, 200, 2)}
        check(len(kap) == 1, 'tail lemma constant %d/%d' % (a, b))
        kappa = kap.pop()
        for idx, r in enumerate(sample):
            wit['attempted'] += 1
            t = r - kappa
            if t % (2 ** s):
                wit['bad'].append([a, b, r, 'parity'])
                continue
            t //= 2 ** s
            got = 'skip'
            if idx % 5 == 4:            # every fifth target: head by Theorem A0's own route (residues of L, g(L))
                got = a0_witness(a2, b2, t)
                if isinstance(got, str):
                    wit['A0_branch_over_cap_fell_back_to_direct_head'] += 1
                else:
                    wit['via_A0_head'] += 1
            if isinstance(got, str):
                got = direct_witness(a2, b2, t)
            if isinstance(got, str):
                wit['bad'].append([a, b, r, got])
                continue
            m, _, _, h = got
            n = apply_word(m, tail)
            ok = n % 2 == 1 and n >= 3 and F3(a, b, n) == r and (n - 1) % (1 << L) == low
            if ok:
                wit['ok'] += 1
            else:
                wit['bad'].append([a, b, r])
            wit['max_bits'] = max(wit['max_bits'], n.bit_length())
            wit['max_reps'] = max(wit['max_reps'], h)
    R['T2_recipe_witnesses'] = wit
    check(not wit['bad'], 'T2 recipe witnesses')

    log('extras')
    # ------------------------------------------------------------------ 7. adversarial extras: larger b, exterior slopes
    extra = [(a, b) for b in (96, 100, 112, 120, 128, 136, 160, 192, 200, 224, 240, 255, 256)
             for a in range(b // 2 + 1, b) if gcd(a, b) == 1]
    rng2 = random.Random(11)
    extra = sorted(rng2.sample(extra, 60)) + [(129, 256), (255, 256), (97, 192), (191, 192), (65, 128), (127, 128)]
    exterior = [(1, 3), (2, 5), (1, 2), (3, 7), (5, 11), (1, 4), (3, 8), (5, 12), (4, 3), (5, 4), (7, 6), (9, 8),
                (11, 12 * 0 + 10), (-1, 4), (-3, 8), (0, 1), (1, 1), (3, 2), (13, 8), (7, 16)]
    rows2, meta2 = [], {}
    for (a, b) in extra + exterior:
        s = v2(b)
        M = 2 ** (s + 1)
        inside = HALF < Fraction(a, b) < ONE and (b >> s) > 1
        cs = cstar_from_digits(a, b) if (inside and s >= 1) else ((a - b) % 2 if inside else -1)
        meta2[(a, b)] = (cs, inside)
        rows2.append((a, b, M, cs, 0, 0, 200, 1000))
    _, cen2 = H.census(22, rows2)
    ext = {'interior_large_b': 0, 'interior_ok': 0, 'exterior': 0, 'exterior_ok': 0, 'bad': []}
    for (a, b) in extra + exterior:
        cs, inside = meta2[(a, b)]
        s = v2(b)
        dat = cen2[(a, b)]['classes']
        two = sorted(c for c, v in dat.items() if v and v['min'] < -X and v['max'] > X)
        x = Fraction(a, b)
        if HALF < x < ONE:
            ext['interior_large_b'] += 1
            ok = two == ([cs] if (b >> s) > 1 else [])
            if s >= 1:
                _, _, cls = classify(a, b)
                for c, (kind, bd, key) in cls.items():
                    v = dat.get(c)
                    if kind == 'FULL':
                        ok = ok and c == cs
                    elif kind == 'POS':
                        ok = ok and (v is None or v['min'] >= bd)
                    elif kind == 'NEG':
                        ok = ok and (v is None or v['max'] <= bd)
                    else:
                        ok = False
                plus, minus, cplus, cminus = t3_pred(a, b)
                if (b >> s) > 1:
                    outside = [v for c, v in dat.items() if c != cs and v]
                    omax = max(v['max'] for v in outside)
                    omin = min(v['min'] for v in outside)
                    ok = ok and ((omax == cplus) if plus else (omax > X)) and ((omin == cminus) if minus else (omin < -X))
            ext['interior_ok'] += ok
        else:
            ext['exterior'] += 1
            allv = [v for v in dat.values() if v]
            lo, hi = min(v['min'] for v in allv), max(v['max'] for v in allv)
            if x >= 1:
                ok = lo == 3 * (a - b) and not two
            else:
                ok = hi == 3 * (a - b) and not two
            ext['exterior_ok'] += ok
        if not ok:
            ext['bad'].append([a, b])
    R['census_extra_slopes_n_below_2^22'] = ext
    check(not ext['bad'], 'extra slopes')

    log('examples')
    # ------------------------------------------------------------------ 8. worked examples of PROOF.md / README / THM
    exm = {}
    for (a, b, cls_claim, word_claim) in [(5, 6, 3, 'ZZ'), (7, 12, 5, 'ZOOZ'), (11, 12, 7, 'ZZZ')]:
        s, ds, L, tail, low, a2, b2 = t2_data(a, b)
        exm['%d/%d' % (a, b)] = {'C*': cstar_from_digits(a, b), 'suffix': 'Z' + tail, 'claim': [cls_claim, word_claim],
                                 'ok': cstar_from_digits(a, b) == cls_claim and 'Z' + tail == word_claim}
        check(exm['%d/%d' % (a, b)]['ok'], 'worked example %d/%d' % (a, b))
    # S(7,10) on (-15, oo) = {3 mod 4};  S(5,6) on (-oo,-3) = {3 mod 4};  S(5,6) on (-oo,0) = {r<0, 3 mod 4} u {-3};
    # S(7,12) on (0, oo) = {5 mod 8}: the "outside C*" parts from the census, the C* parts are T1/T2
    d710 = cen[(7, 10)]['classes']
    d56 = cen[(5, 6)]['classes']
    d712 = cen[(7, 12)]['classes']
    exm['S(7,10)_outside_C*_max'] = max(v['max'] for c, v in d710.items() if c != 3 and v)
    exm['S(5,6)_outside_C*_min'] = min(v['min'] for c, v in d56.items() if c != 3 and v)
    exm['S(7,12)_outside_C*_max'] = max(v['max'] for c, v in d712.items() if c != 5 and v)
    check(exm['S(7,10)_outside_C*_max'] == -15 and exm['S(5,6)_outside_C*_min'] == -3
          and exm['S(7,12)_outside_C*_max'] < 0, 'README half-line examples')
    R['worked_examples'] = exm

    log('T4')
    # ------------------------------------------------------------------ 9. T4: symbolic identity, bank, gcd classes, literal C2 graphs
    import sympy as sp
    Dv, Ev, Uv, dv, D, U = sp.symbols('Dv Ev Uv dv D U')
    v = (Dv, Ev, Uv, dv)
    Oq = (D, 0, 2 * U, 1)
    Jq = (D - 1 + U, D - 1, D + U, D)
    diff = sp.expand(M_form(v, Jq) - 2 * N_form(v, Oq))
    al_s = Dv + 4 * Uv - 2 * dv
    claim = sp.expand(al_s * D - (al_s + Ev - dv) * U + (-Dv + Uv + dv + 2 * Ev))
    R['T4_symbolic_identity'] = bool(sp.simplify(diff - claim) == 0)
    check(R['T4_symbolic_identity'], 'T4 symbolic')
    # J on a strong unit source q = (D,1,U,1): BILINEAR (9) says e_a = D - E, R_a = D + U, H = D - E + U, d = D
    jq_literal_bad = 0
    for n in [3, 5, 7, 9, 11, 13, 19, 27, 37, 45, 63]:
        out, root = tour_from_word(word_of(n) + 'J')
        p = H.graphs([(out, root)])[0]
        jq_literal_bad += p[:5] != (G20[n] - 1 + n, G20[n] - 1, G20[n] + n, G20[n], 1)
    R['T4_Jq_literal_profile_bad'] = jq_literal_bad
    check(jq_literal_bad == 0, 'Jq literal profile')

    # bank of physical sources by O/Z/J words from t (U <= 300); every short word also counted literally
    bank = {}
    frontier = [('', (1, 0, 2, 1))]
    while frontier:
        w, p = frontier.pop()
        if p in bank and len(bank[p]) <= len(w):
            continue
        bank[p] = w
        for ch, f in (('O', O_), ('Z', Z_), ('J', J_)):
            q = f(p)
            if q[2] <= 300:
                frontier.append((w + ch, q))
    lit_items, lit_keys = [], []
    for p, w in sorted(bank.items()):
        if len(w) <= 9:
            lit_items.append(tour_from_word(w))
            lit_keys.append(p)
    lit = H.graphs(lit_items)
    lit_bad = sum(1 for p, l in zip(lit_keys, lit) if l[:4] != p or l[4] != 1)
    bank_stats = {'sources': len(bank), 'literal_checked': len(lit_keys), 'literal_bad': lit_bad,
                  'alpha_even': 0, 'E_gt_d': 0, 'cone_fail': 0, 'slope_claim_fail': 0, 'gamma_gt_1': 0,
                  'gamma_values': {}, 'U_plus_E_minus_D_even': 0}
    for p in bank:
        bank_stats['U_plus_E_minus_D_even'] += (p[2] + p[1] - p[0]) % 2 == 0
    check(lit_bad == 0, 'bank literal profiles')
    gam_examples = []
    for p in sorted(bank):
        Dp, Ep, Up, dp = p
        al, be, c = t4_coeffs(p)
        bank_stats['alpha_even'] += al % 2 == 0
        bank_stats['E_gt_d'] += Ep > dp
        bank_stats['cone_fail'] += not (1 <= dp <= Dp <= Up and 0 <= Ep <= Dp)
        if Ep < dp:
            gm = gcd(al, dp - Ep)
            x = Fraction(be, al)
            bank_stats['slope_claim_fail'] += not (HALF < x < ONE and x.denominator % 2 == 1 and x.denominator == al // gm)
            if gm > 1:
                bank_stats['gamma_gt_1'] += 1
                bank_stats['gamma_values'][str(gm)] = bank_stats['gamma_values'].get(str(gm), 0) + 1
                gam_examples.append((len(bank[p]), p, bank[p], gm))
        else:
            bank_stats['slope_claim_fail'] += be != al
    # every bank source with E < d: offsets over q = (g(n),1,n,1), odd n < 2^14, lie in the class U-D+E+2d mod 2 gamma,
    # and N(v,Oq), M(v,Jq) are odd with M - 2N equal to the T4 formula
    n14 = np.arange(3, 1 << 14, 2, dtype=np.int64)
    g14 = GARR[3:1 << 14:2]
    cls_bad, form_bad, pairs = 0, 0, 0
    for p in sorted(bank):
        al, be, c = t4_coeffs(p)
        offs = al * g14 - be * n14 + c
        Np = (p[0] - p[1] + p[2]) * g14 + 2 * (p[0] + 3 * p[2] - p[3]) * n14 - (p[2] + p[1])          # N(v, Oq)
        Mp = ((p[0] + 2 * p[2]) * (g14 - 1 + n14) - (p[2] + p[3]) * (g14 - 1)
              + (2 * p[0] - p[1] + 6 * p[2] - p[3]) * (g14 + n14) - (p[2] + p[1]) * g14)                 # M(v, Jq)
        pairs += len(n14)
        form_bad += int(np.sum((Mp - 2 * Np) != offs)) + int(np.sum(Np % 2 == 0)) + int(np.sum(Mp % 2 == 0))
        if p[1] < p[3]:
            gm = gcd(al, p[3] - p[1])
            cls_bad += int(np.sum(offs % (2 * gm) != (p[2] - p[0] + p[1] + 2 * p[3]) % (2 * gm)))
    bank_stats['offset_pairs_checked'] = pairs
    bank_stats['offset_formula_or_parity_bad'] = form_bad
    bank_stats['offset_explicit_class_bad'] = cls_bad
    R['T4_bank'] = bank_stats
    check(bank_stats['alpha_even'] == 0 and bank_stats['E_gt_d'] == 0 and bank_stats['cone_fail'] == 0
          and bank_stats['slope_claim_fail'] == 0 and form_bad == 0 and cls_bad == 0
          and bank_stats['U_plus_E_minus_D_even'] == 0, 'T4 bank claims')

    # offsets over q = (g(n),1,n,1), odd n < 2^20: one class mod 2 gamma; window coverage of that class
    nq = np.arange(3, 1 << 20, 2, dtype=np.int64)
    gq = GARR[3:1 << 20:2]
    off_rep = []
    gam_examples.sort()
    test_sources = [((1, 0, 2, 1), 't'), ((1, 0, 4, 1), 'O'), ((3, 1, 3, 1), 'J'), ((15, 6, 18, 9), 'JJZJ')]
    for (_, p, w, gm) in gam_examples[:4]:
        test_sources.append((p, w))
    for p, w in test_sources:
        al, be, c = t4_coeffs(p)
        offs = al * gq - be * nq + c
        Ep, dp = p[1], p[3]
        rec = {'v': list(p), 'word': w, 'alpha': al, 'beta': be, 'c': c}
        if Ep < dp:
            gm = gcd(al, dp - Ep)
            mod = 2 * gm
            resid = sorted(set((offs % mod).tolist()))
            win = set(offs[(offs >= -199) & (offs <= 199)].tolist())
            cls_targets = [r for r in range(-199, 200) if r % mod == resid[0]]
            rec.update({'gamma': gm, 'residues_mod_2gamma': resid, 'window_hit': sum(r in win for r in cls_targets),
                        'window_class_size': len(cls_targets), 'odd_offsets_minus1_plus1_present': [(-1 in win), (1 in win)],
                        'explicit_class_U-D+E+2d_mod_2gamma': (p[2] - p[0] + p[1] + 2 * p[3]) % mod})
            check(len(resid) == 1, 'T4 offsets in one class %s' % (p,))
            check(resid[0] == (p[2] - p[0] + p[1] + 2 * p[3]) % mod, 'T4 explicit class %s' % (p,))
        else:
            rec.update({'E_eq_d_max_offset': int(offs.max()), 'c': c})
            check(int(offs.max()) <= c, 'T4 E = d bound %s' % (p,))
        off_rep.append(rec)
    R['T4_offset_classes'] = off_rep

    # the claimed example "any J-born source over a strong source" (d_v - E_v = 1): J over a strong source has d - E = E
    jb = []
    for word in ['JJZ', 'JJZZ', 'JZJZ', 'OJZ', 'JJJZ']:
        out, root = tour_from_word(word)
        strong = strongly_connected(out)
        p = run_word(word)
        q = J_(p)
        outJ, rootJ = tour_from_word(word + 'J')
        litp, litq = H.graphs([(out, root), (outJ, rootJ)])
        al, be, c = t4_coeffs(q)
        jb.append({'strong_source_word': word, 'profile': list(p), 'literal_profile_DEUds': list(litp[:5]),
                   'strong': strong, 'J_profile': list(q), 'J_literal_profile_DEUds': list(litq[:5]),
                   'd_minus_E_of_J': q[3] - q[1], 'gamma': gcd(al, q[3] - q[1])})
        check(tuple(litp[:4]) == p and tuple(litq[:4]) == q and litp[4] == 1 and litq[4] == 1, 'J-born literal %s' % word)
    R['T4_J_born_over_strong_sources'] = jb

    # literal C2 parent/child graphs built from BILINEAR_C2's text, counted by the subset DP
    lit_pairs = []
    lefts = [('', 'v=t'), ('O', 'O t'), ('OJO', 'OJO'), ('JO', 'J O'), ('J', 'J t (E = d)'), ('JJZJ', 'gamma = 3 source')]
    for lw, label in lefts:
        V, c0 = tour_from_word(lw)
        pv = run_word(lw)
        al, be, cc = t4_coeffs(pv)
        unit_gap = pv[1] < pv[3] and gcd(al, pv[3] - pv[1]) == 1
        chosen = {}
        for n in range(3, 1 << 12, 2):
            if len(V) + (n - 1).bit_length() + 4 > 20:          # child order <= 20 (subset DP)
                break
            dl = al * G20[n] - be * n + cc
            if dl in chosen:
                continue
            if unit_gap and dl not in (-1, 1):
                continue
            chosen[dl] = n
            if len(chosen) >= 2:
                break
        for dl, n in sorted(chosen.items()):
            Qo, ro = tour_from_word(word_of(n) + 'O')
            Qj, rj = tour_from_word(word_of(n) + 'J')
            TN, aN, cN = pair_N(V, c0)
            TM, uM, cM = pair_M(V, c0)
            parent = graft(TN, aN, cN, converse(Qo), ro)
            child = graft(TM, uM, cM, converse(Qj), rj)
            okt = is_tournament(parent) and is_tournament(child)
            (Pp, Pc) = H.graphs([(parent, 0), (child, 0)])
            q = (G20[n], 1, n, 1)
            nN = N_form(pv, O_(q))
            mM = M_form(pv, J_(q))
            rec = {'left_word': lw, 'label': label, 'left_profile': list(pv), 'q_n': n, 'offset_formula': dl,
                   'parent_order': len(parent), 'child_order': len(child), 'parent_H': Pp[0], 'child_H': Pc[0],
                   'parent_C': Pp[5], 'child_C': Pc[5], 'N_formula': nN, 'M_formula': mM}
            rec['ok'] = (okt and Pp[0] == nN and Pc[0] == mM and Pp[5] == 2 and Pc[5] == 2 and Pc[0] - 2 * Pp[0] == dl)
            lit_pairs.append(rec)
            check(rec['ok'], 'literal C2 pair %s n=%d' % (lw, n))
    R['T4_literal_C2_pairs'] = lit_pairs

    log('T5')
    # ------------------------------------------------------------------ 10. T5 clock uniqueness (and its scope)
    def word_changes(a, b, k, heads):
        """For all 2^k words (A = sum e_i 2^i, e_i = 1 for Z): m . w = 2^k m - rev_k(A); returns change matrix."""
        A = np.arange(2 ** k, dtype=np.int64)
        rev = np.zeros_like(A)
        for i in range(k):
            rev |= ((A >> i) & 1) << (k - 1 - i)
        rows = []
        for m in heads:
            nn = (2 ** k) * m - rev
            rows.append(a * nn - b * GARR[nn] - (a * m - b * G20[m]))
        return A, np.array(rows)

    t5 = {'interior_slopes': 0, 'words': 0, 'bad': [], 'outside_interval_source_independent_words': []}
    for (a, b) in interior_slopes(3, 31 if not quick else 13):
        t5['interior_slopes'] += 1
        Wc = clock(a, b)[0] if b % 2 else None
        for k in range(1, 15 if not quick else 9):
            A, ch = word_changes(a, b, k, (3, 5, 7, 9, 11))
            indep = np.all(ch == ch[0], axis=0)
            pw = np.zeros(len(A), dtype=bool)
            if Wc is not None and k % len(Wc) == 0:
                Aw = sum(1 << i for i, cch in enumerate(Wc * (k // len(Wc))) if cch == 'Z')
                pw[Aw] = True
            t5['words'] += len(A)
            badA = np.nonzero(indep != pw)[0]
            for x in badA[:5]:
                t5['bad'].append([a, b, k, int(x)])
    for (a, b) in [(1, 3), (2, 5), (1, 5), (3, 7), (1, 7), (2, 9), (4, 9), (0, 1), (1, 1), (4, 3), (6, 5)]:
        found = []
        for k in range(1, 11):
            A, ch = word_changes(a, b, k, (3, 5, 7, 9, 11, 13))
            for x in np.nonzero(np.all(ch == ch[0], axis=0))[0]:
                w = ''.join('Z' if (int(x) >> i) & 1 else 'O' for i in range(k))
                found.append((w, int(ch[0][x])))
        t5['outside_interval_source_independent_words'].append({'slope': '%d/%d' % (a, b), 'count_len_le_10': len(found),
                                                                'shortest': found[:2]})
    R['T5'] = t5
    check(not t5['bad'], 'T5 interior')

    log('P5')
    # ------------------------------------------------------------------ 11. P5 densities at 2^22 versus 2^26 (stability)
    p5slopes = [(3, 4), (5, 6), (7, 10), (9, 10), (7, 12), (11, 12), (9, 14), (13, 14), (5, 8), (7, 8),
                (13, 24), (17, 24), (19, 24), (23, 24)]
    rows5 = []
    for (a, b) in p5slopes:
        s = v2(b)
        rows5.append((a, b, 2 ** (s + 1), -1, 0, 0, 200, 20000))
    dens = {}
    boundary = {}
    for lg in ((22, 26) if not quick else (20, 22)):
        _, c5 = H.census(lg, rows5)
        for (a, b) in p5slopes:
            _, M, cls = classify(a, b)
            for c, (kind, bd, key) in cls.items():
                if kind not in ('POS', 'NEG'):
                    continue
                boundary['%d/%d:%d' % (a, b, c)] = any(("'1'" in kk or "'1/2'" in kk) for kk in key)
                v = c5[(a, b)]['classes'].get(c)
                if not v:
                    continue
                if kind == 'POS':
                    hit = v['hit_pos']
                    tot = sum(1 for r in range(500, 20001) if r % M == c)
                else:
                    hit = v['hit_neg']
                    tot = sum(1 for r in range(-20000, -499) if r % M == c)
                dens.setdefault('%d/%d:%d' % (a, b, c), {})[str(lg)] = [hit, tot, round(hit / tot, 6)]
    changed = {k: v for k, v in dens.items() if len(set(x[0] for x in v.values())) > 1}
    big50 = {k: v for k, v in dens.items() if v[max(v)][0] >= 50}
    dec = {k: v for k, v in big50.items() if not boundary[k]}
    R['P5_density_hits_n_2^22_vs_2^26'] = {
        'classes_with_data': len(dens), 'classes_ge_50_hits': len(big50),
        'without_slope_1_or_half_piece': len(dec), 'with_slope_1_or_half_piece': len(big50) - len(dec),
        'changed_between_ranges': changed,
        'changed_without_boundary_piece': sorted(k for k in changed if not boundary[k]),
        'decided_fraction_range_2^22': [min(v[min(v)][2] for v in dec.values()), max(v[min(v)][2] for v in dec.values())],
        'boundary_fraction_range_2^22': [min((v[min(v)][2] for k, v in big50.items() if boundary[k]), default=None),
                                         max((v[min(v)][2] for k, v in big50.items() if boundary[k]), default=None)],
        'boundary_fraction_range_largest_n': [min((v[max(v)][2] for k, v in big50.items() if boundary[k]), default=None),
                                              max((v[max(v)][2] for k, v in big50.items() if boundary[k]), default=None)]}

    # exterior slopes (every class one-side bounded): window fractions on the unbounded side at 2^22 and 2^26, and the
    # boundary pieces S(1,1) = {n - g(n)}, S(1,2) = {n - 2 g(n)}: every value of modulus <= X already at n <= 4X + 8
    ext_sl = [(2, 1), (3, 1), (3, 2), (5, 2), (4, 3), (5, 3), (5, 4), (1, 1), (0, 1), (1, 3), (1, 4), (2, 5), (1, 2),
              (3, 8), (9, 8)]
    rowsx = [(a, b, 2 ** (v2(b) + 1), -1, 0, 0, 200, 20000) for (a, b) in ext_sl]
    exd = {}
    for lg in ((22, 26) if not quick else (20, 22)):
        _, cx = H.census(lg, rowsx)
        for (a, b) in ext_sl:
            M = 2 ** (v2(b) + 1)
            for c, v in sorted(cx[(a, b)]['classes'].items()):
                if not v:
                    continue
                if a >= b:
                    hit, tot = v['hit_pos'], sum(1 for r in range(500, 20001) if r % M == c)
                else:
                    hit, tot = v['hit_neg'], sum(1 for r in range(-20000, -499) if r % M == c)
                exd.setdefault('%d/%d:%d' % (a, b, c), {})[str(lg)] = round(hit / tot, 6)
    Xw = 20000
    nn = np.arange(3, 4 * Xw + 9, 2, dtype=np.int64)
    s11 = np.unique((nn - GARR[nn])[((nn - GARR[nn]) >= 500) & ((nn - GARR[nn]) <= Xw)])
    v12 = nn - 2 * GARR[nn]
    s12 = np.unique(v12[(v12 <= -500) & (v12 >= -Xw)])
    big = max(exd['1/1:0'])
    R['exterior_and_boundary_densities'] = {
        'fractions': exd,
        'changed_between_ranges': sorted(k for k, v in exd.items() if len(set(v.values())) > 1),
        'S(1,1)_window_hits_n_le_4X+8': int(len(s11)), 'S(1,1)_fraction_n_le_4X+8': round(len(s11) / len(range(500, Xw + 1, 2)), 6),
        'S(1,2)_window_hits_n_le_4X+8': int(len(s12)),
        'S(1,2)_fraction_n_le_4X+8': round(len(s12) / (len(range(-Xw, -499)) // 2), 6)}
    check(abs(R['exterior_and_boundary_densities']['S(1,1)_fraction_n_le_4X+8'] - exd['1/1:0'][big]) < 1e-6,
          'slope-1 boundary piece decided by n <= 4X + 8')

    log('P6')
    # ------------------------------------------------------------------ 12. P6 witness lengths (own DFS, depth 32)
    if not quick:
        p6 = {}
        for (a, b) in ((3, 5), (8, 9)):
            best = H.wl(a, b, 32, 1024)
            par = (a - b) % 2
            rowsj = {}
            for j in range(1, 11):
                w = 2 ** j
                tg = [r for r in range(-w, w + 1) if r % 2 == par]
                found = [best[r][0] for r in tg if r in best]
                rowsj[j] = {'targets': len(tg), 'missing_at_depth_32': len(tg) - len(found),
                            'maxlen': max(found) if found else None}
            p6['%d/%d' % (a, b)] = rowsj
        R['P6_own_dfs_depth_32'] = p6
        check(p6['8/9'][4]['maxlen'] == 21 and p6['8/9'][4]['missing_at_depth_32'] == 0, 'P6 8/9 j=4')
        # P6 (b) at j = 9 allows 3*9 + 8 = 35 letters, so the packet's depth-34 search cannot decide it: depth 35
        best = H.wl(8, 9, 35, 512)
        tg = [r for r in range(-512, 513) if r % 2 == 1]
        miss = [r for r in tg if r not in best]
        R['P6_8_9_j9_depth_35'] = {'targets': len(tg), 'missing_at_depth_35': miss,
                                   'maxlen_found': max(best[r][0] for r in tg if r in best),
                                   'bound_3j_plus_8': 35, 'clause_b_fails_at_j9': len(miss) > 0}

    R['failures'] = FAILS
    R['verdict_of_checks'] = 'all checks passed' if not FAILS else 'FAILURES PRESENT'
    text = json.dumps(R, sort_keys=True, indent=1, default=str) + '\n'
    outp.write_text(text)
    print(text[:3000])
    print('elapsed %.1f s' % (time.time() - t_start), file=sys.stderr)


if __name__ == '__main__':
    main()
