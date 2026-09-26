#!/usr/bin/env python3
"""collatz_dipspectrum_20260926_audit.py

Independent adversarial audit of THM-4487 (the dip spectrum of 3n+-1), i.e. of
05-knowledge/results/collatz_dipspectrum_20260926_entropy_curve.md.  Written blind to
collatz_dipspectrum_20260926.py; the only data taken from the author's lane are the
published counts of collatz_dipspectrum_20260926.out, used for comparison in section C.

Definitions (as in the note):
  b = +-1, T_b(x) = x/2 (x even), (3x+b)/2 (x odd), alpha = log_2 3, h = binary entropy,
  D_b(X, gamma) = #{ n in [1, X] : T_b^i(n) >= n^gamma for all 0 <= i <= floor(log_2 n) },
  F_b(X, theta) = { y in [1, X] : T_b^i(y) >= y X^(-theta) for all 0 <= i <= floor(log_2 X) }.

Sections:
  A. constants and the Corollary 3 identities (mpmath, 30 digits)
  B. cycle lemma: exhaustive check on all binary words of length <= 12 (steps +1 even,
     1 - alpha odd), and the count 'at least C(t,o)/t good words' up to t = 16
  C. exact counts D_b(2^t, gamma), both sheets, t = 10..max_t (vectorised by dyadic block),
     compared with the note's table; log D/log X and four-doubling slopes
  D. upper-bound mechanism (1.2): every counted n >= n_0(gamma) has o_t(n) >= rho t - 2
  E. lower-bound construction (1.3, rotation): representatives of rotated words and the
     no-dip property over the whole window, both sheets, gamma in {0.85, 0.91, 0.97, 1}
  E2. auditor's alternative repair for (gamma = 1, minus sheet): tilted cycle lemma
  F. block construction (1.3, gamma = 1, minus sheet): W_L, c_L, concatenations, failures
  G. Proposition 2: brute-force #F_b(2^k, theta) and the (X/2, X] construction of 1.4
  H. the entropy / binomial bounds quoted in 1.1

Usage: python3 collatz_dipspectrum_20260926_audit.py [max_t]      (default 24; numpy needed)
"""
import sys, math, itertools, time
from fractions import Fraction
from math import comb, ceil

try:
    import numpy as np
except ImportError:          # pragma: no cover
    np = None
try:
    from mpmath import mp, mpf, log as mlog, sqrt as msqrt, findroot
    mp.dps = 30
    HAVE_MP = True
except ImportError:          # pragma: no cover
    HAVE_MP = False

ALPHA = math.log2(3.0)
LOG2_32 = math.log2(1.5)
LOG4_3 = ALPHA / 2

def h(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)

def T(x, b):
    return x >> 1 if x % 2 == 0 else (3 * x + b) >> 1

def hdr(s):
    return "\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78

def rho_of(gamma):
    return max(0.5, gamma / ALPHA)

def ceil_safe(x):
    """ceil of a float that is never an integer in our uses (irrational multiples)."""
    c = math.ceil(x)
    if abs(x - round(x)) < 1e-9:
        raise ValueError("ceil at an integer boundary: %r" % x)
    return c

# --------------------------------------------------------------------------------------
# Terras inversion: residue r in [0, 2^t) whose parity word of length t is `word`
# --------------------------------------------------------------------------------------
def terras_rep(word, b):
    t = len(word)
    r = 0
    for k in range(t - 1, -1, -1):
        m = t - k
        mod = 1 << m
        if word[k] == 0:
            r = (2 * r) % mod
        else:
            r = ((2 * r - b) * pow(3, -1, mod)) % mod
    return r

def parity_word(n, b, t):
    w = []
    x = n
    for _ in range(t):
        w.append(x & 1)
        x = T(x, b)
    return w

def orbit_min(n, b, steps):
    x = n
    mn = n
    for _ in range(steps):
        x = T(x, b)
        if x < mn:
            mn = x
    return mn

def partial_sums(word):
    S = [0.0]
    s = 0.0
    for c in word:
        s += 1.0 if c == 0 else 1.0 - ALPHA
        S.append(s)
    return S

# ======================================================================================
def section_A():
    print(hdr("A. constants and the Corollary 3 identities"))
    if not HAVE_MP:
        print("mpmath not available; skipping high-precision checks")
        return
    L2 = lambda x: mlog(x) / mlog(2)
    a = L2(3)
    H = lambda p: -(p * L2(p) + (1 - p) * L2(1 - p))
    p0 = 1 / a
    hstar = H(p0)
    closed = L2(a) - (1 - 1 / a) * L2(a - 1)
    hp = L2(a - 1)                                   # claimed h'(1/alpha) = log2(alpha - 1)
    eps = mpf(10) ** -12
    hp_num = (H(p0 + eps) - H(p0 - eps)) / (2 * eps)
    E1 = hp / a                                      # E'(1) = h'(1/alpha)/alpha
    lam_closed = mlog(1 / (a - 1)) / mlog(3)         # claimed lambda* = log_3(1/log_2(3/2))
    M = lambda l: (mpf(2) ** (-l) + (mpf(3) / 2) ** l) / 2
    dM = lambda l: (-mlog(2) * mpf(2) ** (-l) + mlog(mpf(3) / 2) * (mpf(3) / 2) ** l) / 2
    lam_num = findroot(dM, mpf("0.5"))
    Mmin = M(lam_num)
    rate = -L2(Mmin)
    print("alpha = log_2 3            = %s" % mp.nstr(a, 12))
    print("log_3 2 = 1/alpha          = %s" % mp.nstr(p0, 12))
    print("h* = h(log_3 2)            = %s   (note: 0.949956)" % mp.nstr(hstar, 12))
    print("closed form log2(a)-(1-1/a)log2(a-1) = %s   |diff| = %s" % (mp.nstr(closed, 12), mp.nstr(abs(closed - hstar), 3)))
    print("1 - h*                     = %s   (note: 0.050044)" % mp.nstr(1 - hstar, 12))
    print("1/h*                       = %s   (note: 1.052681)" % mp.nstr(1 / hstar, 12))
    print("log_4 3 = alpha/2          = %s   (note: 0.792481)" % mp.nstr(a / 2, 12))
    print("theta_0 = 1 - alpha/2      = %s   (note: 0.207519)" % mp.nstr(1 - a / 2, 12))
    print("log_2(3/2) = alpha - 1     = %s   (note: 0.584963)" % mp.nstr(a - 1, 12))
    print("h'(1/alpha) = log2(alpha-1)= %s   numeric derivative = %s   |diff| = %s"
          % (mp.nstr(hp, 12), mp.nstr(hp_num, 12), mp.nstr(abs(hp - hp_num), 3)))
    print("   (note now prints -0.77358; the original draft printed -0.77353)")
    print("E'(1) = h'(1/alpha)/alpha  = %s ;  -E'(1) = %s   (note now 0.48807; original draft 0.48796)"
          % (mp.nstr(E1, 12), mp.nstr(-E1, 12)))
    print("lambda* closed form log_3(1/log_2(3/2)) = %s" % mp.nstr(lam_closed, 12))
    print("lambda* numeric minimiser of (2^-l + (3/2)^l)/2 = %s   |diff| = %s"
          % (mp.nstr(lam_num, 12), mp.nstr(abs(lam_num - lam_closed), 3)))
    print("identity lambda* = -E'(1):  |lambda* + E'(1)| = %s" % mp.nstr(abs(lam_closed + E1), 3))
    print("min value M(lambda*)       = %s ;  -log2 M = %s ;  1 - h* = %s ;  |diff| = %s"
          % (mp.nstr(Mmin, 12), mp.nstr(rate, 12), mp.nstr(1 - hstar, 12), mp.nstr(abs(rate - (1 - hstar)), 3)))
    trap = (1 - hstar) + L2((3 + msqrt(13)) / 4)
    print("numerological trap: (1-h*) + log2((3+sqrt13)/4) = %s ; |h'(log_3 2)| = %s  (note: 0.7737 vs 0.7735; correct rounding 0.7736)"
          % (mp.nstr(trap, 8), mp.nstr(abs(hp), 8)))
    ok = (abs(closed - hstar) < mpf(10) ** -25 and abs(hp - hp_num) < mpf(10) ** -10
          and abs(lam_num - lam_closed) < mpf(10) ** -20 and abs(lam_closed + E1) < mpf(10) ** -25
          and abs(rate - (1 - hstar)) < mpf(10) ** -25)
    print("all Corollary 3 identities exact to working precision: %s" % ok)
    # the note's rounded decimals
    vals = {"0.949956": hstar, "0.050044": 1 - hstar, "0.792481": a / 2, "0.207519": 1 - a / 2,
            "1.052681": 1 / hstar, "0.48807": -E1, "-0.77358": hp}
    for k, v in vals.items():
        nd = len(k.split(".")[1])
        print("  note prints %s -> computed %s : %s" % (k, mp.nstr(v, nd + 1), "match" if abs(mpf(k) - v) < mpf(10) ** (-nd) * 0.51 else "MISMATCH"))

# ======================================================================================
def section_B(tmax_rot=12, tmax_count=16):
    print(hdr("B. cycle lemma (steps +1 for even letter, 1 - alpha for odd letter)"))
    steps = (1.0, 1.0 - ALPHA)
    bad = 0
    total = 0
    bad_examples = []
    for t in range(1, tmax_rot + 1):
        for mask in range(1 << t):
            w = [(mask >> i) & 1 for i in range(t)]
            S = partial_sums(w)
            mx = max(S)
            cands = {S.index(mx), t - S[::-1].index(mx)}       # first and last maximal index
            cap = max(0.0, S[t])
            for p in cands:
                w2 = w[p:] + w[:p]
                s = 0.0
                ok = True
                for c in w2:
                    s += steps[c]
                    if s > cap + 1e-9:
                        ok = False
                        break
                total += 1
                if not ok:
                    bad += 1
                    if len(bad_examples) < 5:
                        bad_examples.append((t, w, p))
    print("rotation after a maximal partial sum (first and last maximal index) gives all partial sums <= max(0, S_t):")
    print("  words of length 1..%d: %d rotations checked, %d failures %s"
          % (tmax_rot, total, bad, bad_examples if bad_examples else ""))
    # the counting consequence: for every (t, o), #good words >= C(t,o)/t
    print("count of words with o odd letters and all partial sums <= max(0, S_t), vs C(t,o)/t:")
    worst = None
    fails = []
    for t in range(1, tmax_count + 1):
        good = [0] * (t + 1)
        for mask in range(1 << t):
            o = bin(mask).count("1")
            cap = max(0.0, t - o * ALPHA)
            s = 0.0
            ok = True
            for i in range(t):
                s += steps[(mask >> i) & 1]
                if s > cap + 1e-9:
                    ok = False
                    break
            if ok:
                good[o] += 1
        for o in range(t + 1):
            need = comb(t, o) / t
            ratio = good[o] / need
            if good[o] + 1e-9 < need:
                fails.append((t, o, good[o], need))
            if t >= 4 and 1 <= o <= t - 1 and (worst is None or ratio < worst[0]):
                worst = (ratio, t, o, good[o], comb(t, o))
        if t in (8, 12, 16):
            print("  t=%2d: good counts by o = %s" % (t, good))
            print("        C(t,o)/t        = %s" % [round(comb(t, o) / t, 2) for o in range(t + 1)])
    print("  failures of good(t,o) >= C(t,o)/t for t <= %d: %s" % (tmax_count, fails if fails else "none"))
    print("  smallest ratio good/(C(t,o)/t): %.4f at t=%d, o=%d (good=%d, C=%d)" % worst)

# ======================================================================================
# The note's table labels its first column "0.792"; its counts turn out to be those of gamma = 0.7925 (slightly
# above Korec's log_4 3 = 0.792481), so that column is compared against 0.7925 below and log_4 3 is reported as an
# extra column.  The note's script also does not count n = 1 (which is in D_b for every gamma: the window has
# length 0), so its counts are ours minus 1.
NOTE_GAMMAS = [0.7925, 0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 1.0]
AUDIT_GAMMAS = [0.792, LOG4_3, 0.7925, 0.85, 0.91, 0.97, 1.0]
ALL_GAMMAS = sorted(set(NOTE_GAMMAS) | set(AUDIT_GAMMAS))
EXACT_GAMMA = {g: Fraction(str(g)) for g in ALL_GAMMAS if g != LOG4_3}
NOTE_OFFSET = 1   # our counts include n = 1
# published counts, collatz_dipspectrum_20260926.out (rows X = 2^t, columns NOTE_GAMMAS)
NOTE_COUNTS = {
    1: {10: [414, 344, 298, 250, 177, 145, 113, 89],
        12: [1659, 1367, 1161, 896, 670, 521, 350, 281],
        14: [6633, 5463, 4622, 3485, 2719, 1724, 1195, 874],
        16: [26700, 21709, 17971, 13358, 10361, 6479, 4394, 2903],
        18: [107492, 87504, 68659, 54136, 39110, 24655, 16056, 9245],
        20: [432355, 353655, 271543, 215550, 149967, 95636, 56262, 31730],
        22: [1739485, 1421172, 1076683, 823855, 545789, 368356, 196095, 105669],
        24: [6991365, 5701356, 4302718, 3188718, 2166439, 1376462, 675317, 367698]},
    -1: {10: [405, 346, 291, 247, 175, 144, 112, 89],
         12: [1646, 1374, 1154, 891, 668, 518, 349, 281],
         14: [6624, 5463, 4615, 3477, 2715, 1718, 1197, 874],
         16: [26691, 21702, 17975, 13345, 10364, 6474, 4396, 2903],
         18: [107475, 87494, 68611, 54133, 39112, 24646, 16053, 9245],
         20: [432293, 353635, 271447, 215557, 149954, 95586, 56252, 31730],
         22: [1739409, 1421041, 1076614, 823811, 545761, 368320, 196085, 105669],
         24: [6991232, 5701160, 4302522, 3188730, 2166364, 1376381, 675303, 367698]},
}

def block_iter(s, b, chunk=1 << 20):
    """n in [2^s, 2^(s+1)): yields (n, mn, odd): mn = min_{0<=i<=s} T_b^i(n), odd = #odd letters in the word of length s."""
    lo, hi = 1 << s, 1 << (s + 1)
    for start in range(lo, hi, chunk):
        n = np.arange(start, min(start + chunk, hi), dtype=np.int64)
        x = n.copy()
        mn = n.copy()
        odd = np.zeros_like(n)
        for _ in range(s):
            par = x & 1
            x = np.where(par == 1, (3 * x + b) >> 1, x >> 1)
            np.minimum(mn, x, out=mn)
            odd += par
        yield n, mn, odd

def masks_for_block(n, mn, gammas):
    ln = np.log2(n.astype(np.float64))
    lm = np.log2(mn.astype(np.float64))
    r = lm / ln
    res = {}
    for g in gammas:
        if g == 1.0:
            res[g] = mn >= n
            continue
        mask = r >= g - 1e-9
        amb = np.nonzero((r >= g - 1e-9) & (r < g + 1e-9))[0]
        if len(amb):
            if g in EXACT_GAMMA:
                fg = EXACT_GAMMA[g]
                p, q = fg.numerator, fg.denominator
                for idx in amb:
                    mask[idx] = int(mn[idx]) ** q >= int(n[idx]) ** p
            else:               # gamma = log_4 3: mn >= n^gamma  <=>  ln(mn) ln 4 >= ln 3 ln n (50 digits)
                with mp.workdps(50):
                    for idx in amb:
                        mask[idx] = mlog(int(mn[idx])) * mlog(4) >= mlog(3) * mlog(int(n[idx]))
        res[g] = mask
    return res

def section_CD(max_t):
    print(hdr("C. exact counts D_b(2^t, gamma) (window floor(log_2 n), n <= 2^t), both sheets"))
    if np is None:
        print("numpy not available; skipping")
        return None
    t0 = time.time()
    results = {}
    viol = {}
    for b in (1, -1):
        blockcount = {g: {} for g in ALL_GAMMAS}
        viol[b] = {g: [0, 0, 0] for g in ALL_GAMMAS}     # [#counted n >= n_0, violations of o >= rho t - 2, violations of o alpha >= t gamma + gamma - 2]
        for s in range(1, max_t):
            cnt = {g: 0 for g in ALL_GAMMAS}
            for n, mn, odd in block_iter(s, b):
                masks = masks_for_block(n, mn, ALL_GAMMAS)
                for g in ALL_GAMMAS:
                    m = masks[g]
                    cnt[g] += int(m.sum())
                    # section D bookkeeping
                    n0 = 2.0 ** (1.0 / (g - LOG2_32))
                    rho = g / ALPHA
                    sel = m & (n >= n0)
                    viol[b][g][0] += int(sel.sum())
                    viol[b][g][1] += int((sel & (odd < rho * s - 2 - 1e-9)).sum())
                    viol[b][g][2] += int((sel & (odd * ALPHA < s * g + g - 2 - 1e-9)).sum())
            for g in ALL_GAMMAS:
                blockcount[g][s] = cnt[g]
        # n = 1 always counts (window of length 0); n = 2^t never (T^t(2^t) = 1 < 2^(t gamma))
        assert all(orbit_min(1 << t, b, t) == 1 for t in range(1, 8))
        D = {g: {} for g in ALL_GAMMAS}
        for g in ALL_GAMMAS:
            acc = 1
            for t in range(1, max_t + 1):
                if t >= 2:
                    acc += blockcount[g][t - 1]
                D[g][t] = acc
        results[b] = D
    print("(elapsed %.1f s)" % (time.time() - t0))
    for b in (1, -1):
        D = results[b]
        print("\n-- sheet 3n%s1 --  (columns: gamma; the column 0.792481 is gamma = log_4 3; counts include n = 1)" % ("+" if b == 1 else "-"))
        print("gamma:          " + "".join("%10.6f" % g for g in ALL_GAMMAS))
        print("h(rho):         " + "".join("%10.4f" % h(rho_of(g)) for g in ALL_GAMMAS))
        for t in range(10, max_t + 1):
            if t % 2:
                continue
            print("X=2^%2d D:       " % t + "".join("%10d" % D[g][t] for g in ALL_GAMMAS))
        for t in range(10, max_t + 1):
            if t % 2:
                continue
            print("X=2^%2d logD/logX" % t + "".join("%10.4f" % (math.log2(D[g][t]) / t) for g in ALL_GAMMAS))
        for t in range(14, max_t + 1):
            if t % 2:
                continue
            print("slope 2^%2d->2^%2d" % (t - 4, t) + "".join("%10.4f" % ((math.log2(D[g][t]) - math.log2(D[g][t - 4])) / 4) for g in ALL_GAMMAS))
        # comparison with the note (its first column read as gamma = log_4 3; its counts exclude n = 1)
        mism = []
        ncmp = 0
        for t, row in NOTE_COUNTS[b].items():
            if t > max_t:
                continue
            for g, v in zip(NOTE_GAMMAS, row):
                ncmp += 1
                if D[g][t] - NOTE_OFFSET != v:
                    mism.append((t, round(g, 6), D[g][t], v))
        print("comparison with the note's published counts (t <= %d, %d entries), ours - 1 vs note: %s"
              % (max_t, ncmp, "ALL MATCH" if not mism else "MISMATCHES (t, gamma, ours, note) %s" % mism))
        if max_t >= 20:
            print("  D_%s(2^20, 1) = %d (note: 31730 + the n = 1 that the note's script does not count)" % ("+" if b == 1 else "-", D[1.0][20]))
    d = max(abs(results[1][g][t] - results[-1][g][t]) for g in ALL_GAMMAS for t in range(10, max_t + 1))
    print("\nmax |D_+ - D_-| over t in [10,%d] and all gammas: %d ; at gamma=1: %d"
          % (max_t, d, max(abs(results[1][1.0][t] - results[-1][1.0][t]) for t in range(10, max_t + 1))))
    print("per-gamma |D_+ - D_-| at 2^%d: %s" % (max_t, {g: results[1][g][max_t] - results[-1][g][max_t] for g in ALL_GAMMAS}))

    print(hdr("D. upper-bound mechanism of 1.2: counted n >= n_0(gamma) must have o_t(n) >= rho t - 2"))
    print("n_0(gamma) = 2^(1/(gamma - log2(3/2))) is where (3/2)^t <= n^gamma/2 first holds; rho = gamma/alpha")
    for b in (1, -1):
        for g in ALL_GAMMAS:
            tot, v1, v2 = viol[b][g]
            print("  b=%2d gamma=%.4f rho=%.4f n_0=%8.1f : counted n >= n_0: %9d ; violations of o >= rho t - 2: %d ; of o alpha >= t gamma + gamma - 2: %d"
                  % (b, g, g / ALPHA, 2.0 ** (1.0 / (g - LOG2_32)), tot, v1, v2))
    return results

# ======================================================================================
def good_words(t, o, cap_fn):
    """words of length t with o odd letters whose partial sums are all <= cap (cap = cap_fn(S_t))."""
    out = []
    for pos in itertools.combinations(range(t), o):
        w = [0] * t
        for p in pos:
            w[p] = 1
        S = partial_sums(w)
        cap = cap_fn(S[t])
        if max(S[1:]) <= cap + 1e-9:
            out.append(w)
    return out

def section_E():
    print(hdr("E. lower-bound construction of 1.3 (rotation): representatives of the good words"))
    print("o = ceil(rho t) + 1, rho = max(1/2, gamma/alpha); good = all partial sums <= max(0, S_t);")
    print("representative n in [2^t, 2^(t+1)); failure = some i <= t with T_b^i(n) < n^gamma (exact comparison)")
    mind = min(abs(i - o * ALPHA) for i in range(1, 21) for o in range(0, i + 1))
    print("min_{1<=i<=20, 0<=o<=i} |i - o alpha| = %.6f  (Diophantine margin available at these sizes)" % mind)
    for g in (0.85, 0.91, 0.97, 1.0):
        rho = rho_of(g)
        fg = Fraction(str(g))
        p, q = fg.numerator, fg.denominator
        for t in (8, 12, 16, 20):
            o = ceil_safe(rho * t) + 1
            S_t = t - o * ALPHA
            words = good_words(t, o, lambda st: max(0.0, st))
            need = comb(t, o) / t
            line = "  gamma=%.2f t=%2d o=%2d S_t=%7.3f C(t,o)=%7d good=%7d (>= C/t=%9.1f: %s)" % (
                g, t, o, S_t, comb(t, o), len(words), need, len(words) >= need)
            for b in (1, -1):
                fails = 0
                worst = None
                for w in words:
                    r = terras_rep(w, b)
                    n = r + (1 << t)
                    assert parity_word(n, b, t) == w
                    mn = orbit_min(n, b, t)
                    ok = (mn ** q >= n ** p)
                    if not ok:
                        fails += 1
                        if worst is None:
                            worst = (n, mn)
                line += " | b=%2d fails=%d%s" % (b, fails, "" if worst is None else " e.g. n=%d min=%d" % worst)
            if g < 1:
                line += " | proof needs t(1-gamma) > 2 alpha (t > %.0f) or n >= 3^(1/(1-gamma)) = 2^%.1f" % (2 * ALPHA / (1 - g), math.log2(3) / (1 - g))
            print(line)
    print("note: at gamma = 1 on the minus sheet the rotation argument alone does not PROVE T^i(n) >= n (carry negative,")
    print("      S_i may be within o(1) of 0); the counts above show whether it FAILS at these sizes.")

def section_E2():
    print(hdr("E2. optional strengthening (auditor): tilted cycle lemma at gamma = 1 on the minus sheet"))
    print("steps x_j + alpha/t, o = ceil(t/alpha) + 2 (so S_t <= -2 alpha and the tilted total is <= -alpha < 0); rotate after")
    print("the maximal tilted partial sum; then S'_i <= -i alpha/t for 1 <= i <= t, so 2^(-S'_i) >= 3^(i/t) >= 1 + (i/t) ln 3 and")
    print("T_-^i(n) >= n + n (i/t) ln 3 - (3/2)^i > n for n >= 2^t (margin inequality 2^t (i/t) ln 3 >= (3/2)^i, all 1<=i<=t):")
    print("this would give c X^(h*)/log^(3/2) X on the minus sheet too, with the same C(t,o)/t count.")
    for t in (8, 10, 12, 14, 16, 18, 20):
        o = ceil_safe(t / ALPHA) + 2
        if o > t:
            continue
        c = ALPHA
        good = set()
        for pos in itertools.combinations(range(t), o):
            w = [0] * t
            for p in pos:
                w[p] = 1
            S = partial_sums(w)
            U = [S[i] + i * c / t for i in range(t + 1)]
            p = U.index(max(U))
            w2 = w[p:] + w[:p]
            S2 = partial_sums(w2)
            assert all(S2[i] <= -i * c / t + 1e-9 for i in range(1, t + 1)), (t, w, p)
            good.add(tuple(w2))
        need = comb(t, o) / t
        fails = 0
        for w in good:
            n = terras_rep(list(w), -1) + (1 << t)
            if orbit_min(n, -1, t) < n:
                fails += 1
        marg = all((1 << t) * (i / t) * math.log(3) >= 1.5 ** i for i in range(1, t + 1))
        print("  t=%2d o=%2d C(t,o)=%7d distinct rotated words=%6d (>= C/t=%8.1f: %s) S'_i <= -i alpha/t: asserted; margin inequality: %s; minus-sheet no-dip failures: %d"
              % (t, o, comb(t, o), len(good), need, len(good) >= need, marg, fails))

# ======================================================================================
def section_F():
    print(hdr("F. block construction of 1.3 (gamma = 1, minus sheet)"))
    print("W_L = words of length L, o_L = ceil(L/alpha) + 1 odd letters, all partial sums <= 0;")
    print("c_L(W) = -max over W_L and 1<=i<=L of S_i ; c_L(theory) = min{o alpha - i > 0 : 1<=i<=L, 0<=o<=i};")
    print("n_2(L) = (3/2)^L/(2^c_L - 1) is where n 2^c_L - (3/2)^L >= n starts to hold (first block).")
    WL = {}
    for L in range(3, 15):
        oL = ceil_safe(L / ALPHA) + 1
        words = good_words(L, oL, lambda st: 0.0)
        cW = -max(max(partial_sums(w)[1:]) for w in words) if words else float("nan")
        cth = min(o * ALPHA - i for i in range(1, L + 1) for o in range(0, i + 1) if o * ALPHA - i > 0)
        n2 = 1.5 ** L / (2 ** cW - 1) if words else float("nan")
        WL[L] = (oL, words, cW)
        print("  L=%2d o_L=%2d S_L=%7.3f C(L,o_L)=%5d |W_L|=%5d (>= C/L=%7.2f: %s) c_L(W)=%.4f c_L(theory)=%.4f n_2(L)=2^%.1f h(o_L/L)=%.4f"
              % (L, oL, L - oL * ALPHA, comb(L, oL), len(words), comb(L, oL) / L, len(words) >= comb(L, oL) / L,
                 cW, cth, math.log2(n2) if n2 > 0 else float("nan"), h(oL / L)))
    print("concatenations of j blocks (t = jL), minus sheet, representative n in [2^t, 2^(t+1)):")
    for L, j in ((4, 3), (4, 5), (5, 4), (6, 2), (6, 3), (7, 2), (8, 2), (9, 2), (10, 2)):
        t = L * j
        oL, words, cW = WL[L]
        combos = list(itertools.product(words, repeat=j))
        cap = 300000
        combos = combos[:cap]
        fails = 0
        ps_fail = 0
        for tpl in combos:
            w = [c for blk in tpl for c in blk]
            S = partial_sums(w)
            # partial-sum claims: block 1: <= -c_L ; block r >= 2: <= -(r-1) alpha
            for i in range(1, t + 1):
                r = (i - 1) // L + 1
                bound = -cW if r == 1 else -(r - 1) * ALPHA
                if S[i] > bound + 1e-9:
                    ps_fail += 1
                    break
            n = terras_rep(w, -1) + (1 << t)
            assert parity_word(n, -1, t) == w
            if orbit_min(n, -1, t) < n:
                fails += 1
        n2 = 1.5 ** L / (2 ** cW - 1)
        print("  L=%2d j=%d t=%2d: %6d words; partial-sum claim failures %d; no-dip failures (T_-^i(n) < n) %d; 2^t >= n_2(L): %s"
              % (L, j, t, len(combos), ps_fail, fails, (1 << t) >= n2))
    print("exponent of the block lower bound, per letter: e_L = log2(C(L,o_L)/L)/L, vs h(o_L/L) and h* = %.6f" % h(1 / ALPHA))
    for L in (10, 20, 50, 100, 200, 500, 1000, 5000):
        oL = ceil_safe(L / ALPHA) + 1
        eL = (math.log2(comb(L, oL)) - math.log2(L)) / L
        print("  L=%5d o_L/L=%.5f h(o_L/L)=%.6f e_L=%.6f h*-e_L=%.6f" % (L, oL / L, h(oL / L), eL, h(1 / ALPHA) - eL))

# ======================================================================================
def brute_F(k, theta_frac, b):
    """#F_b(2^k, theta) by direct enumeration (numpy), exact at ties."""
    X = 1 << k
    tn, td = theta_frac.numerator, theta_frac.denominator
    total = 0
    chunk = 1 << 20
    for start in range(1, X + 1, chunk):
        y = np.arange(start, min(start + chunk, X + 1), dtype=np.int64)
        x = y.copy()
        mn = y.copy()
        for _ in range(k):
            par = x & 1
            x = np.where(par == 1, (3 * x + b) >> 1, x >> 1)
            np.minimum(mn, x, out=mn)
        # condition mn >= y 2^(-k theta)  <=>  log2 mn - log2 y + k theta >= 0
        d = np.log2(mn.astype(np.float64)) - np.log2(y.astype(np.float64)) + k * float(theta_frac)
        mask = d >= -1e-9
        amb = np.nonzero((d >= -1e-9) & (d < 1e-9))[0]
        for idx in amb:       # exact: mn^td * 2^(k tn) >= y^td
            mask[idx] = int(mn[idx]) ** td * (1 << (k * tn)) >= int(y[idx]) ** td
        total += int(mask.sum())
    return total

PRIOR_F = {  # collatz_thin_20260925_audit.out (independent audit of THM-4476), theta -> k -> (#F_{+1}, #F_{-1})
    "0.03": {12: (274, 273), 14: (1206, 1202), 16: (3989, 3987), 18: (12339, 12337), 20: (46612, 46611)},
    "0.10": {12: (736, 734), 14: (2774, 2771), 16: (9657, 9654), 18: (38146, 38142), 20: (157652, 157648)},
}

def section_G():
    print(hdr("G. Proposition 2: #F_b(X, theta), X = 2^k, window k, both sheets"))
    if np is None:
        print("numpy not available; skipping")
        return
    # (i) the upper bound is the THM-4476 lemma only if X^(log2(3/2)+theta) is dominated by X^h(rho)
    worst = min(h((1 - th) / ALPHA) - (LOG2_32 + th) for th in [i * (1 - ALPHA / 2) / 1000 for i in range(1, 1000)])
    print("h((1-theta)/alpha) - (log2(3/2) + theta) on a 999-point grid of (0, theta_0): min = %.4f (> 0: the lemma's first term is dominated)" % worst)
    for ths in ("0.03", "0.10"):
        th = Fraction(ths)
        rho = (1 - float(th)) / ALPHA
        print("theta=%s rho=%.5f h(rho)=%.5f" % (ths, rho, h(rho)))
        for k in (12, 14, 16, 18, 20):
            X = 1 << k
            Fp = brute_F(k, th, 1)
            Fm = brute_F(k, th, -1)
            prior = PRIOR_F[ths].get(k)
            # construction of 1.4: words of length k-1, o = ceil(rho(k-1)) + 2, all partial sums <= max(0, S_{k-1}),
            # representative y in (X/2, X]; check T_b^i(y) >= y X^(-theta) for 0 <= i <= k (exact)
            o = ceil_safe(rho * (k - 1)) + 2
            words = good_words(k - 1, o, lambda st: max(0.0, st))
            tn, td = th.numerator, th.denominator
            res = []
            for b in (1, -1):
                fails = 0
                for w in words:
                    r = terras_rep(w, b)
                    y = (1 << k) if r == 0 else (1 << (k - 1)) + r
                    assert (1 << (k - 1)) < y <= X
                    assert parity_word(y, b, k - 1) == w
                    mn = orbit_min(y, b, k)
                    if not (mn ** td * (1 << (k * tn)) >= y ** td):
                        fails += 1
                res.append(fails)
            print("  k=%2d #F_+=%8d #F_-=%8d %s | logF/logX(+)=%.4f | construction: o=%2d C(k-1,o)=%6d good=%6d (>= C/(k-1)=%8.1f: %s) fails(+)=%d fails(-)=%d"
                  % (k, Fp, Fm, ("prior audit %s: %s" % (prior, "match" if prior == (Fp, Fm) else "MISMATCH")) if prior else "",
                     math.log2(Fp) / k, o, comb(k - 1, o), len(words), comb(k - 1, o) / (k - 1), len(words) >= comb(k - 1, o) / (k - 1), res[0], res[1]))

# ======================================================================================
def section_H():
    print(hdr("H. entropy and binomial bounds quoted in 1.1 - 1.3"))
    rhos = sorted(set([g / ALPHA for g in NOTE_GAMMAS if g / ALPHA > 0.5] + [1 / ALPHA, 0.9, 0.95]))
    print("(i) C(t, ceil(rho t)) >= 2^(t h(rho))/(t+1), t <= 2000 (literal 1.1 claim; the exponent is evaluated at rho, not at ceil(rho t)/t)")
    for rho in rhos:
        fails = [t for t in range(1, 2001) if math.log2(comb(t, ceil(rho * t))) < t * h(rho) - math.log2(t + 1) - 1e-12]
        print("   rho=%.5f: failures %s" % (rho, (fails[:10] + ["..."] if len(fails) > 10 else fails) if fails else "none"))
    print("(ii) sum_{o >= rho t - 2} C(t,o) <= (t+1) (rho/(1-rho))^2 2^(t h(rho)) for t in [t_0, 1000], t_0 = min{t : rho t - 2 > t/2}")
    for rho in rhos:
        t0 = ceil(2 / (rho - 0.5)) + 1
        fails = []
        worst = 0.0
        for t in range(t0, 1001):
            omin = ceil(rho * t - 2 - 1e-12)
            s = sum(comb(t, o) for o in range(max(0, omin), t + 1))
            bound = (t + 1) * (rho / (1 - rho)) ** 2 * 2 ** (t * h(rho))
            ratio = s / bound
            worst = max(worst, ratio)
            if ratio > 1 + 1e-12:
                fails.append(t)
        print("   rho=%.5f t_0=%3d: max ratio sum/bound = %.4f ; failures %s" % (rho, t0, worst, fails if fails else "none"))
    print("(iii) ratio C(t, ceil(rho t)+1)/C(t, ceil(rho t)) for t in [20, 2000] (bounded below for rho < 1):")
    for rho in rhos:
        rs = [comb(t, ceil(rho * t) + 1) / comb(t, ceil(rho * t)) for t in range(20, 2001)]
        print("   rho=%.5f: min %.4f max %.4f limit (1-rho)/rho = %.4f" % (rho, min(rs), max(rs), (1 - rho) / rho))
    print("(iv) log-power bookkeeping: with C(t,o) >= c 2^(t h)/sqrt t the top block gives D >= c X^h/log^(3/2) X (note claims /log^3: valid, weaker);")
    print("     the upper bound sum_{t <= log2 X} (t+1) 2^(t h) <= C X^h log X (note claims log^2: valid, weaker).")
    rho = 1 / ALPHA
    for t in (100, 1000, 5000):
        o = ceil(rho * t) + 1
        print("   t=%5d: log2(C(t,o)/t) - t h(rho) = %.3f ; -(3/2) log2 t = %.3f ; -3 log2 t = %.3f" % (t, math.log2(comb(t, o)) - math.log2(t) - t * h(rho), -1.5 * math.log2(t), -3 * math.log2(t)))

# ======================================================================================
def main():
    max_t = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    print("collatz_dipspectrum_20260926_audit.py  max_t=%d  numpy=%s mpmath=%s" % (max_t, np is not None, HAVE_MP))
    section_A()
    section_B()
    section_CD(max_t)
    section_E()
    section_E2()
    section_F()
    section_G()
    section_H()
    print("\ndone.")

if __name__ == "__main__":
    main()
