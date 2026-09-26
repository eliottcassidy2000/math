#!/usr/bin/env python3
"""collatz_landing_20260926_audit.py -- adversarial audit of
  05-knowledge/results/collatz_landing_20260926_multiplicity_reassessment.md (sections 0, 1, 3, 7, 8) and
  05-knowledge/results/collatz_thin_20260926_little_o_thin_divergence.md section 3b (discrepancy corollary),
with independent re-derivations of the probe scripts
  collatz_landing_20260926_probe.py / _probe2.py / collatz_zeckendorf_20260926_interlock_probe.py.

Every partial-sum comparison is exact (integer comparison of 3^o 2^W against 2^j; 2^(-3.6) = 2^(-18/5) is
handled by fifth powers; 2^(-1.6) = 2^(-8/5) likewise). The output is written by the caller to
05-knowledge/results/collatz_landing_20260926_audit.out.

Checks:
  [1] section 0: a*(beta) = beta lambda*/h* - 3/2 (bootstrap with c_1 h* = beta + eta), arithmetic.
  [2] section 1: hover-then-drop and climb-then-drop residue classes (Terras), exact landing multiplicities at
      L = 24, theta = 0.15 (theta L = 3.6); exhaustive over all [-2,0]-hover words of length 12; maximum
      multiplicity over hover words of length 18-20 in [-3,0]/[-4,0]; the 1-bit-band bound (<= 2k/3 + O(1));
      single-step factor; distance of a dipper to its landing point; reproduction of probe (B) and probe2.
  [3] section 3: strip entropies h_W by an independent exact DP (m = 150 -> 300 and 300 -> 600), W = 1,
      thresholds (1-h*)/(1-h_W).
  [4] section 8: carry sign and bound (exhaustive over words of length <= 12), the peak ratio bound, and the
      leader/peak word checks on long segments (2^40-1, random 40-bit odd, 5n+1 orbit of 7 to 2^60) at
      L = 30, 40.
  [5] section 3b: Proposition 6 identity (exact rationals) on real orbits, sign convention, constants.
  [6] provenance: sha256 (raw bytes and LF-normalised) of the probe scripts/outputs and notes; THM-4499's
      recorded hashes; statistical significance of the interlock probe's numbers; hand check of the mutual
      information on n = 2..13.
Usage: python3 collatz_landing_20260926_audit.py
"""
import hashlib
import importlib.util
import math
import os
import random
import subprocess
import sys
from collections import Counter
from fractions import Fraction

sys.stdout.reconfigure(newline='\n')

ALPHA = math.log2(3)
RHO = math.log(2) / math.log(3)
H = -(RHO * math.log2(RHO) + (1 - RHO) * math.log2(1 - RHO))
LAM = math.log2(RHO / (1 - RHO)) / ALPHA
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

FAIL = []
NCHECK = [0]
REFUTED = []
NCLAIM = [0]


def check(name, ok, detail=''):
    """internal consistency / reproduction check: a FAIL here would be a bug in this audit or in the probes."""
    NCHECK[0] += 1
    tag = 'PASS' if ok else 'FAIL'
    if not ok:
        FAIL.append(name)
    print('  [%s] %s%s' % (tag, name, ('  -- ' + detail) if detail else ''))


def claim(name, holds, detail=''):
    """a statement made in the audited notes: CLAIM OK / CLAIM FALSE."""
    NCLAIM[0] += 1
    tag = 'CLAIM OK' if holds else 'CLAIM FALSE'
    if not holds:
        REFUTED.append(name)
    print('  [%s] %s%s' % (tag, name, ('  -- ' + detail) if detail else ''))


def T(n, q=3, b=1):
    return n // 2 if n % 2 == 0 else (q * n + b) // 2


def orbit_until(n, q=3, b=1, stop_at_one=True, max_steps=10 ** 6, bound=None):
    """orbit with distinct terms: stops before the first repeated term (the 1,2,1,2 tail) or above `bound`."""
    ys = [n]
    seen = {n}
    for _ in range(max_steps):
        n = T(n, q, b)
        if n in seen:
            break
        if bound is not None and n > bound:
            break
        ys.append(n)
        seen.add(n)
    return ys


def terras_residue(word, q=3, b=1):
    """least r >= 0 with parity word `word` under T_(q,b) (bit-by-bit; T^j(r) mod 2 depends on r mod 2^(j+1))."""
    r = 0
    for j in range(len(word)):
        x = r
        for _ in range(j):
            x = T(x, q, b)
        if x % 2 != word[j]:
            r += 1 << j
    return r


def parity_word(n, k, q=3, b=1):
    w = []
    for _ in range(k):
        w.append(n % 2)
        n = T(n, q, b)
    return w


def landing_exact(ys, k, tL_num, tL_den, X):
    """dippers/landing points for window k, threshold factor 2^(-tL) with tL = tL_num/tL_den (exact: compare
    y_(i+s)^den 2^num < y_i^den), X = size cut (None: no cut). Returns dict landing -> list of dippers, #ND, total."""
    mult = {}
    nd = 0
    total = 0
    for i, y in enumerate(ys):
        if (X is not None and y > X) or i + k >= len(ys):
            continue
        total += 1
        rhs = y ** tL_den
        land = None
        for s in range(1, k + 1):
            if (ys[i + s] ** tL_den) << tL_num < rhs:
                land = i + s
                break
        if land is None:
            nd += 1
        else:
            mult.setdefault(land, []).append(i)
    return mult, nd, total


def partial_sums_ok_band(word, W):
    """all S_j = o_j alpha - j in [-W, 0] for 1 <= j <= |word| (exact)."""
    o = 0
    for j, w in enumerate(word, 1):
        o += w
        if not (3 ** o <= 2 ** j and (3 ** o) << W >= 2 ** j):
            return False
    return True


def band_words(m, W):
    """all words of length m with all partial sums in [-W, 0] (exact DFS)."""
    out = []
    def rec(prefix, o, j):
        if j == m:
            out.append(list(prefix)); return
        for w in (0, 1):
            o2 = o + w
            if 3 ** o2 <= 2 ** (j + 1) and (3 ** o2) << W >= 2 ** (j + 1):
                prefix.append(w); rec(prefix, o2, j + 1); prefix.pop()
    rec([], 0, 0)
    return out


def S(o, j):
    return o * ALPHA - j


# ----------------------------------------------------------------------------------------------------
print('=' * 100)
print('AUDIT collatz_landing_20260926 (reassessment note sections 0/1/3/7/8; THM-4499 note section 3b)')
print('h* = %.7f  lambda* = %.6f  lambda*/h* = %.6f  1/h* = %.6f' % (H, LAM, LAM / H, 1 / H))
print('=' * 100)

# [1] --------------------------------------------------------------------------------------------------
print('\n[1] Section 0: a*(beta) = beta lambda*/h* - 3/2')
for beta, claimed in ((1.0, -0.986), (0.5, -1.243), (0.0, -1.5)):
    a = beta * LAM / H - 1.5
    check('a*(%.1f) = %.6f (note: %.3f)' % (beta, a, claimed), abs(a - claimed) < 6e-4)
# bookkeeping of the modified bootstrap: (D) term M(L) N(X L^(-c1)) with M(L) = L^beta needs beta - c1 h* = -eta
for beta in (1.0, 0.5, 0.0):
    eta = 0.01
    c1 = (beta + eta) / H
    s = min(0.1, eta * math.log(2) / c1)               # THM-4499's choice: e^(s theta_X L) = L^(s c1/ln 2) <= L^eta
    nd_exp = LAM * c1 - 1.5 + s * c1 / math.log(2)   # exponent of L in the (ND) term at theta_X = c1 log2 L / L
    check('beta=%.1f: L^beta * L^(-c1 h*) = L^(-eta) and (ND) exponent <= a*(beta) + eta(lambda*/h* + 1)'
          % beta, abs((beta - c1 * H) + eta) < 1e-12 and nd_exp <= beta * LAM / H - 1.5 + eta * (LAM / H + 1) + 1e-12,
          'ND exponent %.6f, a*(beta) + 1.514 eta = %.6f' % (nd_exp, beta * LAM / H - 1.5 + eta * (LAM / H + 1)))
print('  Verdict: the only change to THM-4499 section 3 is c_1 h* = beta + eta (theta_X = c_1 log_2 L/L); the (D) term')
print('  L^beta N(X L^(-c_1)) <= 2^|a| K X^(h*) L^(a-eta), the (ND) term X^(h*) L^(lambda* c_1 - 3/2 + eta); a*(beta) as stated.')
print('  HYP-9161\'s average form suffices: #(D) = sum of multiplicities <= C L^beta * #landing points <= C L^beta N(X^(1-theta)).')
claim('section 0: a*(beta) = beta lambda*/h* - 3/2; a*(1) = -0.986, a*(1/2) = -1.243, a*(0) = -3/2; bootstrap changed only by c_1 h* = beta + eta', True)

# [2] --------------------------------------------------------------------------------------------------
print('\n[2] Section 1: residue classes and the landing multiplicity (L = 24, theta = 0.15, theta L = 3.6 = 18/5)')
L = 24; X = 1 << L; NUM, DEN = 18, 5      # 2^(-3.6) = 2^(-18/5): y' < y 2^(-3.6) iff y'^5 2^18 < y^5
# single-step factor
mn, mx = 10.0, 0.0
for y in range(1, 1 << 16):
    r = T(y) / y
    mn, mx = min(mn, r), max(mx, r)
check('single step: T(y)/y in [1/2, 2] for 1 <= y < 2^16 (min %.4f, max %.4f; odd: 3/2 + 1/(2y))' % (mn, mx), mn == 0.5 and mx == 2.0)
print('  hence y_(i+s) >= y_i 2^(-s); a dipper (y_(i+s) < y_i 2^(-theta L)) has s > theta L: distance >= floor(theta L) + 1 = %d here.' % (int(3.6) + 1))

hover_words = band_words(12, 2)
print('  words of length 12 with all partial sums in [-2, 0]: %d' % len(hover_words))
best = (0, None)
spread_stats = []
for hw in hover_words:
    for D in range(4, 13):
        word = hw + [0] * D
        k = len(word)
        r = terras_residue(word)
        assert parity_word(r, k) == word
        # largest representative <= X (keeps carries negligible, all hover values <= y_0 <= X)
        y0 = r + ((X - r) >> k << k)
        if y0 > X:
            y0 -= 1 << k
        assert 0 < y0 <= X and y0 < 1 << (2 * L + 8)
        ys = orbit_until(y0, max_steps=3 * L, stop_at_one=False)
        if len(ys) < 2 * L + 1:
            ys = [y0]; n = y0
            for _ in range(2 * L):
                n = T(n); ys.append(n)
        mult, nd, total = landing_exact(ys, L, NUM, DEN, X)
        # where do the 12 hover indices land?
        where = {}
        for lp, dips in mult.items():
            for i in dips:
                if i < 12:
                    where[i] = lp
        lps = sorted(set(where.values()))
        m_max = max((len(v) for v in mult.values()), default=0)
        spread_stats.append((len(lps), m_max))
        if m_max > best[0]:
            best = (m_max, (hw, D, y0, {lp: sorted(v) for lp, v in mult.items()}, where))
m_max, (hw, D, y0, mult, where) = best
print('  best over all [-2,0]-hover words (m = 12) and drop lengths D = 4..12: max multiplicity of one landing point = %d' % m_max)
print('    word = %s + 0^%d, y_0 = %d (< 2^%d), landing points of the hover indices: %s'
      % (''.join(map(str, hw)), D, y0, 2 * L + 8, {lp: [i for i in range(12) if where.get(i) == lp] for lp in sorted(set(where.values()))}))
print('    hover indices are split over %d..%d consecutive landing points (min..max over the words); max multiplicity per word: %s'
      % (min(s[0] for s in spread_stats), max(s[0] for s in spread_stats), sorted(set(s[1] for s in spread_stats))))
claim('section 1, hover-then-drop (W = 2, m = 12): "that landing point has multiplicity m + O(theta L)", i.e. some landing point receives >= m = 12 dippers', m_max >= 12,
      'the maximum is %d (7 hover indices + 2 drop indices); the m hover indices are shared by the halvings of the drop, one 1-bit sub-band each' % m_max)

# climb then drop: word 1^12 0^D, tested at L = 34 (so that the whole climb is <= X) with the same theta L = 3.6
L2 = 34; X2 = 1 << L2
for D in (11, 12):
    word = [1] * 12 + [0] * D
    r = terras_residue(word)
    assert parity_word(r, len(word)) == word
    ys = [r]; n = r
    for _ in range(2 * L2):
        n = T(n); ys.append(n)
    assert max(ys[:13]) <= X2
    mult, nd, total = landing_exact(ys, L2, NUM, DEN, X2)
    where = {i: lp for lp, dips in mult.items() for i in dips if i < 12}
    counts = Counter(where.values())
    print('  climb-then-drop 1^12 0^%d at L = %d, theta L = 3.6: y_0 = %d; climb indices -> landing points %s'
          % (D, L2, r, dict(sorted(counts.items()))))
    claim('section 1, climb-then-drop 1^12 0^%d: "every index of the climb lands at the same point"' % D, len(counts) == 1,
          'the 12 climb indices land at %d distinct points, at most %d per point (thresholds spaced 0.585 bits, the drop crosses 1 bit per step)'
          % (len(counts), max(counts.values())))

# the true maximum at L = 24 over hover words of length 18-20 (bands [-3,0], [-4,0]) followed by halvings
print('  maximum landing multiplicity at L = 24, theta L = 3.6, over words = hover(m, W) + 0^(24-m):')
gbest = (0, None)
for (m, W) in ((20, 3), (18, 4), (16, 4), (20, 4)):
    words = band_words(m, W)
    bm = (0, None)
    for hw in words:
        word = hw + [0] * (L - m)
        r = terras_residue(word)
        y0 = r + ((X - r) >> L << L)
        if y0 > X:
            y0 -= 1 << L
        ys = [y0]; n = y0
        for _ in range(2 * L):
            n = T(n); ys.append(n)
        mult, nd, total = landing_exact(ys, L, NUM, DEN, X)
        mm = max((len(v) for v in mult.values()), default=0)
        if mm > bm[0]:
            bm = (mm, (hw, y0, {lp: sorted(v) for lp, v in mult.items() if len(v) == mm}))
    print('    hover length %2d in [-%d,0] (%6d words): max multiplicity %2d = %.2f L; best word %s, y_0 = %d, landing %s'
          % (m, W, len(words), bm[0], bm[0] / L, ''.join(map(str, bm[1][0])), bm[1][1], bm[1][2]))
    if bm[0] > gbest[0]:
        gbest = bm
# the 1-bit band bound: dippers of j have y_i in (2^(tL) y_j, 2^(tL+1) y_j] (half-open: three values a < b < c with c < 2a is impossible
# for consecutive orbit terms, since a halving leaves any (x, 2x] and two odd steps multiply by (9a+5)/(4a) > 2)
viol = 0
for y0 in list(range(1, 3000)) + [gbest[1][1], 2 ** 40 - 1]:
    ys = [y0]; n = y0
    for _ in range(60):
        n = T(n); ys.append(n)
    for i in range(len(ys) - 2):
        a, b, c = ys[i], ys[i + 1], ys[i + 2]
        if max(a, b, c) < 2 * min(a, b, c):
            viol += 1
check('no three consecutive orbit values lie in one half-open 1-bit band (x, 2x] (3002 starts, 60 steps each)', viol == 0)
# and the band statement itself on the searched optimum: every dipper of the best landing point has y_i in (2^3.6 y_j, 2^4.6 y_j]
lp, dips = next(iter(gbest[1][2].items()))
ys = [gbest[1][1]]; n = ys[0]
for _ in range(2 * L):
    n = T(n); ys.append(n)
check('all %d dippers of the best landing point j=%d lie in the 1-bit band (2^3.6 y_j, 2^4.6 y_j] and y_j = y_(j-1)/2' % (len(dips), lp),
      ys[lp] * 2 == ys[lp - 1] and all((ys[i] ** 5 > (ys[lp] ** 5) << 18) and (ys[i] ** 5 <= (ys[lp] ** 5) << 23) for i in dips))
print('  => multiplicity of a landing point j <= #{i in [j-k, j): y_i in (2^(theta L) y_j, 2^(theta L + 1) y_j]} <= 2 ceil(k/3) <= 2k/3 + 4/3')
print('     (y_j < y_(j-1) forces y_j = y_(j-1)/2, so y_i 2^(-theta L) in (y_j, 2 y_j]; at most two of three consecutive indices are in the band).')
print('     So "multiplicity ~ L" holds with a constant in [%.2f, 0.67]: the observed maximum at L = 24 is %.2f L; THM-4476\'s k is improvable to 2k/3 + O(1), no more.' % (gbest[0] / L, gbest[0] / L))
claim('section 1: "multiplicity ~ L is realised by integers / the pigeonhole is tight on segments" (qualitatively, up to the constant 2/3)', True,
      'max %d = %.2f L at L = 24 from a 20-step hover in [-3, 0]; upper bound 2k/3 + 4/3' % (gbest[0], gbest[0] / L))
claim('section 1: a single step multiplies by a factor in [1/2, 3/2 + o(1)] (b = 1: 3/2 + 1/(2y) <= 2), so a dipper at distance s has s > theta L', True)

# reproduction of probe (B) (L = 20, theta = 0.05, 0.1, 0.2; theta L = 1, 2, 4 bits, exact)
print('  reproduction of probe (B) [collatz_landing_20260926_probe.out], L = 20 (independent code, exact thresholds):')


def orbit_padded(n, steps, q=3, b=1):
    ys = [n]
    for _ in range(steps):
        n = T(n, q, b); ys.append(n)
    return ys


def synthetic_hover_as_probe(L, K):
    word = []; Ssum = 0.0
    for j in range(K):
        if Ssum + (ALPHA - 1) <= 0:
            word.append(1); Ssum += ALPHA - 1
        else:
            word.append(0); Ssum -= 1
    word += [0] * (L // 4) + [1] * (L // 4)
    r = terras_residue(word)
    return r + (1 << len(word)), word


syn, synword = synthetic_hover_as_probe(20, 15)
o = 0; smin = 0.0
for j, w in enumerate(synword[:15], 1):
    o += w; smin = min(smin, S(o, j))
print('    probe.py synthetic_hover(20, 15): docstring says partial sums in [-1, 0]; actual minimum partial sum of the 15-step hover = %.3f (band width %.2f bits)' % (smin, -smin))
probeB_expected = {
    (0.05, '5n+1'): (53, 35, 18, 8, [5, 4, 3, 2, 1, 1, 1, 1]), (0.05, '27'): (181, 140, 41, 19, [5, 5, 3, 3, 3, 3, 2, 2]),
    (0.05, '2^18-1'): (299, 241, 58, 27, [7, 6, 4, 4, 3, 3, 3, 3]), (0.05, '2^19-1'): (528, 468, 60, 30, [5, 5, 4, 4, 4, 3, 3, 3]),
    (0.05, 'syn'): (120, 45, 75, 37, [5, 5, 5, 5, 5, 3, 3, 3]),
    (0.1, '5n+1'): (53, 46, 7, 3, [4, 2, 1]), (0.1, '27'): (181, 147, 34, 15, [5, 4, 4, 3, 3, 3, 2, 2]),
    (0.1, '2^18-1'): (299, 242, 57, 22, [7, 6, 5, 5, 5, 3, 3, 2]), (0.1, '2^19-1'): (528, 469, 59, 25, [5, 4, 4, 4, 4, 4, 4, 3]),
    (0.1, 'syn'): (120, 56, 64, 28, [8, 5, 5, 5, 4, 4, 3, 3]),
    (0.2, '5n+1'): (53, 53, 0, 0, []), (0.2, '27'): (181, 152, 29, 11, [5, 5, 5, 4, 2, 2, 2, 1]),
    (0.2, '2^18-1'): (299, 249, 50, 18, [7, 7, 6, 4, 3, 3, 3, 3]), (0.2, '2^19-1'): (528, 477, 51, 19, [5, 5, 4, 4, 4, 4, 4, 4]),
    (0.2, 'syn'): (120, 64, 56, 22, [7, 5, 5, 5, 5, 5, 4, 2]),
}
cases = [('5n+1', orbit_padded(7, 2000, q=5)), ('27', orbit_padded(27, 200)), ('2^18-1', orbit_padded(2 ** 18 - 1, 400)),
         ('2^19-1', orbit_padded(2 ** 19 - 1, 600)), ('syn', orbit_padded(syn, 200))]
allok = True
for theta in (0.05, 0.1, 0.2):
    tl = int(round(theta * 20))
    for name, ys in cases:
        mult, nd, total = landing_exact(ys, 20, tl, 1, 1 << 20)
        ms = sorted((len(v) for v in mult.values()), reverse=True)
        got = (total, nd, sum(ms), len(ms), ms[:8])
        exp = probeB_expected[(theta, name)]
        ok = got == exp
        allok &= ok
        if not ok:
            print('    MISMATCH theta=%.2f %s: got %s expected %s' % (theta, name, got, exp))
check('probe (B): all 15 rows (total, ND, D, landing, top-8 multiplicities) reproduced exactly', allok)

# reproduction of probe2 (same seed and call order; their float threshold replicated)
print('  reproduction of probe2 [collatz_landing_20260926_probe2.out] (seed 20260926, same call order, float threshold as in the script):')


def orbit_p2(n, steps, q=3):
    ys = [n]
    for _ in range(steps):
        n = T(n, q, 1); ys.append(n)
        if n == 1:
            break
    return ys


def landing_float(ys, L, theta):
    X = 2 ** L; thr = 2.0 ** (-theta * L); mult = {}; nd = 0; total = 0
    for i, y in enumerate(ys):
        if y > X or i + L >= len(ys):
            continue
        total += 1; land = None
        for s in range(1, L + 1):
            if ys[i + s] < y * thr:
                land = i + s; break
        if land is None:
            nd += 1
        else:
            mult.setdefault(land, []).append(i)
    return mult, nd, total


p2_expected = """20 43 6 37 18 4|20 51 30 21 7 6|20 114 61 53 18 6|20 121 74 47 18 5|20 53 53 0 0 0|20 94 77 17 9 3|30 134 58 76 26 7|30 97 26 71 23 9|30 117 55 62 27 9|30 115 54 61 28 7|30 130 125 5 3 2|30 80 75 5 2 4|40 184 47 137 31 11|40 70 0 70 26 9|40 137 28 109 37 11|40 137 25 112 33 8|40 181 177 4 3 2|40 100 96 4 3 2|60 287 52 235 65 16|60 303 37 266 66 21|60 163 0 163 52 12|60 326 43 283 65 14|60 508 443 65 20 7|60 230 198 32 14 4|80 294 1 293 85 18|80 303 1 302 93 11|80 427 0 427 90 16|80 353 4 349 86 18|80 563 468 95 25 11|80 191 188 3 1 3""".split('|')
random.seed(20260926)
rows = []
climb_info = []
for Lp in (20, 30, 40, 60, 80):
    theta = 1.05 * math.log2(Lp) / Lp
    segs = [('2^L-1', orbit_p2(2 ** Lp - 1, 40 * Lp)), ('2^(L-1)+1', orbit_p2(2 ** (Lp - 1) + 1, 40 * Lp)),
            ('rnd~2^L', orbit_p2(random.randrange(2 ** (Lp - 1), 2 ** Lp) | 1, 40 * Lp)),
            ('rnd~2^(L-3)', orbit_p2(random.randrange(2 ** (Lp - 4), 2 ** (Lp - 3)) | 1, 40 * Lp)),
            ('5n+1 from 7', orbit_p2(7, 40 * Lp, q=5)),
            ('5n+1 rnd', orbit_p2(random.randrange(2 ** (Lp // 2 - 1), 2 ** (Lp // 2)) | 1, 40 * Lp, q=5))]
    for name, ys in segs:
        mult, nd, total = landing_float(ys, Lp, theta)
        ms = [len(v) for v in mult.values()]
        rows.append('%d %d %d %d %d %d' % (Lp, total, nd, sum(ms), len(ms), max(ms) if ms else 0))
        if name in ('2^L-1', '2^(L-1)+1') and Lp in (40, 60) and ms:
            lp = max(mult, key=lambda j: len(mult[j]))
            dips = sorted(mult[lp])
            first_below = next(i for i, y in enumerate(ys) if y <= 2 ** Lp)
            climb_info.append((Lp, name, lp, dips, first_below, ys[lp], [round(math.log2(ys[i] / ys[lp]), 2) for i in dips]))
check('probe2: all 30 rows (total, ND, D, landing, maxmult) reproduced exactly', rows == p2_expected,
      '' if rows == p2_expected else 'first mismatch: %s' % next((a, b) for a, b in zip(rows, p2_expected) if a != b))
print('  where is the maximal-multiplicity landing point on the "climb" starts? (note attributes max ~0.25 L to the built-in climb)')
for Lp, name, lp, dips, fb, ylp, logs in climb_info:
    print('    L=%d %-10s: first index with y <= X is %d (the climb of 2^L-1 leaves [1, X] at step 1); max landing point j=%d, its %d dippers at indices %s,'
          % (Lp, name, fb, lp, len(dips), dips))
    print('               log_2(y_i / y_j) of the dippers: %s (all inside one 1-bit band (theta L, theta L + 1] = (%.2f, %.2f], as the band bound says)'
          % (logs, 1.05 * math.log2(Lp), 1.05 * math.log2(Lp) + 1))
    claim('probe2 reading, L=%d %s: "the climbs built into the starts produce one landing point with ~0.25 L dippers" (its dippers would have to be the initial climb)' % (Lp, name),
          min(dips) <= 1, 'the maximal landing point j=%d is far from the start; its dippers (indices >= %d) are a later hover in one 1-bit band; random starts show the same max/L' % (lp, min(dips)))

# [3] --------------------------------------------------------------------------------------------------
print('\n[3] Section 3: strip entropies by an independent exact DP')
pow3 = [3 ** i for i in range(700)]
pow2 = [2 ** i for i in range(700)]


def strip_count(m, W):
    cur = {0: 1}
    for j in range(1, m + 1):
        nxt = {}
        for o, c in cur.items():
            for st in (0, 1):
                o2 = o + st
                if pow3[o2] <= pow2[j] and (pow3[o2] << W) >= pow2[j]:
                    nxt[o2] = nxt.get(o2, 0) + c
        cur = nxt
    return sum(cur.values())


print('  W = 1 counts for m = 1..6: %s' % [strip_count(m, 1) for m in range(1, 7)])
check('W = 1 admits no word of length >= 3 (steps -1 and +0.585 span 1.585 > 1)', all(strip_count(m, 1) == 0 for m in range(3, 8)))
table = {2: 0.2600, 3: 0.6165, 4: 0.7520, 6: 0.8542, 8: 0.8933, 12: 0.9244, 16: 0.9347}
c300 = {2: 79.00, 3: 184.87, 4: 224.88, 6: 256.04, 8: 268.03, 12: 277.81, 16: 282.08}
print('  W   h_W(150->300)  note   log2 count(300)  note   | h_W(300->600)  (1-h*)/(1-h_W)  note threshold')
thr_note = {4: 0.20, 8: 0.47, 16: 0.77}
allok = True
for W in (2, 3, 4, 6, 8, 12, 16):
    c1, c2, c3 = strip_count(150, W), strip_count(300, W), strip_count(600, W)
    hW = (math.log2(c2) - math.log2(c1)) / 150
    hW2 = (math.log2(c3) - math.log2(c2)) / 300
    thr = (1 - H) / (1 - hW)
    ok = abs(hW - table[W]) < 6e-5 and abs(math.log2(c2) - c300[W]) < 6e-3
    allok &= ok
    print('  %2d   %.4f        %.4f  %8.2f        %8.2f | %.4f         %.3f          %s' % (W, hW, table[W], math.log2(c2), c300[W], hW2, thr, thr_note.get(W, '')))
check('h_W table (150 -> 300) and the log2 counts at m = 300 reproduced to 4 decimals', allok)
check('thresholds (1-h*)/(1-h_W): 0.20 (W=4), 0.47 (W=8), 0.77 (W=16)',
      all(abs((1 - H) / (1 - table[W]) - thr_note[W]) < 6e-3 for W in (4, 8, 16)))
drift = max(abs((math.log2(strip_count(600, W)) - math.log2(strip_count(300, W))) / 300 - table[W]) for W in table)
claim('section 3: h_W table (exact DP, tie rule S_j = -j in [-W,0] iff j <= W), W = 1 empty for m >= 3, thresholds 0.20/0.47/0.77 L', True,
      'the 4-decimal values are 150->300 growth rates; the 300->600 rates differ by up to %.4f, so "PROVED numerically" means +-0.003' % drift)
print('  the "one-sided 0.94996" entry is the limit h*, not a 150->300 growth rate (which would be h* - 1.5*log2(2)/150 = %.4f): cosmetic.' % (H - 1.5 / 150))
print('  consequence check: #words with a [-W,0]-hover of length m <= L is 2^(h_W m + o(m)); each is one class mod 2^m with <= X/2^m + 1 elements in [1, X];')
print('  so the count is <= X 2^(-(1-h_W) m + o(m)) for m <= L, and it is <= X^(h* + o(1)) once m >= (1-h*)/(1-h_W) L. Elementary; the exponent o(1) is omitted in the note.')

# [4] --------------------------------------------------------------------------------------------------
print('\n[4] Section 8 (leaders and peaks)')
# carries: exhaustive, words of length <= 12, b in {1, -1, 3, -5}: beta has the sign of b and |beta_s| <= |b|((3/2)^s - 1)
ok = True
for b in (1, -1, 3, -5):
    for k in range(1, 13):
        for wbits in range(1 << k):
            beta = Fraction(0); M = Fraction(1)
            for j in range(k):
                if (wbits >> j) & 1:
                    beta = (3 * beta + b) / 2; M *= Fraction(3, 2)
                else:
                    beta = beta / 2; M /= 2
            if not (beta * b >= 0 and abs(beta) <= abs(b) * (Fraction(3, 2) ** k - 1)):
                ok = False
check('carry beta_s has the sign of b and |beta_s| <= |b|((3/2)^s - 1) (exhaustive, s <= 12, b in {1,-1,3,-5})', ok)
print('  (i) leaders: y_(i+s) = M_s y_i + beta_s > y_i with |beta_s| <= |b|(3/2)^k <= Y_0/2 <= y_i/2 gives M_s > 1 - |beta_s|/y_i >= 1/2: direction right.')
claim('section 8 (i): leaders, M_s > 1 - |beta_s|/y_i >= 1/2, all partial sums > -1, counted by M_k(1)', True)
print('  (ii) peaks. Version at commit 3987cf42b (the text this audit was commissioned on): "M_k > M_s (1 - o(1)), i.e. the word ends within')
print('       0.01 bits of its strict maximum": NOT what y_i >= Y_0 gives (see the witness below); the count then used M_k(1.6), inconsistent with 0.01.')
print('       Version at commit ab9ef18b4 (revised during this audit, now in the worktree): "M_k > M_s - (|beta_k| + |beta_s|)/z >= M_s - 1" and')
print('       "(s = 0) M_k z + beta_k > z, so M_k > 1/2", then the two cases M_s >= 3/2 / M_s < 3/2 give M_k/M_s > 1/3. The step ">= M_s - 1"')
print('       needs z = y_(i-k) >= |beta_k| + |beta_s|, and "M_k > 1/2" by that route needs z >= 2|beta_k|; the hypothesis is y_i >= Y_0, which')
print('       bounds nothing about z (z < y_i is all that is known). Neither inequality is derived: a gap in the proof as written.')
print('       Correct argument needing only y_i >= Y_0: M_k z = y_i - beta_k >= y_i - Y_0/2 >= y_i/2 and M_s z = y_(i-k+s) - beta_s < y_i + Y_0/2 <= (3/2) y_i,')
print('       so M_k/M_s > 1/3 = 0.3333 > 2^(-1.6) = %.4f, i.e. S_k - S_s > -log_2 3 = -1.585 > -1.6 for every 0 <= s < k (s = 0 gives M_k > 1/2).' % 2 ** -1.6)
print('       With the carry sign: b > 0 gives M_s z <= y_(i-k+s) < y_i and M_k/M_s > 1/2 (S_k - S_s > -1); b < 0 gives the ratio > 2/3.')
print('       z cancels; no condition on z or on the intermediate values is needed. The count by M_k(1.6) stands (1.6 is the right constant: 1/3 > 2^-1.6).')
print('       (Alternative repair of the z-route: peaks with y_(i-k) < Y_0 number at most Y_0 since i -> y_(i-k) is injective; assume z >= Y_0 for the rest.)')
# is the revised chain at least true? search all peaks at index k from starts z <= 2^15, k = 8..16, with y_k >= Y_0 (b = 1)
tot = small_z = violM = 0
for k in range(8, 17):
    Y0k = 2 * Fraction(3, 2) ** k
    for z in range(1, 1 << 15):
        ys = [z]; n = z
        for _ in range(k):
            n = T(n); ys.append(n)
        if len(set(ys)) < len(ys) or ys[k] < Y0k or ys[k] <= max(ys[:k]):
            continue
        tot += 1
        beta = [Fraction(0)]; M = [Fraction(1)]
        for j in range(k):
            if ys[j] % 2:
                beta.append((3 * beta[-1] + 1) / 2); M.append(M[-1] * Fraction(3, 2))
            else:
                beta.append(beta[-1] / 2); M.append(M[-1] / 2)
        if any(z < abs(beta[k]) + abs(beta[s]) for s in range(k)):
            small_z += 1
        if any(M[k] <= M[s] - 1 for s in range(k)):
            violM += 1
print('       search: %d peaks at index k (k = 8..16, z < 2^15, y_k >= Y_0): peaks with z < |beta_k| + |beta_s| for some s: %d; violations of M_k > M_s - 1: %d'
      % (tot, small_z, violM))
print('       (so the revised chain is not refuted by examples; it is unproved as written. Note for b = 1 it can be proved differently: a violation')
print('       M_k <= M_s - 1 with y_k >= Y_0 forces z > (3/2)^(k-s) while the peak condition forces the suffix carry gamma > z, but gamma <= (3/2)^(k-s) - 1.)')
claim('section 8 (ii): the displayed justification of M_k/M_s > 1/3 is complete as written (either version)', False,
      '3987cf42b: "(1 - o(1))/0.01 bits" false at the threshold; ab9ef18b4: divides by z without a lower bound on z')
claim('section 8: Proposition statement (leaders and peaks <= X are O(X^(h*) L^(-3/2)) + 2|b| X^0.585 + k)', True, 'with the y_i-comparison closing the peak case')
# a numerical witness that the ratio can be far from 1 at the threshold: all-odd word, b = 1, y_i = Y_0
k = 20
Y0 = 2 * Fraction(3, 2) ** k
beta = Fraction(0); M = Fraction(1); betas = [beta]; Ms = [M]
for j in range(k):
    beta = (3 * beta + 1) / 2; M *= Fraction(3, 2); betas.append(beta); Ms.append(M)
# y_i = M_k z + beta_k = Y0 -> z = (Y0 - beta_k)/M_k ; y_(i-1) = M_(k-1) z + beta_(k-1)
z = (Y0 - betas[k]) / Ms[k]
ratio = float((Y0 - betas[k]) / (Ms[k - 1] * z + betas[k - 1] - betas[k - 1]))   # = M_k/M_(k-1)... trivially 3/2 here
ratio2 = float((Y0 - betas[k]) / (Y0))     # (y_i - beta_k)/y_i at the threshold
print('       witness: all-odd word, b = 1, y_i = Y_0 = 2(3/2)^k: (y_i - beta_k)/y_i = %.4f (= 1/2 + o(1)), so M_k z is only y_i/2: "1 - o(1)" fails at the stated threshold.' % ratio2)

# (iv) numerical verification on long segments. A convergent 3n+1 segment ending at 1 has a single leader (its last term), so the
# test is run in the form the proof consumes: i is a leader of SOME sub-segment containing its forward window iff y_(i+s) > y_i for
# 1 <= s <= k, and i is a peak of some sub-segment containing its backward window iff y_(i-s) < y_i for 1 <= s <= k. This covers
# every leader/peak of every sub-segment (the segment leaders are the suffix minima; the Proposition's proof uses only the window).
print('  (iv) leaders/peaks on long segments, k = L, X = 2^L (window form: y_(i+s) > y_i resp. y_(i-s) < y_i for 1 <= s <= k, which')
print('       is what the proof uses and covers all records of all sub-segments); exact checks: leaders 2 q^(o_s) > 2^s; peaks q^(5(o_k-o_s)) 2^8 > 2^(5(k-s)).')
random.seed(20260926)
rnd40 = random.randrange(2 ** 39, 2 ** 40) | 1
rnd80 = random.randrange(2 ** 79, 2 ** 80) | 1
segments = [('3n+1 orbit of 2^40-1', orbit_until(2 ** 40 - 1), 3, (30, 40, 64)),
            ('3n+1 orbit of random 40-bit odd %d' % rnd40, orbit_until(rnd40), 3, (30, 40, 64)),
            ('3n+1 orbit of random 80-bit odd %d' % rnd80, orbit_until(rnd80), 3, (30, 40, 64)),
            ('3n+1 orbit of 27', orbit_until(27), 3, (12, 20)),
            ('5n+1 orbit of 7 to 2^60', orbit_until(7, q=5, bound=1 << 60), 5, (30, 40))]
for name, ys, q, Ls in segments:
    check('%s: %d terms, all distinct, max = 2^%.1f' % (name, len(ys), math.log2(max(ys))), len(set(ys)) == len(ys))
    n_ = len(ys)
    lq = math.log2(q)
    for Lc in Ls:
        Xc = 1 << Lc; kc = Lc
        Y0 = 2 * (Fraction(q, 2) ** kc)
        lead_tot = lead_big = lead_viol = lead_viol_big = 0; lead_min = 99.0; lead_min_big = 99.0; lead_viol_list = []
        for i in range(n_ - kc):
            if ys[i] > Xc or not all(ys[i + s] > ys[i] for s in range(1, kc + 1)):
                continue
            lead_tot += 1
            big = ys[i] >= Y0
            lead_big += big
            w = parity_word(ys[i], kc, q); o = 0; smin = 99.0; bad = False
            for s, x in enumerate(w, 1):
                o += x
                smin = min(smin, o * lq - s)
                if not (2 * q ** o > 2 ** s):
                    bad = True
            lead_min = min(lead_min, smin)
            if big:
                lead_min_big = min(lead_min_big, smin)
            if bad:
                lead_viol += 1; lead_viol_big += big; lead_viol_list.append((i, ys[i], round(math.log2(float(Y0)), 1)))
        peak_tot = peak_big = peak_viol = peak_viol_big = 0; peak_min = 99.0; peak_min_big = 99.0; peak_viol_list = []
        for i in range(kc, n_):
            if ys[i] > Xc or not all(ys[i - s] < ys[i] for s in range(1, kc + 1)):
                continue
            peak_tot += 1
            big = ys[i] >= Y0
            peak_big += big
            w = parity_word(ys[i - kc], kc, q)
            os_ = [0]
            for x in w:
                os_.append(os_[-1] + x)
            dmin = min(os_[kc] * lq - kc - (os_[s] * lq - s) for s in range(kc))
            peak_min = min(peak_min, dmin)
            if big:
                peak_min_big = min(peak_min_big, dmin)
            if not all(q ** (5 * (os_[kc] - os_[s])) * 2 ** 8 > 2 ** (5 * (kc - s)) for s in range(kc)):
                peak_viol += 1; peak_viol_big += big; peak_viol_list.append((i, ys[i], round(math.log2(float(Y0)), 1)))
        fmt = lambda v: ('%.3f' % v) if v < 99 else 'n/a'
        print('    %s, L=%d: Y_0 = 2(%d/2)^k = 2^%.1f (%s X)' % (name, Lc, q, math.log2(Y0), '>' if Y0 > Xc else '<='))
        print('        leaders <= X with k successors: %d (>= Y_0: %d); min partial sum: all %s, >= Y_0 %s; violations of "> -1": %d (of which >= Y_0: %d) %s'
              % (lead_tot, lead_big, fmt(lead_min), fmt(lead_min_big), lead_viol, lead_viol_big, lead_viol_list[:6]))
        print('        peaks <= X with i >= k: %d (>= Y_0: %d); min (S_k - max_s S_s): all %s, >= Y_0 %s; violations of "> -1.6": %d (of which >= Y_0: %d) %s'
              % (peak_tot, peak_big, fmt(peak_min), fmt(peak_min_big), peak_viol, peak_viol_big, peak_viol_list[:6]))
        check('%s L=%d: no violation among records with y_i >= Y_0 (the Proposition\'s hypothesis)' % (name, Lc), lead_viol_big == 0 and peak_viol_big == 0)
        check('%s L=%d: no violation at all (leaders > -1: %d; peaks > -1.6: %d)' % (name, Lc, lead_viol, peak_viol), lead_viol == 0 and peak_viol == 0)
print('  note: for 5n+1 the Proposition\'s threshold Y_0 = 2(5/2)^k exceeds X, so it asserts nothing there; the words still pass because the actual')
print('  carries are bounded by the orbit\'s reciprocal sum (Proposition 6: y_(i+s) = M_s y_i prod(1 + 1/(5 y_j))), not by the worst case (5/2)^s.')
print('  (iii) counting: each k-word is one class mod 2^k, 2^k > X/2 so <= 2 representatives in [1, X]; M_k(1), M_k(1.6) <= D_s 2^(hk) k^(-3/2) 2^(1.6 lambda*) e^(1.6 s)')
print('  (Lemma M); peaks with i < k number <= k; leaders with fewer than k successors (finite segment) number <= k -- the note omits this (E)-type term for leaders.')
claim('section 8 (iii): the counting (two representatives per class, Lemma M at y = 1 and 1.6, the i < k peaks, O(X^(h*) L^(-3/2)))', True,
      'cosmetic: the leaders of a finite segment with fewer than k successors (at most k) are not mentioned')
claim('section 8 (iv): every leader k-word has all partial sums > -1 and every peak backward k-word ends within 1.6 bits of its maximum on the tested segments', True)

# [5] --------------------------------------------------------------------------------------------------
print('\n[5] Section 3b of the THM-4499 note (discrepancy corollary)')
a_star = LAM / H - 1.5
print('  1/h* = %.4f (note 1.0527), a*/h* = %.4f (note -1.038): displayed coefficient of log_2 log_2 L is a/h* = (a* + eps)/h* = -1.0382 + eps/h*.' % (1 / H, a_star / H))
check('1/h* = 1.0527 and a*/h* = -1.038 to the displayed precision', abs(1 / H - 1.0527) < 6e-5 and abs(a_star / H + 1.038) < 6e-4)


def odd_iterates(n, q, b, count):
    ms = [n]; ds = [0]
    d = 0
    for _ in range(count):
        x = q * n + b; v = 0
        while x % 2 == 0:
            x //= 2; v += 1
        d += v; n = x; ms.append(n); ds.append(d)
    return ms, ds


for name, n0, q, b, cnt in (('3n+1, m_0 = 27', 27, 3, 1, 41), ('3n+1, m_0 = 2^40-1', 2 ** 40 - 1, 3, 1, 150), ('5n+1, m_0 = 7', 7, 5, 1, 80)):
    ms, ds = odd_iterates(n0, q, b, cnt)
    ok = True
    P = Fraction(1); logP_max = 0.0
    for l in range(len(ms)):
        # m_l 2^(d_l) = q^l m_0 + b B_l  and  m_l 2^(d_l)/q^l = m_0 prod_(j<l)(1 + b/(q m_j))
        if Fraction(ms[l] * 2 ** ds[l], q ** l) != n0 * P:
            ok = False
        logP_max = max(logP_max, math.log2(P))
        P *= 1 + Fraction(b, q * ms[l])
    Delta = [ds[l] - l * math.log2(q) for l in range(len(ms))]
    lmax = max(range(len(ms)), key=lambda l: ms[l]); lmin = min(range(len(ms)), key=lambda l: Delta[l])
    lhs = math.log2(max(ms)); rhs = math.log2(n0) - min(Delta)
    check('%s: Proposition 6 identity m_l 2^(d_l)/q^l = m_0 prod_(j<l)(1 + b/(q m_j)) exact for l <= %d' % (name, cnt), ok)
    print('    argmax m_l = %d, argmin Delta_l = %d; log_2 max m_l = %.4f, log_2 m_0 - min Delta_l = %.4f, difference %.4f in [0, log_2 P_max = %.4f]'
          % (lmax, lmin, lhs, rhs, lhs - rhs, logP_max))
    check('%s: log_2 max m_l - (log_2 m_0 - min Delta_l) in [0, max log_2 P] (sign convention m_l = m_0 2^(-Delta_l) P_l)' % name,
          -1e-9 <= lhs - rhs <= logP_max + 1e-9)
print('  (i) L+1 distinct odd iterates <= X = max m_l: CONFIRMED (odd iterates are T-orbit terms; distinct since not eventually periodic).')
print('  (ii) log_2 X = log_2 m_0 - min Delta_l + O(1): CONFIRMED; 1 <= P_l <= exp((1/3) sum 1/m_j) <= e^(K/3) by THM-4476 Cor. 6 (sum over all orbit terms).')
print('  (iii) from L+1 <= K X^(h*) Lam^a (Lam = log_2 X, a < 0): h* Lam >= log_2 L + |a| log_2 Lam - O(1); dropping |a| log_2 Lam >= 0 gives')
print('       Lam >= (1/h*) log_2 L - O(1), so log_2 Lam >= log_2 log_2 L - O(1); re-inserting: -min Delta_l >= (1/h*) log_2 L + (|a|/h*) log_2 log_2 L - O(1),')
print('       i.e. min Delta_l <= -(1/h*) log_2 L + (a/h*) log_2 log_2 L + O(1). Only the LOWER bound on Lam is needed (a < 0): CONFIRMED.')
print('       Caveat: with a = a* + eps the coefficient is -(1.038 - eps/h*); the displayed "-1.038 log_2 log_2 L" holds for every coefficient')
print('       below 1.038, not with 1.038 itself (that would need a = a*, which THM-4499 does not give).')
print('  (iv) THM-4476 Cor. 3 (no m_j <= K j^a, a < 1/h*) follows: max_(l<=L) m_l >= c L^(1/h*) (log_2 L)^(1.038 - eps): consistent and stronger.')
claim('section 3b: min_(l<=L) Delta_l <= -(1/h*) log_2 L + (a/h*) log_2 log_2 L + O(1) for every a > a* (statement and proof)', True)
claim('section 3b displayed instantiation "with a = a* + eps: ... - 1.038 log_2 log_2 L + O(1)" (also in THM-4499\'s status and the commit message)', False,
      'the proved coefficient is -(1.038 - eps/h*); exactly 1.038 would need a = a*, which THM-4499 does not give')

# [6] --------------------------------------------------------------------------------------------------
print('\n[6] Provenance and the interlock probe')
files = ['04-computation/experiments/collatz_landing_20260926_probe.py', '04-computation/experiments/collatz_landing_20260926_probe.out',
         '04-computation/experiments/collatz_landing_20260926_probe2.py', '04-computation/experiments/collatz_landing_20260926_probe2.out',
         '04-computation/experiments/collatz_zeckendorf_20260926_interlock_probe.py', '04-computation/experiments/collatz_zeckendorf_20260926_interlock_probe.out',
         '05-knowledge/results/collatz_landing_20260926_multiplicity_reassessment.md', '05-knowledge/results/collatz_thin_20260926_little_o_thin_divergence.md',
         '04-computation/experiments/collatz_thin_20260926_movingbarrier.py', '04-computation/experiments/collatz_thin_20260926_movingbarrier.out',
         '04-computation/experiments/collatz_thin_20260926_movingbarrier_audit.py', '05-knowledge/results/collatz_thin_20260926_movingbarrier_audit.out']
head = subprocess.run(['git', 'log', '-3', '--format=%h %ci %s'], cwd=ROOT, capture_output=True, text=True).stdout.strip().splitlines()
print('  worktree HEAD at audit time: %s' % head[0][:110])
print('  previous: %s' % head[1][:110])
print('  NOTE: the reassessment note was revised during this audit (ab9ef18b4, section 8 peak argument); hashes below are of the revised text.')
print('  sha256 (raw bytes as checked out | LF-normalised | git HEAD blob):')
lfhash = {}
for f in files:
    p = os.path.join(ROOT, f)
    raw = open(p, 'rb').read()
    lf = raw.replace(b'\r\n', b'\n')
    lfhash[f] = hashlib.sha256(lf).hexdigest()
    try:
        blob = subprocess.run(['git', 'show', 'HEAD:' + f], cwd=ROOT, capture_output=True).stdout
        gh = hashlib.sha256(blob).hexdigest()[:16]
    except Exception:
        gh = 'n/a'
    print('    %s\n      raw %s  LF %s  HEAD %s  (CRLF lines in checkout: %d)' % (f, hashlib.sha256(raw).hexdigest(), lfhash[f], gh, raw.count(b'\r\n')))
rec = {'04-computation/experiments/collatz_thin_20260926_movingbarrier.py': 'da2d12874ab368dbff10136fea23e780687633e27b92c38eaf3d3963df05722e',
       '04-computation/experiments/collatz_thin_20260926_movingbarrier.out': '0b2c11a427dd218cf4edef6a21d0823aa9ad9c41303dbe91383380c8584878b6',
       '04-computation/experiments/collatz_thin_20260926_movingbarrier_audit.py': '9af70c922047ce0a5dd2bb8aeadcc471cebb5da95802a00a7ccad6acd7bbbd08',
       '05-knowledge/results/collatz_thin_20260926_movingbarrier_audit.out': '645ef0d39fab9b85441bf05650ba0ce788baa46a0c07152d0390290fcd7cbee3'}
for f, hsh in rec.items():
    check('THM-4499 recorded sha256 matches LF bytes of %s' % os.path.basename(f), lfhash[f] == hsh)
print('  the reassessment note and the probe scripts carry NO hashes (scripts cited by name only).')
claim('THM-4499 canon file: recorded script/output/audit hashes match the LF bytes in the worktree', all(lfhash[f] == h for f, h in rec.items()))
note_txt = open(os.path.join(ROOT, '05-knowledge/results/collatz_landing_20260926_multiplicity_reassessment.md'), encoding='utf-8').read()
claim('section 7: the Zeckendorf conclusion is qualified by "At this resolution"', 'At this resolution' in note_txt,
      'stated in section 7; the commit message a810b9bb4 ("no interlock with the Collatz step") and the status line of HYP-9161 carry no such qualifier')
claim('status lines: section 1 "(PROVED)" label', False, 'the climb-then-drop bullet is false and the hover bullet overstates by a factor ~W+1; the qualitative claim survives')
claim('status lines: section 8 "(PROVED)" label', False, 'statement true, proof text incomplete in both versions (see [4](ii)); a corrected proof is supplied above')
claim('status lines: THM-4499 status "tight on residue classes (hover/climb then drop)"', False, 'climb-then-drop gives multiplicity <= 2; hover gives ~ (2/3) L at best')

# interlock: chi-square of the printed transition matrices (counts reconstructed from row-normalised entries and row sizes)
def chi2_sf(x, dof):   # even dof
    m = dof // 2
    return math.exp(-x / 2) * sum((x / 2) ** i / math.factorial(i) for i in range(m))


mats = {'cZ_diag under T': ([[0.335, 0.328, 0.337], [0.330, 0.339, 0.330], [0.332, 0.332, 0.336]], [55691, 55507, 55467]),
        'cZ_diag(n) -> cZ_diag(3n+1), odd n': ([[0.320, 0.355, 0.325], [0.341, 0.315, 0.344], [0.344, 0.329, 0.327]], [11005, 11265, 11062])}
for name, (M, ns) in mats.items():
    chi = 0.0; maxdev = 0.0
    for row, n in zip(M, ns):
        for p in row:
            e = n / 3; chi += (p * n - e) ** 2 / e; maxdev = max(maxdev, abs(p - 1 / 3))
    print('  %s: max |entry - 1/3| = %.3f (%.1f points), chi-square vs uniform rows = %.1f on 6 dof, p = %.1e' % (name, maxdev, 100 * maxdev, chi, chi2_sf(chi, 6)))
    if '3n+1' in name:
        claim('section 7: transition matrices "uniform to within 2 percent" (3n+1 matrix)', maxdev <= 0.02,
              '%.1f points off 1/3 and highly significant (p ~ 1e-8); negligible in magnitude but not noise' % (100 * maxdev))

# mutual information: hand check on n = 2..13 (cZ_len vs parity) and the probe's own function on the same data
hand_zeck = {2: [2], 3: [3], 4: [3, 1], 5: [5], 6: [5, 1], 7: [5, 2], 8: [8], 9: [8, 1], 10: [8, 2], 11: [8, 3], 12: [8, 3, 1], 13: [13]}
hand_table = {(1, 0): 2, (1, 1): 3, (2, 0): 3, (2, 1): 3, (0, 0): 1}     # (len mod 3, parity) counts, by hand (12 = 8+3+1 has 3 terms, 3 mod 3 = 0)
hand_I = (2 / 12) * math.log2((2 / 12) / ((5 / 12) * (6 / 12))) + (3 / 12) * math.log2((3 / 12) / ((5 / 12) * (6 / 12))) + 0 + 0 + (1 / 12) * math.log2((1 / 12) / ((1 / 12) * (6 / 12)))
spec = importlib.util.spec_from_file_location('interlock', os.path.join(ROOT, '04-computation/experiments/collatz_zeckendorf_20260926_interlock_probe.py'))
ilk = importlib.util.module_from_spec(spec); spec.loader.exec_module(ilk)
F = ilk.fib_upto(200000)
zk_ok = all(sorted(sum(F[i - 2] for i in []) for _ in []) == [] for _ in [0])  # placeholder
zk_ok = True
for n, rep in hand_zeck.items():
    idx = ilk.zeck(n, F)
    if sorted(F[i - 2] for i in idx) != sorted(rep):
        zk_ok = False
check('probe zeck(n) agrees with the hand Zeckendorf representations of n = 2..13', zk_ok)
pairs = [(len(hand_zeck[n]) % 3, n % 2) for n in range(2, 14)]
check('hand contingency table (cZ_len, parity) for n = 2..13 equals the probe-derived one', Counter(pairs) == Counter(hand_table))
I_probe, Hy = ilk.mutual_info(pairs)
check('probe mutual_info on n = 2..13 = hand value %.6f bits (probe %.6f, H(parity) = %.3f)' % (hand_I, I_probe, Hy), abs(I_probe - hand_I) < 1e-9 and abs(Hy - 1.0) < 1e-9)
# independent recomputation of the largest reported MI (cZ_len vs v2(3n+1) mod 3, odd n <= 200000) and its noise floor
N = 200000
cnt = Counter(); cx = Counter(); cy = Counter(); tot = 0
for n in range(3, N + 1, 2):
    z = ilk.zeck(n, F); c = len(z) % 3
    x = 3 * n + 1; v = 0
    while x % 2 == 0:
        x //= 2; v += 1
    t = v % 3
    cnt[(c, t)] += 1; cx[c] += 1; cy[t] += 1; tot += 1
I = sum((k / tot) * math.log2(k * tot / (cx[c] * cy[t])) for (c, t), k in cnt.items())
chi = 2 * tot * math.log(2) * I
print('  independent I(cZ_len; v2(3n+1) mod 3) over odd n <= %d: %.5f bits (probe 0.00027); null expectation (r-1)(c-1)/(2N ln2) = %.1e bits;'
      % (N, I, 4 / (2 * tot * math.log(2))))
print('  G-statistic 2N ln2 I = %.1f on 4 dof, p = %.1e: ~9x the noise floor and statistically significant, though 0.02%% of H(target).' % (chi, chi2_sf(chi, 4)))
check('reported 0.00027 reproduced', abs(I - 0.00027) < 2e-5)
claim('section 7: "all mutual informations are at the noise floor" (cZ_len vs v2(3n+1) mod 3: %.1fx the null expectation, p = %.0e)' % (I / (4 / (2 * tot * math.log(2))), chi2_sf(chi, 4)),
      chi2_sf(chi, 4) > 0.01, 'false for this entry; the note\'s "at this resolution" qualifier (stated in section 7) is what saves the conclusion')

print('\n' + '=' * 100)
print('internal/reproduction checks: %d, failures: %d (must be 0)' % (NCHECK[0], len(FAIL)))
for f in FAIL:
    print('  FAILED: ' + f)
print('note claims tested: %d, refuted: %d' % (NCLAIM[0], len(REFUTED)))
for f in REFUTED:
    print('  CLAIM FALSE: ' + f)
print('=' * 100)
