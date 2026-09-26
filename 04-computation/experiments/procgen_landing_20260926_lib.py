#!/usr/bin/env python3
"""procgen_landing_20260926_lib.py -- landing multiplicity of the THM-4476 / THM-4499 recursion.

Lane `landing`, session collatz-procgen-20260922, 2026-09-26.

Setting (THM-4476, section 1.5). T_b(x) = x/2 (x even), (3x+b)/2 (x odd), b odd. For a
segment y_0, y_1, ... of distinct positive integers, a scale X, the window k = floor(log2 X)
and a depth D > 0 (D = theta log2 X), an index i with y_i <= X and i + k <= M (M = last
index) is a DIPPER if y_(i+s) < y_i 2^(-D) for some 1 <= s <= k; its LANDING index is i + s
for the least such s.  The landing multiplicity of j is the number of dippers landing at j.

Everything here is exact integer arithmetic except where a float fast path is used; every
float decision closer than 1e-7 (in log2 units) to a tie is re-decided exactly.
"""
import math
from fractions import Fraction

ALPHA = math.log2(3.0)                      # log2 3
RHO = math.log(2.0) / math.log(3.0)         # log_3 2
H_STAR = -(RHO * math.log2(RHO) + (1 - RHO) * math.log2(1 - RHO))
LAMBDA_STAR = math.log2(RHO / (1 - RHO)) / ALPHA


def a_star(mu):
    """Exponent of THM-4499's bootstrap when the (averaged) multiplicity is O(L^mu)."""
    return mu * LAMBDA_STAR / H_STAR - 1.5


def T(y, b=1):
    return y >> 1 if y % 2 == 0 else (3 * y + b) // 2


def segment(y, n, b=1):
    out = [y]
    for _ in range(n):
        out.append(T(out[-1], b))
    return out


def orbit_to_one(y, b=1, cap=10 ** 6):
    out = [y]
    while out[-1] != 1 and len(out) <= cap:
        out.append(T(out[-1], b))
    return out


def as_fraction(D):
    """Depths are dyadic rationals with denominator <= 64 (callers use depth())."""
    F = D if isinstance(D, Fraction) else Fraction(D)
    if F.denominator > 64:
        raise ValueError("depth must have denominator <= 64; use depth()")
    return F


def depth(x):
    """Round a real depth theta*log2(X) to the nearest multiple of 1/64 (exact thereafter)."""
    return Fraction(round(64 * x), 64)


def below_exact(u, v, D):
    """Exact test u < v 2^(-D) for positive integers u, v and rational D = p/q >= 0."""
    F = as_fraction(D)
    p, q = F.numerator, F.denominator
    if q == 1:
        return (u << p) < v
    return (pow(u, q) << p) < pow(v, q)


def below(u, v, D):
    """u < v 2^(-D), float fast path with an exact fallback within 1e-7 of a tie."""
    d = math.log2(u) + float(D) - math.log2(v)
    if d < -1e-7:
        return True
    if d > 1e-7:
        return False
    return below_exact(u, v, D)


def landing_map(ys, X, D, k=None, logs=None):
    """Dippers and landing indices of the segment ys at scale X and depth D (window k).

    Returns (land, nd, e): land = {j: [dippers in increasing order]}, nd = number of
    no-dip indices with y_i <= X, e = number of end indices (i + k > M) with y_i <= X.
    """
    if k is None:
        k = X.bit_length() - 1              # floor(log2 X) for integer X >= 1
    M = len(ys) - 1
    Df = float(D)
    L2 = logs if logs is not None else [math.log2(v) for v in ys]
    land, nd, e = {}, 0, 0
    for i, v in enumerate(ys):
        if v > X:
            continue
        if i + k > M:
            e += 1
            continue
        thr = L2[i] - Df
        hit = None
        for s in range(1, k + 1):
            d = L2[i + s] - thr
            if d < -1e-7:
                hit = i + s
                break
            if d < 1e-7 and below_exact(ys[i + s], v, D):
                hit = i + s
                break
        if hit is None:
            nd += 1
        else:
            land.setdefault(hit, []).append(i)
    return land, nd, e


def ub_mult(k, D, b):
    """Proposition U (D >= 1): b > 0: m(j) <= ceil((k - D)/alpha) for every landing index;
    b < 0: m(j) <= ceil((k - D + log2(3/2))/alpha) for landing values y_j >= k|b|."""
    x = (k - float(D)) / ALPHA if b > 0 else (k - float(D) + math.log2(1.5)) / ALPHA
    return math.ceil(x - 1e-12)


def check_structure(ys, j, dippers, D, b):
    """Lemma S (shell), Lemma T (time), Lemma O (odd separation) for one landing index.
    Returns a list of violated items (empty if all hold)."""
    bad = []
    F = as_fraction(D)
    yj = ys[j]
    for i in dippers:
        # shell: 2^D y_j < y_i <= 2^(D+1) y_j, i.e. y_j < y_i 2^-D and y_i 2^-D <= 2 y_j
        if not below(yj, ys[i], F):
            bad.append(("shell-low", i))
        if below(2 * yj, ys[i], F):
            bad.append(("shell-high", i))
        if j - i < math.floor(F) + 1:
            bad.append(("time", i))
    for a, c in zip(dippers, dippers[1:]):
        if not any(ys[t] % 2 for t in range(a, c)):
            bad.append(("odd-sep", a, c))
    if len(dippers) >= 2 and ys[j - 1] != 2 * yj:
        bad.append(("halving-into-landing", j))
    return bad


def depth_intervals(ys, j, k):
    """For the landing index j, the integer depths D >= 1 with landing_D(i) = j, for every
    i in [j-k, j-1] (ignoring the X-restriction).  Lemma A says each i has at most one."""
    out = {}
    for i in range(max(0, j - k), j):
        mids = ys[i + 1:j]
        lo_int = []
        for D in range(1, 2 * ys[i].bit_length() + 2):
            if below(ys[j], ys[i], D) and all(not below(m, ys[i], D) for m in mids):
                lo_int.append(D)
        out[i] = lo_int
    return out


# ----------------------------------------------------------------------------------------
# Terras inverse and the Sturmian hover-then-halve construction
# ----------------------------------------------------------------------------------------

def terras_inverse(word, b=1):
    """The residue r in [0, 2^n) whose T_b parity word of length n is `word` (list of 0/1).
    Adding 2^t u to a number with a fixed t-prefix changes T^t by 3^(o_t) u (odd), so the
    t-th letter is fixed by one bit of r."""
    r, v, o = 0, 0, 0
    for t, w in enumerate(word):
        if (v & 1) != w:
            r += 1 << t
            v += 3 ** o
        if w:
            v = (3 * v + b) // 2
            o += 1
        else:
            v //= 2
    return r


def parity_word(y, n, b=1):
    w = []
    for _ in range(n):
        w.append(y & 1)
        y = T(y, b)
    return w


def sturmian_hover(O, a):
    """Word visiting the unit band [-a, 1-a) of the pure walk S_t = o_t alpha - t at every
    odd count o = 0..O, at times n(o) = floor(o alpha + a); returns (word, visit times)."""
    word, times = [], [0]
    n_prev = 0
    for o in range(1, O + 1):
        n_o = math.floor(o * ALPHA + a)
        gap = n_o - n_prev
        word += [1] if gap == 1 else [1, 0]
        times.append(n_o)
        n_prev = n_o
    return word, times


def hostile_word(k, D, a, O=None):
    """Sturmian hover of O+1 band visits followed by floor(D)+1 halvings, of length <= k."""
    c = math.floor(float(D)) + 1
    if O is None:
        O = 0
        while math.floor((O + 1) * ALPHA + a) + c <= k:
            O += 1
    hw, times = sturmian_hover(O, a)
    return hw + [0] * c, times, c


def realize(word, b, lo, hi):
    """Largest integer y in [lo, hi] with T_b-word prefix `word` (None if none)."""
    n = len(word)
    r = terras_inverse(word, b)
    m = 1 << n
    y = hi - ((hi - r) % m)
    return y if y >= lo else None


def dippers_at(ys, j, X, D, k):
    """The dippers of the segment ys that land at index j (scale X, depth D, window k);
    the segment must extend at least to j (dippers i >= j - k with a full window are assumed
    to exist beyond j when j + k > len(ys) - 1 is irrelevant for landing at j)."""
    Df = float(D)
    lj = math.log2(ys[j])
    out = []
    for i in range(max(0, j - k), j):
        if ys[i] > X:
            continue
        li = math.log2(ys[i])
        if not (lj < li - Df - 1e-7 or (abs(lj - li + Df) <= 1e-7 and below_exact(ys[j], ys[i], D))):
            continue
        if any(below(ys[t], ys[i], D) for t in range(i + 1, j)):
            continue
        out.append(i)
    return out


def best_hostile(k, D, b=1, X=None, tries=64):
    """Search offsets a for a Sturmian hover-then-halve integer y <= X with the largest actual
    landing multiplicity at scale X (window k), depth D.  Returns (m, y, a, len(word))."""
    if X is None:
        X = (1 << (k + 1)) - 1
    fracD = float(D) - math.floor(float(D))
    best = (0, None, None, None)
    for O_shift in (0, 1, 2):
        for t in range(tries):
            eps = (t + 0.5) / (tries * 50.0)
            _, _, c = hostile_word(k, D, 0.5)
            # choose a so that frac(O alpha + a) = 1 - frac(D) - eps (shell = band shifted by eps)
            O = 0
            while math.floor((O + 1) * ALPHA + 0.999) + c <= k:
                O += 1
            O = max(1, O - O_shift)
            a = (1 - fracD - eps - O * ALPHA) % 1.0
            word, times, c = hostile_word(k, D, a, O)
            if len(word) > k:
                continue
            y = realize(word, b, X >> 3, int(X * 2 ** (a - 1) / 1.01))
            if y is None:
                continue
            ys = segment(y, len(word), b)
            if len(set(ys)) < len(ys):
                continue
            j = len(word)
            m = len(dippers_at(ys, j, X, D, k))
            if m > best[0]:
                best = (m, y, a, len(word))
    return best


# ----------------------------------------------------------------------------------------
# Proposition H: a residue class of heavy dippers with a uniform margin
# ----------------------------------------------------------------------------------------

def dist_int(x):
    return abs(x - round(x))


def margin_word(O, D):
    """Sturmian hover of O+1 visits with offset a chosen so that the shell of the landing
    point coincides with the visit band up to eps; then floor(D)+1 halvings.  Returns
    (word, visit times, safe visits) where a visit o (r = O - o) is 'safe' if r = 0 or
    ||r alpha + D|| >= 1/8 (Proposition H: at least one of any two consecutive r >= 1)."""
    eps = 1e-6
    a = (1 - (float(D) % 1.0) - eps - O * ALPHA) % 1.0
    hw, times = sturmian_hover(O, a)
    c = math.floor(float(D)) + 1
    safe = [o for o in range(O + 1) if o == O or dist_int((O - o) * ALPHA + float(D)) >= 0.125]
    return hw + [0] * c, times, safe


def check_margin_class(O, D, b, X, samples, rng):
    """For random y in the residue class of margin_word(O, D) with 48 b l < y <= X/4, the
    landing index l = len(word) receives every safe visit as a dipper.  Returns the list of
    (y, m, number of safe visits)."""
    word, times, safe = margin_word(O, D)
    l = len(word)
    r = terras_inverse(word, b)
    k = X.bit_length() - 1
    lo, hi = 48 * abs(b) * l + 1, X // 4
    out = []
    for _ in range(samples):
        u = rng.randrange((lo - r) // (1 << l) + 1, (hi - r) // (1 << l) + 1)
        y = r + (u << l)
        if not (lo <= y <= hi):
            continue
        ys = segment(y, l + k + 1, b)
        land, _, _ = landing_map(ys, X, D, k)
        ds = set(land.get(l, []))
        out.append((y, len(ds), len(safe), all(times[o] in ds for o in safe)))
    return out


# ----------------------------------------------------------------------------------------
# Pure walk (2-adic model): positions z0 + o_t alpha - t, exact comparisons
# ----------------------------------------------------------------------------------------

def walk_from_word(word):
    o, t, out = 0, 0, [(0, 0)]
    for w in word:
        o += w
        t += 1
        out.append((o, t))
    return out


def pbelow(p, q, D):
    """Pure walk: position q < position p - D, positions (o, t), value o alpha - t."""
    do, dt = q[0] - p[0], q[1] - p[1]
    x = do * ALPHA - dt + float(D)
    if x < -1e-9:
        return True
    if x > 1e-9:
        return False
    F = as_fraction(D)                          # exact: 3^(do) 2^(D) < 2^(dt)
    lhs = pow(3, 64 * do) if do >= 0 else Fraction(1, pow(3, -64 * do))
    return lhs * Fraction(2) ** (64 * F) < Fraction(2) ** (64 * dt)


def landing_map_walk(pos, z0, logX, D, k):
    """landing_map for the pure walk with starting height z0 (log2), scale 2^logX."""
    M = len(pos) - 1
    land, nd, e = {}, 0, 0
    for i, p in enumerate(pos):
        if z0 + p[0] * ALPHA - p[1] > logX:
            continue
        if i + k > M:
            e += 1
            continue
        hit = None
        for s in range(1, k + 1):
            if pbelow(p, pos[i + s], D):
                hit = i + s
                break
        if hit is None:
            nd += 1
        else:
            land.setdefault(hit, []).append(i)
    return land, nd, e


# ----------------------------------------------------------------------------------------
# Theorem 1 (saturation): psi(X) = X^h* L^(a*(mu)) solves the recursion for every depth
# ----------------------------------------------------------------------------------------

def saturation_margin(mu, L, D, C=1.0):
    """log2 of [C L^mu psi(X 2^-D) + E X^h* L^-3/2 2^(lambda* D)] / psi(X) with K = 1 and
    E = C^(-lambda*/h*); Theorem 1 says it is >= 0 for every D >= 1 (L > D)."""
    a = a_star(mu)
    E = C ** (-LAMBDA_STAR / H_STAR)
    t1 = math.log2(C) + mu * math.log2(L) - D * H_STAR + a * (math.log2(L - D) - math.log2(L))
    t2 = math.log2(E) - 1.5 * math.log2(L) + LAMBDA_STAR * D - a * math.log2(L)
    return max(t1, t2) + math.log2(1 + 2 ** (min(t1, t2) - max(t1, t2)))


# ----------------------------------------------------------------------------------------
# Records (OEIS b-files; cached in scratch, fetched with a generic User-Agent if absent)
# ----------------------------------------------------------------------------------------

UA = "Mozilla/5.0 (research; math-repo)"


def load_bfile(name, cache, fallback_dirs=()):
    import os
    import shutil
    import subprocess
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, name)
    if not os.path.exists(path):
        for d in fallback_dirs:
            src = os.path.join(d, name)
            if os.path.exists(src):
                shutil.copy(src, path)
                break
    if not os.path.exists(path):
        url = "https://oeis.org/A%s/%s" % (name[1:7], name)
        subprocess.run(["curl", "-s", "-f", "-L", "-A", UA, "-o", path, url], check=True, timeout=120)
    vals = []
    for line in open(path):
        line = line.strip()
        if line and not line.startswith("#"):
            vals.append(int(line.split()[1]))
    return path, vals


def aligned_hostile_word(k, D, eps=1e-3):
    """Sturmian hover whose visit band is the landing shell shifted by eps, then floor(D)+1
    halvings, with the largest O such that the word has length <= k.  Returns (word, a, O)."""
    c = math.floor(float(D)) + 1
    fracD = float(D) % 1.0
    O = int(k / ALPHA) + 2
    while O > 0:
        a = (1 - fracD - eps - O * ALPHA) % 1.0
        if math.floor(O * ALPHA + a) + c <= k:
            hw, _ = sturmian_hover(O, a)
            return hw + [0] * c, a, O
        O -= 1
    raise ValueError("window too short")
