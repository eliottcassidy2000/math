#!/usr/bin/env python3
"""procgen_sporadic_20260926_lib.py -- helpers for the lane "sporadic" (free and sporadic cycles).

Session collatz-procgen-20260922, lane "sporadic", 2026-09-26.

Conventions.
  * T = T_{q,d}(y) = y/2 (y even), (q y + d)/2 (y odd), q odd >= 3, d odd, gcd(q, d) = 1.
  * A parity word w of length p with a ones at positions s_0 < ... < s_{a-1} has carry
        c_w = sum_i q^(a-1-i) 2^(s_i),
    gap D = D(p, a) = 2^p - q^a, and 2^p T^p(y) = q^a y + d c_w on the cylinder of w, so the unique
    2-adic point with itinerary w^infinity is  y_w = d c_w / D  (Boehm-Sontacchi; Lagarias 1990).
  * A shape (p, a) is FREE for T_{q,d} when every word of the shape gives an integer y_w.

Everything is exact integer / Fraction arithmetic except the few mpmath evaluations that are used
only to locate candidates (every such candidate is re-checked exactly).
"""
import hashlib
import math
import os
import resource
import subprocess
import sys
import time
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SCRATCH = os.path.join(ROOT, "scratch", "procgen_sporadic")

N_CHECKS = [0]


def check(cond, msg):
    """print a claim only after verifying it; raise on failure."""
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)
    N_CHECKS[0] += 1
    print("  [ok] " + msg, flush=True)


def peak_rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 * 1024) if sys.platform == "darwin" else r / 1024


def child_peak_rss_mb():
    r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return r / (1024 * 1024) if sys.platform == "darwin" else r / 1024


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------
# words, carries, periodic points
# ------------------------------------------------------------------------------------------------
def carry(w, q):
    """c_w for a 0/1 sequence w (Horner form: c <- q c + 2^t at an odd letter t)."""
    c = 0
    for t, b in enumerate(w):
        if b:
            c = q * c + (1 << t)
    return c


def gap(p, a, q):
    return (1 << p) - q ** a


def Tmap(y, q, d):
    return y // 2 if y % 2 == 0 else (q * y + d) // 2


def periodic_point(w, q, d):
    p, a = len(w), sum(w)
    return Fraction(d * carry(w, q), gap(p, a, q))


def parity_word(y, q, d, p):
    out = []
    for _ in range(p):
        out.append(y & 1)
        y = Tmap(y, q, d)
    return tuple(out), y


def words_of_shape(p, a):
    for pos in combinations(range(p), a):
        w = [0] * p
        for s in pos:
            w[s] = 1
        yield tuple(w)


def primitive_period(w):
    p = len(w)
    for k in range(1, p + 1):
        if p % k == 0 and w == w[k:] + w[:k]:
            return k
    return p


def mobius(n):
    res, m, f = 1, n, 2
    while f * f <= m:
        if m % f == 0:
            m //= f
            if m % f == 0:
                return 0
            res = -res
        f += 1
    if m > 1:
        res = -res
    return res


def lyndon_count(p, a):
    """number of aperiodic necklaces (Lyndon words) of length p with a ones."""
    g = math.gcd(p, a)
    tot = sum(mobius(k) * math.comb(p // k, a // k) for k in range(1, g + 1) if g % k == 0)
    assert tot % p == 0
    return tot // p


def necklace_rep(w):
    return min(w[k:] + w[:k] for k in range(len(w)))


def cycle_of(y, q, d, maxlen=10 ** 7):
    """the cycle through y (list, starting at y) if y is periodic, else None."""
    orb = [y]
    x = Tmap(y, q, d)
    while x != y:
        orb.append(x)
        if len(orb) > maxlen:
            return None
        x = Tmap(x, q, d)
    return orb


def cycle_data(orb, q, d):
    """(p, a, least element, word of the least element, D) for a cycle given as a list."""
    p = len(orb)
    a = sum(v & 1 for v in orb)
    m = min(orb)
    i = orb.index(m)
    rot = orb[i:] + orb[:i]
    w = tuple(v & 1 for v in rot)
    return p, a, m, w, gap(p, a, q)


# ------------------------------------------------------------------------------------------------
# Diophantine side: |2^p - q^a| small
# ------------------------------------------------------------------------------------------------
def gersonides_solutions(q, pmax):
    """all (p, a), p >= 1, a >= 1, with |2^p - q^a| = 1 and p <= pmax (exact)."""
    out = []
    for p in range(1, pmax + 1):
        a = 1
        while q ** a <= (1 << p) + 1:
            if abs((1 << p) - q ** a) == 1:
                out.append((p, a))
            a += 1
    return out


def gersonides_predicted(q, pmax):
    """the elementary classification: a = 1 and q = 2^p +- 1, or (q, p, a) = (3, 3, 2)."""
    out = []
    for p in range(1, pmax + 1):
        if q in ((1 << p) - 1, (1 << p) + 1):
            out.append((p, 1))
    if q == 3:
        out.append((3, 2))
    return sorted(out)


ELLISON_S = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 19, 27}


def small_gaps_q3(N):
    """all (p, a) with p, a >= 1 and |2^p - 3^a| <= N, complete by Ellison 1971 (|2^x - 3^y| >
    2^x e^(-x/10) for x not in S): returns the list and the bound x0 beyond which x not in S is
    excluded.  The exceptional x in S are enumerated directly."""
    x0 = 1
    while not ((1 << x0) * math.exp(-x0 / 10) > N and all((1 << x) * math.exp(-x / 10) > N
                                                           for x in range(x0, x0 + 60))):
        x0 += 1
    out = []
    for p in sorted(set(range(1, x0)) | ELLISON_S):
        a = 1
        while 3 ** a <= (1 << p) + N:
            if abs((1 << p) - 3 ** a) <= N:
                out.append((p, a))
            a += 1
    return sorted(out), x0


# ------------------------------------------------------------------------------------------------
# Stern-Brocot path of theta = log_q 2 (exact: a/p < log_q 2  <=>  q^a < 2^p)
# ------------------------------------------------------------------------------------------------
def below_theta(a, p, q):
    """True iff a/p < log_q 2 (p >= 1)."""
    return q ** a < (1 << p)


def sb_path(q, maxden):
    """the Stern-Brocot path (mediants) of log_q 2 with denominators <= maxden, as a list of
    (a, p, side) with side = 'L' (below theta: best lower approximation) or 'U' (above)."""
    ln, ld, rn, rd = 0, 1, 1, 0
    out = []
    while True:
        mn, md = ln + rn, ld + rd
        if md > maxden:
            break
        if below_theta(mn, md, q):
            out.append((mn, md, "L"))
            ln, ld = mn, md
        else:
            out.append((mn, md, "U"))
            rn, rd = mn, md
    return out


def on_sb_path(a, p, q):
    """is the reduced fraction a/p a node of the Stern-Brocot path of log_q 2?  0/1 counts as the
    left boundary (best lower approximation) and 1/1 is the first mediant."""
    g = math.gcd(a, p)
    a0, p0 = a // g, p // g
    if (a0, p0) == (0, 1):
        return True, "L"
    for (x, y, side) in sb_path(q, p0):
        if (x, y) == (a0, p0):
            return True, side
    return False, ("L" if below_theta(a0, p0, q) else "U")


# ------------------------------------------------------------------------------------------------
# the lattice of clocks L_d = {(p, a) : 2^p = 3^a mod d} and its near-critical cone
# ------------------------------------------------------------------------------------------------
def ordmod(b, d):
    if d == 1:
        return 1
    e, x = 1, b % d
    while x != 1:
        x = x * b % d
        e += 1
    return e


def lattice_basis(d, q=3):
    """Gauss-reduced basis of L_d = {(p, a) in Z^2 : 2^p = q^a (mod d)} and its index."""
    e = ordmod(2, d)
    pos = {}
    x = 1
    for p in range(e):
        pos[x] = p
        x = x * 2 % d
    a0, t = 1, q % d
    while t not in pos:
        a0 += 1
        t = t * q % d
    b1, b2 = (e, 0), (pos[t], a0)
    dot = lambda u, v: u[0] * v[0] + u[1] * v[1]
    u, v = b1, b2
    while True:
        if dot(u, u) > dot(v, v):
            u, v = v, u
        m = round(Fraction(dot(u, v), dot(u, u)))
        if m == 0:
            break
        v = (v[0] - m * u[0], v[1] - m * u[1])
    if u[1] < 0 or (u[1] == 0 and u[0] < 0):
        u = (-u[0], -u[1])
    idx = abs(u[0] * v[1] - u[1] * v[0])
    assert idx == e * a0
    return u, v, idx, e, a0


def _frac_log2(num_expr):
    import mpmath
    mpmath.mp.dps = 70
    x = num_expr(mpmath)
    s = mpmath.nstr(x, 66, strip_zeros=False)
    return Fraction(s)


def cone_min_clock(d, X, q=3, pmax=None):
    """least p over the clocks (p, a) of L_d with q^a < 2^p < (q + d/(X+1))^a, i.e. the clocks whose
    perigee window d/(2^(p/a) - q) reaches X + 1.  Every primitive T_d-cycle with period p below this
    value therefore has least element <= X.  Returns ((p, a), e); with pmax given, returns (None, e)
    if there is no such clock with p < pmax.  Exact: both comparisons use 66-digit rationals and are
    asserted to be far (> a 10^-60) from ties; the optimum is re-verified with integers."""
    e = ordmod(2, d)
    tab = {}
    x = 1 % d
    for p in range(e):
        tab[x] = p
        x = x * 2 % d
    th = _frac_log2(lambda M: M.log(q) / M.log(2))
    # a cycle with least element m >= X + 1 needs d/(2^(p/a) - q) >= X + 1, i.e. 2^p < (q + d/(X+1))^a
    r = _frac_log2(lambda M: M.log(q + M.mpf(d) / (X + 1)) / M.log(2))
    thN, thD, rN, rD = th.numerator, th.denominator, r.numerator, r.denominator
    E = 10 ** 60
    best = None
    tq = 1 % d
    a = 0
    while True:
        a += 1
        tq = tq * q % d
        fl, rem = divmod(a * thN, thD)
        assert rem * E > a * thD and (thD - rem) * E > a * thD, "log_2 q tie"
        pl = fl + 1                      # least p with 2^p > q^a
        if best is not None and pl >= best[0]:
            break
        if pmax is not None and pl >= pmax:
            break
        p0 = tab.get(tq)
        if p0 is None:
            continue
        pst = pl + ((p0 - pl) % e)
        diff = pst * rD - a * rN
        assert abs(diff) * E > a * rD, "cone tie"
        if diff < 0 and (best is None or pst < best[0]) and (pmax is None or pst < pmax):
            best = (pst, a)
    if best is None:
        return None, e
    p, a = best
    assert (1 << p) > q ** a and (pow(2, p, d) - pow(q, a, d)) % d == 0
    if p <= 250000:
        assert (1 << p) * (X + 1) ** a < (q * (X + 1) + d) ** a
    return best, e


# ------------------------------------------------------------------------------------------------
# C helpers
# ------------------------------------------------------------------------------------------------
def compile_c(src, exe):
    os.makedirs(os.path.dirname(exe), exist_ok=True)
    cmd = ["cc", "-O2", "-Wall", "-o", exe, src]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("compile failed: " + r.stderr)
    return exe


def run_parallel(cmds, outs, nproc=2):
    """run shell-free commands with at most nproc at once; stdout to the given files."""
    pending = list(zip(cmds, outs))
    running = []
    while pending or running:
        while pending and len(running) < nproc:
            cmd, out = pending.pop(0)
            fo = open(out, "w")
            running.append((subprocess.Popen(cmd, stdout=fo), fo, cmd))
        time.sleep(0.2)
        still = []
        for pr, fo, cmd in running:
            if pr.poll() is None:
                still.append((pr, fo, cmd))
            else:
                fo.close()
                if pr.returncode != 0:
                    raise RuntimeError("command failed: %r" % (cmd,))
        running = still


def parse_sweep(files):
    """parse sweep outputs: per-d summaries {d: [(ylo, yhi, nprim, nnon, nunres, nover), ...]} and
    cycles [(d, m, p, a, prim)], plus the other lines."""
    summ, cyc, other = {}, [], []
    for f in files:
        with open(f) as fh:
            for line in fh:
                t = line.split()
                if not t:
                    continue
                if t[0] == "D":
                    v = tuple(map(int, t[1:]))
                    summ.setdefault(v[0], []).append(v[1:])
                elif t[0] == "c":
                    cyc.append(tuple(map(int, t[1:])))
                else:
                    other.append(line.strip())
    return summ, cyc, other


def parse_traj(files):
    cyc, ent, bad, nstarts = [], [], [], 0
    for f in files:
        with open(f) as fh:
            for line in fh:
                t = line.split()
                if t[0] == "C":
                    cyc.append(tuple(int(v) for v in t[1:]))
                elif t[0] == "E":
                    ent.append(tuple(int(v) for v in t[1:]))
                elif t[0] == "S":
                    nstarts += int(t[1])
                else:
                    bad.append(line.strip())
    return cyc, ent, bad, nstarts
