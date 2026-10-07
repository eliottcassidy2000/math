#!/usr/bin/env python3
"""The excess-rate family F_q: weak resets paid by a 1/q deadline.

A first-descent segment of an odd source x is its orbit under U(n)=oddpart(3n+1)
up to the first value y < x; it has length l (odd steps), total valuation A and
word w. It is q-admissible iff  2^(qA - l) >= 3^(ql), i.e. the excess
A - l log_2 3 is at least l/q. The family F_q consists of the odd sources all of
whose chained segments (down to ROOT) are q-admissible; q(n) is the least such q.

Theorem (deadline). For n in F_q: every segment with start x >= q satisfies
y^(2q) 2^l <= x^(2q), hence sum of those lengths <= 2q log_2 n, and the segments
starting below q are bounded by the finite table of odd sources below q. A
source-only deadline T_q(n) follows conditional on a supplied checked finite seed table, and Codex's backward compiler turns it
into a positive weight floor eta(n, T) = 2/((C+2) binom(C+1, floor((C+1)/2))).

Also: a first-descent word w (non-rising, 2^A > 3^l) descends from x exactly
when x > x_w = c_w/(2^A - 3^l), a POSITIVE rational cycle point; weak resets are
shadows of large positive rational cycles.

Companion note: 05-knowledge/results/collatz_weak_reset_family_20261005.md
Usage: python3 <this file> [--census-bits 20] [--json PATH]
"""
from fractions import Fraction as F
from math import comb, factorial, log, log2, ceil
import argparse
import json

import numpy as np

try:
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap if not a or not callable(a[0]) else a[0]

CHECKS = 0
LOG2_3 = log2(3.0)


def require(cond, witness=None):
    global CHECKS
    CHECKS += 1
    if not cond:
        raise RuntimeError(witness)


def v2(n):
    return (n & -n).bit_length() - 1


def U(n):
    t = 3 * n + 1
    return t >> v2(t)


def root_word(n):
    w = []
    while n != 1:
        t = 3 * n + 1
        a = v2(t)
        w.append(a)
        n = t >> a
    return w


def counters(word):
    if not word:
        return 0, 0
    return len(word) - 1, sum((a - 1) // 2 for a in word)


def wW(L, K):
    return F(2 * factorial(K) * factorial(L + 1), factorial(L + K + 2))


def carry(word):
    l = len(word)
    c = 0
    Apre = 0
    for i, a in enumerate(word):
        c += 3 ** (l - 1 - i) * (1 << Apre)
        Apre += a
    return c, Apre


# ---------------------------------------------------------------------------
# segments, admissibility, q(n)
# ---------------------------------------------------------------------------
def segment(x):
    """First-descent segment of odd x>1: (word, endpoint)."""
    w = []
    v = x
    while True:
        t = 3 * v + 1
        a = v2(t)
        v = t >> a
        w.append(a)
        if v < x:
            return w, v


def chain(n):
    """Chained first-descent segments from n down to ROOT."""
    segs = []
    x = n
    while x != 1:
        w, y = segment(x)
        segs.append((x, w, y))
        x = y
    return segs


def admissible(l, A, q):
    return (1 << (q * A - l)) >= 3 ** (q * l) if q * A - l >= 0 else False


def q_segment(l, A):
    q = 1
    while not admissible(l, A, q):
        q += 1
        if q > 10 ** 6:
            raise RuntimeError("no q")
    return q


def q_source(n):
    return max((q_segment(len(w), sum(w)) for (_, w, _) in chain(n)), default=0)


@njit(cache=True)
def _census(limit):
    """For odd n<limit: q(n) via float excess (verified exactly elsewhere), tau(n),
    longest segment length, and the start of the critical segment."""
    m = (limit + 1) // 2
    qn = np.zeros(m, dtype=np.int64)
    tau = np.zeros(m, dtype=np.int64)
    lmax = np.zeros(m, dtype=np.int64)
    crit = np.zeros(m, dtype=np.int64)
    log23 = np.log2(3.0)
    for i in range(1, m):
        n = 2 * i + 1
        x = n
        qbest = 0
        t = 0
        lbest = 0
        cstart = 0
        while x != 1:
            v = x
            l = 0
            A = 0
            while True:
                tt = 3 * v + 1
                a = 0
                while (tt & 1) == 0:
                    tt >>= 1
                    a += 1
                v = tt
                l += 1
                A += a
                if v < x:
                    break
            t += l
            exc = A - l * log23
            qs = int(np.ceil(l / exc - 1e-9))
            if qs < 1:
                qs = 1
            if qs > qbest:
                qbest = qs
                cstart = x
            if l > lbest:
                lbest = l
            x = v
        qn[i] = qbest
        tau[i] = t
        lmax[i] = lbest
        crit[i] = cstart
    return qn, tau, lmax, crit


# ---------------------------------------------------------------------------
# the deadline, the compiled floor (Codex backward compiler, sections 2-3)
# ---------------------------------------------------------------------------
def deadline(n, q, tau_table_max):
    """Historical numerical formula, conditional on a truthful checked seed maximum.
    Floating log2 is for the finite display, not an all-input exact certificate.
    """
    return int(2 * q * log2(n)) + tau_table_max


def deadline_sharp(n, q, tau_table_max_10q):
    """T'_q(n) = floor(1.051 q log2 n) + max tau over odd sources below 10q:
    conditional on checked seeds; the repaired contraction is 197/207,
    whose reciprocal 207/197 is below 1.051. This display uses floating log2."""
    return int(1.051 * q * log2(n)) + tau_table_max_10q


def log10_frac(fr):
    """log10 of a positive Fraction without float underflow."""
    return (fr.numerator.bit_length() - fr.denominator.bit_length()) * log(2) / log(10) +         log(fr.numerator / (1 << (fr.numerator.bit_length() - 1))) / log(10) -         log(fr.denominator / (1 << (fr.denominator.bit_length() - 1))) / log(10)


def compiled_floor(n, T):
    """eta(n,T) from the backward compiler: 2^A <= n (10/3)^T, N <= C."""
    Q = (n * 10 ** T) // 3 ** T
    a_max = Q.bit_length() - 1
    C = (a_max + T - 1) // 2 - 1
    if C < 1:
        return None, C
    return F(2, (C + 2) * comb(C + 1, (C + 1) // 2)), C


def section_theory(limit_check=1 << 16):
    print("== S1. Segments, positive rational cycle points, admissibility, deadline inequality ==")
    # positive cycle points: a first-descent word descends from x iff x > x_w
    checked = 0
    for n in range(3, limit_check, 2):
        for (x, w, y) in chain(n):
            l = len(w)
            c, A = carry(w)
            require((1 << A) > 3 ** l, ("segment word is non-rising", n, w))
            xw = F(c, (1 << A) - 3 ** l)
            require(xw > 0 and F(x) > xw, ("descent iff x > x_w", n, x, w))
            require(y == (3 ** l * x + c) >> A, ("forward formula", n, w))
            # proper prefixes are rising-or-not-yet-descended: no earlier descent
            checked += 1
            # a cylinder mate below x_w does not descend along w (THM-4512 exceptional member)
            r2 = (c * pow(((1 << A) - 3 ** l) % (1 << (A + 1)), -1, 1 << (A + 1))) % (1 << (A + 1))
            require(x % (1 << (A + 1)) == r2, ("cylinder", n, w))
            if r2 < xw and r2 >= 1:
                v = r2
                ok = True
                for a in w:
                    t = 3 * v + 1
                    if v2(t) != a:
                        ok = False
                        break
                    v = t >> a
                if ok:
                    require(v >= r2, ("below x_w must not descend", n, w, r2))
    print("  %d chained segments below 2^%d: non-rising words, descent iff start > x_w = c_w/(2^A-3^l) > 0,"
          " cylinder residue and forward formula exact" % (checked, limit_check.bit_length() - 1))
    # the deadline inequality: q-admissible segment with start x >= q has y^(2q) 2^l <= x^(2q)
    viol = 0
    tested = 0
    for n in range(3, limit_check, 2):
        for (x, w, y) in chain(n):
            l = len(w)
            A = sum(w)
            q = q_segment(l, A)
            if x >= q:
                require(y ** (2 * q) * (1 << l) <= x ** (2 * q), ("segment inequality", n, x, w, q))
                tested += 1
    print("  segment inequality y^(2q) 2^l <= x^(2q) holds at the least admissible q for all %d segments with x >= q" % tested)
    # examples
    print("  examples (source: chain of (start, word, end) with per-segment least q):")
    for n in (5, 7, 23, 27, 47, 55, 739):
        segs = chain(n)
        desc = "; ".join("%d-%s->%d q%d" % (x, "".join(str(a) if a < 10 else "(%d)" % a for a in w), y,
                                             q_segment(len(w), sum(w))) for (x, w, y) in segs)
        print("   %5d: q(n)=%d, tau=%d : %s" % (n, q_source(n), len(root_word(n)), desc[:150] + ("..." if len(desc) > 150 else "")))
    return


def section_floors(tau_table):
    print("== S2. Source-only deadlines and compiled floors for members ==")
    rows = []
    for n in (5, 7, 23, 27, 47, 55, 739, 1161):
        q = q_source(n)
        tmax = max((int(tau_table[(y - 1) // 2]) for y in range(1, q, 2)), default=0)
        tmax10 = max((int(tau_table[(y - 1) // 2]) for y in range(1, 10 * q, 2)), default=0)
        T = deadline(n, q, tmax)
        Ts = deadline_sharp(n, q, tmax10)
        eta, C = compiled_floor(n, Ts)
        L, K = counters(root_word(n))
        actual = wW(L, K)
        tau_n = len(root_word(n))
        require(eta is not None and eta <= actual, ("floor below actual", n, eta, actual))
        require(tau_n <= Ts <= T or tau_n <= T, ("deadline", n, T, Ts))
        # designated leaf: n itself if 3 | n, else rho(n); selector degree with 8(16/25)^d <= eta_leaf/2
        if n % 3 != 0:
            a0 = {1: 6, 2: 5, 4: 4, 5: 1, 7: 2, 8: 3}[n % 9]
            z = ((1 << a0) * n - 1) // 3
            eta_z, Cz = compiled_floor(z, Ts + 1)
        else:
            z = n
            eta_z, Cz = eta, C
        Lz, Kz = counters(root_word(z))
        require(eta_z <= wW(Lz, Kz), ("leaf floor", z))
        d = int(ceil((log(16) / log(10) - log10_frac(eta_z)) / (log(25 / 16) / log(10))))
        rows.append((n, q, T, Ts, C, eta, actual, z, eta_z, d, tau_n))
        print("   n=%5d q=%4d tau=%3d T_q=%5d T'_q=%5d C=%5d log10 eta=%9.1f log10 W=%7.1f leaf=%6d log10 eta_leaf=%9.1f degree=%d"
              % (n, q, tau_n, T, Ts, C, log10_frac(eta), log10_frac(actual), z, log10_frac(eta_z), d))
    return rows


# ---------------------------------------------------------------------------
# run-block tree of Codex (strong guard a >= r+2), containment in F_3
# ---------------------------------------------------------------------------
def run_block_member(n):
    """Codex run-block recognizer: blocks (1^r, a) with a >= r+2, unit parents, down to 5."""
    x = n
    while x != 5:
        if x < 5 or x % 3 == 0:
            return False
        r = v2(x + 1) - 1
        v = x
        for _ in range(r):
            t = 3 * v + 1
            if v2(t) != 1:
                return False
            v = t >> 1
        t = 3 * v + 1
        a = v2(t)
        y = t >> a
        if a < r + 2 or y % 3 == 0 or y >= x or y < 5:
            return False
        x = y
    return True


def section_census(bits, qn, tau, lmax, crit):
    print("== S3. Census of q(n) below 2^%d ==" % bits)
    n = 2 * np.arange(len(qn)) + 1
    sel = n > 1
    qs = qn[sel]
    total = int(sel.sum())
    thresholds = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 100, 150, 200, 300, 500]
    cover = {q: float((qs <= q).mean()) for q in thresholds}
    print("  proportion with q(n) <= q:", " ".join("%d:%.4f" % (q, cover[q]) for q in thresholds))
    imax = int(np.argmax(qn))
    print("  max q(n) = %d at n = %d (critical segment start %d, longest segment %d, tau %d)"
          % (int(qn[imax]), 2 * imax + 1, int(crit[imax]), int(lmax[imax]), int(tau[imax])))
    # exact verification of q at the maximum and at the examples
    for m in (7, 27, 55, 2 * imax + 1):
        qe = q_source(m)
        require(qe == int(qn[(m - 1) // 2]), ("exact q", m, qe, int(qn[(m - 1) // 2])))
    print("  exact integer verification of q(n) at 7, 27, 55 and the maximizer: agrees with the float census")
    # Codex run-block tree members below 2^bits are all in F_3
    members = [int(x) for x in n[sel] if run_block_member(int(x))]
    require(all(int(qn[(m - 1) // 2]) <= 3 for m in members), "run-block tree not inside F_3")
    print("  Codex run-block tree: %d members below 2^%d (incl. 5), all with q(n) <= 3; F_3 has %d members"
          % (len(members), bits, int((qs <= 3).sum())))
    # deadline theorem over the census: tau(n) <= floor(2 q log2 n) + max tau below q
    tau_prefix_max = np.maximum.accumulate(tau)   # max tau over odd y <= 2i+1
    bad = 0
    skipped = 0
    for i in range(1, len(qn)):
        q = int(qn[i])
        nn = 2 * i + 1
        j = (q - 2) // 2
        if j >= len(tau):
            skipped += 1
            continue
        tmax = int(tau_prefix_max[j]) if j >= 0 else 0
        if int(tau[i]) > int(2 * q * log2(nn)) + tmax:
            bad += 1
    require(bad == 0, ("deadline violations", bad))
    print("  deadline below-q table: checked %d sources; skipped %d whose required seed table exceeds this census" %
          (len(qn) - 1 - skipped, skipped))
    bad = 0
    skipped = 0
    worst = 0.0
    for i in range(1, len(qn)):
        q = int(qn[i])
        nn = 2 * i + 1
        j = (10 * q - 2) // 2
        if j >= len(tau):
            skipped += 1
            continue
        tmax = int(tau_prefix_max[j]) if 10 * q >= 3 else 0
        Ts = int(1.051 * q * log2(nn)) + tmax
        if int(tau[i]) > Ts:
            bad += 1
        worst = max(worst, int(tau[i]) / max(Ts, 1))
    require(bad == 0, ("sharp deadline violations", bad))
    print("  sharp deadline below-10q table: checked %d sources; skipped %d whose required table exceeds this census; worst ratio %.3f" %
          (len(qn) - 1 - skipped, skipped, worst))
    # distribution of the critical segment length among hard sources
    hard = qs >= 50
    print("  sources with q(n) >= 50: %.5f of the census; their longest segment has mean length %.1f (overall mean %.2f)"
          % (float(hard.mean()), float(lmax[sel][hard].mean()) if hard.any() else float('nan'), float(lmax[sel].mean())))
    return cover, dict(max_q=int(qn[imax]), argmax=2 * imax + 1, run_block_members=len(members))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--census-bits', type=int, default=20)
    ap.add_argument('--json', type=str, default='')
    args = ap.parse_args()
    section_theory(1 << min(16, args.census_bits))
    qn, tau, lmax, crit = _census(1 << args.census_bits)
    rows = section_floors(tau)
    cover, extra = section_census(args.census_bits, qn, tau, lmax, crit)
    print("== Summary ==")
    print("  checks: %d" % CHECKS)
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(dict(checks=CHECKS, coverage=cover, extra=extra,
                           rows=[dict(n=r[0], q=r[1], T=r[2], T_sharp=r[3], C=r[4], log10_eta=log10_frac(r[5]), log10_W=log10_frac(r[6]), leaf=r[7],
                                      log10_eta_leaf=log10_frac(r[8]), degree=r[9], tau=r[10]) for r in rows],
                           status="PROVED scoped family theorem; FINITE-EXACT checks; census VERIFIED"),
                      fh, indent=1, default=str)
        print("  json written:", args.json)


if __name__ == '__main__':
    main()
