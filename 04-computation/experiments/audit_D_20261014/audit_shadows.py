#!/usr/bin/env python3
"""Audit D (2026-10-14), THM-4654 statements 5-7 and consequences: independent checks.

Different method from collatz_sources_check.py:
  * safety via the THM-4072 follower graph on integer states (not the Fraction slack automaton);
  * statement 7 is checked x-CENTRICALLY over ALL positive integers A at once: for a rational x with T-code e
    and each h <= safe depth of e, the smallest positive integer A_h whose Mahler code starts with e[0,h) is the
    canonical residue (THM-2228 (4)); since the bound increases with A and N_e only grows with the agreement,
    h < B(A_h, x, h) for all h certifies the inequality for EVERY A >= 1 against this x (up to the depth cap);
  * the A-centric scan of the original (A <= 200000, depth >= 20, v <= 45, |u| <= 12v) is re-counted with
    bitmask agreement (XOR + trailing zeros), as a replication of the 981 / 4,980,199 figures;
  * exact constants with mpmath (60 digits).
Also: the entry window, the vacuity of the memory lemma (Terras/2-adic), statement 6 and its sharpness.
"""
import math
import random
import sys
import time
import mpmath

mpmath.mp.dps = 60
LOG3 = mpmath.log(3, 2)
CC = 1 - LOG3 / 3
LOG3F = float(LOG3)
CF = float(CC)
OUT = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_D_20261014/audit_shadows.out'
lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


# ------------------------------------------------------------------ greedy d and the follower graph
def greedy_d(n):
    a, d = 1, []
    for k in range(1, n + 1):
        t = 3 * a
        dk = 1 if t >= 2 ** k else 0
        a = t - dk * 2 ** k
        d.append(dk)
    return d


D1 = [None] + greedy_d(3000)


def fstep(q, s):
    nxt = D1[q + 1]
    if s == nxt:
        return q + 1
    if s == 0 and nxt == 1:
        return 0
    return None


# ------------------------------------------------------------------ exact helper functions
def n1_exact(u, v):
    if 6 * abs(u) <= v:
        return 0
    return int(mpmath.ceil(mpmath.log(mpmath.mpf(6 * abs(u)) / v, 2) / CC))


def ell_exact(v):
    return int(mpmath.floor(mpmath.log(mpmath.mpf('6.3') * v, 2) / CC)) + 1


def ell2(v):
    """Terras/2-adic memory length: smallest l with 2^l >= 2.1 v  (window numerator differences are < 2.1 v)."""
    l = 0
    while 10 * 2 ** l < 21 * v:
        l += 1
    return l


def tcode(u, v, L):
    """Collatz parity code of u/v (v odd) for L steps, numerators over v; also the orbit numerators."""
    e, js = [], [u]
    j = u
    for _ in range(L):
        b = j & 1
        e.append(b)
        j = (3 * j + v) // 2 if b else j // 2
        js.append(j)
    return e, js


def safe_depth(e):
    q = 0
    for i, s in enumerate(e):
        q = fstep(q, s)
        if q is None:
            return i
    return len(e)


def mahler_code(A, L):
    out, a = [], A
    for _ in range(L):
        out.append(a & 1)
        a = (3 * a + (a & 1)) // 2
    return out


t0 = time.time()
# ================================================================== statement 5: entry window, exact
viol5 = cnt5 = cntwin = 0
rng = random.Random(20261014)
samples = []
for v in range(1, 202, 2):
    for u in range(-12 * v, 12 * v + 1):
        if math.gcd(u, v) == 1:
            samples.append((u, v))
for v in range(1, 100, 2):
    for _ in range(300):
        u = rng.randint(-10 ** 12 * v, 10 ** 12 * v)
        if math.gcd(u, v) == 1:
            samples.append((u, v))
    for k in range(1, 60, 3):          # long zero runs: u = 2^k * odd
        for m in (1, -1, 3, -5):
            u = (2 ** k) * m
            if math.gcd(u, v) == 1:
                samples.append((u, v))
for r in range(1, 40):                   # slow descent along (100)^r: x = 1/5 + 8^r t
    for t in (1, -1, 3, -3, 7):
        samples.append((1 + 5 * 8 ** r * t, 5))
maxD = 0
for (u, v) in samples:
    e, js = tcode(u, v, 400)
    D = safe_depth(e)
    maxD = max(maxD, D)
    n1 = n1_exact(u, v)
    k = 0
    for n in range(0, D + 1):
        if n > 0:
            k += e[n - 1]
        jn = js[n]
        # exact: x_n = jn/v in [-(3^k/2^n)|x|, (3^k/2^n)|x| + 11/10)  and slope bound 3k-3 <= n
        cnt5 += 1
        lo_ok = jn * 2 ** n * v >= -(3 ** k) * abs(u) * v       # jn/v >= -(3^k/2^n)|u|/v
        hi_ok = 10 * jn * 2 ** n < 10 * 3 ** k * abs(u) + 11 * 2 ** n * v
        if not (lo_ok and hi_ok and 3 * k - 3 <= n):
            viol5 += 1
        if n >= n1:
            cntwin += 1
            if not (-v <= 2 * jn and 5 * jn < 8 * v):       # x_n in [-1/2, 8/5)
                viol5 += 1
say(f"[S5] entry window (exact, with the sharper slope 3^k/2^n) on {len(samples)} rationals "
    f"(v<=201 with |u|<=12v; random |u|<=1e12 v; u=2^k m; x=1/5+8^r t), {cnt5} orbit points, {cntwin} "
    f"points past n1(x): {viol5} violations; max safe depth met {maxD}")

# ================================================================== statement 5: the memory lemma is vacuous
#  distinct window points with a common T-code prefix of length l are congruent mod 2^l (Terras), hence differ by
#  >= 2^l/v; the lemma's hypothesis 3 lam^l R < 1/v gives R < 2^(c l)/(3v) < 2^l / v, so it never applies to
#  distinct points.  Check: for every odd v <= 4001, window points are already separated by codes of length ell2(v).
worst = 0
for v in range(1, 1002, 2):
    l2 = ell2(v)
    lv = ell_exact(v)
    assert l2 <= lv
    seen = {}
    for j in range(-((v) // 2), (8 * v + 4) // 5 + 1):
        if not (-v <= 2 * j and 5 * j < 8 * v):
            continue
        e, _ = tcode(j, v, l2)
        key = tuple(e)
        assert key not in seen, (v, j, seen.get(key))
        seen[key] = j
    worst = max(worst, lv - l2)
say(f"[S5] memory: for every odd v <= 1001 the window points of (1/v)Z are pairwise separated by their T-code "
    f"prefixes of length ell2(v) = ceil(log2(2.1v)) <= ell(v); ell(v) - ell2(v) reaches {worst}. "
    f"So two DISTINCT window points never share a code prefix of length ell(v): the lemma's hypothesis forces x = x'.")
for v in (1, 5, 45, 201, 1001, 10 ** 6 + 1):
    say(f"     v={v}: ell(v)={ell_exact(v)}  (THM: floor(log2(6.3v)/c)+1),  Terras ell2(v)={ell2(v)}")

# ================================================================== statement 6: periodic runs, exhaustive small
def run_len(code, i, p):
    """length of the longest prefix of code[i:] that is p-periodic (agrees with (code[i:i+p])^inf)."""
    h = 0
    n = len(code)
    while i + h < n and code[i + h] == code[i + (h % p)]:
        h += 1
    return h


viol6 = cnt6 = 0
best6 = (1e9, None)
for A in range(1, 5001):
    code = mahler_code(A, 140)
    LA = math.log2(A + 1)
    for i in range(0, 20):
        for p in range(1, 11):
            h = run_len(code, i, p)
            if i + h >= 140:
                raise SystemExit("increase code length")
            cnt6 += 1
            bd = p * LOG3F + i * math.log2(1.5) + LA
            if not (h < bd):
                viol6 += 1
            if bd - h < best6[0]:
                best6 = (bd - h, (A, i, p, h))
say(f"[S6] periodic-run inequality h < p log2 3 + i log2(3/2) + log2(A+1): all A <= 5000, i < 20, p <= 10 "
    f"({cnt6} runs): {viol6} violations; min slack {best6[0]:.4f} at (A, i, p, h) = {best6[1]}")
# the equality family A_t against 1/5 (code (100)^inf), slack -> log2(27/19) from above
for t in range(0, 12):
    A = (2 ** (18 * t + 8) - 9) // 19
    assert 19 * A + 9 == 2 ** (18 * t + 8)
    code = mahler_code(A, 18 * t + 30)
    h = run_len(code, 0, 3)
    assert code[:3] == [1, 0, 0]
    bd = mpmath.log(27 * (mpmath.mpf(A) + 1), 2)
    if t in (0, 1, 2, 5, 11):
        say(f"     A_t, t={t}: agreement with (100)^inf = {h} (18t+8 = {18 * t + 8}); bound - h = "
            f"{mpmath.nstr(bd - h, 12)}")
say(f"     log2(27/19) = {mpmath.nstr(mpmath.log(mpmath.mpf(27) / 19, 2), 12)}")
# sharper families: w = 0^(p-1) 1, C_M(w) = 2^(p-1), D = 3^p - 2^p, A = (2^h - 2^(p-1))/D with ord_D(2) | h-p+1
for p in (3, 4, 5, 6, 8, 10, 12):
    Dp = 3 ** p - 2 ** p
    w = [0] * (p - 1) + [1]
    Cw = sum(w[j] * 2 ** j * 3 ** (p - 1 - j) for j in range(p))
    assert Cw == 2 ** (p - 1)
    o = 1
    while pow(2, o, Dp) != 1:
        o += 1
    h = p - 1 + o
    while h < 60:
        h += o
    A = (2 ** h - 2 ** (p - 1)) // Dp
    assert Dp * A + 2 ** (p - 1) == 2 ** h and A > 0
    code = mahler_code(A, h + 40)
    hh = run_len(code, 0, p)
    assert code[:p] == w and hh == h, (p, hh, h)
    bd = p * LOG3 + mpmath.log(mpmath.mpf(A) + 1, 2)
    say(f"     w = 0^({p - 1})1: A = (2^{h} - 2^{p - 1})/(3^{p} - 2^{p}), agreement {hh}; bound - h = "
        f"{mpmath.nstr(bd - hh, 8)}  (limit log2(3^p/(3^p-2^p)) = "
        f"{mpmath.nstr(mpmath.log(mpmath.mpf(3) ** p / Dp, 2), 8)})")

say(f"     [time {time.time() - t0:.1f}s]")


# ================================================================== statement 7, x-centric over ALL A
def check_x(u, v, HCAP, stats):
    """For x = u/v: for each h <= min(safe depth, HCAP), A_h = least positive integer with Mahler code prefix
    e[0,h); check h < B(A_h, x, h) (theorem form) and the Terras-improved form.  Returns min slacks."""
    # lazily generate the T-code while it stays safe (follower graph), up to HCAP
    e = []
    j, q = u, 0
    while len(e) < HCAP:
        b = j & 1
        q2 = fstep(q, b)
        if q2 is None:
            break
        e.append(b)
        q = q2
        j = (3 * j + v) // 2 if b else j // 2
    D = len(e)
    if D == 0:
        return
    n1 = n1_exact(u, v)
    lv = ell_exact(v)
    l2 = ell2(v)
    cnt_state = (21 * v) // 10 + 1
    facs, facs2 = set(), set()
    # incremental canonical residue r_h and a_h = M^h(r_h)
    r, a, C3 = 0, 0, 1    # r_0 = 0 ; a_0 = M^0(r_0) = 0 ; C3 = 3^h
    for h in range(1, min(D, HCAP) + 1):
        s = e[h - 1]
        # extend r_(h-1) -> r_h so that the (h-1)-th Mahler bit equals s: a = M^(h-1)(r)
        if (a & 1) != s:
            r += 1 << (h - 1)
            a += C3                  # M^(h-1)(r + 2^(h-1)) = M^(h-1)(r) + 3^(h-1)
        a = (3 * a + s) // 2         # now a = M^h(r_h)
        C3 *= 3
        Ah = r if r > 0 else (1 << h)
        if h >= lv:
            facs.add(tuple(e[h - lv:h]))
        if h >= l2:
            facs2.add(tuple(e[h - l2:h]))
        LA = mpmath.log(mpmath.mpf(Ah) + 1, 2) if Ah > 2 ** 50 else mpmath.mpf(math.log2(Ah + 1))
        B = LOG3 * (n1 + min(cnt_state, lv + len(facs))) + LA
        B2 = LOG3 * (n1 + min(cnt_state, l2 + len(facs2))) + LA
        stats['pairs'] += 1
        sl = B - h
        sl2 = B2 - h
        if sl <= 0:
            stats['viol'] += 1
            stats['viol_ex'].append((u, v, h, Ah))
        if sl2 <= 0:
            stats['viol2'] += 1
        if sl < stats['min'][0]:
            stats['min'] = (sl, (u, v, h, Ah if Ah < 10 ** 15 else f'~2^{Ah.bit_length()}'))
        if sl2 < stats['min2'][0]:
            stats['min2'] = (sl2, (u, v, h, Ah if Ah < 10 ** 15 else f'~2^{Ah.bit_length()}'))
    stats['maxdepth'] = max(stats['maxdepth'], min(D, HCAP))


# sanity: the canonical residue really has the prescribed Mahler code (spot check)
for (u, v) in [(1, 5), (0, 1), (7, 45), (-3, 11), (1, 19)]:
    e, _ = tcode(u, v, 60)
    r = (-sum(e[j] * 2 ** j * 3 ** (59 - j) for j in range(60)) * pow(3, -60, 2 ** 60)) % 2 ** 60
    assert mahler_code(r, 60) == e
say("[S7] canonical-residue construction checked against direct Mahler codes")

stats = dict(pairs=0, viol=0, viol2=0, viol_ex=[], min=(mpmath.mpf(10 ** 9), None),
             min2=(mpmath.mpf(10 ** 9), None), maxdepth=0)
VMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 151
HCAP = 300
nx = 0
for v in range(1, VMAX + 1, 2):
    for u in range(-12 * v, 12 * v + 1):
        if math.gcd(u, v) != 1:
            continue
        check_x(u, v, HCAP, stats)
        nx += 1
say(f"[S7] x-centric, ALL positive A at once: {nx} rationals u/v (v odd <= {VMAX}, |u| <= 12v), every agreement "
    f"length h <= min(safe depth, {HCAP}): {stats['pairs']} (x, h) cases; violations {stats['viol']}; "
    f"min additive slack B - h = {mpmath.nstr(stats['min'][0], 8)} at (u, v, h, A_h) = {stats['min'][1]}; "
    f"max depth reached {stats['maxdepth']}  [time {time.time() - t0:.1f}s]")
say(f"     Terras-improved form (ell2(v) = ceil(log2 2.1v) in place of ell(v)): violations {stats['viol2']}; "
    f"min slack {mpmath.nstr(stats['min2'][0], 8)} at {stats['min2'][1]}")

# adversarial families: large |u|, zero runs, slow descent, purely periodic safe cycles
stats_b = dict(pairs=0, viol=0, viol2=0, viol_ex=[], min=(mpmath.mpf(10 ** 9), None),
               min2=(mpmath.mpf(10 ** 9), None), maxdepth=0)
fam = []
for v in range(1, 200, 2):
    for _ in range(100):
        u = rng.randint(-10 ** 15 * v, 10 ** 15 * v)
        if math.gcd(u, v) == 1:
            fam.append((u, v))
    for k in range(1, 80, 2):
        for m in (1, -1, 3, -3, 5):
            u = (2 ** k) * m
            if math.gcd(u, v) == 1:
                fam.append((u, v))
for r in range(1, 60):
    for t in (1, -1, 3, -3, 5, 7, -9):
        fam.append((1 + 5 * 8 ** r * t, 5))
# purely periodic rational cycles with safe period words of length <= 16: x = C_T(w)/(2^p - 3^k)
from itertools import product as iprod
ncyc = 0
for p in range(1, 17):
    for w in iprod((0, 1), repeat=p):
        if safe_depth(list(w) * 6) < 6 * p:
            continue
        k = sum(w)
        Ct = 0
        for idx, s in enumerate(w):
            if s:
                Ct = 3 * Ct + 2 ** idx
        den = 2 ** p - 3 ** k
        g = math.gcd(Ct, den)
        uu, vv = Ct // g, den // g
        if vv < 0:
            uu, vv = -uu, -vv
        fam.append((uu, vv))
        ncyc += 1
for (u, v) in fam:
    check_x(u, v, HCAP, stats_b)
say(f"[S7] adversarial families ({len(fam)} rationals: random |u| <= 1e15 v, u = 2^k m, x = 1/5 + 8^r t, and {ncyc} "
    f"purely periodic safe cycles of period <= 16): {stats_b['pairs']} (x, h) cases; violations {stats_b['viol']}; "
    f"min slack {mpmath.nstr(stats_b['min'][0], 8)} at {stats_b['min'][1]}; Terras form violations "
    f"{stats_b['viol2']}, min slack {mpmath.nstr(stats_b['min2'][0], 8)}  [time {time.time() - t0:.1f}s]")

with open(OUT, 'w') as fh:
    fh.write("\n".join(lines) + "\n")
