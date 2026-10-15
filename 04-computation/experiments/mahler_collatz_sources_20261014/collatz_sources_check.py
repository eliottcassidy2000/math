#!/usr/bin/env python3
"""Rational shadows of safe codes: exact checks of the logarithmic-memory and complexity-approximation bounds.

For a positive integer A with Mahler code e (M(x) = ceil(3x/2)) and a rational x = u/v (odd v) whose Collatz (T) code
agrees with e on a SAFE common prefix of length h, the theorem (THM in mahler_collatz_sources_20261014.md) gives

  h < log2(3) * (n1(x) + min(floor(2.1 v) + 1, l(v) + N_e(l(v)))) + log2(A + 1),

n1(x) = 0 if |x| <= 1/6 else ceil(log2(6|x|)/c),  c = 1 - log2(3)/3,  l(v) = floor(log2(6.3 v)/c) + 1,
N_e(l) = number of distinct length-l factors of e in the agreement window.

Checks: (a) the slope/offset lemmas on every safe word to length 22; (b) the entry window [-1/2, 1.6) for every
rational of height <= 60 along every safe word to length 18; (c) the bound for every A <= AMAX with safe depth >= 20
against every rational u/v, v odd <= VMAX, |u| <= UMAX; (d) the A_t = (2^(18t+8) - 9)/19 family against x = 1/5.
"""
import math
import sys
from fractions import Fraction as F
from itertools import product

AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
VMAX, UMAX = 45, 120
LOG3 = math.log2(3)
C = 1 - LOG3 / 3


def step(m, e):
    t = F(3, 2) * m - e
    return None if t <= 0 else min(t, F(1))


def mahler_code(A, L):
    out, a = [], A
    for _ in range(L):
        out.append(a & 1)
        a = (3 * a + (a & 1)) // 2
    return out


def safe_depth(code):
    m = F(1)
    for i, e in enumerate(code):
        m = step(m, e)
        if m is None:
            return i
    return len(code)


def collatz_code_rational(u, v, L):
    out, x = [], u                       # x = numerator over the fixed odd denominator v
    for _ in range(L):
        b = x & 1                        # parity of x/v in Z_(2) equals parity of x (v odd)
        out.append(b)
        x = (3 * x + v) // 2 if b else x // 2
    return out


def n1(u, v):
    ax = abs(u) / v
    return 0 if ax <= 1 / 6 else math.ceil(math.log2(6 * ax) / C)


def ell(v):
    return math.floor(math.log2(6.3 * v) / C) + 1


def factors(code, lo, hi, l):
    return len({tuple(code[i:i + l]) for i in range(max(lo, 0), hi - l + 1)})


# (a) slope and offset lemmas
worst_slope, worst_off, nw = 0.0, F(0), 0
for L in range(1, 23):
    for w in product((0, 1), repeat=L):
        m = F(1)
        ok = True
        for e in w:
            m = step(m, e)
            if m is None:
                ok = False
                break
        if not ok:
            continue
        nw += 1
        k = sum(w)
        assert k <= L // 3 + 1
        off = F(0)
        for e in w:
            off = (3 * off + 1) / 2 if e else off / 2
        assert off < F(11, 10)
        worst_off = max(worst_off, off)
        worst_slope = max(worst_slope, 3 ** k / 2 ** L / (3 * (3 ** (1 / 3) / 2) ** L))
print(f"(a) {nw} safe words of length 1..22: ones <= floor(l/3)+1, offset T_w(0) < 11/10 (max {float(worst_off):.6f}), "
      f"slope/(3 lambda^l) <= {worst_slope:.4f}")

# (b) entry window: every rational of height <= 60 along every safe word of length <= 18
viol = 0
cnt = 0
safe18 = [w for w in product((0, 1), repeat=18) if safe_depth(list(w)) == 18]
for v in range(1, 61, 2):
    for u in range(-60, 61):
        if math.gcd(u, v) != 1:
            continue
        code = collatz_code_rational(u, v, 18)
        for L in range(1, 19):
            pass
        # the rational follows its own code; check the window whenever its own prefix is safe
        d = safe_depth(code)
        x = F(u, v)
        lam = 3 ** (1 / 3) / 2
        for n in range(d + 1):
            if n > 0:
                b = code[n - 1]
                x = (3 * x + 1) / 2 if b else x / 2
            lo, hi = -3 * lam ** n * abs(u / v), 3 * lam ** n * abs(u / v) + 1.1
            cnt += 1
            if not (lo - 1e-12 <= float(x) < hi + 1e-12):
                viol += 1
print(f"(b) window x_n in [-3 lam^n |x|, 3 lam^n |x| + 11/10) along safe own-prefixes: {cnt} checks, {viol} violations")

# (c) the bound on actual positive integers
LCODE = 120
found = []
for A in range(1, AMAX + 1):
    code = mahler_code(A, LCODE)
    D = safe_depth(code)
    if D >= 20:
        found.append((A, D, code))
found.sort(key=lambda t: -t[1])
print(f"(c) positive integers A <= {AMAX} with safe depth >= 20: {len(found)}; deepest: "
      + ", ".join(f"A={A} (depth {D})" for A, D, _ in found[:8]))
nchk = nviol = 0
tight = []
for A, D, code in found:
    safe = code[:D]
    for v in range(1, VMAX + 1, 2):
        for u in range(-UMAX * v // 10, UMAX * v // 10 + 1):
            if math.gcd(u, v) != 1:
                continue
            cx = collatz_code_rational(u, v, D)
            h = 0
            while h < D and cx[h] == safe[h]:
                h += 1
            if h == 0:
                continue
            l = ell(v)
            N = factors(safe, n1(u, v), h, l)
            bound_state = LOG3 * (n1(u, v) + math.floor(2.1 * v) + 1) + math.log2(A + 1)
            bound_cplx = LOG3 * (n1(u, v) + l + N) + math.log2(A + 1)
            bound = min(bound_state, bound_cplx)
            nchk += 1
            if h >= bound:
                nviol += 1
            tight.append((h / bound, A, u, v, h, round(bound, 2)))
tight.sort(reverse=True)
print(f"    bound checked on {nchk} (A, u/v) pairs with nonempty safe agreement: {nviol} violations; "
      f"tightest h/bound: {[(round(t[0], 3),) + t[1:] for t in tight[:5]]}")
# the 2-adic separation (9a) of the interface note, for comparison, at the tightest pairs
for r, A, u, v, h, bd in tight[:3]:
    Q = max(abs(u), v)
    print(f"    A={A}, x={u}/{v}: agreement {h}; new bound {bd}; interface bound 14Q+2+bitlength(A+1) = "
          f"{14 * Q + 2 + (A + 1).bit_length()}")

# (d) the family A_t against 1/5: agreement exactly L_t = 18t+8 bits
for t in range(0, 6):
    A = (2 ** (18 * t + 8) - 9) // 19
    code = mahler_code(A, 18 * t + 40)
    cx = collatz_code_rational(1, 5, 18 * t + 40)
    h = 0
    while code[h] == cx[h]:
        h += 1
    D = safe_depth(code)
    print(f"(d) t={t}: A_t={A if A < 10**12 else '2^' + str(18*t+8) + '...'}: agreement with 1/5 = {h} "
          f"(= 18t+8: {h == 18 * t + 8}), safe depth {D}, bound {LOG3 * 3 + math.log2(A + 1):.2f}")
