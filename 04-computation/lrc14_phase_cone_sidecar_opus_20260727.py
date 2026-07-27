#!/usr/bin/env python3
"""THM-2368 Sec 8: exact phase/event sidecar (37) of the BARE drift tensor
H_E on the canonical typed row, and the phase-cone certificate decision.

OBJECT.  Canonical typed row THM-2309 (25) (the THM-2541/THM-2550 row):

    w = (H,q1..q5,c1,c2,c3) = (1,14,27,40,53,66,13,2197,742586),
    owner j = 1 (c1 = 13, depth 1), targets a = c2 = 13^3, c = c3 = 2*13^5.

Grafts (THM-2365 Sec 1 concrete basis / drift companion):

    eta = e_c2 - e_q1,   ell = e_c3 - e_q2,   k_a = q1 = 14,  k_c = q2 = 27,
    k_a = k_c = 1 (mod 13).

THM-2368 (8): w_j = 13 A_j, w_a = 13 A_a, w_c = 13 A_c with

    A_j = 1,   A_a = 169,   A_c = 57122.

Under x = (y+h)/13 (THM-2368 (9)-(10), (29); minus-shift convention of
THM-2365 (5) exactly as in the audited drift companion
04-computation/lrc14_owner_loop_drift_typed_row_opus_20260727.py):

    A_y(r) = g((14 y + r)/13),          C_y(r) = g((27 y + r)/13),
    P_y(h) = gamma((y+h)/13) * prod_{q in {40,53,66}} g(q (y+h)/13),
    Q_{s,t}(y) = (1/13) sum_h A_y(h+s) C_y(h+t) P_y(h)        [k_a=k_c=1 mod 13],
    M_{r,s,t}(y) = d(y) g(169 y - s/13) g(57122 y - t/13) d(57122 y - r/13),
    H_E(r,s,t) = int_T M_{r,s,t}(y) Q_{s,t}(y) dy.

CONVENTION LOG (every choice):
  * g = 1 - d, d(z) = 1_{||z||<1/14}, gamma(z) = 1_{||z||>=1/7} (THM-2368 (2)).
  * The c3-blocker safe factor and the probe share one argument word:
    g(57122 y - t/13) = 1 - d(57122 y - t/13); we write dc[r] for the probe
    word and gc[t] = 1 - dc[t].  This forces the diagonal zero H(t,s,t)=0.
  * All breakpoints lie on the integer grid Z/T_DEN, T_DEN = 297836897838480
    (the drift companion's bare scale); chamber midpoints are evaluated on
    the doubled grid D2 = 2*T_DEN.  All decisions are integer/exact.
  * Only chambers inside supp d(y) = (13/14, 15/14) (wrapped) contribute:
    d(y) is a factor of M.  The full chamber partition (34) is still counted.
  * Walk object per nonzero target character (lam,mu) != (0,0):
      c_i = len_i * nd_i * V_i(lam,mu),
      V_i(lam,mu) = sum_{s,t} ga[s] gc[t] N[s,t] zeta^{-(lam s + mu t)},
      N[s,t] = 13 Q_{s,t},  nd_i = sum_r dc[r]  (r marginalised on the
      common chamber base BEFORE the DFT -- MISTAKE-281 discipline).
    The endpoint is sum_i c_i = 13 * T_DEN * sum_{s,t} S_E(s,t)
    zeta^{-(lam s + mu t)} = 13^4 * T_DEN * B(0,-lam,-mu) in THM-2365 (10)
    notation: these are exactly the 168 deep-colour-a=0 nonzero-target
    modes.  For a circulant tensor H(r,s,t) = G(r-t) every such mode is 0,
    so a nonzero endpoint certifies non-circulant, hence D_{H_E} > 0 by
    THM-2365 (17)-(18).  A cone-aligned walk has nonzero endpoint
    (THM-2368 Sec 8 phase-cone certificate).
  * AMBIGUITY (logged): the task's 'stored exact 13^3 masses' -- the stored
    artifact 05-knowledge/results/lrc14_owner_loop_drift_typed_row_opus_
    20260727.out stores the 169 successor-mass numerators S_E(s,t)*T_DEN,
    not the 13^3 table.  GATE therefore has two parts: (A) the full 13^3
    bare table is INDEPENDENTLY recomputed here with the audited companion
    interval machinery (copied verbatim) and the chamber reconstruction
    must equal it exactly on all 2197 cells; (B) its row sums must equal
    the byte-stored S_E numerators.  (A)+(B) is strictly stronger than the
    stated gate.
  * AMBIGUITY (logged): 'aggregate chamber contributions that are rational
    multiples of the same zeta-power' -- implemented as the natural
    generalisation: aggregate by primitive direction in Z[zeta_13]
    (integer coefficient vector reduced mod Phi_13 to last-coordinate-0
    form, divided by the gcd).  Two contributions are positive-rational
    multiples of one another iff their primitive directions coincide; this
    includes all same-zeta-power aggregations.
  * Embedding: the standard archimedean embedding zeta = exp(2 pi i/13).
    Character (k*lam, k*mu) under the standard embedding equals character
    (lam,mu) under the embedding zeta -> zeta^k, so running all 168
    characters covers every (character, embedding) pair Galois-equivariantly.
  * No floats in decisions: zero tests are algebraic (a vector in Z^13
    maps to 0 in Z[zeta_13] iff its coordinates are all equal); signs of
    provably-nonzero real cyclotomic values are decided by fixed-point
    integer arithmetic at scale 10^140.  Exact rational alternating-series
    intervals certify every stored sine/cosine value to error < 10^6 units.
    Sign evaluation first replaces the input by its conjugate-difference or
    conjugate-sum vector, so its rounding error is at most B*10^6 for the
    same B used in the algebraic-integer norm bound |value| >= 1/B^11.
    Floats appear only in descriptive prints.

SCOPE (MISTAKE-281 / THM-2541 guardrail): the row is TYPED, not an asserted
scalar cover; no physical current or lawful ancestry/semantic intertwiner is
claimed; no scalar row is excluded; no LRC(14) progress beyond what the cited
theorems literally license; full support is never noncancellation.

Script: 04-computation/lrc14_phase_cone_sidecar_opus_20260727.py
Output: 05-knowledge/results/lrc14_phase_cone_sidecar_opus_20260727.out
"""

import os
import sys
import time
from bisect import bisect_left, bisect_right
from fractions import Fraction
from functools import cmp_to_key
from math import gcd

import numpy as np

T_START = time.time()


def log(msg=""):
    print(msg, flush=True)


def elapsed():
    return f"[{time.time() - T_START:8.1f}s]"


def require(cond, msg):
    if not cond:
        raise AssertionError(msg)


# --------------------------------------------------------------------------
# Row data and scales (identical to the audited drift companion)
# --------------------------------------------------------------------------
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
UNIT_IDX = (1, 2, 3, 4, 5)
C3 = W[TB]

LCM_W = 1
for v in W:
    LCM_W = LCM_W * v // gcd(LCM_W, v)
T_DEN = 182 * LCM_W
require(T_DEN == 297836897838480, "T_DEN anchor")
D2 = 2 * T_DEN

L0 = LCM_W // C3
PP_T = T_DEN // (13 * C3)
HPROBE_T = 13 * L0
HSUCC_T = L0
require(PP_T == 14 * L0 and T_DEN == 13 * C3 * PP_T, "comb grid exact")

AJ, AA, AC = 1, 169, 57122
require(13 * AJ == W[OWNER] and 13 * AA == W[TA] and 13 * AC == W[TB],
        "THM-2368 (8) blocker scales")
KA, KC = 14, 27
require(KA % 13 == 1 and KC % 13 == 1, "grafts are 1 mod 13")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DRIFT_OUT = os.path.join(
    REPO, "05-knowledge", "results",
    "lrc14_owner_loop_drift_typed_row_opus_20260727.out")


# --------------------------------------------------------------------------
# Audited interval machinery (verbatim from the drift companion)
# --------------------------------------------------------------------------
def in_comb(i, elli):
    wsp = W[i]
    U = T_DEN // (182 * wsp)
    lo = (-13 - 14 * elli) % 182
    out = []
    for n in range(wsp):
        s = (lo + 182 * n) * U
        t = s + 26 * U
        if t <= T_DEN:
            out.append((s, t))
        else:
            out.append((s, T_DEN))
            out.append((0, t - T_DEN))
    out.sort()
    return out


def subtract_comb(iv, wsp, PD, lo, hi):
    U = T_DEN // (PD * wsp)
    step = PD * U
    lenW = (hi - lo) * U
    base = (lo % PD) * U
    out = []
    ap = out.append
    for A, B in iv:
        p0 = A - lenW + 1
        k0 = -((base - p0) // step)
        p = base + k0 * step
        cur = A
        while p < B:
            we = p + lenW
            if we > cur:
                if p > cur:
                    ap((cur, p))
                cur = we
                if cur >= B:
                    break
            p += step
        if cur < B:
            ap((cur, B))
    return out


def build_set(pattern, ell):
    ins = [i for i, mmode in pattern.items() if mmode == "in"]
    start = min(ins, key=lambda i: W[i])
    iv = in_comb(start, ell[start])
    for i, mmode in pattern.items():
        if mmode == "gout":
            iv = subtract_comb(iv, W[i], 91, -13 - 7 * ell[i], 13 - 7 * ell[i])
    rest = sorted((W[i], i) for i, mmode in pattern.items()
                  if i != start and mmode in ("in", "out"))
    for _, i in rest:
        e = ell[i]
        if pattern[i] == "out":
            iv = subtract_comb(iv, W[i], 182, -13 - 14 * e, 13 - 14 * e)
        else:
            iv = subtract_comb(iv, W[i], 182, 13 - 14 * e, 169 - 14 * e)
    return iv


def check_intervals(iv, top):
    last = -1
    for a, b in iv:
        require(0 <= a < b <= top and a >= last, "interval list corrupt")
        last = b
    return sum(b - a for a, b in iv)


def accumulate_interval(a, b, PP, hp, hs, probe_acc):
    succ = 0
    c = (a - hp) // PP + 1
    cend = (b + hp) // PP
    while c <= cend:
        ctr = c * PP
        lo = ctr - hp
        hi = ctr + hp
        ov = (b if b < hi else hi) - (a if a > lo else lo)
        if ov > 0:
            probe_acc[c % 13] += ov
            lo2 = ctr - hs
            hi2 = ctr + hs
            ov2 = (b if b < hi2 else hi2) - (a if a > lo2 else lo2)
            if ov2 > 0:
                succ += ov2
        c += 1
    return succ


PAT_E = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
         6: "in", 7: "out", 8: "out"}
ZELL = (0,) * 9


def twist_shift_vector(s, t):
    e = [0] * 9
    e[1] = s % 13
    e[2] = t % 13
    e[TA] = (-s) % 13
    e[TB] = (-t) % 13
    return tuple(e)


# --------------------------------------------------------------------------
# Fixed-point trig table (scale 10^140) with validated error budget
# --------------------------------------------------------------------------
PREC = 140
SCALE = 10 ** PREC
ERR_UNITS = 10 ** 6          # conservative per-value fixed-point error bound


def _atan_inv(x):
    total = 0
    p = SCALE // x
    n = 1
    sgn = 1
    xx = x * x
    while p:
        total += sgn * (p // n)
        p //= xx
        n += 2
        sgn = -sgn
    return total


def _fsin(x):
    term = x
    s = x
    n = 1
    sgn = 1
    while term:
        term = term * x // SCALE * x // SCALE // ((n + 1) * (n + 2))
        n += 2
        sgn = -sgn
        s += sgn * term
    return s


def _fcos(x):
    term = SCALE
    c = SCALE
    n = 0
    sgn = 1
    while term:
        term = term * x // SCALE * x // SCALE // ((n + 1) * (n + 2))
        n += 2
        sgn = -sgn
        c += sgn * term
    return c


def build_trig():
    pi_fp = 16 * _atan_inv(5) - 4 * _atan_inv(239)
    SIN = []
    COS = []
    for j in range(13):
        th = 2 * pi_fp * j // 13
        SIN.append(_fsin(th))
        COS.append(_fcos(th))
    tol = 10 ** (PREC - 130)  # 10^10 units, vastly above the true error
    for j in range(13):
        require(abs(SIN[j] * SIN[j] + COS[j] * COS[j] - SCALE * SCALE)
                < 3 * SCALE * tol, "trig table: pythagoras")
        require(abs(SIN[j] + SIN[(13 - j) % 13]) < tol
                if j else SIN[0] == 0, "trig table: sin symmetry")
        require(abs(COS[j] - COS[(13 - j) % 13]) < tol, "trig: cos symmetry")
    require(abs(sum(SIN)) < 13 * tol and abs(sum(COS)) < 13 * tol,
            "trig table: root-of-unity sums")
    return SIN, COS


SIN13, COS13 = build_trig()


CERT_SCALE = 10 ** (PREC + 20)
CERT_UP = 10 ** 20


def _ceildiv(a, b):
    return -((-a) // b)


def _atan_interval_fp(q, n_terms=128):
    """Outward-rounded fixed-point enclosure of atan(1/q)."""
    lo = hi = 0
    for n in range(n_terms):
        den = (2 * n + 1) * q ** (2 * n + 1)
        tl = CERT_SCALE // den
        th = _ceildiv(CERT_SCALE, den)
        if n % 2 == 0:
            lo += tl
            hi += th
        else:
            lo -= th
            hi -= tl
    n = n_terms
    den = (2 * n + 1) * q ** (2 * n + 1)
    th = _ceildiv(CERT_SCALE, den)
    # The alternating remainder has the sign of the next term and magnitude
    # at most that term.  Widen only on that side.
    if n % 2 == 0:
        hi += th
    else:
        lo -= th
    return lo, hi


def _alt_series_interval_fp(x, sine, n_terms=72):
    """Outward-rounded enclosure of sin/cos(x), with x fixed-point exact.

    Here 0 <= x < pi/2.  By the chosen truncation the omitted alternating
    tail is decreasing, so it lies between zero and the next signed term.
    """
    require(0 <= x < 2 * CERT_SCALE, "acute angle out of range")
    den_scale = CERT_SCALE * CERT_SCALE
    if sine:
        tl = th = x
    else:
        tl = th = CERT_SCALE
    lo = tl
    hi = th
    for n in range(1, n_terms):
        den = (2 * n) * (2 * n + 1) if sine else (2 * n - 1) * (2 * n)
        full_den = den_scale * den
        tl = tl * x * x // full_den
        th = _ceildiv(th * x * x, full_den)
        if n % 2 == 0:
            lo += tl
            hi += th
        else:
            lo -= th
            hi -= tl
    n = n_terms
    den = (2 * n) * (2 * n + 1) if sine else (2 * n - 1) * (2 * n)
    full_den = den_scale * den
    th = _ceildiv(th * x * x, full_den)
    if n % 2 == 0:
        hi += th
    else:
        lo -= th
    return lo, hi


def certify_trig_error():
    """Rigorously certify the fixed trig table against exact intervals.

    Machin's formula encloses pi.  Symmetry reduces every 13th-root angle to
    an acute angle, where monotonicity plus alternating Taylor bounds gives a
    rational interval.  The whole interval must lie within ERR_UNITS/SCALE
    of the integer table entry; no numerical transcendental library enters.
    """
    a5_lo, a5_hi = _atan_interval_fp(5)
    a239_lo, a239_hi = _atan_interval_fp(239)
    pi_lo = 16 * a5_lo - 4 * a239_hi
    pi_hi = 16 * a5_hi - 4 * a239_lo
    require(pi_lo < pi_hi, "pi enclosure corrupt")
    eps = ERR_UNITS * CERT_UP

    for j in range(13):
        if j == 0:
            sin_iv = (0, 0)
            cos_iv = (CERT_SCALE, CERT_SCALE)
        else:
            r = min(j, 13 - j)
            if r <= 3:
                alo = (2 * r * pi_lo) // 13
                ahi = _ceildiv(2 * r * pi_hi, 13)
                cos_sgn = 1
            else:
                k = 13 - 2 * r
                alo = (k * pi_lo) // 13
                ahi = _ceildiv(k * pi_hi, 13)
                cos_sgn = -1
            require(0 <= alo <= ahi and 2 * ahi < pi_lo,
                    "acute-angle reduction failed")
            sl, _ = _alt_series_interval_fp(alo, True)
            _, su = _alt_series_interval_fp(ahi, True)
            cl, _ = _alt_series_interval_fp(ahi, False)
            _, cu = _alt_series_interval_fp(alo, False)
            if j > 6:
                sin_iv = (-su, -sl)
            else:
                sin_iv = (sl, su)
            cos_iv = (cl, cu) if cos_sgn > 0 else (-cu, -cl)

        sa = SIN13[j] * CERT_UP
        ca = COS13[j] * CERT_UP
        require(sa - eps <= sin_iv[0] <= sin_iv[1] <= sa + eps,
                f"rigorous sine error certificate failed at j={j}")
        require(ca - eps <= cos_iv[0] <= cos_iv[1] <= ca + eps,
                f"rigorous cosine error certificate failed at j={j}")


certify_trig_error()


# --------------------------------------------------------------------------
# Exact Z[zeta_13] helpers.  Elements are integer 13-vectors of exponent
# coefficients; the evaluation kernel is spanned by the all-ones vector.
# --------------------------------------------------------------------------
def conv13(u, v):
    out = [0] * 13
    for i, ui in enumerate(u):
        if ui:
            for j, vj in enumerate(v):
                if vj:
                    out[(i + j) % 13] += ui * vj
    return out


def conj13(u):
    return [u[(-i) % 13] for i in range(13)]


def elem_is_zero(u):
    f = u[0]
    return all(x == f for x in u)


def _checked_sign(val, B, what):
    # val approximates a nonzero real/imaginary algebraic integer alpha;
    # B bounds every conjugate and the fixed-point rounding error by
    # B*ERR_UNITS.  Norm(alpha) is a nonzero integer of degree at most 12,
    # so |alpha| >= 1/B^11.  The exponent 12 below makes the separation
    # dominate the B-weighted table error.
    require(B > 0, f"{what}: zero coefficient bound")
    require(SCALE > 4 * (B + 1) ** 12 * ERR_UNITS,
            f"{what}: precision insufficient for norm bound")
    bound = B * ERR_UNITS
    require(abs(val) > bound, f"{what}: fixed-point dead zone (impossible)")
    return 1 if val > 0 else -1


def im_sign(u):
    """Exact sign of Im(u) under zeta = exp(2 pi i/13)."""
    w = [u[j] - u[(-j) % 13] for j in range(13)]
    if elem_is_zero(w):
        return 0
    # Im(w)=2 Im(u), and using w makes the certified rounding error depend
    # on exactly B=sum|w_j| rather than on an arbitrary kernel representative.
    val = sum(w[j] * SIN13[j] for j in range(13))
    return _checked_sign(val, sum(abs(x) for x in w), "im_sign")


def re_sign(u):
    """Exact sign of Re(u)."""
    w = [u[j] + u[(-j) % 13] for j in range(13)]
    if elem_is_zero(w):
        return 0
    # Re(w)=2 Re(u), with the same representative-safe error accounting.
    val = sum(w[j] * COS13[j] for j in range(13))
    return _checked_sign(val, sum(abs(x) for x in w), "re_sign")


def cross_sign(u, v):
    """sign of Im(conj(u)*v):  >0  iff v is CCW of u (within pi)."""
    return im_sign(conv13(conj13(u), v))


def dot_sign(u, v):
    return re_sign(conv13(conj13(u), v))


def reduce_prim(vec):
    """Canonical primitive direction of a Z[zeta13] element given by an
    exponent-coefficient vector: subtract coord 12 (mod Phi13 normal form
    with last coordinate 0), divide by gcd.  Returns (primitive, gcd), or
    None if the element is zero.  Two chamber contributions are
    positive-rational multiples of each other iff their reduced primitives
    coincide; the element equals gcd * primitive exactly."""
    b = vec[12]
    red = [x - b for x in vec]
    g = 0
    for x in red:
        g = gcd(g, abs(x))
    if g == 0:
        return None
    return tuple(x // g for x in red), g


def rot13(u, k):
    """u * zeta^k."""
    return [u[(j - k) % 13] for j in range(13)]


# --------------------------------------------------------------------------
# Step 1: rational y-breakpoints of all rooted and blocker words (34)
# --------------------------------------------------------------------------
def gen_breakpoints():
    """Returns (all_points_sorted_unique (np.int64 circle grid T_DEN),
    slow_points set, counts dict).  Families:
      f1 owner blocker d(y);  f2 target-a blocker g(169y - s/13), 13 shifts;
      f3 target-c blocker/probe word g/d(57122y - u/13), 13 shifts (the two
         share every breakpoint);  f4 rooted A_y(r);  f5 rooted C_y(r);
      f6 rooted guard gamma((y+h)/13);  f7 rooted units q in {40,53,66}."""
    counts = {}
    slow = []
    # f1
    f1 = [T_DEN // 14, 13 * T_DEN // 14]
    counts["f1 owner d(y)"] = len(f1)
    slow += f1
    # f2
    K2 = T_DEN // (182 * 169)
    f2 = []
    for s in range(13):
        for n in range(169):
            base = 182 * n + 14 * s
            f2.append((K2 * (base + 13)) % T_DEN)
            f2.append((K2 * (base - 13)) % T_DEN)
    counts["f2 target-a blocker (13 shifts)"] = len(f2)
    slow += f2
    # f3 (fast; generated as numpy)
    K3 = T_DEN // (182 * 57122)
    n_arr = np.arange(57122, dtype=np.int64) * 182
    f3_parts = []
    for u in range(13):
        for sgn in (13, -13):
            f3_parts.append((n_arr + (14 * u + sgn)) * K3 % T_DEN)
    f3 = np.concatenate(f3_parts)
    counts["f3 target-c blocker+probe (13 shifts, shared)"] = int(f3.size)
    # f4
    K4 = T_DEN // 196
    f4 = []
    for r in range(13):
        for n in range(-1, 3):
            for sgn in (13, -13):
                v = 182 * n - 14 * r + sgn
                if 0 <= v < 196:
                    f4.append(K4 * v)
    counts["f4 rooted A_y"] = len(f4)
    slow += f4
    # f5
    K5 = T_DEN // 378
    f5 = []
    for r in range(13):
        for n in range(-1, 4):
            for sgn in (13, -13):
                v = 182 * n - 14 * r + sgn
                if 0 <= v < 378:
                    f5.append(K5 * v)
    counts["f5 rooted C_y"] = len(f5)
    slow += f5
    # f6
    G7 = 13 * T_DEN // 7
    f6 = []
    for h in range(13):
        for n in range(0, 2):
            for sgn in (1, -1):
                v = T_DEN * (13 * n - h) + sgn * G7
                if 0 <= v < T_DEN:
                    f6.append(v)
    counts["f6 rooted guard"] = len(f6)
    slow += f6
    # f7
    f7 = []
    for q in (40, 53, 66):
        K7 = T_DEN // (14 * q)
        for h in range(13):
            for n in range((q * h) // 13 - 1, (q * (h + 1)) // 13 + 2):
                for sgn in (13, -13):
                    v = K7 * (182 * n + sgn) - h * T_DEN
                    if 0 <= v < T_DEN:
                        f7.append(v)
    counts["f7 rooted units 40/53/66"] = len(f7)
    slow += f7

    allpts = np.unique(np.concatenate(
        [f3, np.array(slow, dtype=np.int64)]))
    return allpts, set(slow), counts


# --------------------------------------------------------------------------
# Step 1b: exact word evaluation at chamber midpoints (grid D2)
# --------------------------------------------------------------------------
P13D2 = 13 * D2
TH14 = P13D2 // 14
TH7 = P13D2 // 7
D2_14 = D2 // 14
D2_13 = D2 // 13


def _danger(v, per, th):
    require(v != th and v != per - th, "midpoint hit a word boundary")
    return v < th or v > per - th


def rooted_snapshot(m):
    """(a_bits, c_bits, p_bits) of the rooted words at y = m/D2."""
    a = 0
    c = 0
    p = 0
    for r in range(13):
        if not _danger((KA * m + r * D2) % P13D2, P13D2, TH14):
            a |= 1 << r
        if not _danger((KC * m + r * D2) % P13D2, P13D2, TH14):
            c |= 1 << r
    for h in range(13):
        base = m + h * D2
        if _danger(base % P13D2, P13D2, TH7):
            continue                      # gamma = 0 (near-integer danger)
        ok = True
        for q in (40, 53, 66):
            if _danger((q * base) % P13D2, P13D2, TH14):
                ok = False
                break
        if ok:
            p |= 1 << h
    return (a, c, p)


def ga_bits(m):
    ga = 0
    v0 = (AA * m) % D2
    for s in range(13):
        if not _danger((v0 - s * D2_13) % D2, D2, D2_14):
            ga |= 1 << s
    return ga


def dc_bits(m):
    dc = 0
    v0 = (AC * m) % D2
    for r in range(13):
        if _danger((v0 - r * D2_13) % D2, D2, D2_14):
            dc |= 1 << r
    return dc


def owner_d(m):
    v = m % D2
    return _danger(v, D2, D2_14)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    log("THM-2368 Sec 8 phase/event sidecar + phase-cone certificate,"
        " bare tensor H_E, canonical typed row")
    log("script=04-computation/lrc14_phase_cone_sidecar_opus_20260727.py")
    log("machine=opus  date=2026-07-27")
    log("")
    log("[1] conventions (see module docstring; files override summaries)")
    log(f"    w={W}, grafts k_a=q1=14, k_c=q2=27 (both = 1 mod 13)")
    log(f"    blocker scales A_j={AJ}, A_a={AA}, A_c={AC} (THM-2368 (8))")
    log("    gc[t] = 1 - dc[t]: c3 safe factor and probe share one word")
    log(f"    grids: T_DEN={T_DEN}, midpoints on D2=2*T_DEN")
    log("    fixed cyclotomic sign table: exact rational interval error"
        " certificate PASS")
    log("    walk per character (lam,mu)!=(0,0): the a=0 deep-colour modes"
        " 13*T_DEN*sum_{s,t} S_E zeta^(-lam s - mu t)")
    log("    (circulant tensors kill every such mode => nonzero endpoint"
        " certifies D_(H_E) > 0 via THM-2365 (17)-(18))")
    log("")

    # ------------------------------------------------------------------
    log("[2] chamber partition (THM-2368 (34))")
    allpts, slowset, counts = gen_breakpoints()
    tot_gen = sum(counts.values())
    for k, v in counts.items():
        log(f"    {k}: {v} generated breakpoints")
    n_chambers = int(allpts.size)
    log(f"    generated (with multiplicity): {tot_gen}")
    log(f"    distinct breakpoints = chambers of the circle partition:"
        f" {n_chambers}")
    LJ = 13 * T_DEN // 14
    RJ = 15 * T_DEN // 14
    # build the support sweep list (wrapped interval J = (13/14, 15/14))
    pts_list = allpts.tolist()
    lo_i = bisect_right(pts_list, LJ)
    hi_wrap = T_DEN // 14
    inner = [p for p in pts_list[lo_i:] if p < T_DEN] \
        + [p + T_DEN for p in pts_list[:bisect_left(pts_list, hi_wrap)]]
    sweep = [LJ] + inner + [RJ]
    for i in range(len(sweep) - 1):
        require(sweep[i] < sweep[i + 1], "sweep not strictly increasing")
    n_sup = len(sweep) - 1
    log(f"    owner-blocker support supp d(y) = (13/14, 15/14) (wrapped);"
        f" contributing chambers: {n_sup}")
    log(f"    {elapsed()} partition done")
    log("")

    # ------------------------------------------------------------------
    log("[3] single exact 1-D sweep over the support (integer arithmetic)")
    states = {}
    chamber_seq = []
    rooted_ids = {}
    rooted_list = []
    cur_rid = -1
    cur_ga = -1
    slow_in_J = set()
    for p in slowset:
        slow_in_J.add(p)
        slow_in_J.add(p + T_DEN)
    first = True
    for i in range(n_sup):
        x0 = sweep[i]
        x1 = sweep[i + 1]
        m = x0 + x1
        if first or x0 in slow_in_J:
            require(owner_d(m), "chamber outside owner support in sweep")
            snap = rooted_snapshot(m)
            rid = rooted_ids.get(snap)
            if rid is None:
                rid = len(rooted_list)
                rooted_ids[snap] = rid
                rooted_list.append(snap)
            cur_rid = rid
            cur_ga = ga_bits(m)
            first = False
        key = (cur_rid, cur_ga, dc_bits(m))
        ln = x1 - x0
        if key in states:
            states[key] += ln
        else:
            states[key] = ln
        chamber_seq.append((x0, key))
    require(sum(states.values()) == T_DEN // 7, "support length mismatch")
    log(f"    chambers walked: {n_sup}")
    log(f"    distinct rooted snapshots (A,C,P): {len(rooted_list)}")
    log(f"    distinct sidecar states (rooted, ga, dc): {len(states)}")
    log(f"    aggregation stage 1 (logged): chambers with identical"
        f" (rooted,ga,dc) word data merged: {n_sup} -> {len(states)}")
    log(f"    {elapsed()} sweep done")
    log("")

    # ------------------------------------------------------------------
    log("[4] rooted N tables and per-state constant chamber tables")
    N_tabs = []
    for (a, c, p) in rooted_list:
        av = [(a >> r) & 1 for r in range(13)]
        cv = [(c >> r) & 1 for r in range(13)]
        pv = [(p >> h) & 1 for h in range(13)]
        require(0 < sum(av) < 13 and 0 < sum(cv) < 13, "A/C not proper")
        require(3 <= sum(pv) <= 10, "P support outside THM-2368 (14)")
        N = np.zeros((13, 13), dtype=np.int64)
        for s in range(13):
            for t in range(13):
                N[s, t] = sum(av[(h + s) % 13] * cv[(h + t) % 13] * pv[h]
                              for h in range(13))
        N_tabs.append(N)
    log(f"    all rooted words nonempty proper; |supp P| in [3,10]"
        f" (THM-2368 (12)/(14)): PASS ({len(rooted_list)} snapshots)")

    # per-character scatter matrix: row ci*13+k picks (s,t) with
    # (lam*s + mu*t + k) = 0 mod 13
    CHARS = [(lam, mu) for lam in range(13) for mu in range(13)
             if (lam, mu) != (0, 0)]
    SCAT = np.zeros((168 * 13, 169), dtype=np.int64)
    for ci, (lam, mu) in enumerate(CHARS):
        for j in range(169):
            s, t = divmod(j, 13)
            k = (-(lam * s + mu * t)) % 13
            SCAT[ci * 13 + k, j] = 1

    Hnum = np.zeros((13, 13, 13), dtype=np.int64)   # scale 13*T_DEN
    dirs = [dict() for _ in range(168)]             # primitive dir -> weight
    state_dir_idx = {}                              # key -> per-char prims
    for key, ln in states.items():
        rid, ga, dc = key
        gav = np.array([(ga >> s) & 1 for s in range(13)], dtype=np.int64)
        dcv = np.array([(dc >> r) & 1 for r in range(13)], dtype=np.int64)
        nd = int(dcv.sum())
        Y = gav[:, None] * (1 - dcv)[None, :] * N_tabs[rid]
        w = ln * nd
        for r in range(13):
            if dcv[r]:
                Hnum[r] += ln * Y
        Vall = (SCAT @ Y.reshape(169)).reshape(168, 13)
        prims = []
        for ci in range(168):
            pr = reduce_prim(Vall[ci].tolist())
            if pr is None:
                prims.append(None)
                continue
            prim, g = pr
            prims.append(prim)
            if w > 0:
                d = dirs[ci]
                wg = w * g          # contribution = (len*nd*g) * primitive
                if prim in d:
                    d[prim] += wg
                else:
                    d[prim] = wg
        state_dir_idx[key] = prims
    log(f"    per-state constant tables built; H reconstruction and 168"
        f" walks accumulated")
    log(f"    {elapsed()} state processing done")
    log("")

    # ------------------------------------------------------------------
    log("[5] GATE (MISTAKE-281 common-base discipline)")
    log("    (A) independent recomputation of the full 13^3 bare table with"
        " the audited drift-companion")
    log("        interval machinery (verbatim copy), all 169 twists:")
    E0 = build_set(PAT_E, ZELL)
    lenE0 = check_intervals(E0, T_DEN)
    require(Fraction(lenE0, T_DEN) == Fraction(1882176, 28589561),
            "E_1 measure anchor")
    log(f"        E_1 measure anchor 1882176/28589561: PASS")
    HE_tab = {}
    SE_comp = {}
    for s in range(13):
        for t in range(13):
            iv = build_set(PAT_E, twist_shift_vector(s, t))
            lenE = 0
            HE = [0] * 13
            succ = 0
            for A, B in iv:
                lenE += B - A
                succ += accumulate_interval(A, B, PP_T, HPROBE_T, HSUCC_T,
                                            HE)
            require(sum(HE) == 2 * lenE - succ,
                    f"(19b) identity fails at {(s, t)}")
            for r in range(13):
                HE_tab[(r, s, t)] = HE[r]
            SE_comp[(s, t)] = sum(HE)
    log(f"        (19b) successor identity at all 169 twists: PASS")
    log(f"        {elapsed()} recomputation done")
    gate_a = all(int(Hnum[r, s, t]) == 13 * HE_tab[(r, s, t)]
                 for r in range(13) for s in range(13) for t in range(13))
    log(f"    (A) chamber reconstruction == independent 13^3 table"
        f" (x13 scale), all 2197 cells: {'PASS' if gate_a else 'FAIL'}")
    require(gate_a, "GATE (A) failed -- fix construction before"
                    " interpretation")
    # (B) stored artifact
    with open(DRIFT_OUT, "r", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    idx = lines.index("    S_E(s,t)*T_DEN:")
    SE_stored = {}
    for s in range(13):
        row = [int(x) for x in lines[idx + 1 + s].split()]
        require(len(row) == 13, "stored S_E row malformed")
        for t in range(13):
            SE_stored[(s, t)] = row[t]
    gate_b = all(SE_comp[(s, t)] == SE_stored[(s, t)]
                 and int(Hnum[:, s, t].sum()) == 13 * SE_stored[(s, t)]
                 for s in range(13) for t in range(13))
    log(f"    (B) row sums == byte-stored S_E numerators"
        f" (lrc14_owner_loop_drift_typed_row_opus_20260727.out):"
        f" {'PASS' if gate_b else 'FAIL'}")
    require(gate_b, "GATE (B) failed")
    log("    NOTE (ambiguity logged): the stored artifact holds the 169"
        " successor numerators, not the")
    log("    13^3 table; gate (A) supplies the full-table equality against"
        " an in-script independent")
    log("    reconstruction by the audited machinery, strictly stronger"
        " than the stated gate.")
    diag_ok = all(int(Hnum[t, s, t]) == 0
                  for s in range(13) for t in range(13))
    log(f"    diagonal-plane zero H_E(t,s,t)=0 (all 169): "
        f"{'PASS' if diag_ok else 'FAIL'}")
    require(diag_ok, "diagonal zero failed")
    log("")

    # ------------------------------------------------------------------
    log("[6] per-character endpoint cross-check (walks vs stored S_E DFT)")
    n_nonzero_end = 0
    endpoint_nonzero = [False] * 168
    for ci, (lam, mu) in enumerate(CHARS):
        acc = [0] * 13
        for prim, wt in dirs[ci].items():
            for j in range(13):
                acc[j] += wt * prim[j]
        ref = [0] * 13
        for s in range(13):
            for t in range(13):
                ref[(-(lam * s + mu * t)) % 13] += 13 * SE_stored[(s, t)]
        da = acc[0] - ref[0]
        require(all(acc[j] - ref[j] == da for j in range(13)),
                f"endpoint mismatch at character {(lam, mu)}")
        nz = not elem_is_zero(acc)
        endpoint_nonzero[ci] = nz
        if nz:
            n_nonzero_end += 1
    log("    walk endpoints == 13 * DFT(stored S_E) as elements of"
        " Q(zeta_13), all 168 characters: PASS")
    log(f"    characters with NONZERO endpoint (a=0 deep-colour mode"
        f" alive): {n_nonzero_end} / 168")
    log(f"    {elapsed()} endpoint checks done")
    log("")

    # ------------------------------------------------------------------
    log("[7] exact phase-cone decision per character (no floats in"
        " decisions)")
    log("    aggregation stage 2 (logged): contributions merged by"
        " primitive direction in Z[zeta_13]")
    log("    (positive-rational-multiple classes; includes every"
        " same-zeta-power merge)")

    def sort_dirs(dlist):
        def half(d):
            i = im_sign(d)
            if i > 0:
                return 0
            if i < 0:
                return 1
            return 0 if re_sign(d) > 0 else 1

        halves = {d: half(d) for d in dlist}

        def cmp(u, v):
            hu, hv = halves[u], halves[v]
            if hu != hv:
                return -1 if hu < hv else 1
            s = cross_sign(list(u), list(v))
            if s > 0:
                return -1
            if s < 0:
                return 1
            return 0
        return sorted(dlist, key=cmp_to_key(cmp))

    def cone_test(dlist):
        """Returns (span_class, info).  span_class in
        {'<pi', '=pi-cert', '=pi-fail', '>pi', 'empty'}.
        For '<pi': info = (start_dir, end_dir) of the minimal closed arc
        (the two extreme contribution directions)."""
        if not dlist:
            return "empty", None
        if len(dlist) == 1:
            return "<pi", (dlist[0], dlist[0])
        sd = sort_dirs(dlist)
        m = len(sd)
        big_gap_at = None
        pi_gaps = []
        for i in range(m):
            u = sd[i]
            v = sd[(i + 1) % m]
            cs = cross_sign(list(u), list(v))
            if cs < 0:
                big_gap_at = i
                break
            if cs == 0 and dot_sign(list(u), list(v)) < 0:
                pi_gaps.append(i)
        if big_gap_at is not None:
            i = big_gap_at
            return "<pi", (sd[(i + 1) % m], sd[i])
        if pi_gaps:
            # span == pi: extremes antipodal; certificate iff some
            # direction is NOT collinear with the extreme line
            i = pi_gaps[0]
            u = sd[i]
            if any(cross_sign(list(u), list(d)) != 0 for d in sd):
                return "=pi-cert", (sd[(i + 1) % m], sd[i])
            return "=pi-fail", None
        return ">pi", None

    span_counts = {"<pi": 0, "=pi-cert": 0, "=pi-fail": 0, ">pi": 0,
                   "empty": 0}
    cone_holders = []
    dir_stats = []
    results = {}
    for ci, (lam, mu) in enumerate(CHARS):
        dlist = list(dirs[ci].keys())
        dir_stats.append(len(dlist))
        cls, info = cone_test(dlist)
        span_counts[cls] += 1
        results[(lam, mu)] = (cls, info, len(dlist))
        if cls in ("<pi", "=pi-cert"):
            cone_holders.append((lam, mu))
    log(f"    distinct aggregated directions per character:"
        f" min={min(dir_stats)} max={max(dir_stats)}"
        f" mean={sum(dir_stats)/168:.1f}")
    log(f"    span statistics over the 168 nonzero target characters:")
    log(f"      span < half turn (cone holds):        {span_counts['<pi']}")
    log(f"      span = half turn, certificate holds:  "
        f"{span_counts['=pi-cert']}")
    log(f"      span = half turn, certificate fails:  "
        f"{span_counts['=pi-fail']}")
    log(f"      span > half turn (cone fails):        {span_counts['>pi']}")
    log(f"      no nonzero contributions:             {span_counts['empty']}")
    log(f"    {elapsed()} cone tests done")
    log("")

    # ------------------------------------------------------------------
    log("[8] DECISION")
    if cone_holders:
        log("    PHASE-CONE CERTIFICATE HOLDS for at least one nonzero"
            " target character.")
        log(f"    certified characters ({len(cone_holders)}):"
            f" {cone_holders}")
        # certifying units
        shown = 0
        for (lam, mu) in cone_holders:
            cls, info, ndir = results[(lam, mu)]
            dlist = list(dirs[CHARS.index((lam, mu))].keys())
            unit = None
            for k in range(13):
                for sg in (1, -1):
                    sig = [re_sign([sg * x for x in rot13(list(d), k)])
                           for d in dlist]
                    if all(s >= 0 for s in sig) and any(s > 0 for s in sig):
                        unit = (sg, k)
                        break
                if unit:
                    break
            start, end = info
            if unit:
                sgs = "" if unit[0] == 1 else "-"
                log(f"    character (lam,mu)=({lam},{mu}): certifying unit"
                    f" u = {sgs}zeta^{unit[1]} (exact 26th-root sector);"
                    f" {ndir} directions;")
            else:
                log(f"    character (lam,mu)=({lam},{mu}): no 26th-root"
                    f" unit certifies; exact sector = closed arc from")
            log(f"        extreme directions {start} .. {end}"
                f" (reduced primitive Z[zeta_13] vectors; arc"
                f" {'< pi' if cls == '<pi' else '= pi'})")
            shown += 1
            if shown >= 12:
                log(f"    ... ({len(cone_holders) - shown} further certified"
                    f" characters suppressed in listing)")
                break
        log("")
        log("    CONSEQUENCE: for each certified character the walk endpoint"
            " is a nonzero a=0 deep-colour")
        log("    nonzero-target mode of H_E, so H_E is NOT circulant and"
            " THM-2365 (17)-(18) give D_(H_E) > 0.")
        log("    This is the first MECHANISM-LEVEL positive-drift"
            " certificate: the (37) sidecar phases of the")
        log("    certified character(s) fit one open half-plane, the anchor"
            " the uniform THM-2368 (38) argument needs.")
        log("    (Consistent with THM-2550's certified D_(H_E) ="
            " 9497453727128823229273/7328279741345331835978099188.)")
    else:
        log("    PHASE-CONE CERTIFICATE FAILS for all 168 nonzero target"
            " characters.")
        log("    Combined with THM-2550's certified D_(H_E) > 0: the"
            " canonical bare drift is POSITIVE but NOT")
        log("    half-plane-explainable at the a=0 target-character level:"
            " every nonzero mode has")
        log("    non-semicircular direction support but a nonzero endpoint."
            "  No winding claim is made.")
        # First cone-breaking event for the first character with nonzero
        # endpoint.  This tracks feasibility of the accumulated distinct
        # direction set; it is not a chronological winding-number test.
        ci0 = next(i for i in range(168) if endpoint_nonzero[i])
        lam, mu = CHARS[ci0]
        log(f"    first cone-breaking phase event (Sec 8), character"
            f" (lam,mu)=({lam},{mu}):")
        seen = []
        seen_set = set()
        event = None
        for (x0, key) in chamber_seq:
            prim = state_dir_idx[key][ci0]
            if prim is None or prim in seen_set:
                continue
            seen.append(prim)
            seen_set.add(prim)
            cls, _ = cone_test(seen)
            if cls in ("=pi-fail", ">pi"):
                event = (x0, prim)
                break
        if event is None:
            log("      (walk stays cone-feasible chamber-by-chamber; the"
                " failure is at the closed-arc boundary)")
        else:
            x0, prim = event
            yfrac = Fraction(x0 % T_DEN, T_DEN)
            edge_gap = abs(yfrac - Fraction(13, 14))
            log(f"      chamber left endpoint y = {yfrac} (grid"
                f" {x0 % T_DEN}/T_DEN)")
            log(f"      exact gap from owner-window edge 13/14 = {edge_gap}")
            log(f"      first direction making the accumulated set fail the"
                f" phase-cone certificate: {prim}")
    log("")

    # ------------------------------------------------------------------
    log("[9] scope (MISTAKE-281 / THM-2541 guardrails, unchanged)")
    log("    * The row is the canonical TYPED row THM-2309 (25): typed, NOT"
        " an asserted scalar cover.")
    log("    * No physical current is claimed; no row is excluded; the"
        " scalar ledger stays 165;")
    log("      LRC(14) remains OPEN.  Full support is never read as"
        " noncancellation.")
    log("    * The certificate (when it holds) licenses exactly: this row's"
        " bare tensor is non-circulant")
    log("      with a mechanism-level (half-plane) explanation of one (or"
        " more) surviving modes -- the")
    log("      anchor named by THM-2368 (38); nothing more.")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
