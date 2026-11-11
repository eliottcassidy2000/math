#!/usr/bin/env python3
"""S91 bridges (3): the F1 data of the owner's numbers 2, 3, 5, 7, 13, 19, 27, 1039 under the local-density theory
of this session (CONTEXT_S91.md), computed with own code.

Notation: psi, q = h - psi, x = 2 psi + 1 = y - 2q, y = 2h + 1; F1: N = x1 x2 + 2^s y1 y2 - [s = 0], ordered triples;
sheet 0: s = 0 (distinguished class N = -1 mod p), sheet 1: s >= 1 (distinguished class N = 0 mod p).
A state h is p-divisible when p | gcd(x(h), y(h)) = gcd(2h + 1, q(h)).

Sections.
  A  r(N), its split by sheets and the one-apex words, at 2, 3, 5, 7, 13, 19, 21, 27, 43, 67, 163, 1039 (complete search).
  B  p-divisible states below 2^26 for the moduli 3, 5, 7, 13, 19, 27, 1039: first state, count, count * p^2 / 2^26;
     diagonal states x = y = 2^(j+1) + 1 at h = 2^j;  p has a diagonal divisible state iff ord_p(2) is even.
  C  rates max_k rho(T_k)/2 of the twisted digit operators: own construction, compared with ../twisted_spectrum.out
     for m = 3, 5, 7, 13, 19, 27; for m = 1039 (not in that table) every character class is bounded by
     ||T^512||^(1/512)/2 (power iteration) and the largest candidates are computed by dense eigenvalues.
     Floating point: status VERIFIED, not exact.
  D  first-order prediction E_pred = p D/(1 - D), D = 2f - 1/p^2, on odd N in [2^24, 2^27): own implementation of the
     weights W(h1, s); compared with ../prereg_window_27_24.json for p = 3, 5, 7, 13, 19; NEW for p = 1039
     (written to b3_prediction_1039.txt BEFORE the window is enumerated by this script).
  E  own enumeration of the window: measured excess at 3, 5, 7, 13, 19 (compared with the session's
     window_split_prereg_27_24.out) and at 1039 (tested against D with the session's control rule
     |E - E_pred| <= 0.008 + 4 se, se = sqrt(0.75 p / #odd N)).
  F  limiting singular series S_1(N) = (8/pi^2) prod_{odd p | N} (1 + p/(p^2 - 1)), S_0(N) = S_1(N + 1) (odd part),
     at the eight numbers and at 7, 21, 43, 67, 163; and the misses of F1 among odd N < 4096 against S.

Usage: python3 b3_f1_data.py [K LO]      (window [2^LO, 2^K), default 27 24; use 21 18 for a quick run: then the
                                          comparisons with the session's files for [2^24, 2^27) are skipped)
"""
import json
import math
import os
import re
import sys

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""), flush=True)


K = int(sys.argv[1]) if len(sys.argv) > 2 else 27
LO = int(sys.argv[2]) if len(sys.argv) > 2 else 24
FULL = (K, LO) == (27, 24)
A, B = 1 << LO, 1 << K
MODULI = np.array([3, 5, 7, 13, 19, 27, 1039], dtype=np.int64)
OWNER = [2, 3, 5, 7, 13, 19, 27, 1039]


@njit(cache=False)
def build_q(size):
    q = np.zeros(size, np.int32)
    for h in range(1, size):
        half = h >> 1
        q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
    return q


QSIZE = max(B // 2 + 2, (1 << 26) + 2)
q = build_q(QSIZE)
check("0  q table: q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2); spot values x(10) = x(13) = 19, x(17) = x(23) = 27, x(256) = y(256) = 513",
      all(2 * h + 1 - 2 * int(q[h]) == v for h, v in ((10, 19), (13, 19), (17, 27), (23, 27), (256, 513))))

# ---------------------------------------------------------------- A. r(N)
X = 4096
st = [(2 * h + 1 - 2 * int(q[h]), 2 * h + 1, h) for h in range(X // 2 + 2) if 2 * h + 1 <= X + 1]
r0 = [0] * (X + 1)
r1 = [0] * (X + 1)
for x1, y1, h1 in st:
    for x2, y2, h2 in st:
        yy = y1 * y2
        if yy > X + 1:
            break
        base = x1 * x2
        s = 0
        while True:
            N = base + (yy << s) - (1 if s == 0 else 0)
            if N > X:
                break
            if s == 0:
                r0[N] += 1
            else:
                r1[N] += 1
            s += 1
print("A. N: r = r0 + r1 (sheet s = 0, sheet s >= 1); one-apex words h with x(h) = N")
table_r = {}
for N in [2, 3, 5, 7, 13, 19, 21, 27, 43, 67, 163, 1039]:
    words = [h for _x, _y, h in st if _x == N]
    table_r[N] = (r0[N] + r1[N], r0[N], r1[N])
    print("      N = %4d: r = %d = %d + %d; words %s" % (N, r0[N] + r1[N], r0[N], r1[N], words))
check("A1 r = 0 at 2 (even), 7, 13, 21, 43, 67, 163, 1039; r(3) = 1, r(5) = 3, r(19) = 4, r(27) = 5",
      all(table_r[N][0] == 0 for N in (2, 7, 13, 21, 43, 67, 163, 1039)) and table_r[3][0] == 1 and table_r[5][0] == 3
      and table_r[19][0] == 4 and table_r[27][0] == 5)
check("A2 every F1 value is odd (r(N) = 0 for all even N <= 4096): the 2-adic factor is the parity alone",
      all(r0[N] + r1[N] == 0 for N in range(0, X + 1, 2)))


# ---------------------------------------------------------------- B. divisible states
@njit(cache=False)
def divisible_scan(q, moduli, limit):
    nm = moduli.shape[0]
    first = np.full(nm, -1, np.int64)
    count = np.zeros(nm, np.int64)
    for h in range(limit):
        y = 2 * h + 1
        qq = q[h]
        for j in range(nm):
            p = moduli[j]
            if y % p == 0 and qq % p == 0:
                count[j] += 1
                if first[j] < 0:
                    first[j] = h
    return first, count


LIM = 1 << 26
first, count = divisible_scan(q, MODULI, LIM)


def order2(p):
    k, v = 1, 2 % p
    while v != 1:
        v = v * 2 % p
        k += 1
    return k


print("B. modulus: first divisible state h (x, y), count below 2^26, count * m^2 / 2^26, ord_m(2), least j with m | 2^j + 1")
first_states = {}
for j, m in enumerate(MODULI):
    m = int(m)
    h = int(first[j])
    xx, yy = 2 * h + 1 - 2 * int(q[h]), 2 * h + 1
    od = order2(m)
    diag = [jj for jj in range(1, 2 * od + 1) if (2 ** jj + 1) % m == 0]
    first_states[m] = h
    print("      m = %4d: h = %7d (x, y) = (%d, %d); count %8d; ratio %.4f; ord = %d; m | 2^j + 1 first at j = %s"
          % (m, h, xx, yy, int(count[j]), int(count[j]) * m * m / LIM, od, diag[0] if diag else "never"))
check("B1 first divisible states: 3: h = 1, 5: 2, 7: 38, 13: 32, 19: 218 (as in ../divisible_states.out)",
      [first_states[m] for m in (3, 5, 7, 13, 19)] == [1, 2, 38, 32, 218])
check("B2 x = y at h = 2^j (x = y = 2^(j+1) + 1) for j < 26, and h = 256 gives x = y = 513 = 27 * 19",
      all(2 * (1 << j) + 1 - 2 * int(q[1 << j]) == 2 * (1 << j) + 1 and int(q[1 << j]) == 0 for j in range(26)) and 513 == 27 * 19)
check("B3 ord_7(2) = 3 and ord_1039(2) = 519 are odd (no m with p | 2^m + 1: no diagonal divisible state); "
      "ord_3, ord_5, ord_13, ord_19, ord_27 = 2, 4, 12, 18, 18 are even",
      order2(7) == 3 and order2(1039) == 519 and [order2(m) for m in (3, 5, 13, 19, 27)] == [2, 4, 12, 18, 18]
      and all((2 ** j + 1) % 7 for j in range(1, 7)) and all(pow(2, j, 1039) != 1038 for j in range(1, 1039)))
ratio1039 = int(count[6]) * 1039 * 1039 / LIM
check("B4 density of 1039-divisible states below 2^26 is 1/p^2 within 4 Poisson standard deviations",
      abs(int(count[6]) - LIM / 1039 ** 2) <= 4 * math.sqrt(LIM / 1039 ** 2), "count %d, expected %.1f" % (int(count[6]), LIM / 1039 ** 2))


# ---------------------------------------------------------------- C. rates
def setup(m, j):
    """T (character k = 2j): (eta, beta) -> (2 eta, 0) weight 1; -> (2 eta + 1, 1) weight zeta^(j (eta + beta)); index 2 eta + beta."""
    eta = np.repeat(np.arange(m), 2)
    beta = np.tile(np.arange(2), m)
    c0 = 2 * ((2 * eta) % m)
    c1 = 2 * ((2 * eta + 1) % m) + 1
    w = np.exp(2j * np.pi * ((j * (eta + beta)) % m) / m)
    inv2 = pow(2, -1, m)
    etap = np.arange(m)
    rr0 = 2 * ((etap * inv2) % m)
    rr1 = 2 * (((etap - 1) * inv2) % m)
    return c0, c1, w, rr0, rr1


def dense(m, j):
    c0, c1, w, _, _ = setup(m, j)
    T = np.zeros((2 * m, 2 * m), dtype=np.complex128)
    idx = np.arange(2 * m)
    T[idx, c0] += 1.0
    T[idx, c1] += w
    return T


def rate_dense(m, j):
    return float(np.max(np.abs(np.linalg.eigvals(dense(m, j))))) / 2.0


def rate_bound(m, j, n=512, steps=12, seed=1):
    """||T^n||_2^(1/n) / 2 by power iteration on (T^n)^* T^n (an upper bound of the rate up to the convergence of the iteration)."""
    c0, c1, w, rr0, rr1 = setup(m, j)
    wc = np.conj(w)
    rng = np.random.default_rng(seed)
    x = rng.normal(size=2 * m) + 1j * rng.normal(size=2 * m)
    x /= np.linalg.norm(x)
    val = 0.0
    for _ in range(steps):
        v = x
        logs = 0.0
        for t in range(n):
            v = v[c0] + w * v[c1]
            if t % 16 == 15:
                sc = np.linalg.norm(v)
                v = v / sc
                logs += math.log(sc)
        sc = np.linalg.norm(v)
        v = v / sc
        logs += math.log(sc)
        val = logs
        for t in range(n):
            u = np.empty(2 * m, dtype=np.complex128)
            u[0::2] = v[rr0] + v[rr0 + 1]
            u[1::2] = wc[rr1] * v[rr1] + wc[rr1 + 1] * v[rr1 + 1]
            v = u
            if t % 16 == 15:
                v = v / np.linalg.norm(v)
        x = v / np.linalg.norm(v)
    return math.exp(val / n) / 2.0


# identity S_n(k, l) = sum_{h < 2^n} zeta^(k q(h) + l h) = (T_k^n v_l)[(0, 0)] at m = 1039, n = 16
m = 1039
nb = 16
hh = np.arange(1 << nb, dtype=np.int64)
qq = q[:1 << nb].astype(np.int64)
ok_id = True
for j, l in ((1, 0), (1, 5), (7, 1038), (300, 77), (519, 2)):
    k = (2 * j) % m
    brute = np.sum(np.exp(2j * np.pi * ((k * qq + l * hh) % m) / m))
    c0, c1, w, _, _ = setup(m, j)
    v = np.exp(2j * np.pi * ((l * np.repeat(np.arange(m), 2)) % m) / m)
    for _ in range(nb):
        v = v[c0] + w * v[c1]
    ok_id &= abs(v[0] - brute) < 1e-6 * (1 << nb)
check("C1 own twisted operator at m = 1039 reproduces the exponential sums S_16(k, l) (5 pairs)", ok_id)

table = {}
path = os.path.join(PARENT, "twisted_spectrum.out")
for line in open(path):
    mm = re.match(r"\s*(\d+)\s+([0-9.]+)\s+(\d+)\s+([0-9.]+)", line)
    if mm:
        table[int(mm.group(1))] = (float(mm.group(2)), int(mm.group(3)))
print("C. rates max_k rho(T_k)/2 (own dense eigenvalues; table = ../twisted_spectrum.out)")
ok_tab = True
rates_small = {}
for m in (3, 5, 7, 13, 19, 27):
    rr = {j: rate_dense(m, j) for j in range(1, m)}
    best = max(rr, key=rr.get)
    rates_small[m] = rr[best]
    ok_tab &= abs(rr[best] - table[m][0]) < 2e-6
    print("      m = %2d: own %.6f (at k = +-%d), table %.6f (k = %d)" % (m, rr[best], (2 * best) % m, table[m][0], table[m][1]))
check("C2 own rates agree with ../twisted_spectrum.out at m = 3, 5, 7, 13, 19, 27 (0.866025, 0.741760, 0.932039, 0.827243, 0.805253, 0.866025)", ok_tab)
# validation of the norm bound at m = 127 (all classes): bound >= dense - 1e-3 and the same maximiser
m = 127
dn = {j: rate_dense(m, j) for j in range(1, (m - 1) // 2 + 1)}
bd = {j: rate_bound(m, j) for j in dn}
check("C3 validation at m = 127: norm bound within [-0.001, +0.012] of the dense rate for all 63 classes, same maximiser, max 0.875417",
      all(-1e-3 <= bd[j] - dn[j] <= 0.012 for j in dn) and max(dn, key=dn.get) == max(bd, key=bd.get)
      and abs(max(dn.values()) - table[127][0]) < 2e-6,
      "max deviation %.5f, min %.5f" % (max(bd[j] - dn[j] for j in dn), min(bd[j] - dn[j] for j in dn)))
m = 1039
bounds = {j: rate_bound(m, j) for j in range(1, (m - 1) // 2 + 1)}
top = sorted(bounds, key=bounds.get, reverse=True)[:16]
dense_top = {j: rate_dense(m, j) for j in top}
jbest = max(dense_top, key=dense_top.get)
RATE1039 = dense_top[jbest]
rest = max(bounds[j] for j in bounds if j not in dense_top)
print("      m = 1039: 519 classes +-j (k = 2j); largest norm bounds:", [(j, round(bounds[j], 5)) for j in top[:8]])
print("      m = 1039: dense rates of the 16 largest:", [(j, round(dense_top[j], 6)) for j in sorted(dense_top, key=dense_top.get, reverse=True)[:8]])
print("      m = 1039: max_k rho(T_k)/2 = %.6f at k = +-%d (j = %d); largest bound among the other 503 classes: %.5f; median bound %.5f"
      % (RATE1039, (2 * jbest) % m, jbest, rest, float(np.median(list(bounds.values())))))
check("C4 m = 1039: the dense maximum exceeds the norm bound of every class not computed densely", RATE1039 > rest,
      "%.6f > %.6f" % (RATE1039, rest))
check("C5 m = 1039: rate < 1 (strict cancellation) and all dense values below their norm bounds + 0.001",
      RATE1039 < 1 and all(dense_top[j] <= bounds[j] + 1e-3 for j in dense_top))


# ---------------------------------------------------------------- D. prediction on the window
@njit(cache=False)
def w_closed(x1, Y, A, B):
    f0 = 0.5 * x1 + Y
    f1 = 1.0 * x1 + Y
    w = 0.0
    if f0 < A:
        top = f1 if f1 < A else A
        w += (B - A) * math.log(top / f0)
    lo = f0 if f0 > A else A
    hi = f1 if f1 < B else B
    if lo < hi:
        w += B * math.log(hi / lo)
    return w / x1


@njit(parallel=True, cache=False)
def predicted(q, A, B, K, moduli, nchunks):
    nm = moduli.shape[0]
    hmax = B // 2 + 1
    T = np.zeros((nchunks, 2))
    F = np.zeros((nchunks, nm, 2))
    step = hmax // nchunks + 1
    for c in prange(nchunks):
        for h1 in range(c * step, min(hmax, (c + 1) * step)):
            y1 = 2 * h1 + 1
            x1 = y1 - 2 * q[h1]
            for s in range(K):
                Y = y1 << s
                if Y >= B:
                    break
                w = w_closed(float(x1), float(Y), float(A), float(B))
                sheet = 0 if s == 0 else 1
                T[c, sheet] += w
                for j in range(nm):
                    p = moduli[j]
                    if y1 % p == 0 and q[h1] % p == 0:
                        F[c, j, sheet] += w
    return T.sum(axis=0), F.sum(axis=0)


Tp, Fp = predicted(q, A, B, K, MODULI, 48)
nodd = (B - A) // 2
pred = {}
for j, m in enumerate(MODULI):
    for sheet in (0, 1):
        f = Fp[j, sheet] / Tp[sheet]
        D = 2 * f - 1.0 / (int(m) ** 2)
        pred[(int(m), sheet)] = (float(f), float(int(m) * D / (1 - D)))
print("D. window odd N in [2^%d, 2^%d): predicted triples per odd N: %.4f (s = 0), %.4f (s >= 1)" % (LO, K, Tp[0] / nodd, Tp[1] / nodd))
for m in (3, 5, 7, 13, 19, 1039):
    print("      p = %4d: E_pred = %+.5f (s = 0), %+.5f (s >= 1);  f = %.3e, %.3e" % (m, pred[(m, 0)][1], pred[(m, 1)][1], pred[(m, 0)][0], pred[(m, 1)][0]))
se1039 = math.sqrt(0.75 * 1039 / nodd)
with open(os.path.join(HERE, "b3_prediction_1039_%d_%d.txt" % (K, LO)), "w") as fh:
    fh.write("written before the enumeration of odd N in [2^%d, 2^%d) by b3_f1_data.py\n" % (LO, K))
    fh.write("p = 1039: first divisible state h = %d; E_pred (s = 0) = %+.6f; E_pred (s >= 1) = %+.6f\n"
             % (first_states[1039], pred[(1039, 0)][1], pred[(1039, 1)][1]))
    fh.write("test: |E - E_pred| <= 0.008 + 4 se, se = sqrt(0.75 p / #odd N) = %.5f\n" % se1039)
if FULL:
    spec = json.load(open(os.path.join(PARENT, "prereg_window_27_24.json")))
    okp = True
    for m in (3, 5, 7, 13, 19):
        for sheet in (0, 1):
            ref = spec["predictions"]["%d:%d" % (m, sheet)]
            okp &= abs(ref["E"] - pred[(m, sheet)][1]) < 2e-6 and abs(ref["f"] - pred[(m, sheet)][0]) < 1e-8
    check("D1 own E_pred and f agree with ../prereg_window_27_24.json at p = 3, 5, 7, 13, 19 (both sheets)", okp)
    check("D2 1039 is not among the 46 primes of the pre-registration (its prediction is new here)",
          "1039:0" not in spec["predictions"])
check("D3 prediction for 1039: |E_pred| < 0.01 on both sheets (no early divisible state: a control prime)",
      abs(pred[(1039, 0)][1]) < 0.01 and abs(pred[(1039, 1)][1]) < 0.01)


# ---------------------------------------------------------------- E. enumeration of the window
@njit(parallel=True, cache=False)
def enumerate_classes(q, A, B, K, moduli, nchunks):
    nm = moduli.shape[0]
    pmax = 0
    for j in range(nm):
        if moduli[j] > pmax:
            pmax = moduli[j]
    cls = np.zeros((nchunks, nm, 2, pmax), np.int64)
    tot = np.zeros((nchunks, 2), np.int64)
    hmax = B // 2 + 1
    H = q.shape[0]
    for c in prange(nchunks):
        for h1 in range(c, hmax, nchunks):
            y1 = 2 * h1 + 1
            x1 = y1 - 2 * q[h1]
            for s in range(K):
                Y = y1 << s
                if Y >= B:
                    break
                d0 = 1 if s == 0 else 0
                sheet = 0 if s == 0 else 1
                # x2 <= y2 gives N <= (x1 + Y) y2 - d0;  2 x2 >= y2 + 1 gives N >= x1 (y2 + 1)/2 + Y y2 - d0
                y2lo = (A + d0) // (x1 + Y)
                h2lo = (y2lo - 1) // 2 - 1
                if h2lo < 0:
                    h2lo = 0
                y2hi = (2 * (B + d0) - x1) // (x1 + 2 * Y) + 1
                h2hi = (y2hi - 1) // 2 + 1
                if h2hi > H - 1:
                    h2hi = H - 1
                for h2 in range(h2lo, h2hi + 1):
                    y2 = 2 * h2 + 1
                    x2 = y2 - 2 * q[h2]
                    N = x1 * x2 + Y * y2 - d0
                    if N >= A and N < B:
                        tot[c, sheet] += 1
                        for j in range(nm):
                            cls[c, j, sheet, N % moduli[j]] += 1
    return cls.sum(axis=0), tot.sum(axis=0)


cls, tot = enumerate_classes(q, A, B, K, MODULI, 96)
print("E. measured triples per odd N: %.4f (s = 0), %.4f (s >= 1), total %.4f" % (tot[0] / nodd, tot[1] / nodd, tot.sum() / nodd))
check("E1 total mass within 0.2 % of the prediction, each sheet within 2 %",
      abs(tot.sum() / Tp.sum() - 1) < 0.002 and abs(tot[0] / Tp[0] - 1) < 0.02 and abs(tot[1] / Tp[1] - 1) < 0.02)


def nclass(p, c):
    """number of odd N in [A, B) with N = c mod p (p odd)."""
    cc = c if c % 2 == 1 else c + p            # the odd representative modulo 2p
    return (B - 1 - cc) // (2 * p) - (A - 1 - cc) // (2 * p)


measured = {}
for j, m in enumerate(MODULI):
    m = int(m)
    for sheet in (0, 1):
        cs = (m - 1) if sheet == 0 else 0
        means = np.array([cls[j, sheet, c] / nclass(m, c) for c in range(m)])
        others = np.delete(means, cs)
        mo = others.mean()
        measured[(m, sheet)] = float((means[cs] - mo) / mo)
for m in (3, 5, 7, 13, 19, 27, 1039):
    print("      p = %4d: E measured = %+.5f (s = 0), %+.5f (s >= 1)" % (m, measured[(m, 0)], measured[(m, 1)])
          + ("" if m == 27 else "   predicted %+.5f, %+.5f" % (pred[(m, 0)][1], pred[(m, 1)][1])))
check("E2 class counts sum to the totals for every modulus", all(int(cls[j, s].sum()) == int(tot[s]) for j in range(len(MODULI)) for s in (0, 1)))
if FULL:
    ref = {}
    for line in open(os.path.join(PARENT, "window_split_prereg_27_24.out")):
        mm = re.match(r"\s*(\d+)\s+(s=0|s>=1)\s+([+-][0-9.]+)\s+([+-][0-9.]+)", line)
        if mm:
            ref[(int(mm.group(1)), 0 if mm.group(2) == "s=0" else 1)] = float(mm.group(3))
    okm = all(abs(ref[(m, s)] - measured[(m, s)]) < 2e-5 for m in (3, 5, 7, 13, 19) for s in (0, 1))
    check("E3 own measured excess equals the session's window_split_prereg_27_24.out at p = 3, 5, 7, 13, 19 (5 decimals)", okm,
          str({k_: round(v, 5) for k_, v in measured.items() if k_[0] in (3, 5, 7, 13, 19)}))
tol = 0.008 + 4 * se1039
check("E4 p = 1039: |E - E_pred| <= 0.008 + 4 se = %.4f on both sheets (control behaviour, as predicted)" % tol,
      abs(measured[(1039, 0)] - pred[(1039, 0)][1]) <= tol and abs(measured[(1039, 1)] - pred[(1039, 1)][1]) <= tol,
      "E = %+.5f, %+.5f; E_pred = %+.5f, %+.5f; se = %.5f" % (measured[(1039, 0)], measured[(1039, 1)], pred[(1039, 0)][1], pred[(1039, 1)][1], se1039))
if FULL:
    check("E5 slow prime 7: measured excess exceeds the first-order prediction by more than 0.04 on both sheets (second-order term)",
          measured[(7, 0)] - pred[(7, 0)][1] > 0.04 and measured[(7, 1)] - pred[(7, 1)][1] > 0.04)

# ---------------------------------------------------------------- F. singular series
C0 = 8 / math.pi ** 2


def odd_primes(n):
    out = []
    while n % 2 == 0:
        n //= 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        out.append(n)
    return out


def S1(N):
    v = C0
    for p in odd_primes(N):
        v *= 1 + p / (p * p - 1)
    return v


W0, W1 = 0.43186124, 0.56813876
print("F. limiting singular series: N, S_1(N) [sheet s >= 1], S_0(N) = S_1(N + 1) [sheet s = 0], S = w_0 S_0 + w_1 S_1, r(N)")
for N in [2, 3, 5, 7, 13, 19, 27, 1039, 21, 43, 67, 163]:
    if N % 2 == 0:
        print("      N = %4d: even, the 2-adic factor is 0 (every F1 value is odd)" % N)
        continue
    print("      N = %4d: S_1 = %.5f  S_0 = %.5f  S = %.5f   r = %d" % (N, S1(N), S1(N + 1), W0 * S1(N + 1) + W1 * S1(N), r0[N] + r1[N]))
prod = 1.0
for p in range(3, 200000, 2):
    if all(p % d for d in range(3, int(p ** 0.5) + 1, 2)):
        prod *= 1 - 1 / p ** 2
check("F1 prod over odd p of (1 - 1/p^2) = 8/pi^2 = 0.810569 (partial product to 2*10^5)", abs(prod - C0) < 1e-5, "%.6f" % prod)
vals = {N: (S1(N), S1(N + 1)) for N in (3, 5, 7, 13, 19, 27, 1039, 21, 43, 67, 163)}
expect = {3: (1.11453, 0.81057), 5: (0.97944, 1.11453), 7: (0.92878, 0.81057), 13: (0.87329, 0.92878), 19: (0.85335, 0.97944),
          27: (1.11453, 0.92878), 1039: (0.81135, 1.05523), 21: (1.27707, 0.88487), 43: (0.82943, 0.88487),
          67: (0.82267, 0.85842), 163: (0.81554, 0.83035)}
check("F2 the table of S_1, S_0 (hand values to 5 decimals)", all(abs(vals[N][0] - expect[N][0]) < 1e-5 and abs(vals[N][1] - expect[N][1]) < 1e-5 for N in expect))
check("F3 no local obstruction at the F1 misses 7, 13, 21, 43, 67, 163, 1039: S_0, S_1 >= 8/pi^2 > 0.81 always",
      all(min(vals[N]) >= C0 - 1e-12 for N in vals))
odd = list(range(1025, 4096, 2))
Stot = {N: W0 * S1(N + 1) + W1 * S1(N) for N in odd}
miss = [N for N in odd if r0[N] + r1[N] == 0]
hit = [N for N in odd if r0[N] + r1[N] > 0]
mS_miss = sum(Stot[N] for N in miss) / len(miss)
mS_all = sum(Stot[N] for N in odd) / len(odd)
rbar = sum(r0[N] + r1[N] for N in odd) / len(odd)
lowS = sorted(odd, key=lambda N: Stot[N])[:len(odd) // 4]
highS = sorted(odd, key=lambda N: Stot[N])[-(len(odd) // 4):]
mr_low = sum(r0[N] + r1[N] for N in lowS) / len(lowS)
mr_high = sum(r0[N] + r1[N] for N in highS) / len(highS)
print("      odd N in [1025, 4096): %d values, %d misses (%.1f %%); mean r = %.3f; mean S over all = %.4f, over the misses = %.4f"
      % (len(odd), len(miss), 100 * len(miss) / len(odd), rbar, mS_all, mS_miss))
print("      mean r on the quarter with the smallest S: %.3f; on the quarter with the largest S: %.3f; exp(-mean r) = %.4f"
      % (mr_low, mr_high, math.exp(-rbar)))
check("F4 misses have smaller limiting singular series on average, and r increases with S (direction of HYP-9657 (B))",
      mS_miss < mS_all and mr_low < mr_high)
check("F5 the miss rate is far above the Poisson value exp(-mean r): misses at this height are not a local-density effect",
      len(miss) / len(odd) > 10 * math.exp(-rbar))

if all(RESULTS):
    print("ALL CHECKS PASSED (%d checks)" % len(RESULTS))
    sys.exit(0)
print("CHECK FAILED (%d of %d)" % (RESULTS.count(False), len(RESULTS)))
sys.exit(1)
