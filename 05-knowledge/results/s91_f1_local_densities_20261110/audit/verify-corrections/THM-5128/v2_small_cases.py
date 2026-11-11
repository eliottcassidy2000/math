#!/usr/bin/env python3
"""Verifier check 2 (THM-5128, corrected text): the new or rewritten proof steps, tested on small cases.

Written from the statements of the theorem file only (no session or audit code imported).

 (1) definitions, (B3), (F);
 (2) (D3) 'uniform states modulo p^k': exact enumeration against (1 - p^-2) sigma_1(p^v)/p^v, the value at
     c = 0 (mod p^k), the mean over v >= 1 and over v >= 0;
 (3) (E) unitarity for odd k (and failure for even k), the bound 2^((n+j+1)/2), the cell-sum bound used in
     the by-value proof, the quadratic form of q mod 2;
 (4) (E) the matrix M: M M^T = 2I, M^6 = 8J with J^2 = -I, M^12 = -64 I, characteristic polynomial, the
     twelve initial values of S_n, D_n against brute force, zero sets, corank of Q_n;
 (5) (E) by value on a small census: odd classes modulo 2^j, boxes, 'larger word odd', the closed forms of
     the limits 0.5296 / 0.5120, kappa_0, kappa_1;
 (6) (E) by-value proof: cut-cell count 6 * 2^(n-i) and the block bound O(2^(3n/4)) for fixed smaller states;
 (7) (A): identity (1), bound (2), the constant p^((r-1)/2) for c != 0 (and its comparison with (2) at p = 2),
     the equality case;
 (8) (C)/(D3) on the census: confinement of divisible smaller states (primes, 15, 9), first-order law at 9.
"""
import itertools
import math
import random
import sys
from fractions import Fraction

import numpy as np

random.seed(5128)
ok = True


def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name)
    ok = ok and bool(cond)


# ---------------------------------------------------------------- (1) definitions
LOGH = 20
H = 1 << LOGH
psi = np.zeros(H, dtype=np.int64)
for n in range(1, LOGH + 1):
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    par = h >> 1
    psi[lo:hi] = np.where(h & 1 == 0, psi[par] + par, psi[par] + (par >> 1) + 1)
hh = np.arange(H, dtype=np.int64)
q = hh - psi
x = 2 * psi + 1
y = 2 * hh + 1

print("== (1) definitions, (B3), (F) ==")
check("q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2), h < 2^19",
      np.array_equal(q[0:H:2], q[:H // 2]) and np.array_equal(q[1:H:2], q[:H // 2] + (hh[:H // 2] + 1) // 2))
check("0 <= q <= h/2, so y/2 < x <= y", bool(np.all(q >= 0) and np.all(2 * q <= hh) and np.all(2 * x > y) and np.all(x <= y)))
eq = np.nonzero(x == y)[0]
check("x = y exactly at h = 0 and h = 2^j (h < 2^20), with x = y = 2^(j+1) + 1",
      list(eq) == [0] + [1 << j for j in range(LOGH)] and all(x[1 << j] == (1 << (j + 1)) + 1 for j in range(LOGH)))
# (B3): 0 < rho - (1 - phi/2) <= 3n/(4(2h+1)), equality only at h = 1  (exact integers: 2^n * S, S = 2(2h+1)(rho-1+phi/2))
worst = Fraction(0)
allpos = True
eq_only_1 = True
for h in range(1, 1 << 14):
    n = h.bit_length()
    rev = int(format(h, "b")[::-1], 2)          # phi = rev / 2^n
    S = Fraction(-4 * int(q[h])) + Fraction((2 * h + 1) * rev, 1 << n)   # = 2 (2h+1)(rho - 1 + phi/2)
    if S <= 0:
        allpos = False
    r = S / (2 * n)                              # (2h+1)(rho - 1 + phi/2)/n
    if r > Fraction(3, 4):
        eq_only_1 = False
    if r == Fraction(3, 4) and h != 1:
        eq_only_1 = False
    worst = max(worst, r)
check("(B3) one-signed: 0 < rho - (1 - phi/2) for 1 <= h < 2^14", allpos)
check("(B3) <= 3n/(4(2h+1)) with equality only at h = 1 (max ratio %s)" % worst, eq_only_1 and worst == Fraction(3, 4))
g = np.gcd(x, y)
first = {}
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 61, 127, 257, 65537):
    idx = np.nonzero(g % p == 0)[0]
    first[p] = int(idx[0]) if len(idx) else None
print("   first p-divisible states:", first)
check("first p-divisible states quoted in (F)",
      [first[p] for p in (3, 5, 17, 11, 13, 43, 257, 65537, 7, 19, 23, 29, 31, 37, 41, 61, 127)] ==
      [1, 2, 8, 16, 32, 64, 128, 32768, 38, 218, 103, 1435, 8664, 1350, 266, 5886, 109283])
check("13: divisible states 32, then 721", list(np.nonzero(g % 13 == 0)[0][:2]) == [32, 721])


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


small_m = []
for p in range(3, 300, 2):
    if not is_prime(p):
        continue
    m = next((mm for mm in range(1, p) if pow(2, mm, p) == p - 1), None)
    if m is not None and (1 << (m - 1)) < p * p:
        small_m.append(p)
check("primes < 300 with 2^(m-1) < p^2 (least m): %s" % small_m, small_m == [3, 5, 11, 13, 17, 19, 41, 43, 113, 241, 257])
check("29, 61: first divisible state beyond p^2", first[29] > 29 ** 2 and first[61] > 61 ** 2)

# ---------------------------------------------------------------- (2) uniform states modulo p^k
print("== (2) (D3) uniform states modulo p^k ==")


def vp(c, p, k):
    if c == 0:
        return k
    v = 0
    while c % p == 0:
        c //= p
        v += 1
    return v


for p, k, u in ((3, 2, 1), (3, 2, 2), (5, 2, 3), (3, 3, 2), (7, 2, 3)):
    m = p ** k
    cnt = [0] * m
    # count over (a1, b1, a2, b2) in (Z/m)^4 of a1 a2 + u b1 b2 = c, by convolution of the two product laws
    prod = [0] * m
    for a in range(m):
        for b in range(m):
            prod[a * b % m] += 1
    uprod = [0] * m
    for c in range(m):
        uprod[u * c % m] += prod[c]
    for c1 in range(m):
        if prod[c1]:
            for c2 in range(m):
                cnt[(c1 + c2) % m] += prod[c1] * uprod[c2]
    good = True
    for c in range(m):
        v = vp(c, p, k)
        P = Fraction(cnt[c], m ** 4)
        if v < k:
            sig = sum(Fraction(1, p ** i) for i in range(v + 1))        # sigma_1(p^v)/p^v
            want = (1 - Fraction(1, p * p)) * sig
            want2 = (1 + Fraction(1, p)) * (1 - Fraction(1, p ** (v + 1)))
            good = good and m * P == want == want2
        else:
            want0 = (1 + Fraction(1, p)) * (1 - Fraction(1, p ** k)) + Fraction(1, p ** k)
            good = good and m * P == want0
    check("p^k = %d, u = %d: p^k P(c) = (1 - p^-2) sigma_1(p^v)/p^v = (1 + 1/p)(1 - p^(-v-1)) for v < k" % (m, u), good)
    cls = sum(Fraction(cnt[c], m ** 4) for c in range(0, m, p)) * p
    check("   class 0 (mod p) has factor 1 - 1/p^2 + 1/p = %s" % cls, cls == 1 - Fraction(1, p * p) + Fraction(1, p))
for p in (3, 5, 7, 11):
    f = lambda v: (1 + Fraction(1, p)) * (1 - Fraction(1, p ** (v + 1)))
    cf = 1 - Fraction(1, p * p) + Fraction(1, p)
    check("p = %d: factor at v = 0 is 1 - 1/p^2; class factor - factor(v = 1) = p^-3; factor(v >= 2) > class factor" % p,
          f(0) == 1 - Fraction(1, p * p) and cf - f(1) == Fraction(1, p ** 3) and all(f(v) > cf for v in range(2, 9)))
    V = 400
    mean1 = float(sum((1 - Fraction(1, p)) * Fraction(1, p ** (v - 1)) * f(v) for v in range(1, V)))
    mean0 = float(sum((1 - Fraction(1, p)) * Fraction(1, p ** v) * f(v) for v in range(0, V)))
    check("p = %d: mean over v >= 1 = class factor, mean over v >= 0 = 1" % p,
          abs(mean1 - float(cf)) < 1e-12 and abs(mean0 - 1) < 1e-12)
print("   limits at 9: v = 0, 1, c = 0 (mod 9): %.4f %.4f %.4f" % (8 / 9, 32 / 27, 35 / 27))

# (8/pi^2) sigma_1(n)/n has mean 1 over odd N on both sheets (numerical, N < 2^20) and minimum 8/pi^2
LIM = 1 << 20
sig = np.zeros(LIM + 2, dtype=np.float64)       # sigma_1(odd part)/odd part
oddpart = np.arange(LIM + 2, dtype=np.int64)
oddpart[0] = 1
while True:
    ev = (oddpart % 2 == 0)
    if not ev.any():
        break
    oddpart[ev] //= 2
for d in range(1, LIM + 2, 2):
    sig[d::d] += 1.0 / d                         # sum over odd d | n of 1/d  = sigma_1(odd part)/odd part
N_odd = np.arange(1, LIM, 2)
m1 = (8 / math.pi ** 2) * sig[N_odd].mean()
m0 = (8 / math.pi ** 2) * sig[N_odd + 1].mean()
check("(8/pi^2) sigma_1(n)/n: mean over odd N < 2^20 = %.5f (sheet 1), %.5f (sheet 0); minimum 8/pi^2" % (m1, m0),
      abs(m1 - 1) < 2e-3 and abs(m0 - 1) < 2e-3 and sig[N_odd].min() == 1.0)

# ---------------------------------------------------------------- (3) the prime 2: unitarity, bounds, cells
print("== (3) (E) unitarity, bound, cell sums ==")


def Tk2(j, k):
    """twisted matrices of the modulus 2^j: states w mod 2^(j+1); returns (T0, T1)."""
    n_st = 1 << (j + 1)
    T0 = np.zeros((n_st, n_st), dtype=complex)
    T1 = np.zeros((n_st, n_st), dtype=complex)
    for w in range(n_st):
        T0[w, (2 * w) % n_st] = 1
        c = (w + (w & 1)) // 2
        T1[w, (2 * w + 1) % n_st] = np.exp(2j * np.pi * k * c / (1 << j))
    return T0, T1


uni = True
even_fail = True
for j in range(1, 7):
    for k in range(1, 1 << j):
        T0, T1 = Tk2(j, k)
        T = T0 + T1
        is2I = np.allclose(T @ T.conj().T, 2 * np.eye(1 << (j + 1)))
        if k % 2 == 1:
            uni = uni and is2I
        else:
            even_fail = even_fail and (not is2I)
check("T_k T_k^* = 2I for every odd k, j = 1..6", uni)
check("T_k T_k^* != 2I for every even k != 0, j = 1..6", even_fail)
bound_ok = True
for j in range(1, 6):
    for k in range(1, 1 << j, 2):
        for n in range(1, 17):
            s = np.exp(2j * np.pi * k * (q[: 1 << n] % (1 << j)) / (1 << j)).sum()
            bound_ok = bound_ok and abs(s) <= 2 ** ((n + j + 1) / 2) + 1e-9
check("|sum_{h < 2^n} e(k q(h)/2^j)| <= 2^((n+j+1)/2), odd k, j <= 5, n <= 16", bound_ok)
bits = lambda h, i: (h >> i) & 1
qf = all((sum(bits(h, i) * (bits(h, i + 1) + bits(h, i + 2)) for i in range(18)) - int(q[h])) % 2 == 0 for h in range(1 << 16))
check("q(h) = sum_i h_i (h_(i+1) + h_(i+2)) (mod 2), h < 2^16", qf)
check("ceil(H/2) = H_0 + H_1 (mod 2), H < 2^12", all(((H_ + 1) // 2 - (H_ & 1) - ((H_ >> 1) & 1)) % 2 == 0 for H_ in range(1 << 12)))
# cell sums: cells of B_n with i top and i bottom digits; |sum over cell of e(k q/2^j')| <= 2^(r/2) 2^((j'+1)/2)
cell_ok = True
worst_cell = 0.0
for trial in range(400):
    n = random.randint(6, 16)
    i = random.randint(1, n // 2)
    r = n - 2 * i
    jp = random.randint(1, 4)
    k = random.randrange(1, 1 << jp, 2)
    top = random.randrange(1 << (i - 1), 1 << i)
    bot = random.randrange(0, 1 << i)
    free = np.arange(1 << r, dtype=np.int64)
    cell = (top << (n - i)) + (free << i) + bot
    s = np.exp(2j * np.pi * k * (q[cell] % (1 << jp)) / (1 << jp)).sum()
    ratio = abs(s) / (2 ** (r / 2) * 2 ** ((jp + 1) / 2))
    worst_cell = max(worst_cell, ratio)
    cell_ok = cell_ok and ratio <= 1 + 1e-9
check("cell sums of the by-value proof: |sum| <= 2^(r/2) 2^((j'+1)/2) on 400 random cells (worst ratio %.3f)" % worst_cell, cell_ok)

# ---------------------------------------------------------------- (4) the matrix M, period 12, corank
print("== (4) (E) mod 4: M, period 12, corank ==")
M = np.array([[1, 0, 1, 0], [1, 0, -1, 0], [0, 1, 0, -1], [0, 1, 0, 1]], dtype=object)
I4 = np.array([[int(a == b) for b in range(4)] for a in range(4)], dtype=object)


def mpow(A, e):
    R = I4
    for _ in range(e):
        R = R.dot(A)
    return R


check("M M^T = 2I", (M.dot(M.T) == 2 * I4).all())
M6 = mpow(M, 6)
J = M6 // 8
check("M^6 = 8J with J^2 = -I", (M6 == 8 * J).all() and (J.dot(J) == -I4).all())
check("M^12 = -64 I", (mpow(M, 12) == -64 * I4).all())
cp = np.round(np.poly(np.array(M, dtype=float))).astype(int).tolist()
check("characteristic polynomial lam^4 - 2 lam^3 + 2 lam^2 - 4 lam + 4: %s" % cp, cp == [1, -2, 2, -4, 4])
ev = np.linalg.eigvals(np.array(M, dtype=float))
args = sorted(round((np.angle(e) / (2 * np.pi)) * 24) for e in ev)
check("eigenvalues sqrt2 e(+-1/24), sqrt2 e(+-7/24)", np.allclose(abs(ev), math.sqrt(2)) and args == [-7, -1, 1, 7])
s3 = math.sqrt(3)
check("factorisation (lam^2 - (1+sqrt3) lam + 2)(lam^2 + (sqrt3-1) lam + 2)",
      np.allclose(np.polymul([1, -(1 + s3), 2], [1, s3 - 1, 2]), [1, -2, 2, -4, 4]))
# transfer against brute force
e00 = np.array([1, 0, 0, 0], dtype=object)
one = np.array([1, 1, 1, 1], dtype=object)
vv = np.array([1, 1, -1, -1], dtype=object)      # v(a, b) = (-1)^a in the order (0,0), (0,1), (1,0), (1,1)
S_tr, D_tr = [], []
P = I4
for n in range(0, 49):
    S_tr.append(int(e00.dot(P).dot(one)))
    D_tr.append(int(e00.dot(P).dot(vv)))
    P = P.dot(M)
S_bf = [int(((-1) ** (q[: 1 << n] % 2)).sum()) for n in range(0, 19)]
D_bf = [int((x[: 1 << n] % 4 == 1).sum() - (x[: 1 << n] % 4 == 3).sum()) for n in range(0, 19)]
check("S_n, D_n by transfer = brute force for n <= 18", S_tr[:19] == S_bf and D_tr[:19] == D_bf)
check("(D_0..D_11) = (1, 0, 2, 4, 0, 0, 8, 0, -16, 0, 0, -64)", D_tr[:12] == [1, 0, 2, 4, 0, 0, 8, 0, -16, 0, 0, -64])
check("(S_0..S_11) = (1, 2, 2, 0, 0, 0, -8, -16, -16, -32, -64, -64)", S_tr[:12] == [1, 2, 2, 0, 0, 0, -8, -16, -16, -32, -64, -64])
check("S_(n+12) = -64 S_n, D_(n+12) = -64 D_n for n <= 36",
      all(S_tr[n + 12] == -64 * S_tr[n] and D_tr[n + 12] == -64 * D_tr[n] for n in range(37)))
check("D_n = 0 exactly for n = 1, 4, 5, 7, 9, 10 (mod 12); S_n = 0 exactly for n = 3, 4, 5 (mod 12) (n <= 48)",
      all((D_tr[n] == 0) == (n % 12 in (1, 4, 5, 7, 9, 10)) and (S_tr[n] == 0) == (n % 12 in (3, 4, 5)) for n in range(49)))
check("D_18 = -512, D_20 = 1024, D_27 = 16384", (D_tr[18], D_tr[20], D_tr[27]) == (-512, 1024, 16384))


def rank_f2(rows, n):
    rows = list(rows)
    rk = 0
    for col in range(n):
        piv = next((i for i in range(rk, len(rows)) if (rows[i] >> col) & 1), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and (rows[i] >> col) & 1:
                rows[i] ^= rows[rk]
        rk += 1
    return rk


cor = []
for n in range(1, 61):
    rows = []
    for i in range(n):
        r_ = 0
        for d in (1, 2):
            if i - d >= 0:
                r_ |= 1 << (i - d)
            if i + d < n:
                r_ |= 1 << (i + d)
        rows.append(r_)
    cor.append(n - rank_f2(rows, n))
check("corank of Q_n = 1, 0, 1, 2, 1, 0 for n = 1, 2, 3, 4, 5, 0 (mod 6), n <= 60",
      all(cor[n - 1] == [0, 1, 0, 1, 2, 1][n % 6] for n in range(1, 61)))
check("S_n, D_n are 0 or +-2^(n - rank/2), n <= 48",
      all(v == 0 or abs(v) == 2 ** (n - (n - [0, 1, 0, 1, 2, 1][n % 6]) // 2)
          for n in range(1, 49) for v in (S_tr[n], D_tr[n])))

# ---------------------------------------------------------------- (5) census by value
print("== (5) census by value (X = 2^18 and 2^20) ==")


def census(LOGX):
    """ordered triples (h1, h2, s) with N <= 2^LOGX; returns per sheet arrays (N, hsmall, hlarge) with weight
    (2 for h1 != h2, 1 for h1 = h2)."""
    X = 1 << LOGX
    out = {0: [], 1: []}
    s = 0
    while (1 << s) <= X:
        lim = X + (1 if s == 0 else 0)            # M <= lim
        h1 = 0
        while True:
            y1 = 2 * h1 + 1
            if (y1 * y1 << s) >= lim + 1 and (int(x[h1]) * int(x[h1]) + (y1 * y1 << s)) > lim:
                break
            h2max = min(H, (lim // (y1 << s) + 1) // 2 + 1)
            if h2max <= h1:
                break
            Mv = int(x[h1]) * x[h1:h2max] + (y1 << s) * y[h1:h2max]
            sel = np.nonzero(Mv <= lim)[0]
            if len(sel):
                Nv = Mv[sel] - (1 if s == 0 else 0)
                h2 = sel + h1
                wt = np.where(h2 == h1, 1, 2)
                out[0 if s == 0 else 1].append((Nv, np.full(len(sel), h1), h2, wt))
            h1 += 1
        s += 1
    res = {}
    for sh in (0, 1):
        res[sh] = tuple(np.concatenate([t[i] for t in out[sh]]) for i in range(4))
    return res


cen = {}
for LOGX, want in ((18, (540939, 671731)), (20, (2398262, 2995785))):
    cen[LOGX] = census(LOGX)
    T = tuple(int(cen[LOGX][sh][3].sum()) for sh in (0, 1))
    check("X = 2^%d: #T_0, #T_1 = %s (session and audit: %s)" % (LOGX, T, want), T == want)

for LOGX in (18, 20):
    for sh in (0, 1):
        Nv, hs, hl, wt = cen[LOGX][sh]
        T = wt.sum()
        check("   every N odd (sheet %d)" % sh, bool(np.all(Nv % 2 == 1)))
        devs = []
        for j in (2, 3, 4, 5):
            cnt = np.bincount(Nv % (1 << j), weights=wt, minlength=1 << j)
            odd = cnt[1::2]
            devs.append(float(np.abs(odd / (T / (1 << (j - 1))) - 1).max()))
        lodd = float((wt * (hl & 1)).sum() / T)
        print("   X = 2^%d sheet %d: max rel. dev. of the odd classes mod 4, 8, 16, 32: %s; larger word odd: %.5f"
              % (LOGX, sh, " ".join("%.5f" % d for d in devs), lodd))
        # y of the larger state by value is NOT uniform mod 8
        cy = np.bincount((2 * hl + 1) % 8, weights=wt, minlength=8)[1::2] / T
        if LOGX == 20:
            check("   by value the larger state is not equidistributed mod 8 (shares of y = 1, 3, 5, 7: %s)"
                  % " ".join("%.4f" % c for c in cy), float(np.abs(cy - 0.25).max()) > 0.01)
# audit values at k = 18 and 20 (a6): sheet 0 larger word odd 0.53516, 0.53458; sheet 1 0.51522, 0.51485
for LOGX, w0, w1 in ((18, 0.53516, 0.51522), (20, 0.53458, 0.51485)):
    got = []
    for sh in (0, 1):
        Nv, hs, hl, wt = cen[LOGX][sh]
        st = hs != hl                             # the audit's share is over the triples with h1 != h2
        got.append(float((wt[st] * (hl[st] & 1)).sum() / wt[st].sum()))
    check("   larger word odd at 2^%d (triples with h1 != h2): %.5f, %.5f (audit a6: %.5f, %.5f)" % (LOGX, got[0], got[1], w0, w1),
          abs(got[0] - w0) < 6e-6 and abs(got[1] - w1) < 6e-6)

try:
    import mpmath as mp
    mp.mp.dps = 30
    Li2 = lambda z: mp.polylog(2, z)
    l = lambda jj: -Li2(-mp.mpf(2) ** (-jj))
    k0 = l(0) - 2 * l(1) + l(2)
    k1 = l(1) - l(2)
    check("kappa_0 = %.10f, kappa_1 = %.10f" % (k0, k1), abs(k0 - mp.mpf("0.1615389173")) < 1e-10 and abs(k1 - mp.mpf("0.2125139092")) < 1e-10)
    check("kappa_0 = pi^2/12 + 2 Li2(-1/2) - Li2(-1/4), kappa_1 = Li2(-1/4) - Li2(-1/2)",
          abs(k0 - (mp.pi ** 2 / 12 + 2 * Li2(-0.5) - Li2(-0.25))) < 1e-25 and abs(k1 - (Li2(-0.25) - Li2(-0.5))) < 1e-25)
    gbar = 2 * (k0 + k1)
    check("g_bar = %.10f; shares %.8f, %.8f; word pairs %.8f, %.8f" % (gbar, k0 / (k0 + k1), k1 / (k0 + k1),
                                                                    k0 / (mp.pi ** 2 / 12), 1 - k0 / (mp.pi ** 2 / 12)),
          abs(gbar - mp.mpf("0.7481056530")) < 1e-10 and abs(k0 / (k0 + k1) - mp.mpf("0.43186124")) < 1e-8
          and abs(k0 / (mp.pi ** 2 / 12) - mp.mpf("0.19640777")) < 1e-8)
    # F_s^odd
    for s in range(0, 4):
        Fodd = lambda r_: (1 / r_) * mp.log((3 * r_ / 4 + 2 ** s) / (r_ / 2 + 2 ** s))
        num = mp.quad(Fodd, [mp.mpf(1) / 2, 1])
        closed = Li2(-mp.mpf(2) ** (-s - 1)) - Li2(-mp.mpf(2) ** (-s - 2)) + Li2(-3 * mp.mpf(2) ** (-s - 3)) - Li2(-3 * mp.mpf(2) ** (-s - 2))
        Fs = lambda r_: (1 / r_) * mp.log((r_ + 2 ** s) / (r_ / 2 + 2 ** s))
        full = mp.quad(Fs, [mp.mpf(1) / 2, 1])
        check("   s = %d: int F_s^odd = closed form (%.10f); int F_s = l_s - 2 l_(s+1) + l_(s+2)" % (s, closed),
              abs(num - closed) < 1e-20 and abs(full - (l(s) - 2 * l(s + 1) + l(s + 2))) < 1e-20)
    lim0 = (Li2(-0.5) - Li2(-0.25) + Li2(-mp.mpf(3) / 8) - Li2(-0.75)) / k0
    lim1 = (Li2(-0.25) - Li2(-mp.mpf(3) / 8)) / k1
    tele = mp.nsum(lambda s: Li2(-mp.mpf(2) ** (-s - 1)) - Li2(-mp.mpf(2) ** (-s - 2)) + Li2(-3 * mp.mpf(2) ** (-s - 3))
                   - Li2(-3 * mp.mpf(2) ** (-s - 2)), [1, mp.inf])
    check("limits of 'larger word odd': %.6f (sheet 0), %.6f (sheet 1); quoted 0.5296 (0.52956 kappa_0), 0.5120 (0.51198 kappa_1)"
          % (lim0, lim1), abs(lim0 - mp.mpf("0.529558")) < 1e-6 and abs(lim1 - mp.mpf("0.511980")) < 1e-6
          and abs(tele - (Li2(-0.25) - Li2(-mp.mpf(3) / 8))) < 1e-15)
    HAVE_MP = True
except ImportError:
    print("   mpmath not available: closed forms not checked")
    HAVE_MP = False

# boxes: every odd class mod 2^j receives 2^(2n+1-j) + O(2^(3n/2)) triples of a box h1, h2 < 2^n
box_worst = 0.0
for n in (6, 8, 10):
    for s in (0, 1, 3):
        xs, ys = x[: 1 << n], y[: 1 << n]
        Nbox = (np.outer(xs, xs) + (np.outer(ys, ys) << s) - (1 if s == 0 else 0))
        for j in (2, 3, 4):
            cnt = np.bincount((Nbox % (1 << j)).ravel(), minlength=1 << j)
            assert cnt[0::2].sum() == 0
            dev = np.abs(cnt[1::2] - 2.0 ** (2 * n + 1 - j)).max() / 2 ** (1.5 * n)
            box_worst = max(box_worst, dev)
check("boxes: |count - 2^(2n+1-j)| <= C 2^(3n/2), n = 6, 8, 10, j = 2, 3, 4 (worst ratio %.3f)" % box_worst, box_worst < 2.0)

# ---------------------------------------------------------------- (6) by-value proof: cut cells and block bound
print("== (6) (E) by-value proof: cut cells, block bound ==")
cut_ok = True
worst_cut = 0.0
blk_worst = 0.0
for trial in range(300):
    n = random.randint(8, 18)
    i = -(-n // 4)
    h1 = random.choice([0, 1, 2, 3, 5, 8, 21, 38, 100, 257])
    s = random.randint(0, 3)
    a, b = int(x[h1]), (2 * h1 + 1) << s
    lo, hi = 1 << (n - 1), 1 << n
    # X chosen so that the hyperbola cuts the block
    Xr = int((a * 0.75 + b) * (2 * random.uniform(lo, hi) + 1))
    hb = np.arange(lo, hi, dtype=np.int64)
    inside = (a * x[lo:hi] + b * y[lo:hi] <= Xr) & (hb > h1)
    r = n - 2 * i
    top = hb >> (n - i)
    bot = hb & ((1 << i) - 1)
    cid = (top - (1 << (i - 1))) * (1 << i) + bot
    ncell = 1 << (2 * i - 1)
    per = np.bincount(cid, weights=inside, minlength=ncell)
    cut = (per > 0) & (per < (1 << r))
    pts = int(cut.sum()) * (1 << r)
    worst_cut = max(worst_cut, pts / (6 * 2 ** (n - i)))
    cut_ok = cut_ok and pts <= 6 * 2 ** (n - i)
    tot = int(inside.sum())
    if tot:
        Nv = a * x[lo:hi][inside] + b * y[lo:hi][inside] - (1 if s == 0 else 0)
        for j in (2, 3, 4):
            cnt = np.bincount(Nv % (1 << j), minlength=1 << j)[1::2]
            blk_worst = max(blk_worst, float(np.abs(cnt - tot / (1 << (j - 1))).max()) / 2 ** (0.75 * n))
check("cut cells (i = ceil(n/4) top and bottom digits) hold at most 6 * 2^(n-i) points on 300 random regions (worst %.3f of the bound)"
      % worst_cut, cut_ok)
check("block bound: |#{N = c mod 2^j} - 2^(1-j) total| <= C 2^(3n/4), j = 2, 3, 4 (worst ratio %.3f)" % blk_worst, blk_worst < 2.0)

# ---------------------------------------------------------------- (7) (A)
print("== (7) (A): identity, bounds, equality ==")


def rand_dist(size, kind):
    if kind == 0:
        w = [random.randint(1, 20) for _ in range(size)]
    elif kind == 1:
        w = [random.randint(0, 3) for _ in range(size)]
        if sum(w) == 0:
            w[random.randrange(size)] = 1
    else:
        w = [1] * size
        w[0] = random.randint(1, 40)
    tot = sum(w)
    return [Fraction(a, tot) for a in w]


A_ok = True
c_ne0_ok = True
ratio_ne0 = 0.0
for p in (2, 3, 5, 7):
    for r in (1, 2, 3):
        if p ** r > 130:
            continue
        vecs = list(itertools.product(range(p), repeat=r))
        for trial in range(6):
            while True:
                Bm = [[random.randrange(p) for _ in range(r)] for _ in range(r)]
                det = round(np.linalg.det(np.array(Bm, dtype=float))) % p
                if det:
                    break
            D1 = rand_dist(len(vecs), trial % 3)
            D2 = rand_dist(len(vecs), (trial + 1) % 3)
            Bval = lambda u, v: sum(u[a_] * Bm[a_][b_] * v[b_] for a_ in range(r) for b_ in range(r)) % p
            P = [Fraction(0)] * p
            for iu, u in enumerate(vecs):
                for iv, v in enumerate(vecs):
                    P[Bval(u, v)] += D1[iu] * D2[iv]
            pi1, pi2 = D1[0], D2[0]
            mm = pi1 + pi2 - Fraction(1, p ** r)
            Rh = [P[c] - Fraction(1, p) - mm * ((1 if c == 0 else 0) - Fraction(1, p)) for c in range(p)]
            n1 = sum((d - Fraction(1, p ** r)) ** 2 for d in D1)
            n2 = sum((d - Fraction(1, p ** r)) ** 2 for d in D2)
            A_ok = A_ok and sum(Rh) == 0
            for c in range(p):
                A_ok = A_ok and Rh[c] ** 2 <= Fraction((p - 1) ** 2, p * p) * p ** r * n1 * n2
                if c != 0:
                    c_ne0_ok = c_ne0_ok and Rh[c] ** 2 <= p ** (r - 1) * n1 * n2
                    if n1 * n2 > 0:
                        ratio_ne0 = max(ratio_ne0, math.sqrt(float(Rh[c] ** 2 / (p ** (r - 1) * n1 * n2))))
            # excess of class 0 over the mean of the others
            if p > 1:
                exc = P[0] - sum(P[1:]) / (p - 1)
                A_ok = A_ok and exc == mm + Fraction(p, p - 1) * Rh[0]
check("(1), sum Rhat = 0, bound (2), 'excess over the mean of the others = pi1 + pi2 - p^-r + (p/(p-1)) Rhat(0)'", A_ok)
check("c != 0: |Rhat(c)| <= p^((r-1)/2) |D1 - U| |D2 - U| (largest ratio seen %.3f)" % ratio_ne0, c_ne0_ok)
for p in (2, 3, 5, 7):
    c2 = (p - 1) / p * p ** 0.5       # constant of (2) / p^((r-1)/2)
    print("   p = %d: constant of (2) divided by p^((r-1)/2) = (p-1)/sqrt(p) = %.4f  (%s)"
          % (p, c2, "the c != 0 constant is an improvement" if c2 > 1 else "NOT an improvement: (2) is already smaller"))
# equality case at c = 0, r = 2
eq_ok = True
for p in (3, 5, 7):
    vecs = list(itertools.product(range(p), repeat=2))
    for s in (0, 1):
        u_ = pow(2, s, p)
        Bv = lambda a_, b_: (a_[0] * b_[0] + u_ * a_[1] * b_[1]) % p
        L = [(t_ % p, 0) for t_ in range(1, p)]
        Lp = [(t_ % p, t_ % p) for t_ in range(1, p)]
        perp = lambda line: [v for v in vecs if v != (0, 0) and all(Bv(u, v) == 0 for u in line)]
        t = Fraction(1, p * p)
        U = Fraction(1, p * p)
        D1 = {v: U for v in vecs}
        D2 = {v: U for v in vecs}
        for v in L:
            D1[v] += t
        for v in Lp:
            D1[v] -= t
        for v in perp(L):
            D2[v] += t
        for v in perp(Lp):
            D2[v] -= t
        P0 = sum(D1[a_] * D2[b_] for a_ in vecs for b_ in vecs if Bv(a_, b_) == 0)
        Rh0 = P0 - Fraction(1, p) - (D1[(0, 0)] + D2[(0, 0)] - U) * (1 - Fraction(1, p))
        n1 = sum((d - U) ** 2 for d in D1.values())
        n2 = sum((d - U) ** 2 for d in D2.values())
        eq_ok = eq_ok and Rh0 == 2 * t * t * (p - 1) ** 2 and Rh0 ** 2 == (p - 1) ** 2 * n1 * n2
check("equality in (2) at c = 0, r = 2 (F1 forms, p = 3, 5, 7): Rhat(0) = 2 t^2 (p-1)^2 = (p-1)|D1-U||D2-U|", eq_ok)

# ---------------------------------------------------------------- (8) (C), (D3) on the census
print("== (8) (C)/(D3) on the census ==")
for LOGX in (18, 20):
    for sh in (0, 1):
        Nv, hs, hl, wt = cen[LOGX][sh]
        Mv = Nv + (1 if sh == 0 else 0)
        gs = g[hs]
        strict = hs != hl
        viol = 0
        for p in (3, 5, 7, 11, 13, 17, 31):
            viol += int(((gs % p == 0) & (Mv % p != 0)).sum())
        viol += int(((np.gcd(gs, 15) == 3) & (Mv % 3 != 0)).sum()) + int(((np.gcd(gs, 15) == 5) & (Mv % 5 != 0)).sum())
        viol += int(((gs % 9 == 0) & (Mv % 9 != 0)).sum())
        check("X = 2^%d sheet %d: a divisible smaller state forces the distinguished class (p, 15, 9): %d violations"
              % (LOGX, sh, viol), viol == 0)
        T = float(wt.sum())
        if LOGX == 18:
            # D of (C2): smaller state divisible, h1 != h2 (the audit's a3_analysis_28.out); the session's
            # first_order_exact.py also counts the diagonal h1 = h2 in F_small
            for p, want_s, want_a in ((3, (0.02852, 0.04830), (0.02854, 0.04839)), (5, (0.00501, 0.00553), (0.00507, 0.00572)),
                                      (7, (0.15987, 0.10714), (0.15991, 0.10715))):
                cnt = np.bincount(Mv % p, weights=wt, minlength=p)
                Rs = []
                for incl_diag in (False, True):
                    D = float((wt * ((gs % p == 0) & (strict | incl_diag))).sum())
                    pred = np.full(p, (T - D) / p)
                    pred[0] += D
                    Rs.append(p * np.abs(cnt - pred).max() / T)
                check("   (C2) p = %d: p max|R_small| = %.5f with h1 != h2 (audit a3: %.5f), %.5f with the diagonal (session: %.5f)"
                      % (p, Rs[0], want_a[sh], Rs[1], want_s[sh]), abs(Rs[0] - want_a[sh]) < 6e-6 and abs(Rs[1] - want_s[sh]) < 6e-6)
            # modulus 9: first-order law with the shares of smaller states with 3 || g, 9 | g
            F1_ = float((wt * ((gs % 3 == 0) & (gs % 9 != 0) & strict)).sum()) / T
            F2_ = float((wt * ((gs % 9 == 0) & strict)).sum()) / T
            F0_ = 1 - F1_ - F2_
            cnt = np.bincount(Mv % 9, weights=wt, minlength=9) / T
            meas = (9 * cnt[[1, 2, 4, 5, 7, 8]].mean(), 9 * cnt[[3, 6]].mean(), 9 * cnt[0])
            law = (F0_, F0_ + 3 * F1_, F0_ + 3 * F1_ + 9 * F2_)
            aud = {0: ((0.8279, 0.8285), (1.2365, 1.2455), (1.5596, 1.5380)),
                   1: ((0.8104, 0.8149), (1.2568, 1.2625), (1.6238, 1.5854))}[sh]
            print("   modulo 9: 9 * share for v = 0, 1, >= 2: measured %.4f %.4f %.4f; first-order law %.4f %.4f %.4f"
                  % (meas + law))
            check("   modulo 9 at 2^18: agrees with the audit's a4_moduli_26.out row k = 18",
                  all(abs(meas[i_] - aud[i_][0]) < 6e-5 and abs(law[i_] - aud[i_][1]) < 6e-5 for i_ in range(3)))

print("ALL PASS" if ok else "SOME FAIL")
sys.exit(0 if ok else 1)
