#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A5: THM-5128 (E), the prime 2 (own code; exact where stated).

(1) T_k on Z/2^(j+1) for j = 1..6: T_k T_k^* = 2I for every odd k (exact Gaussian-rational check through the phase
    exponents, not floating point); for even k != 0 it fails; e_0^T T_k^n 1 = sum_{h < 2^n} e(k q(h)/2^j) (float,
    against brute force); |sum| <= 2^((n + j + 1)/2) for n <= 22 and whether the bound is attained.
(2) q(h) = sum_i h_i (h_(i+1) + h_(i+2)) mod 2 for all h < 2^22 (exact).
(3) Joint law of (x, y) mod 2^j on h < 2^n: largest deviation from 4^-(j-1) on pairs of odd residues, scaled by 2^(n/2).
(4) Band matrix B_n over F_2 (ones at distances 1, 2): corank for n <= 120, against 1,0,1,2,1,0 for n = 1..5,0 mod 6.
    S_n, D_n exactly for n <= 96 by a digit recursion in python integers (checked against brute force for n <= 22);
    S_n, D_n in {0, +-2^(n - rank/2)}; D_18, D_20, D_27; zero set of D_n modulo 12.
"""
import sys
import numpy as np
from a_lib import psi_table

ok = True
NB = 22
psi = psi_table(1 << NB)
hh = np.arange(1 << NB, dtype=np.int64)
q = hh - psi
x = 2 * psi + 1
y = 2 * hh + 1

# ---------------------------------------------------------------- (1)
print("(1) twisted matrices modulo 2^j")
for j in range(1, 7):
    S = 1 << (j + 1)                 # state space Z/2^(j+1)
    mod = 1 << j
    c = np.array([(w + (w & 1)) // 2 for w in range(S)])
    # exact check of T T^* = 2 I through exponents: rows w, w' share a column pair iff 2w = 2w' mod S
    all_odd_ok = True
    even_fail = 0
    even_tot = 0
    for k in range(1, mod):
        good = True
        for w in range(S):
            for w2 in range(S):
                if (2 * w) % S != (2 * w2) % S:
                    continue          # disjoint supports: inner product 0
                # inner product = 1 + e(k (c(w) - c(w2))/2^j)
                e = (k * (c[w] - c[w2])) % mod
                if w == w2:
                    if e != 0:
                        good = False  # norm^2 must be 2
                else:
                    if 2 * e != mod:  # need e(.) = -1, i.e. exponent = 1/2
                        good = False
        if k % 2 == 1:
            all_odd_ok &= good
        else:
            even_tot += 1
            even_fail += (not good)
    # floating check of the sum formula and the bound
    worst_ratio = 0.0
    worst_err = 0.0
    att = None
    for k in range(1, mod, 2):
        T = np.zeros((S, S), dtype=complex)
        for w in range(S):
            T[w, (2 * w) % S] += 1.0
            T[w, (2 * w + 1) % S] += np.exp(2j * np.pi * k * c[w] / mod)
        v = np.ones(S, dtype=complex)
        ph = np.exp(2j * np.pi * k * (q % mod) / mod)
        cs = np.cumsum(ph)
        for n in range(1, NB + 1):
            v = T @ v
            brute = cs[(1 << n) - 1]
            worst_err = max(worst_err, abs(v[0] - brute))
            r = abs(brute) / 2 ** ((n + j + 1) / 2)
            if r > worst_ratio:
                worst_ratio = r
                att = (k, n)
    ok &= all_odd_ok and worst_ratio <= 1 + 1e-9 and worst_err < 1e-6
    print("  j = %d: T_k T_k^* = 2I for all odd k (exact): %s; even k != 0 failing it: %d of %d; "
          "max |e_0 T^n 1 - brute| = %.1e (n <= %d); max |sum| / 2^((n+j+1)/2) = %.4f at (k, n) = %s"
          % (j, all_odd_ok, even_fail, even_tot, worst_err, NB, worst_ratio, att))

# ---------------------------------------------------------------- (2)
bits = [(hh >> i) & 1 for i in range(NB + 2)]
form = np.zeros(1 << NB, dtype=np.int64)
for i in range(NB):
    form += bits[i] * (bits[i + 1] + bits[i + 2])
good2 = bool(np.all((form - q) % 2 == 0))
ok &= good2
print("(2) q(h) = sum_i h_i (h_(i+1) + h_(i+2)) (mod 2) for all h < 2^%d: %s" % (NB, good2))

# ---------------------------------------------------------------- (3)
print("(3) joint law of (x, y) mod 2^j over h < 2^n (index order): max |count - 2^n/4^(j-1)| / 2^(n/2)")
for j in range(2, 7):
    m = 1 << j
    row = []
    for n in (12, 16, 20, 22):
        code = ((x[: 1 << n] % m) >> 1) * (m >> 1) + ((y[: 1 << n] % m) >> 1)
        cnt = np.bincount(code, minlength=(m >> 1) ** 2)
        row.append("n=%d: %.3f" % (n, np.abs(cnt - 2.0 ** n / (m >> 1) ** 2).max() / 2 ** (n / 2)))
    print("  j = %d: " % j + "  ".join(row))


# ---------------------------------------------------------------- (4)
def corank(n):
    rows = []
    for i in range(n):
        r = 0
        for d in (1, 2):
            if i - d >= 0:
                r |= 1 << (i - d)
            if i + d < n:
                r |= 1 << (i + d)
        rows.append(r)
    rank = 0
    for bit in range(n):
        piv = None
        for idx in range(rank, n):
            if (rows[idx] >> bit) & 1:
                piv = idx
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for idx in range(n):
            if idx != rank and (rows[idx] >> bit) & 1:
                rows[idx] ^= rows[rank]
        rank += 1
    return n - rank


pattern = {1: 1, 2: 0, 3: 1, 4: 2, 5: 1, 0: 0}
cor = [corank(n) for n in range(1, 121)]
good4 = all(cor[n - 1] == pattern[n % 6] for n in range(1, 121))
ok &= good4
print("(4) corank of B_n equals 1,0,1,2,1,0 for n = 1,2,3,4,5,0 (mod 6), all n <= 120:", good4)


def SD(nmax):
    """Exact S_n, D_n: read digits from the top; state = (previous digit a = h_0, digit before b = h_1) and the sign.
    Appending digit d: q -> q + d * (h_0 + h_1) (mod 2)."""
    # vec[(a, b)] = signed count
    vec = {(0, 0): 1}
    S, D = [None], [None]
    for n in range(1, nmax + 1):
        new = {}
        for (a, b), val in vec.items():
            for d in (0, 1):
                sgn = -1 if (d and ((a + b) & 1)) else 1
                key = (d, a)
                new[key] = new.get(key, 0) + sgn * val
        vec = new
        S.append(sum(vec.values()))
        D.append(sum(val * (-1 if a else 1) for (a, b), val in vec.items()))
    return S, D


S, D = SD(96)
sg = np.where(q % 2 == 0, 1, -1)
sgd = np.where((hh + q) % 2 == 0, 1, -1)
x4 = np.where(x % 4 == 1, 1, -1)
good_b = True
for n in range(1, NB + 1):
    good_b &= int(sg[: 1 << n].sum()) == S[n] and int(sgd[: 1 << n].sum()) == D[n] and int(x4[: 1 << n].sum()) == D[n]
ok &= good_b
print("    digit recursion for S_n, D_n equals brute force (and D_n = #{x = 1 mod 4} - #{x = 3 mod 4}) for n <= %d: %s"
      % (NB, good_b))
good_m = True
for n in range(1, 97):
    r = n - cor[n - 1]
    for val in (S[n], D[n]):
        if val != 0 and abs(val) != 2 ** (n - r // 2):
            good_m = False
ok &= good_m
print("    S_n, D_n in {0, +-2^(n - rank/2)} for n <= 96:", good_m)
print("    D_18 = %d, D_20 = %d, D_27 = %d   (file: -512, 1024, 16384)" % (D[18], D[20], D[27]))
ok &= (D[18], D[20], D[27]) == (-512, 1024, 16384)
zD = sorted({n % 12 for n in range(1, 97) if D[n] == 0})
nzD = sorted({n % 12 for n in range(1, 97) if D[n] != 0})
zS = sorted({n % 12 for n in range(1, 97) if S[n] == 0})
print("    residues n mod 12 with D_n = 0 (n <= 96): %s; with D_n != 0: %s   (file: zero at 1, 4, 5, 7, 9, 10)" % (zD, nzD))
print("    residues n mod 12 with S_n = 0 (n <= 96): %s" % zS)
ok &= zD == [1, 4, 5, 7, 9, 10] and not (set(zD) & set(nzD))
print("    signs of D_n for n = 2, 3, 6, 8, 11, 12, 14, 15, 18, 20, 23, 24, ...: "
      + " ".join("%d:%s" % (n, "+" if D[n] > 0 else "-") for n in range(1, 61) if D[n] != 0))
print("    |S_n| = 2^((n+2)/2) (the bound of (E) with j = 1 attained) at n =",
      [n for n in range(1, 97) if abs(S[n]) == 2 ** ((n + 2) // 2) and n % 2 == 0])
print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
sys.exit(0 if ok else 1)
