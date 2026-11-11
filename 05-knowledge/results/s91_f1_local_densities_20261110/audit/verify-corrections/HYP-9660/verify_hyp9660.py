#!/usr/bin/env python3
"""Verifier of the corrections to HYP-9660 (2026-11-11): every number and every elementary proof step quoted in the
file is re-derived from the named receipts or from the definitions.  Read-only on the session's files.

Sections
  1. rates_exact/rates_by_modulus.txt against the statements of (B) and against the binary scans
  2. the complement symmetry, the divisor lemma and the S82 dictionary a = -k/2 (own matrices, small moduli)
  3. the formulas and census numbers of (A)
  4. the arithmetic of (D) and (E)
"""
import cmath
import math
import os
import sys

import numpy as np

EXP = "C:/tmp/wt-s91/04-computation/experiments/s91_f1_local_densities_20261110"
OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"{'PASS' if cond else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))


# ---------------------------------------------------------------------------------------------------------------
print("== 1. the table of rates ==")
rate = {}
with open(os.path.join(EXP, "rates_exact/rates_by_modulus.txt")) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        m, r, k, sg = line.split()
        rate[int(m)] = float(r)
mods = sorted(rate)
check("262 moduli: every odd m in [3, 521], and 1023, 1143",
      len(mods) == 262 and mods[:260] == list(range(3, 522, 2)) and mods[260:] == [1023, 1143])
mx = max(rate.values())
check("no rate above 0.932039", abs(mx - 0.932039) < 5e-7, f"max {mx:.6f}")
at_max = [m for m in mods if abs(rate[m] - mx) < 5e-7]
check("the maximum is attained exactly at the multiples of 7 (all of them, and nothing else)",
      at_max == [m for m in mods if m % 7 == 0], f"{len(at_max)} moduli")
p7 = sorted(((rate[m], m) for m in mods if m % 7), reverse=True)
tiers = {}
for r, m in p7:
    tiers.setdefault(round(r, 6), []).append(m)
top = sorted(tiers.items(), reverse=True)[:4]
print("   top tiers prime to 7:", [(r, sorted(ms)) for r, ms in top])
check("prime to 7: 0.913738 at 93, 279, 465 (and 1023); 0.900219 at 31, 155, 341, 403; 0.900128 at 381 (and 1143)",
      sorted(top[0][1]) == [93, 279, 465, 1023] and abs(top[0][0] - 0.913738) < 1e-6
      and sorted(top[1][1]) == [31, 155, 341, 403] and abs(top[1][0] - 0.900219) < 1e-6
      and sorted(top[2][1]) == [381, 1143] and abs(top[2][0] - 0.900128) < 1e-6)
quoted = {3: 0.8660, 5: 0.7418, 25: 0.8153, 15: 0.8922, 45: 0.8922, 75: 0.8922, 51: 0.8690, 65: 0.8481, 85: 0.8503,
          93: 0.9137, 31: 0.9002, 381: 0.9001, 127: 0.8754, 7: 0.932, 11: 0.859, 257: 0.855, 73: 0.852, 241: 0.851}
bad = {m: (rate[m], q) for m, q in quoted.items() if abs(round(rate[m], len(str(q).split('.')[1])) - q) > 1e-9}
check("every rate quoted in (B) is the table value at the quoted precision", not bad, str(bad))
check("63 has the rate of 7; 1023 that of 93; 1143 that of 381",
      abs(rate[63] - rate[7]) < 5e-7 and abs(rate[1023] - rate[93]) < 5e-7 and abs(rate[1143] - rate[381]) < 5e-7)


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def prime_power_divisors(m):
    out = []
    for d in range(3, m + 1, 2):
        if m % d == 0:
            ps = [p for p in range(3, d + 1, 2) if is_prime(p) and d % p == 0]
            if len(ps) == 1:
                out.append(d)
    return out


exc = [m for m in mods if rate[m] > max(rate[d] for d in prime_power_divisors(m) if d in rate) + 5e-5
       and len(prime_power_divisors(m)) and not is_prime(m)]
check("exceptions to the worst-prime-power rule up to 101 are exactly 15, 45, 51, 65, 75, 85, 93",
      [m for m in exc if m <= 101] == [15, 45, 51, 65, 75, 85, 93])
print("   all exceptions in the table:", exc)
prim = [m for m in exc if rate[m] > max(rate[d] for d in range(3, m, 2) if m % d == 0) + 5e-5]
print("   exceptions not inherited from a proper divisor:", [(m, round(rate[m], 6)) for m in prim])
above = [m for m in exc if m > 101]
inher = [m for m in above if m not in prim]
print("   above 101:", len(above), "exceptions;", len(inher), "inherit the rate of a proper divisor:", inher)
check("above 101 the table has 26 more exceptions; 171, 205, 219, 381, 395, 445 exceed every proper divisor "
      "(0.8776, 0.8180, 0.8786, 0.9001, 0.7715, 0.7826), the other 20 inherit",
      len(above) == 26 and [m for m in prim if m > 101] == [171, 205, 219, 381, 395, 445] and len(inher) == 20
      and [round(rate[m], 4) for m in (171, 205, 219, 381, 395, 445)] == [0.8776, 0.8180, 0.8786, 0.9001, 0.7715, 0.7826]
      and all(any(abs(rate[m] - rate[d]) < 5e-7 for d in range(3, m, 2) if m % d == 0) for m in inher))
check("no multiple of 7 is an exception (each has exactly the rate of 7)", not [m for m in exc if m % 7 == 0])
check("381 is an exception, and 45, 75 inherit the rate of 15", 381 in exc and rate[45] == rate[15] == rate[75])
mono = all(rate[m] >= rate[d] - 5e-7 for m in mods for d in range(3, m, 2) if m % d == 0)
check("rate(m) >= rate(d) for every divisor d > 1 (whole table)", mono)
primes = sorted(((rate[m], m) for m in mods if is_prime(m)), reverse=True)
print("   slowest primes:", [(m, round(r, 6)) for r, m in primes[:9]])
check("slowest primes of the scan: 7, 31, 127, 3, 11, 257, 73, 241", [m for r, m in primes[:8]] == [7, 31, 127, 3, 11, 257, 73, 241])
mn = min(rate.values())
check("5 is the fastest modulus of the scan (0.742)", abs(rate[5] - mn) < 1e-12 and round(mn, 3) == 0.742, f"min {mn:.6f}")
oth = [rate[p] for p in mods if is_prime(p) and p <= 89 and p not in (3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 73)]
check("D5 band: the other primes up to 89 lie in [0.7691 (79), 0.8060 (37)]",
      abs(min(oth) - rate[79]) < 1e-12 and abs(max(oth) - rate[37]) < 1e-12 and round(rate[79], 4) == 0.7691
      and round(rate[37], 4) == 0.8060)


def order2(p):
    j, t = 1, 2 % p
    while t != 1:
        t, j = (2 * t) % p, j + 1
    return j


check("orders of 2: 7, 31, 127, 73 -> 3, 5, 7, 9; 5 -> 4; 11 -> 10; 257 -> 16",
      [order2(p) for p in (7, 31, 127, 73, 5, 11, 257)] == [3, 5, 7, 9, 4, 10, 16])
check("1023 = 2^10 - 1 = 3 * 11 * 31, 1143 = 9 * 127", 1023 == 2 ** 10 - 1 == 3 * 11 * 31 and 1143 == 9 * 127)

# the binary scans behind the table
tab = {}
ndense = nops = 0
werr = 0.0
dense_all_ok = True
for f in ("scan_003_101.npz", "scan_103_521.npz", "scan_1023.npz", "scan_1143.npz"):
    d = np.load(os.path.join(EXP, "rates_exact", f))
    for m in np.unique(d["m"]):
        sel = d["m"] == m
        r = d["rate"][sel]
        tab[int(m)] = float(np.nanmax(r))
        if f == "scan_103_521.npz":
            nops += int(sel.sum())
            ndense += int(np.isfinite(r).sum())
            fin = np.isfinite(r)
            werr = max(werr, float(np.abs(r[fin] - d["est"][sel][fin]).max()))
            if m <= 201:
                dense_all_ok &= bool(np.isfinite(r).all())
            else:
                e = d["est"][sel]
                need = (e >= e.max() - 0.012) | (e >= 0.835)
                dense_all_ok &= bool(np.isfinite(r[need]).all())
        if f == "scan_003_101.npz":
            dense_all_ok &= bool(np.isfinite(r).all())
check("the table equals the per-modulus maxima of the four binary scans",
      sorted(tab) == mods and max(abs(tab[m] - rate[m]) for m in mods) < 6e-7)
check("[103, 521]: 65,310 operators, 8,342 dense computations, largest |estimate - dense| = 0.002",
      nops == 65310 and ndense == 8342 and round(werr, 3) == 0.002, f"{nops}, {ndense}, {werr:.5f}")
check("dense eigenvalues of every operator for m <= 201; for m >= 203 of every estimate within 0.012 of the largest or >= 0.835",
      dense_all_ok)

# ---------------------------------------------------------------------------------------------------------------
print("== 2. operators: complement symmetry, divisor lemma, S82 dictionary ==")


def T_full(m, k):
    """THM-5128 (B): (T u)(eta, beta) = u(2 eta, 0) + w^(eta + beta) u(2 eta + 1, 1), w = zeta^(k/2)."""
    inv2 = pow(2, -1, m)
    w = cmath.exp(2j * math.pi * ((k * inv2) % m) / m)
    M = np.zeros((2 * m, 2 * m), dtype=complex)
    for eta in range(m):
        for beta in (0, 1):
            i = 2 * eta + beta
            M[i, 2 * ((2 * eta) % m) + 0] += 1
            M[i, 2 * ((2 * eta + 1) % m) + 1] += w ** (eta + beta)
    return M, w


def C_op(m, w):
    """(C u)(eta, beta) = w^(beta - eta) u(-1 - eta, 1 - beta)."""
    M = np.zeros((2 * m, 2 * m), dtype=complex)
    for eta in range(m):
        for beta in (0, 1):
            M[2 * eta + beta, 2 * ((-1 - eta) % m) + (1 - beta)] += w ** (beta - eta)
    return M


def T_red(m, k, sign):
    """f(eta) -> f(2 eta) +- w^(-(eta + 1)) f(-2 eta - 2)."""
    inv2 = pow(2, -1, m)
    w = cmath.exp(2j * math.pi * ((k * inv2) % m) / m)
    M = np.zeros((m, m), dtype=complex)
    for eta in range(m):
        M[eta, (2 * eta) % m] += 1
        M[eta, (-2 * eta - 2) % m] += sign * w ** (-(eta + 1))
    return M


def T_s82(q, a):
    """chanpak_mass_exploratory_20261106/twisted_scan.py, copied from its definition."""
    n = 2 * q
    M = np.zeros((n, n), dtype=complex)
    for P in range(q):
        for l in (0, 1):
            i = 2 * P + l
            for b in (0, 1):
                add = (P + (1 if l == 1 else 0)) if b == 0 else 0
                j = 2 * ((2 * P + b) % q) + b
                M[j, i] += cmath.exp(2j * math.pi * a * add / q)
    return M


def specmatch(A, B, tol=1e-7):
    """multiset equality of two spectra (greedy)."""
    a = list(np.linalg.eigvals(A))
    b = list(np.linalg.eigvals(B))
    if len(a) != len(b):
        return False
    for z in a:
        j = min(range(len(b)), key=lambda t: abs(b[t] - z))
        if abs(b[j] - z) > tol:
            return False
        b.pop(j)
    return True


comm = c2 = split = conj = True
rad = {}
for m in (3, 5, 7, 9, 15, 21, 31, 45, 51, 63, 65, 93):
    best = 0.0
    for k in range(1, m):
        T, w = T_full(m, k)
        C = C_op(m, w)
        comm &= np.abs(T @ C - C @ T).max() < 1e-9
        c2 &= np.abs(C @ C - cmath.exp(2j * math.pi * k / m) * np.eye(2 * m)).max() < 1e-9
        if m <= 21:
            both = np.block([[T_red(m, k, +1), np.zeros((m, m))], [np.zeros((m, m)), T_red(m, k, -1)]])
            split &= specmatch(T, both, 1e-5)
        r = float(np.abs(np.linalg.eigvals(T)).max()) / 2
        r2 = max(float(np.abs(np.linalg.eigvals(T_red(m, k, s))).max()) for s in (1, -1)) / 2
        split &= abs(r - r2) < 1e-8
        rk = float(np.abs(np.linalg.eigvals(T_full(m, m - k)[0])).max()) / 2
        conj &= abs(r - rk) < 1e-8
        best = max(best, r)
    rad[m] = best
check("T_k commutes with C, and C^2 = zeta^k (12 moduli, every k)", comm and c2)
check("spectrum of T_k = union of the spectra of the two reduced m x m operators (multisets for m <= 21, radii for all)", split)
check("k and m - k have the same spectral radius", conj)
check("own full operators reproduce the table at 3, 5, 7, 9, 15, 21, 31, 45, 51, 63, 65, 93",
      all(abs(rad[m] - rate[m]) < 6e-7 for m in rad), str({m: round(r, 6) for m, r in rad.items()}))

# eigenvector relation u(eta, 1) = +- w^(-eta) u(-1 - eta, 0) on the eigenspaces C = +- w
m, k = 15, 2
T, w = T_full(m, k)
C = C_op(m, w)
vals, vecs = np.linalg.eig(C)
okrel = True
for j in range(2 * m):
    s = +1 if abs(vals[j] - w) < 1e-8 else (-1 if abs(vals[j] + w) < 1e-8 else 0)
    okrel &= s != 0
    u = vecs[:, j]
    for eta in range(m):
        okrel &= abs(u[2 * eta + 1] - s * w ** (-eta) * u[2 * ((-1 - eta) % m)]) < 1e-8
check("C has the eigenvalues +- w only, and there u(eta, 1) = +- w^(-eta) u(-1 - eta, 0) (m = 15, k = 2)", okrel)

# divisor lemma: for k = (m/d) k', the functions of (eta mod d, beta) are invariant and T_k acts as T_k' of d
okdiv = True
for m, d in ((15, 3), (15, 5), (21, 7), (45, 15), (63, 7), (63, 9), (93, 31)):
    for kp in range(1, d):
        k = (m // d) * kp
        T, _ = T_full(m, k)
        Td, _ = T_full(d, kp)
        E = np.zeros((2 * m, 2 * d))                        # extension of functions of (eta mod d, beta)
        for eta in range(m):
            for beta in (0, 1):
                E[2 * eta + beta, 2 * (eta % d) + beta] = 1
        okdiv &= np.abs(T @ E - E @ Td).max() < 1e-9
check("divisor lemma T_k E = E T_k' for k = (m/d) k' (7 pairs (m, d), every k')", okdiv)

# S82 dictionary: the operator of twisted_scan.py at a/q has the spectral radius of T_k with a = -k/2 (mod q)
okdict = True
for q in (7, 15, 31, 63, 93):
    for k in range(1, q):
        a = (-k * pow(2, -1, q)) % q
        r1 = float(np.abs(np.linalg.eigvals(T_s82(q, a))).max())
        r2 = float(np.abs(np.linalg.eigvals(T_full(q, k)[0])).max())
        okdict &= abs(r1 - r2) < 1e-7
check("S82 dictionary: rho(T_(q, a)) = rho(T_k) for a = -k/2 (mod q), q = 7, 15, 31, 63, 93, every k", okdict)
r63 = max(float(np.abs(np.linalg.eigvals(T_full(63, k)[0])).max()) / 2 for k in range(1, 63) if math.gcd(k, 63) == 1)
r93b = float(np.abs(np.linalg.eigvals(T_s82(93, 8))).max()) / 2
check("63: 0.8895 over the characters prime to 63 (S82's 1/63), 0.932039 over all k; 8/93 has 0.8927 (S82 output)",
      round(r63, 4) == 0.8895 and round(r93b, 4) == 0.8927, f"{r63:.6f}, {r93b:.6f}")

# ---------------------------------------------------------------------------------------------------------------
print("== 3. (A) ==")
import mpmath as mp

mp.mp.dps = 30
l = lambda j: -mp.polylog(2, -mp.mpf(2) ** (-j))
kap0 = l(0) - 2 * l(1) + l(2)
w0 = kap0 / (mp.pi ** 2 / 12)
w1 = 1 - w0
sig = 2 * w1
check("w_0 = 0.19640777, w_1 = 0.80359223, sigma = 2 w_1 = 1.60718", abs(w0 - mp.mpf("0.19640777")) < 5e-9
      and abs(w1 - mp.mpf("0.80359223")) < 5e-9 and round(float(sig), 5) == 1.60718, f"{float(sig):.8f}")


def ratio(p, s):
    return (1 + (p - 1) * (s * p - 1) / p ** 3) / (1 - (s * p - 1) / p ** 3)


P8 = (3, 5, 7, 11, 13, 17, 23, 31)
conj_q = (1.495, 1.298, 1.216, 1.140, 1.119, 1.092, 1.068, 1.051)
lett_q = (1.447, 1.274, 1.199, 1.130, 1.110, 1.085, 1.064, 1.047)
check("conjectured limits (sigma = 2 w_1) as quoted", all(round(float(ratio(p, sig)), 3) == c for p, c in zip(P8, conj_q)),
      str([round(float(ratio(p, sig)), 4) for p in P8]))
check("letter-count limits (sigma = 3/2) as quoted", all(round(float(ratio(p, mp.mpf(3) / 2)), 3) == c for p, c in zip(P8, lett_q)),
      str([round(float(ratio(p, mp.mpf(3) / 2)), 4) for p in P8]))
check("sigma = 3/2 is THM-4916: L(0) = 1 + (3p-2)(p-1)/(2p^3), L(r) = 1 - (3p-2)/(2p^3)",
      all(abs((1 + (p - 1) * (1.5 * p - 1) / p ** 3) - (1 + (3 * p - 2) * (p - 1) / (2 * p ** 3))) < 1e-15 for p in P8))
d = np.load(os.path.join(EXP, "census_k2/census_k2_23_24.npz"))
meas = []
for i, p in enumerate(d["primes"]):
    L = d["classes"][i][:int(p)]
    meas.append(float(L[0] / np.mean(L[1:])))
meas_q = (1.857, 1.614, 1.240, 1.374, 1.171, 1.534, 1.124, 1.102)
check("measured ratios on [2^23, 2^24) from the stored census (ties at 3 and 7 resolved)",
      list(map(int, d["primes"])) == list(P8) and all(round(a, 3) == b for a, b in zip(meas, meas_q)),
      str([round(a, 6) for a in meas]))
later = {}
for wnd in ("24_25", "25_26"):
    dd = np.load(os.path.join(EXP, f"census_k2/census_k2_{wnd}.npz"))
    later[wnd] = [float(dd["classes"][i][0] / np.mean(dd["classes"][i][1:int(p)])) for i, p in enumerate(dd["primes"])]
q2526 = (1.842, 1.594, 1.257, 1.361, 1.169, 1.501, 1.130, 1.100)
check("ratios on [2^25, 2^26) as quoted (census_k2_25_26.out)", all(round(a, 3) == b for a, b in zip(later["25_26"], q2526)),
      str([round(a, 6) for a in later["25_26"]]))
seq = list(zip(meas, later["24_25"], later["25_26"]))
fall = [p for p, s in zip(P8, seq) if s[0] > s[1] > s[2]]
rise = [p for p, s in zip(P8, seq) if s[0] < s[1] < s[2]]
check("from [2^23, 2^24) to [2^25, 2^26) the ratios fall at 3, 5, 11, 13, 17, 31 and rise at 7 and 23",
      fall == [3, 5, 11, 13, 17, 31] and rise == [7, 23])
check("on both later windows every ratio still exceeds both candidate limits",
      all(a > float(ratio(p, sig)) and a > float(ratio(p, mp.mpf(3) / 2)) for w in later.values() for a, p in zip(w, P8)))
import re as _re

kl = {}
for wnd in ("12_13", "23_24", "24_25", "25_26"):
    with open(os.path.join(EXP, f"census_k2/census_k2_{wnd}.out")) as fh:
        mm = _re.search(r"kappa = ([0-9.]+), lambda = ([0-9.]+)", fh.read())
    kl[wnd] = (float(mm.group(1)), float(mm.group(2)))
    dd = np.load(os.path.join(EXP, f"census_k2/census_k2_{wnd}.npz"))
    c2 = dd["cnt2"].astype(float)                       # [k1, l2] counts of the two-change words
    assert abs(c2[1, :].sum() / c2.sum() - kl[wnd][0]) < 6e-6 and abs(c2[:, 1].sum() / c2.sum() - kl[wnd][1]) < 6e-6
print("   kappa, lambda by window:", kl)
check("cascade shares kappa, lambda: 0.8041, 0.7792 on [2^23, 2^24); 0.8040, 0.7813 on [2^25, 2^26); 0.7535 on [2^12, 2^13)",
      [round(x, 4) for x in kl["23_24"]] == [0.8041, 0.7792] and [round(x, 4) for x in kl["25_26"]] == [0.8040, 0.7813]
      and round(kl["12_13"][1], 4) == 0.7535)
check("kappa is within 0.001 of w_1; lambda lies between 3/4 and w_1 and rises by about 0.001 per octave",
      all(abs(kl[w][0] - float(w1)) < 1e-3 for w in ("23_24", "24_25", "25_26"))
      and all(0.75 < kl[w][1] < float(w1) for w in kl)
      and round(kl["24_25"][1] - kl["23_24"][1], 3) == 0.001 and round(kl["25_26"][1] - kl["24_25"][1], 3) == 0.001)
exc_c = [a - float(ratio(p, sig)) for a, p in zip(meas, P8)]
exc_l = [a - float(ratio(p, mp.mpf(3) / 2)) for a, p in zip(meas, P8)]
print("   excess over the conjectured limit:", [round(e, 3) for e in exc_c])
check("the measured ratio exceeds both candidate limits at all eight primes", min(exc_c) > 0 and min(exc_l) > 0)
check("excess 0.36, 0.32, 0.23, 0.44 at 3, 5, 11, 17 and 0.02, 0.05, 0.06, 0.05 at 7, 13, 23, 31",
      [round(exc_c[i], 2) for i in (0, 1, 3, 5)] == [0.36, 0.32, 0.23, 0.44]
      and [round(exc_c[i], 2) for i in (2, 4, 6, 7)] == [0.02, 0.05, 0.06, 0.05])
gap = [float(ratio(p, sig) - ratio(p, mp.mpf(3) / 2)) for p in P8]
check("difference of the two candidates: 0.05 at 3, 0.02 at 5 and 7, at most 0.01 from 11 on",
      round(gap[0], 2) == 0.05 and round(gap[1], 2) == 0.02 and round(gap[2], 2) == 0.02 and max(gap[3:]) <= 0.0105,
      str([round(g, 4) for g in gap]))
tabc = {3: (0.104425, 0.074074), 5: (0.063260, 0.032000), 7: (0.037449, 0.017493), 11: (0.025099, 0.007513),
        13: (0.018131, 0.005462), 17: (0.017302, 0.003257), 23: (0.008874, 0.001808), 31: (0.006763, 0.001007)}
rr = [a / b for a, b in tabc.values()]
check("P(S2, -C1) is 1.4 to 6.7 times (p-1)/p^3 (cascade table of census_k2_23_24.out)",
      round(min(rr), 1) == 1.4 and round(max(rr), 1) == 6.7 and all(abs(b - (p - 1) / p ** 3) < 6e-7 for p, (a, b) in tabc.items()),
      str([round(x, 2) for x in rr]))
check("both candidates give P(S2) = share/p^2 + (p-1)/p^3 (columns 'letter c', 'value c' at p = 3: 0.15741, 0.16336)",
      round(0.75 / 9 + 2 / 27, 5) == 0.15741 and round(float(w1) / 9 + 2 / 27, 5) == 0.16336)


def xy(h):
    """x(h) = 2 psi(h) + 1, y(h) = 2h + 1 (THM-5128, Setting)."""
    def psi(h):
        if h == 0:
            return 0
        return psi(h // 2) + (h // 2 if h % 2 == 0 else (h // 2) // 2 + 1)
    return 2 * psi(h) + 1, 2 * h + 1


bal = {p: next(h for h in range(1, 10 ** 5) if xy(h)[0] == xy(h)[1] and xy(h)[0] % p == 0) for p in (3, 5, 11, 13, 17)}
check("balanced divisible states: 3, 5, 11, 17 at h = 1, 2, 16, 8 (x = y = 3, 5, 33, 17); 13 first at h = 32 (65)",
      bal == {3: 1, 5: 2, 11: 16, 13: 32, 17: 8} and [xy(h)[0] for h in (1, 2, 16, 8, 32)] == [3, 5, 33, 17, 65])
check("7, 23, 31 have odd order of 2 (no balanced divisible state)", all(order2(p) % 2 for p in (7, 23, 31)))

# ---------------------------------------------------------------------------------------------------------------
print("== 4. (D), (E) ==")
f3 = {12: (1.0102, 0.9973, 0.9925), 13: (1.0068, 1.0009, 0.9923), 14: (1.0062, 1.0020, 0.9918),
      15: (1.0052, 1.0030, 0.9918), 16: (1.0048, 1.0025, 0.9927)}
defi = [round(100 * (1 - v[2]), 2) for v in f3.values()]
check("deficit of class 2 between 0.73 % and 0.82 %", min(defi) == 0.73 and max(defi) == 0.82, str(defi))
check("class 2 is the largest deviation for n >= 13 and class 0 at n = 12",
      all(abs(v[2] - 1) > max(abs(v[0] - 1), abs(v[1] - 1)) for n, v in f3.items() if n >= 13)
      and abs(f3[12][0] - 1) > abs(f3[12][2] - 1))
dif = [round(f3[n][1] - f3[n][2], 3) for n in (13, 14, 15, 16)]
check("classes 1 and 2 differ by 0.009 to 0.011 for n = 13..16", min(dif) == 0.009 and max(dif) == 0.011, str(dif))
check("largest deviation at n = 16 is 0.73 % (title: up to 0.7 %)", round(100 * max(abs(x - 1) for x in f3[16]), 2) == 0.73)
check("rank-r seam with uniform states: class factors 1 + p^(1-r) - p^(-r) and 1 - p^(-r)",
      all(abs(p * (1 / p + p ** (-r) * (1 - 1 / p)) - (1 + p ** (1 - r) - p ** (-r))) < 1e-15
          and abs(p * (1 / p - p ** (-r) / p) - (1 - p ** (-r))) < 1e-15 for p in (3, 5, 7) for r in (1, 2, 3, 4)))
e8 = 8 / math.pi ** 2
prod = 1.0
for p in range(3, 2 * 10 ** 6, 2):
    if is_prime(p) if p < 10 ** 4 else all(p % q for q in range(3, int(p ** 0.5) + 1, 2)):
        prod *= 1 - 1 / p ** 2
check("prod over odd p of (1 - 1/p^2) = 8/pi^2 = 0.8106", abs(prod - e8) < 1e-6 and round(e8, 4) == 0.8106, f"{prod:.7f}")
d3 = (1 - 1.03549 / 9) / (1 - 1 / 9) - 1
d7 = (1 - 0.77657 / 49) / (1 - 1 / 49) - 1
check("3- and 7-divisible deviations below 2^26 change the coprime proportion by about -0.4 % and +0.5 %, and cancel",
      round(100 * d3, 1) == -0.4 and round(100 * d7, 1) == 0.5 and abs((1 + d3) * (1 + d7) - 1) < 5e-4,
      f"{100 * d3:.2f} %, {100 * d7:.2f} %, product {100 * ((1 + d3) * (1 + d7) - 1):+.3f} %")

print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
sys.exit(0 if OK else 1)
