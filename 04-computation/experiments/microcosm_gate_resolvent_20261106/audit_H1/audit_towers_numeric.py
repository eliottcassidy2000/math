# audit_H1: independent numerical checks of THM-4900 (2)-(3): one-vertex tower spectrum (cot^2 formula, k = 1..14), gate-tower limit (2+s)/z and its rate, D-tower and E-tower limits, positivity of F_T(z) for z > 0.
import math
import random
from fractions import Fraction
import numpy as np
from audit_tournament_lib import rand_tour, dominance, gate, pad, double, F_at, phi

rng = random.Random(7)
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


# ---- one-vertex spectrum ----
for k in range(1, 15):
    X = [[0]]
    for _ in range(k):
        X = gate(X)
    M = np.array(dominance(X), dtype=float)
    ev = np.linalg.eigvals(M)
    check(np.max(np.abs(ev.real)) < 1e-8, "eigenvalues not imaginary")
    im = np.sort(ev.imag)
    pred = [0.0, math.sqrt(3), -math.sqrt(3)]
    for j in range(0, k - 1):
        c = 1 / math.tan((2 * j + 1) * math.pi / (2 * (2 * k - 1)))
        r = math.sqrt(4 + c * c)
        pred += [r, -r]
    pred = np.sort(np.array(pred))
    err = np.max(np.abs(im - pred))
    check(len(pred) == 2 * k + 1 and err < 1e-7, f"spectrum k={k} err={err}")
    nonc = [abs(x) for x in im if abs(abs(x) - math.sqrt(3)) > 1e-9 and abs(x) > 1e-9]
    check(all(x > 2 for x in nonc), "non-core eigenvalue modulus <= 2")
    if k in (2, 3, 6, 14):
        print(f"k={k:2d} (order {2*k+1}): |eig| = {np.round(np.unique(np.round(np.abs(im), 10)), 4).tolist()}  max err vs cot^2 formula {err:.1e}")
print("(2) spectrum of Gamma^k(K1) = {0, +-i sqrt3} u {+-i sqrt(4 + cot^2((2j+1)pi/(2(2k-1)))), j < k-1}: k = 1..14")

# ---- gate tower limits ----
print("(2) gate towers: F_(Gamma^k T)(z) - (2+s)/z, exact rational F, cores of orders 1..6")
for zf in (Fraction(1, 2), Fraction(1), Fraction(2), Fraction(3), Fraction(4)):
    zz = float(zf)
    sv = math.sqrt(zz * zz + 4)
    lim = (2 + sv) / zz
    ratio = ((sv - 1) / (sv + 1)) ** 2
    worst_rate = 0.0
    for trial in range(4):
        A = rand_tour(rng.randint(1, 6), rng)
        X = A
        errs = []
        for k in range(0, 11):
            if k > 0:
                X = gate(X)
            Fv = F_at(X, zf)
            errs.append(abs(float(Fv) - lim))
        # observed contraction per gate
        # the per-gate contraction tends to ((s-1)/(s+1))^2 with relative correction O(ratio^k): use k = 8 -> 9 -> 10
        obs = [errs[i + 1] / errs[i] for i in range(8, 10) if errs[i] > 1e-14]
        if obs:
            worst_rate = max(worst_rate, max(abs(o / ratio - 1) for o in obs))
        check(errs[-1] < 10 * errs[0] * ratio ** 10 + 1e-12, f"rate z={zz}")
    print(f"   z = {zz}: limit {lim:.10f}; predicted contraction ((s-1)/(s+1))^2 = {ratio:.6f}; "
          f"max relative deviation of observed contraction (k = 8..10) = {worst_rate:.2e}")
    check(worst_rate < 0.02, f"contraction mismatch z={zz}")

# ---- positivity of F_T(z) on z > 0 ----
neg = 0
for trial in range(300):
    A = rand_tour(rng.randint(1, 9), rng)
    zf = Fraction(rng.randint(1, 400), rng.randint(1, 100))
    if F_at(A, zf) <= 0:
        neg += 1
check(neg == 0, "F_T(z) <= 0 for some z > 0")
print(f"(2) F_T(z) > 0 at 300 random (T, z > 0): violations {neg}")

# ---- D-tower ----
print("(3) D-tower: F_(D^k T)(z) vs (z/w_k) F_T(w_k), w_k^2 - 1 = (z^2-1)/2^k, and limit z phi(T)")
for zz in (0.5, 2.0, 3.0):
    A = rand_tour(4, rng)
    ph = float(phi(A))
    X = A
    for k in range(1, 5):
        X = double(X)
        M = np.array(dominance(X), dtype=float)
        n = len(X)
        Fd = float(np.ones(n) @ np.linalg.solve(zz * np.eye(n) - M, np.ones(n)))
        wk = math.sqrt(1 + (zz * zz - 1) / 2 ** k)
        M0 = np.array(dominance(A), dtype=float)
        F0 = float(np.ones(4) @ np.linalg.solve(wk * np.eye(4) - M0, np.ones(4)))
        check(abs(Fd - zz / wk * F0) < 1e-9, f"D tower formula z={zz} k={k}")
    print(f"   z = {zz}: F(D^4 T) = {Fd:.6f}, (z/w_4)F_T(w_4) = {zz / wk * F0:.6f}, z*phi(T) = {zz * ph:.6f}")
    # limit via the (proved) formula for large k
    wbig = math.sqrt(1 + (zz * zz - 1) / 2 ** 40)
    F0b = float(np.ones(4) @ np.linalg.solve(wbig * np.eye(4) - M0, np.ones(4)))
    check(abs(zz / wbig * F0b - zz * ph) < 1e-8, "D tower limit")

# ---- E-tower ----
for zz in (0.5, 2.0, 3.0):
    A = rand_tour(5, rng)
    X = A
    for k in range(1, 13):
        X = pad(X)
    n = len(X)
    M = np.array(dominance(X), dtype=float)
    Fe = float(np.ones(n) @ np.linalg.solve(zz * np.eye(n) - M, np.ones(n)))
    check(abs(Fe - 1) < 1e-4, f"E tower z={zz}")
    print(f"(3) E-tower z = {zz}: F(E^12 T) = {Fe:.8f}")

print("FAILURES:", len(fails))
