#!/usr/bin/env python3
"""Audit E, task B: step-by-step numerical checks of THM-4606 Proof 2 (escape lemma).

 (B1) P(tau > 2m) = C(2m,m)/4^m exactly (m <= 30, by path enumeration DP) and C(2m,m)/4^m >= 1/(2 sqrt m) (m <= 10^6).
 (B2) Exact law of V_n = #visits to 0 among S_0..S_{n-1} (SRW from 0) for n <= 400 via the renewal convolution;
      check P(V_n >= m) <= exp(-(m-1)/(2 sqrt n)) for every m, n.
 (B3) ln(1-x) >= -2x on [0, 1/2].
 (B4) The constants: eps = Lambda/(3 ln p), delta_1(n_0), delta_2(C'), and the size of A_p(1/2) (log10), p = 5, 7, 9, 13.
 (B5) Along simulated chains driven by ACTUAL orbit parities: xi_n fair and serially uncorrelated (departure and
      non-departure steps separately), ln lambda_n >= xi_n - ln p * 1{dep}, D_n <= V_n, and at every step with
      |f_s| >= 4: |f_{s+1}| >= lambda_s|f_s| - 1 and ln|f_{s+1}| >= ln|f_s| + ln lambda_s - 2/(lambda_s|f_s|).
 (B6) Empirical escape: P(absorbed within 3000 steps) and P(|f_n| >= 4 exp(Lambda n/3) for all n <= 3000) from (0, e0)
      for growing e0 (the lemma predicts -> 0 and -> 1 respectively, at a rate it does not quantify)."""
import math, random
from fractions import Fraction as Fr

rnd = random.Random(77)

# ---------------- B1 ----------------
def b1():
    # exact P(tau > 2m) by DP over SRW paths avoiding 0 after time 0
    ok = True
    for m in range(1, 31):
        n = 2 * m
        dist = {1: Fr(1, 2), -1: Fr(1, 2)}
        for t in range(1, n):
            nd = {}
            for x, pr in dist.items():
                for d in (1, -1):
                    y = x + d
                    if y == 0: continue
                    nd[y] = nd.get(y, 0) + pr / 2
            dist = nd
        surv = sum(dist.values())
        cm = Fr(math.comb(2 * m, m), 4 ** m)
        if surv != cm: ok = False; print("B1 identity FAIL", m)
    # inequality for m up to 10^6 via the product (2k-1)/(2k)
    r = 1.0; worst = 1e9; worst_m = None
    for m in range(1, 10**6 + 1):
        r *= (2 * m - 1) / (2 * m)
        ratio = r * 2 * math.sqrt(m)
        if ratio < worst: worst, worst_m = ratio, m
    print(f"(B1) P(tau>2m) = C(2m,m)/4^m exact for m<=30: {'PASS' if ok else 'FAIL'}; "
          f"min over m<=1e6 of [C(2m,m)/4^m]/[1/(2 sqrt m)] = {worst:.6f} at m={worst_m} (>=1 needed; -> 2/sqrt(pi)=1.128)")
    return ok and worst >= 1

# ---------------- B2 ----------------
def b2(nmax=400):
    # first-return law f(2k) = C(2k,k)/((2k-1) 4^k)
    f = [0.0] * (nmax + 1)
    for k in range(1, nmax // 2 + 1):
        f[2 * k] = math.comb(2 * k, k) / ((2 * k - 1) * 4 ** k)
    # g_j(t) = P(R_j = t), R_j = time of j-th return
    worst = 0.0; worst_at = None
    g = [0.0] * (nmax + 1); g[0] = 1.0
    cdf_by_j = [None]
    j = 0
    # P(V_n >= m) = P(R_{m-1} <= n-1)
    cdfs = []
    cdfs.append([1.0] * (nmax + 1))          # j = 0: R_0 = 0 <= anything
    while j < nmax // 2 + 1:
        j += 1
        ng = [0.0] * (nmax + 1)
        for t in range(0, nmax + 1, 2):
            if g[t] == 0.0: continue
            for s in range(2, nmax + 1 - t, 2):
                ng[t + s] += g[t] * f[s]
        g = ng
        c = []; acc = 0.0
        for t in range(nmax + 1):
            acc += g[t]; c.append(acc)
        cdfs.append(c)
        if acc < 1e-300: break
    for n in range(1, nmax + 1):
        for m in range(1, len(cdfs) + 1):
            pv = cdfs[m - 1][n - 1] if m - 1 < len(cdfs) else 0.0
            bound = math.exp(-(m - 1) / (2 * math.sqrt(n)))
            if pv > 0 and pv / bound > worst: worst, worst_at = pv / bound, (n, m)
    print(f"(B2) exact SRW visit law, n <= {nmax}: max over (n,m) of P(V_n >= m) / exp(-(m-1)/(2 sqrt n)) = {worst:.4f} at (n,m)={worst_at}"
          f" ({'PASS' if worst <= 1 + 1e-12 else 'FAIL'}: must be <= 1)")
    return worst <= 1 + 1e-12

# ---------------- B3 ----------------
def b3():
    ok = all(math.log(1 - x) >= -2 * x - 1e-15 for x in [i / 100000 * 0.5 for i in range(100001)])
    print(f"(B3) ln(1-x) >= -2x on [0,1/2] (grid of 1e5 points; also concavity argument): {'PASS' if ok else 'FAIL'}")
    return ok

# ---------------- B4 ----------------
def delta1(eps, n0):
    # sum_{n > n0} exp(-eps sqrt(n)/2): tail by integral bound int_{n0}^inf exp(-c sqrt x) dx = (2/c^2)(1+c sqrt n0) e^{-c sqrt n0}
    c = eps / 2
    s = math.sqrt(n0)
    return (2 / c**2) * (1 + c * s) * math.exp(-c * s)

def delta2(L, lnp, C):
    tot = 0.0
    n = 1
    while True:
        term = math.exp(-2 * (L * n / 3 + C) ** 2 / (n * lnp ** 2))
        tot += term
        if n > 10 and term < 1e-30 * max(tot, 1e-300) and L * n / 3 > C: break
        n += 1
        if n > 10**8: break
    return tot

def b4():
    for p in (5, 7, 9, 13):
        L = 0.5 * math.log(p / 4); lnp = math.log(p); eps = L / (3 * lnp)
        # smallest n0 (on a coarse doubling grid then bisection) with delta1 <= 1/4
        lo, hi = 1, 2
        while delta1(eps, hi) > 0.25: hi *= 2
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if delta1(eps, mid) > 0.25: lo = mid
            else: hi = mid
        n0 = hi
        C = 1.0
        while delta2(L, lnp, C) > 0.25: C *= 1.5
        A_log10 = (math.log(4) + C + n0 * lnp + 1 + math.log(max(1, 4 / (1 - math.exp(-L / 3))))) / math.log(10)
        print(f"(B4) p={p:2d}: Lambda={L:.5f}, eps={eps:.5f}, n0(delta1<=1/4)={n0:,}, C'(delta2<=1/4)~{C:.1f};"
              f"  A_p(1/2) = 10^{A_log10:,.0f}  (finite; purely qualitative)")

# ---------------- chain (integer coordinates; f = E / p^|k|) ----------------
def step(p, k, E, b):
    s = E & 1
    if s == 0:
        if b == 0: return k, E >> 1
        if k >= 0: return k, (p * E + 1 - p ** k) >> 1
        return k, (p * E + p ** (-k) - 1) >> 1
    if b == 0:
        if k >= 0: return k + 1, (p * E + 1) >> 1
        return k + 1, (E + p ** (-k - 1)) >> 1
    if k >= 1: return k - 1, (E - p ** (k - 1)) >> 1
    return k - 1, (p * E - 1) >> 1

def lam_xi(p, k, s, b):
    """(lambda, xi) per Proof 1/2."""
    if s == 0:
        lam = 0.5 if b == 0 else p / 2
        return lam, math.log(lam)
    if k == 0:
        return 0.5, (math.log(0.5) if b == 0 else math.log(p / 2))
    if k >= 1:
        lam = 0.5 if b == 0 else p / 2
    else:
        lam = p / 2 if b == 0 else 0.5
    return lam, math.log(lam)

def Tint(p, x):
    return x >> 1 if x % 2 == 0 else (p * x + 1) >> 1

def b5(npaths=600, nsteps=1500):
    ok = True
    stats = {}
    viol_ineq = viol_D = viol_lam = 0
    nineq = 0
    for p in (5, 7, 9):
        hi_dep = hi_non = n_dep = n_non = 0
        pairs = [0, 0, 0, 0]
        for _ in range(npaths):
            e0 = rnd.choice([1, 2, 3, 5, 10, 100, 10**4])
            k, E = 0, e0
            v = rnd.getrandbits(nsteps + 200)
            D = 0; visits = 0; flips = 0
            prev = None
            for n in range(nsteps):
                b = v & 1; v = Tint(p, v)
                s = E & 1
                f = Fr(E, p ** abs(k))
                lam, xi = lam_xi(p, k, s, b)
                dep = (s == 1 and k == 0)
                if s == 1:
                    flips += 1
                    if k == 0: visits += 1
                if dep: D += 1
                if D > visits: viol_D += 1
                if math.log(lam) < xi - math.log(p) * dep - 1e-12: viol_lam += 1
                hi = (xi > 0)
                if dep: n_dep += 1; hi_dep += hi
                else: n_non += 1; hi_non += hi
                if prev is not None: pairs[2 * prev + hi] += 1
                prev = hi
                k2, E2 = step(p, k, E, b)
                f2 = Fr(E2, p ** abs(k2))
                if abs(f) >= 4:
                    nineq += 1
                    lhs = abs(f2)
                    lamF = Fr(1, 2) if lam == 0.5 else Fr(p, 2)
                    if lhs < lamF * abs(f) - 1:          # exact rational check of |f'| >= lambda|f| - 1
                        viol_ineq += 1
                    if math.log(float(lhs)) < math.log(float(abs(f))) + math.log(lam) - 2 / (lam * float(abs(f))) - 1e-12:
                        viol_ineq += 1
                k, E = k2, E2
                if k == 0 and E == 0: break
                if abs(E) > p ** abs(k) * 10**40: break
        tot = sum(pairs)
        corr = (pairs[0] + pairs[3] - pairs[1] - pairs[2]) / max(tot, 1)
        print(f"(B5) p={p}: P(xi = ln(p/2)) at departures {hi_dep/max(n_dep,1):.4f} (n={n_dep}),"
              f" off departures {hi_non/max(n_non,1):.4f} (n={n_non}); lag-1 corr of xi {corr:+.4f}")
    print(f"(B5) violations: D_n > V_n: {viol_D}; ln lambda < xi - ln p 1{{dep}}: {viol_lam};"
          f" induction step inequalities at |f|>=4: {viol_ineq} of {nineq} steps checked")
    return viol_D == viol_lam == viol_ineq == 0

def b6(nsamp=2000, nsteps=3000):
    for p in (5, 7, 9):
        L = 0.5 * math.log(p / 4)
        row = []
        for e0 in (1, 4, 16, 256, 4096, 2**20, 2**40):
            absorbed = 0; envelope_ok = 0
            for _ in range(nsamp):
                k, E = 0, e0
                good = True; ab = False
                for n in range(nsteps):
                    f = abs(E) / p ** abs(k)
                    if f < 4 * math.exp(L * n / 3): good = False
                    k, E = step(p, k, E, rnd.getrandbits(1))
                    if k == 0 and E == 0: ab = True; break
                    if abs(E) > p ** abs(k) * 10**200: break    # far beyond any return (and envelope)
                absorbed += ab; envelope_ok += good and not ab
            row.append(f"e0={e0}: absorbed {absorbed/nsamp:.4f}, envelope held {envelope_ok/nsamp:.3f}")
        print(f"(B6) p={p}: " + "; ".join(row))

if __name__ == '__main__':
    ok = b1() & b2() & b3()
    b4()
    ok &= b5()
    b6()
    print("ALL B CHECKS PASS" if ok else "SOME B CHECK FAILED")
