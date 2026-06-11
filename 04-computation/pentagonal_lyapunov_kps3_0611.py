#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-485 lab: the pentagonal Lyapunov gap and the sign-rigidity of Euler's
recurrence.  kind-pasteur-2026-06-11-S3 (HYP-2416/2417/2418).

Euler: p(n) = sum_{k>=1} eps_k [ p(n - g_k) + p(n - gbar_k) ],
  g_k = k(3k-1)/2, gbar_k = k(3k+1)/2, eps_k = (-1)^{k+1}, p(0) = 1.
Generating function: 1/F(x), F(x) = 1 - sum eps_k (x^{g_k} + x^{gbar_k})
                                  = prod_{n>=1} (1 - x^n)   [pentagonal thm].
Euler's signs: subexponential growth (rate 0).  This lab measures what every
other sign regime does:

  A. all-plus rate:  a_n = sum_k [a(n-g_k)+a(n-gbar_k)]: rate = -log x*,
     where sum_k (x*^{g_k} + x*^{gbar_k}) = 1 (the positive fixed point).
  B. fresh-iid-signs Lyapunov constant: signs sigma_k(n) iid uniform ±1,
     resampled for every (n,k) use: gamma = lim E[log|a_n|]/n  (renormalized
     accumulation; several independent runs; report mean ± spread).
  C. k-FIXED random signs: eps: pairs -> ±1 drawn once; growth rate is then
     deterministic per eps = -log(radius of convergence of 1/F_eps); measured
     empirically as slope of log|a_n|.
  D. RIGIDITY SWEEP (HYP-2417): all 2^K sign patterns on the first K pairs
     (signs eps_k for k > K continued by Euler's alternation so that F_eps is
     a bounded perturbation of the eta product): measure empirical rates;
     count how many achieve rate ~ 0.  Prediction: only Euler's (and possibly
     a symmetry partner).
  E. mod-2 floor (HYP-2418): all sign regimes agree with p(n) mod 2 —
     verified directly (signs invisible over F_2).

Honest numerics: integer arithmetic where possible (C: exact big ints for
fixed signs; rates from log|a_n|/n over windows); for B, float with explicit
renormalization (track log-scale offset; MISTAKE-067 guard: no overflow
artifacts).  All randomness seeded.
"""

import math, random, time
from fractions import Fraction


def pent_pairs(N):
    """All (g_k, gbar_k, k) with g_k <= N."""
    out = []
    k = 1
    while k * (3 * k - 1) // 2 <= N:
        out.append((k * (3 * k - 1) // 2, k * (3 * k + 1) // 2, k))
        k += 1
    return out


# ---------------------------------------------------------------- A. all-plus

def allplus_rate(tol=1e-12):
    """Solve sum_k (x^{g_k} + x^{gbar_k}) = 1 for x in (0,1) by bisection."""
    def f(x):
        s = 0.0
        k = 1
        while True:
            g = k * (3 * k - 1) // 2
            gb = k * (3 * k + 1) // 2
            t = x ** g + x ** gb
            s += t
            if t < 1e-18 or k > 10000:
                break
            k += 1
        return s - 1.0
    lo, hi = 1e-9, 1.0 - 1e-9
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    xstar = (lo + hi) / 2
    return xstar, -math.log(xstar)


def allplus_empirical(N=4000):
    """Exact integer all-plus recurrence; slope of log a_n."""
    a = [0] * (N + 1)
    a[0] = 1
    pp = pent_pairs(N)
    for n in range(1, N + 1):
        s = 0
        for g, gb, k in pp:
            if g > n:
                break
            s += a[n - g]
            if gb <= n:
                s += a[n - gb]
        a[n] = s
    l1 = math.log(a[N])
    l0 = math.log(a[N // 2])
    return (l1 - l0) / (N - N // 2)


# ----------------------------------------------- B. fresh-iid Lyapunov constant

def fresh_lyapunov(N=20000, runs=8, seed=20260611):
    """Fresh iid sign per (n,k) use; renormalized float accumulation.
    Returns per-run gamma estimates (slope of E log|a_n| over the second half)."""
    rates = []
    for r in range(runs):
        rng = random.Random(seed + r)
        vals = [0.0] * (N + 1)   # renormalized values
        logsc = [0.0] * (N + 1)  # log-scale offsets: a_n = vals[n] * exp(logsc[n])
        vals[0] = 1.0
        pp = pent_pairs(N)
        loga = [0.0] * (N + 1)
        for n in range(1, N + 1):
            # accumulate in the scale of the largest contributing term
            mx = -1e18
            terms = []
            for g, gb, k in pp:
                if g > n:
                    break
                s1 = 1 if rng.random() < 0.5 else -1
                terms.append((s1, n - g))
                if logsc[n - g] > mx and vals[n - g] != 0.0:
                    mx = logsc[n - g]
                if gb <= n:
                    s2 = 1 if rng.random() < 0.5 else -1
                    terms.append((s2, n - gb))
                    if logsc[n - gb] > mx and vals[n - gb] != 0.0:
                        mx = logsc[n - gb]
            if mx < -1e17:
                mx = 0.0
            tot = 0.0
            for s, idx in terms:
                if vals[idx] != 0.0:
                    tot += s * vals[idx] * math.exp(logsc[idx] - mx)
            if tot == 0.0:
                vals[n] = 0.0
                logsc[n] = mx
                loga[n] = -1e18
            else:
                logsc[n] = mx + math.log(abs(tot))
                vals[n] = 1.0 if tot > 0 else -1.0
                loga[n] = logsc[n]
        # slope over second half (robust: least squares on (n, loga))
        xs = list(range(N // 2, N + 1))
        ys = [loga[n] for n in xs]
        npts = len(xs)
        mx_ = sum(xs) / npts
        my_ = sum(ys) / npts
        slope = sum((x - mx_) * (y - my_) for x, y in zip(xs, ys)) / sum((x - mx_) ** 2 for x in xs)
        rates.append(slope)
    return rates


# --------------------------------------------- C/D. fixed-sign patterns, exact

def fixed_signs_rate(eps_fn, N=3000):
    """Exact integer recurrence with k-dependent fixed signs eps_fn(k) in {+1,-1};
    returns empirical rate (least-squares slope of log|a_n| on the second half,
    skipping zeros) and whether a_n matches p(n) when eps = Euler."""
    a = [0] * (N + 1)
    a[0] = 1
    pp = pent_pairs(N)
    for n in range(1, N + 1):
        s = 0
        for g, gb, k in pp:
            if g > n:
                break
            e = eps_fn(k)
            s += e * a[n - g]
            if gb <= n:
                s += e * a[n - gb]
        a[n] = s
    xs, ys = [], []
    for n in range(N // 2, N + 1):
        if a[n] != 0:
            xs.append(n)
            ys.append(math.log(abs(a[n])))
    if len(xs) < 10:
        return 0.0, 0.0, a
    # two-regressor least squares: log|a_n| ~ alpha*n + beta*sqrt(n) + c
    # subexponential (Euler-like) <=> alpha ~ 0, beta ~ pi*sqrt(2/3) = 2.5651
    import statistics
    n_ = len(xs)
    X1 = xs
    X2 = [math.sqrt(x) for x in xs]
    m1, m2, my = sum(X1) / n_, sum(X2) / n_, sum(ys) / n_
    s11 = sum((x - m1) ** 2 for x in X1)
    s22 = sum((x - m2) ** 2 for x in X2)
    s12 = sum((x - m1) * (z - m2) for x, z in zip(X1, X2))
    s1y = sum((x - m1) * (y - my) for x, y in zip(X1, ys))
    s2y = sum((z - m2) * (y - my) for z, y in zip(X2, ys))
    det = s11 * s22 - s12 * s12
    alpha = (s1y * s22 - s2y * s12) / det
    beta = (s2y * s11 - s1y * s12) / det
    return alpha, beta, a


def euler_eps(k):
    return 1 if k % 2 == 1 else -1


def rigidity_sweep(K=6, N=2000):
    """All 2^K sign choices on k <= K, Euler continuation beyond; rates."""
    results = []
    for mask in range(1 << K):
        def eps(k, mask=mask):
            if k <= K:
                return 1 if (mask >> (k - 1)) & 1 else -1
            return euler_eps(k)
        alpha, beta, _ = fixed_signs_rate(eps, N)
        results.append((mask, alpha, beta))
    return results


# ----------------------------------------------------------- E. the mod-2 floor

def mod2_floor(N=3000):
    """Verify: recurrence mod 2 is sign-free and reproduces p(n) mod 2."""
    a = [0] * (N + 1)
    a[0] = 1
    pp = pent_pairs(N)
    for n in range(1, N + 1):
        s = 0
        for g, gb, k in pp:
            if g > n:
                break
            s ^= a[n - g]
            if gb <= n:
                s ^= a[n - gb]
        a[n] = s
    # independent p(n) mod 2 via exact Euler recurrence
    p = [0] * (N + 1)
    p[0] = 1
    for n in range(1, N + 1):
        s = 0
        for g, gb, k in pp:
            if g > n:
                break
            e = euler_eps(k)
            s += e * p[n - g]
            if gb <= n:
                s += e * p[n - gb]
        p[n] = s
    mism = sum(1 for n in range(N + 1) if a[n] != p[n] % 2)
    ones = sum(a)
    return mism, ones / (N + 1)


def main():
    t0 = time.time()
    print("=== A. all-plus fixed point and rate ===", flush=True)
    xstar, rate = allplus_rate()
    print(f"   x* = {xstar:.12f}, gamma_+ = -log x* = {rate:.12f}", flush=True)
    emp = allplus_empirical(4000)
    print(f"   empirical all-plus slope (exact ints, n in [2000,4000]): {emp:.12f}"
          f"   (match: {abs(emp-rate) < 1e-3})", flush=True)

    print("\n=== B. fresh-iid random-sign Lyapunov constant gamma_pent ===", flush=True)
    rates = fresh_lyapunov(N=20000, runs=8)
    m = sum(rates) / len(rates)
    sd = (sum((r - m) ** 2 for r in rates) / (len(rates) - 1)) ** 0.5
    print(f"   per-run gammas: {[f'{r:.6f}' for r in rates]}", flush=True)
    print(f"   gamma_pent = {m:.6f} +- {sd:.6f}  (vs gamma_+ = {rate:.6f}, Euler = 0)", flush=True)

    print("\n=== C. Euler control + simple variants (exact ints, N=3000) ===", flush=True)
    a_e, b_e, arr_e = fixed_signs_rate(euler_eps, 3000)
    print(f"   Euler signs: alpha = {a_e:.6f}, beta = {b_e:.4f} (theory: alpha=0, beta=pi*sqrt(2/3)={3.14159265*math.sqrt(2/3):.4f})", flush=True)
    ok_p = all(arr_e[n] >= 0 for n in range(3001))
    print(f"   Euler a_n = p(n) >= 0 check: {ok_p}", flush=True)
    a_f, b_f, _ = fixed_signs_rate(lambda k: -euler_eps(k), 3000)
    print(f"   global flip: alpha = {a_f:.6f}, beta = {b_f:.4f}", flush=True)

    print("\n=== D. rigidity sweep: all 2^6 sign patterns on k<=6, Euler beyond ===", flush=True)
    res = rigidity_sweep(K=6, N=2000)
    res.sort(key=lambda t: t[1])
    euler_mask = sum(1 << (k - 1) for k in range(1, 7) if euler_eps(k) == 1)
    print(f"   Euler mask = {euler_mask:06b}", flush=True)
    zeroish = [(m, al, be) for m, al, be in res if al < 0.005]
    print(f"   patterns with alpha < 0.005 (subexponential): {len(zeroish)} of 64:", flush=True)
    for m_, al_, be_ in zeroish:
        tag = " <-- EULER" if m_ == euler_mask else ""
        print(f"      mask {m_:06b}  alpha {al_:.6f}  beta {be_:.4f}{tag}", flush=True)
    print(f"   alpha range over the sweep: [{res[0][1]:.6f}, {res[-1][1]:.6f}]", flush=True)

    print("\n=== E. mod-2 floor (signs invisible over F_2) ===", flush=True)
    mism, dens = mod2_floor(3000)
    print(f"   recurrence-mod-2 vs p(n) mod 2: {mism} mismatches over n <= 3000; "
          f"odd density {dens:.4f} (Parkin-Shanks ~ 0.5)", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
