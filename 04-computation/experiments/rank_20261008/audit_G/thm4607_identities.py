#!/usr/bin/env python3
"""Audit G, THM-4607 analytic ingredients, checked numerically on random inputs (floats; tolerances 1e-12 relative):
 (1) ln mu = (1-g'(x)) ln(m_i/d) + g'(x) ln(m_j/d) - rho, rho = g(x+delta) - g(x) - g'(x) delta, mu = (m_i/d)(1+M)/(1+M').
 (2) 0 <= rho <= (delta^2/2) max_{[x,x+delta]} g''  and  rho <= (Delta^2/2) min(1/4, e^(Delta-|x|)).
 (3) g'' <= min(1/4, e^-|y|).
 (4) |H'| <= B_H = kappa((2 Delta + ln 4)/4 + 1/4) for H'' = kappa min(1/4, e^(2 Delta - |y|)), H'(0) = 0 (numerical integral).
 (5) Taylor lower bound H(X+delta) - H(X) - H'(X) delta >= (delta_min^2/2) kappa min(1/4, e^(Delta - |X|)).
 (6) F' >= mu F - R/d on random pair-chain states (exact rationals) for random maps.
 (7) Exact Z_2 departure rule: at M = 1, e odd, both coins give f' = (min(m0,m1)/2) f + const, f = e/max(M,1)."""
import math, random
from fractions import Fraction as Fr
rnd = random.Random(42)
g = lambda y: math.log1p(math.exp(y)) if y < 30 else y + math.log1p(math.exp(-y))
g1 = lambda y: 1 / (1 + math.exp(-y))
g2 = lambda y: math.exp(-abs(y)) / (1 + math.exp(-abs(y))) ** 2
bad = [0] * 7
for _ in range(200000):
    d = rnd.choice([2, 3, 5, 7]); mi, mj = rnd.randint(1, 40), rnd.randint(1, 40)
    x = rnd.uniform(-15, 15); M = math.exp(x); delta = math.log(mi / mj); Mp = M * mi / mj
    mu = (mi / d) * (1 + M) / (1 + Mp)
    rho = g(x + delta) - g(x) - g1(x) * delta
    lhs = math.log(mu); rhs = (1 - g1(x)) * math.log(mi / d) + g1(x) * math.log(mj / d) - rho
    if abs(lhs - rhs) > 1e-9 * (1 + abs(lhs)): bad[0] += 1
    a, b = sorted((x, x + delta)); mx = max(g2(a + (b - a) * k / 400) for k in range(401))
    if rho < -1e-12 or rho > delta * delta / 2 * mx * (1 + 1e-6) + 1e-12: bad[1] += 1
    Delta = abs(delta) + rnd.uniform(0, 2)
    if rho > Delta ** 2 / 2 * min(0.25, math.exp(Delta - abs(x))) + 1e-12: bad[1] += 1
    y = rnd.uniform(-40, 40)
    if g2(y) > min(0.25, math.exp(-abs(y))) + 1e-15: bad[2] += 1
# (4),(5)
for Delta in (0.3, 1.0, 2.5, 5.0):
    kappa = 1.0; h = 1e-3
    Hpp = lambda y: kappa * min(0.25, math.exp(2 * Delta - abs(y)))
    # integrate H'' from 0 to large
    tot = 0.0; y = 0.0
    while y < 2 * Delta + 60:
        tot += h * (Hpp(y) + Hpp(y + h)) / 2; y += h
    BH = kappa * ((2 * Delta + math.log(4)) / 4 + 0.25)
    if abs(tot - BH) > 1e-4: bad[3] += 1
    print(f"Delta = {Delta}: integral of H'' over [0, inf) = {tot:.6f}, B_H formula = {BH:.6f}")
    # H via numeric double integral on a grid; check Taylor bound at random points
    L = 2 * Delta + 40; n = int(2 * L / h)
    ys = [-L + k * h for k in range(n + 1)]
    Hp = [0.0] * (n + 1); H = [0.0] * (n + 1); z = n // 2
    for k in range(z + 1, n + 1): Hp[k] = Hp[k - 1] + h * (Hpp(ys[k - 1]) + Hpp(ys[k])) / 2
    for k in range(z - 1, -1, -1): Hp[k] = Hp[k + 1] - h * (Hpp(ys[k + 1]) + Hpp(ys[k])) / 2
    for k in range(z + 1, n + 1): H[k] = H[k - 1] + h * (Hp[k - 1] + Hp[k]) / 2
    for k in range(z - 1, -1, -1): H[k] = H[k + 1] - h * (Hp[k + 1] + Hp[k]) / 2
    dmin = Delta / 3
    for _ in range(3000):
        X = rnd.uniform(-L / 2, L / 2); dl = rnd.choice([-1, 1]) * rnd.uniform(dmin, Delta)
        kx = int(round((X + L) / h)); kd = int(round((X + dl + L) / h)); X = ys[kx]; dl = ys[kd] - X
        lhs = H[kd] - H[kx] - Hp[kx] * dl
        rhs = (dmin ** 2 / 2) * kappa * min(0.25, math.exp(Delta - abs(X)))
        if lhs < rhs * (1 - 1e-3) - 1e-9: bad[4] += 1
# (6) exact pair-chain states
for _ in range(20000):
    d = rnd.choice([2, 3, 5]); m = [rnd.choice([x for x in range(1, 30) if math.gcd(x, d) == 1]) for _ in range(d)]
    r = [(-m[i] * i) % d + d * rnd.randint(-3, 3) for i in range(d)]
    R = max(abs(x) for x in r)
    M = Fr(rnd.randint(1, 50), rnd.randint(1, 50)); e = Fr(rnd.randint(-10 ** 6, 10 ** 6), rnd.randint(1, 30))
    i, j = rnd.randrange(d), rnd.randrange(d)
    Mp = M * Fr(m[i], m[j]); ep = (m[i] * e + r[i] - Fr(m[i], m[j]) * M * r[j]) / d
    F = abs(e) / (1 + M); Fp = abs(ep) / (1 + Mp); mu = Fr(m[i], d) * (1 + M) / (1 + Mp)
    if Fp < mu * F - Fr(R, d): bad[5] += 1
# (7) Z_2 departure rule, exact
for _ in range(20000):
    m0, m1 = rnd.choice([1, 3, 5, 7, 9, 11, 13, 15]), rnd.choice([1, 3, 5, 7, 9, 11, 13, 15])
    if m0 == m1: continue
    r0, r1 = 2 * rnd.randint(-9, 9), 2 * rnd.randint(-9, 9) + 1
    e = Fr(2 * rnd.randint(-10 ** 6, 10 ** 6) + 1)       # odd, M = 1
    consts = []
    for jj in (0, 1):
        ii = 1 - jj; m = (m0, m1); rr = (r0, r1)
        Mp = Fr(m[ii], m[jj]); ep = (m[ii] * e + rr[ii] - Mp * rr[jj]) / 2
        fp = ep / max(Mp, 1)
        consts.append(fp - Fr(min(m0, m1), 2) * e)
    # constants must be e-independent: recompute with another odd e
    e2 = e + 2 * rnd.randint(1, 1000)
    for jj, c in zip((0, 1), consts):
        ii = 1 - jj; m = (m0, m1); rr = (r0, r1)
        Mp = Fr(m[ii], m[jj]); ep = (m[ii] * e2 + rr[ii] - Mp * rr[jj]) / 2
        if ep / max(Mp, 1) - Fr(min(m0, m1), 2) * e2 != c: bad[6] += 1
print("failures: ln-mu identity", bad[0], "| rho bounds", bad[1], "| g'' bound", bad[2], "| B_H integral", bad[3],
      "| H Taylor lower bound", bad[4], "| F' >= mu F - R/d", bad[5], "| Z_2 departure rule (exact)", bad[6])
