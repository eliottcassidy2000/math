#!/usr/bin/env python3
"""General 2-adic Matthews-Watts map T(x) = (m_0 x + r_0)/2 (x even), (m_1 x + r_1)/2 (x odd); m_0, m_1 odd,
r_0 even, r_1 odd.  Pair chain u = M v + e, M = q^k with q = m_1/m_0.  Normalized offset f = e / max(M, 1).
Claim (checked exactly): every step has f' = lambda f + eta with |eta| <= (|r_0| + |r_1|)/2 * C, lambda = v's multiplier
m_beta/2 when M' > 1, u's multiplier when M' <= 1; off departures (flips at k = 0) the two coins give the two
multipliers m_0/2, m_1/2 (one each); at a departure both coins give min(m_0, m_1)/2 (the smaller multiplier; = m_0/2 when m_0 < m_1)."""
import random
from fractions import Fraction as Fr
def par(z): return (z.numerator * pow(z.denominator, -1, 2)) % 2
def step(m0, m1, r0, r1, M, e, beta):
    sig = par(e); m = (m0, m1); r = (r0, r1)
    j = beta; i = beta ^ sig              # v's branch j, u's branch i (u = M v + e, M a 2-adic unit)
    Mn = M * Fr(m[i], m[j])
    en = (m[i] * e + r[i] - Fr(m[i], m[j]) * M * r[j]) / 2
    return Mn, en, i, j
rnd = random.Random(11); bad = 0; checks = 0; dep_ok = 0; dep_n = 0
for _ in range(40000):
    m0 = rnd.choice([1, 3, 5, 7, 9, 11, 13, 15]); m1 = rnd.choice([1, 3, 5, 7, 9, 11, 13, 15])
    if m0 == m1: continue
    r0 = 2 * rnd.randint(-6, 6); r1 = 2 * rnd.randint(-6, 6) + 1
    q = Fr(m1, m0); k = rnd.randint(-5, 5); M = q ** k
    den = (m0 * m1) ** max(0, abs(k)) * rnd.choice([1, 3, 5])
    den = den if den % 2 else den + 0
    e = Fr(rnd.randint(-10**6, 10**6), 1 if rnd.random() < 0.5 else (m1 if k < 0 else 1) ** abs(k))
    f = e / max(M, 1)
    lams = []
    for beta in (0, 1):
        Mn, en, i, j = step(m0, m1, r0, r1, M, e, beta)
        fn = en / max(Mn, 1)
        use = (m0, m1)[j] if Mn > 1 else (m0, m1)[i]
        lam = Fr(use, 2) * (max(M, 1) / max(Mn, 1)) * (Mn if (Mn > 1) != (M > 1) and False else 1)
        # predicted multiplier on f: v's multiplier if M' > 1 (and M > 1), u's if both <= 1; departures: m0/2
        eta = fn - Fr(use, 2) * f * (1 if (M > 1) == (Mn > 1) else (max(M, 1) * (Fr(m1, m0) if Mn > M else Fr(m0, m1)) / max(Mn, 1)) * 1)
        lams.append((use, i, j, Mn, M))
        checks += 1
    # departure check: at k = 0 with e odd, both coins should give m0/2 as the effective multiplier of f
    if k == 0 and par(e) == 1:
        dep_n += 1
        outs = []
        for beta in (0, 1):
            Mn, en, i, j = step(m0, m1, r0, r1, M, e, beta)
            fn = en / max(Mn, 1)
            outs.append(fn - Fr(min(m0, m1), 2) * f)
        R = Fr(abs(r0) + abs(r1), 2) * max(1, Fr(max(m0, m1), min(m0, m1)))
        if all(abs(o) <= R for o in outs): dep_ok += 1
    elif par(e) == 0 or k != 0:
        # off departures: the two coins give the two multipliers m0/2 and m1/2 (one each), additive error bounded
        mults = []
        for beta in (0, 1):
            Mn, en, i, j = step(m0, m1, r0, r1, M, e, beta)
            fn = en / max(Mn, 1)
            best = min((abs(fn - Fr(mm, 2) * f), mm) for mm in (m0, m1))
            mults.append(best[1])
        R = Fr(abs(r0) + abs(r1), 2) * max(1, Fr(max(m0, m1), min(m0, m1)))
        okm = sorted(mults) == sorted([m0, m1])
        errs = []
        for beta, mm in zip((0, 1), mults):
            Mn, en, i, j = step(m0, m1, r0, r1, M, e, beta)
            errs.append(abs(en / max(Mn, 1) - Fr(mm, 2) * f) <= R)
        if not (okm and all(errs)) and abs(f) > 100 * R: bad += 1
print(f"states checked {checks // 2}: off-departure states with |f| > 100 R violating 'two coins = two multipliers, bounded error': {bad}")
print(f"departure states: {dep_n}, all with both coins giving min(m0,m1)/2 up to bounded error: {dep_ok == dep_n} ({dep_ok}/{dep_n})")
