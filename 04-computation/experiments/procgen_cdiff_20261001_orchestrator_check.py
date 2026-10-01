#!/usr/bin/env python3
"""procgen_cdiff_20261001_orchestrator_check.py -- the orchestrator's independent audit of the cdiff lane (2026-10-01).

Written from the statements in the lane note; the lane's code was not read.

S(A) = oddpart(3A+1) (sign kept for negative A), K(A) = (A - S(A))/2, M_v = 2^v - 3, N(n) = #{v >= 3 : M_v | n}.
Checks:
  A. THM-4527 (S15) and Prop 2.6: over all odd A in Z, every d in Z is a drop exactly 2 + N(6d+1) times (|d| <= 1500);
     over positive A: d < 0 once, d >= 0 exactly 1 + N(6d+1) times.
  B. Densities delta_j = P(N = j-1) computed exactly (independent route: prime factorisation of M_v for v <= 64, connected
     components of the shared-prime graph, inclusion-exclusion inside components, convolution across components), compared
     with the lane's 18-digit values; and a sieve over [0, 2*10^7].
  C. Records: m(d_k) = k+1 for the lane's record values d_k (k <= 12, exact big integers); minimality for k <= 5 by sieve;
     the proved bounds sqrt(2 log2 X) - 2.5 <= max m <= 1 + log2(6X+1)/2 on the sieve range.
  D. Summatory identities: sum_{N <= 4^n} K_N = (3*16^n - 4^n - 2)/6 (n <= 9); sum_{A odd < 2^n} K(A) = -J_{n-1} (n <= 22).
  E. k-step identity (2^p - 3^k) Syr^k(A) = 2*3^k*D + c(v); two copies at k = 2 (D <= -1: drops of -16D-5, -16D-7) and the
     law mu_2 = 2 + [5 | D] on D < 0; P(mu_2 = 0) and P(mu_3 = 0) on finite ranges against 0.4066556 / 0.1000254.
  F. Theorem 4.5: A -> (K(A), K(S(A))) injective on odd A < 10^6, for 3x+1 and for 3x-1.
  G. Theorem 5.2: the qx+1 drop map hits every d >= 0 iff q = 2^a - 1, every d < 0 iff q = 2^a + 1 (odd q <= 33, |d| <= 400);
     5x+1 misses about 59.27% of d >= 0.
  H. Mod 3^k: S = (6K+1)(2^v-3)^{-1} mod 3^k and Syr^k(A) = 2^{-p} c(v) mod 3^k (k <= 5, random A).
"""
import math
import random
import sys
import time
from fractions import Fraction
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np
import sympy

OKS = []


def ok(cond, msg):
    OKS.append(bool(cond))
    print(('[OK] ' if cond else '[FAIL] ') + msg, flush=True)


def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n


def v2(n):
    n = abs(n)
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


def S(A, s=1):
    return oddpart(3 * A + s)


def K(A, s=1):
    return (A - S(A, s)) // 2


def Nfun(n, vmax=None):
    n = abs(n)
    if vmax is None:
        vmax = n.bit_length() + 3
    return sum(1 for v in range(3, vmax + 1) if n % (2 ** v - 3) == 0)


def check_A():
    D = 1500
    B = 8 * D + 50
    cnt = Counter()
    pos = Counter()
    for A in range(-B, B + 1):
        if A % 2 == 0:
            continue
        k = K(A)
        cnt[k] += 1
        if A > 0:
            pos[k] += 1
    bad = [d for d in range(-D, D + 1) if cnt[d] != 2 + Nfun(6 * d + 1)]
    ok(not bad, f'A: over odd A in Z (|A| <= {B}), every |d| <= {D} is a drop exactly 2 + N(6d+1) times (bad: {bad[:5]})')
    badp = [d for d in range(-D, D + 1) if pos[d] != (1 if d < 0 else 1 + Nfun(6 * d + 1))]
    ok(not badp, f'A: over odd A > 0: d < 0 once, d >= 0 exactly 1 + N(6d+1) times (|d| <= {D})')


def check_B():
    """exact densities by conditioning on the valuations at the shared primes (independent of the lane's route)"""
    import mpmath
    mpmath.mp.dps = 50
    V = 64
    Ms = {v: 2 ** v - 3 for v in range(3, V + 1)}
    fac = {v: sympy.factorint(Ms[v]) for v in Ms}
    owner = defaultdict(list)
    for v, f in fac.items():
        for q in f:
            owner[q].append(v)
    shared = sorted(q for q, vs in owner.items() if len(vs) > 1)
    req = {v: {q: e for q, e in fac[v].items() if q in shared} for v in Ms}
    priv = {v: Ms[v] // math.prod(q ** e for q, e in req[v].items()) for v in Ms}
    Emax = {q: max(req[v].get(q, 0) for v in Ms) for q in shared}
    # joint law of capped valuations a_q in {0..Emax_q} (state Emax_q means >= Emax_q)
    states = [dict()]
    probs = [mpmath.mpf(1)]
    for q in shared:
        ns, npb = [], []
        for st, pb in zip(states, probs):
            for a in range(Emax[q] + 1):
                pa = (mpmath.mpf(q) ** (-a)) * ((1 - mpmath.mpf(1) / q) if a < Emax[q] else 1)
                d = dict(st); d[q] = a
                ns.append(d); npb.append(pb * pa)
        states, probs = ns, npb
    dist = [mpmath.mpf(0)] * (V + 2)
    for st, pb in zip(states, probs):
        poly = [mpmath.mpf(1)]
        for v in Ms:
            if all(st[q] >= e for q, e in req[v].items()):
                qv = mpmath.mpf(1) / priv[v]
                newp = [mpmath.mpf(0)] * (len(poly) + 1)
                for i, a in enumerate(poly):
                    newp[i] += a * (1 - qv)
                    newp[i + 1] += a * qv
                poly = newp
        for i, a in enumerate(poly):
            dist[i] += pb * a
    lane = {1: '0.6961749034401325952912', 2: '0.266238036677738373892', 3: '0.03538955883931588242586'}
    good = True
    for j in (1, 2, 3):
        diff = abs(dist[j - 1] - mpmath.mpf(lane[j]))
        if diff > mpmath.mpf('1e-17'):
            good = False
        print(f'    delta_{j} (orchestrator, v <= {V}) = {mpmath.nstr(dist[j - 1], 20)}   lane {lane[j]}   |diff| = {mpmath.nstr(diff, 3)}')
    ok(good and len(shared) == 10,
       f'B: exact densities delta_1..3 agree with the lane to < 1e-17 (v <= 64, tail < 2^-63); {len(shared)} shared primes {shared}; '
       f'{len(states)} valuation states')
    total = [float(x) for x in dist[:4]]
    # sieve check
    X = 2 * 10 ** 7
    m = np.ones(X + 1, dtype=np.int8)
    v = 3
    while 2 ** v - 3 <= 6 * X + 1:
        M = 2 ** v - 3
        r = (-pow(6, -1, M)) % M  # 6d + 1 = 0 mod M
        m[r::M] += 1
        v += 1
    freq = np.bincount(m, minlength=6) / (X + 1)
    sv = [float(freq[j]) for j in (1, 2, 3)]
    ok(all(abs(sv[j - 1] - total[j - 1]) < 2e-5 for j in (1, 2, 3)),
       f'B: sieve on [0, 2e7]: delta_1..3 = {[round(x, 7) for x in sv]} (exact {[round(total[j], 7) for j in range(3)]})')
    return m


def check_C(m):
    recs = [2, 24, 314, 7854, 479104, 121213354, 49576261854, 50617363353104, 115043252496711229,
            117459160799142164979, 480760345150888881259729, 2423512899905630850430294729]
    good = all(1 + Nfun(6 * d + 1) == k + 2 for k, d in enumerate(recs))
    ok(good, 'C: m(d_k) = k+1 exactly for the lane\'s 12 record values (exact big-integer divisibility)')
    firsts = {}
    for k in range(2, 7):
        idx = np.nonzero(m >= k)[0]
        firsts[k] = int(idx[0]) if len(idx) else None
    ok([firsts[k] for k in range(2, 7)] == [2, 24, 314, 7854, 479104],
       f'C: first d with m(d) >= 2..6 on the sieve range: {[firsts[k] for k in range(2, 7)]}')
    X = len(m) - 1
    mx = int(m.max())
    ok(math.sqrt(2 * math.log2(X)) - 2.5 <= mx <= 1 + math.log2(6 * X + 1) / 2,
       f'C: max m(d) on [0, {X}] = {mx} lies in [{math.sqrt(2*math.log2(X))-2.5:.2f}, {1+math.log2(6*X+1)/2:.2f}]')


def check_D():
    good = True
    for n in range(1, 10):
        s = sum(K(4 * N - 3) for N in range(1, 4 ** n + 1))
        if s != (3 * 16 ** n - 4 ** n - 2) // 6 or (3 * 16 ** n - 4 ** n - 2) % 6:
            good = False
    ok(good, 'D: sum_{N <= 4^n} K_N = (3*16^n - 4^n - 2)/6 for n <= 9')
    J = [0, 1]
    for i in range(2, 30):
        J.append(J[-1] + 2 * J[-2])
    good = True
    acc = 0
    A = 1
    for n in range(1, 23):
        while A < 2 ** n:
            acc += K(A)
            A += 2
        if acc != -J[n - 1]:
            good = False
    ok(good, 'D: sum_{A odd < 2^n} K(A) = -J_{n-1} (Jacobsthal) for n <= 22')


def syr_k(A, k, s=1):
    vs = []
    x = A
    for _ in range(k):
        y = 3 * x + s
        v = v2(y)
        vs.append(v)
        x = y // 2 ** v
    return x, vs


def c_of(vs, k):
    # Syr^k(A) = (3^k A + c)/2^p with c = sum_i 3^{k-1-i} 2^{v_1+...+v_i}
    c = 0
    acc = 0
    for i in range(k):
        c += 3 ** (k - 1 - i) * 2 ** acc
        acc += vs[i]
    return c


def check_E():
    rng = random.Random(1)
    good = True
    for _ in range(3000):
        A = 2 * rng.randint(0, 10 ** 9) + 1
        k = rng.randint(1, 6)
        y, vs = syr_k(A, k)
        p = sum(vs)
        D = (A - y)
        assert D % 2 == 0
        D //= 2
        if (2 ** p - 3 ** k) * y != 2 * 3 ** k * D + c_of(vs, k):
            good = False
    ok(good, 'E: k-step identity (2^p - 3^k) Syr^k(A) = 2*3^k*D + c(v) on 3000 random (A, k)')
    good = all(syr_k(-16 * D - 5, 2)[0] == (-16 * D - 5) - 2 * D and syr_k(-16 * D - 7, 2)[0] == (-16 * D - 7) - 2 * D
               for D in range(-3000, 0))
    ok(good, 'E: every D in [-3000, -1] is the 2-step drop of -16D-5 and of -16D-7')
    # mu_2 law on D < 0: 2 + [5 | D]
    Dmax = 2000
    cnt = Counter()
    for A in range(1, 16 * Dmax + 200, 2):
        y, _ = syr_k(A, 2)
        D = (A - y) // 2
        if -Dmax <= D < 0:
            cnt[D] += 1
    ok(all(cnt[D] == 2 + (1 if D % 5 == 0 else 0) for D in range(-Dmax, 0)), f'E: mu_2(D) = 2 + [5 | D] for D in [-{Dmax}, -1]')
    X = 300000
    hit2 = np.zeros(X + 1, dtype=bool)
    for A in range(1, 5 * X + 20, 2):
        y, _ = syr_k(A, 2)
        D = (A - y) // 2
        if 0 <= D <= X:
            hit2[D] = True
    p0 = 1 - hit2.mean()
    ok(abs(p0 - 0.4066556) < 2e-3, f'E: P(mu_2 = 0) on [0, {X}] = {p0:.6f} (lane: 0.4066556)')
    X3 = 100000
    hit3 = np.zeros(X3 + 1, dtype=bool)
    for A in range(1, 13 * X3 + 40, 2):
        y, _ = syr_k(A, 3)
        D = (A - y) // 2
        if 0 <= D <= X3:
            hit3[D] = True
    p3 = 1 - hit3.mean()
    ok(abs(p3 - 0.1000254) < 3e-3, f'E: P(mu_3 = 0) on [0, {X3}] = {p3:.6f} (lane: 0.1000254)')


def check_F():
    for s in (1, -1):
        seen = {}
        dup = 0
        for A in range(1, 10 ** 6, 2):
            a = S(A, s)
            key = (K(A, s), (a - S(a, s)) // 2)
            if key in seen:
                dup += 1
            seen[key] = A
        ok(dup == 0, f'F: (K(A), K(S(A))) injective on odd A < 10^6 for 3x{"+" if s > 0 else "-"}1')


def check_G():
    res = {}
    D = 400
    for q in range(1, 34, 2):
        B = 2 * (q + 1) * D + 20
        hit = set()
        for A in range(1, B, 2):
            T = oddpart(q * A + 1)
            hit.add((A - T) // 2)
        pos_all = all(d in hit for d in range(0, D + 1))
        neg_all = all(d in hit for d in range(-D, 0))
        pred_pos = ((q + 1) & q) == 0           # q = 2^a - 1
        pred_neg = q > 1 and ((q - 1) & (q - 2)) == 0  # q = 2^a + 1 (q >= 3)
        res[q] = (pos_all == pred_pos, neg_all == pred_neg, pos_all, neg_all)
    good = all(a and b for a, b, _, _ in res.values())
    both = [q for q, r in res.items() if r[2] and r[3]]
    ok(good and both == [3], f'G: hits all d >= 0 iff q = 2^a - 1, all d < 0 iff q = 2^a + 1 (odd q <= 33); all of Z only for q = {both}')
    q = 5
    Dm = 20000
    hit = set()
    for A in range(1, 2 * (q + 1) * Dm + 20, 2):
        T = oddpart(q * A + 1)
        d = (A - T) // 2
        if 0 <= d <= Dm:
            hit.add(d)
    miss = 1 - len(hit) / (Dm + 1)
    ok(abs(miss - 0.5927) < 5e-3, f'G: 5x+1 misses {miss:.4f} of d in [0, {Dm}] (lane: 0.5927)')


def check_H():
    rng = random.Random(9)
    good = True
    for _ in range(4000):
        A = 2 * rng.randint(0, 10 ** 12) + 1
        k = rng.randint(1, 5)
        mod = 3 ** k
        v = v2(3 * A + 1)
        s = S(A)
        kk = K(A)
        if s % mod != ((6 * kk + 1) * pow(2 ** v - 3, -1, mod)) % mod:
            good = False
        y, vs = syr_k(A, k)
        p = sum(vs)
        if y % mod != (pow(2, -p, mod) * c_of(vs, k)) % mod:
            good = False
    ok(good, 'H: S = (6K+1)(2^v-3)^{-1} mod 3^k and Syr^k(A) = 2^{-p} c(v) mod 3^k (4000 random cases, k <= 5)')


def main():
    t0 = time.time()
    print('==== A ====', flush=True); check_A()
    print('==== B ====', flush=True); m = check_B()
    print('==== C ====', flush=True); check_C(m)
    print('==== D ====', flush=True); check_D()
    print('==== E ====', flush=True); check_E()
    print('==== F ====', flush=True); check_F()
    print('==== G ====', flush=True); check_G()
    print('==== H ====', flush=True); check_H()
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
