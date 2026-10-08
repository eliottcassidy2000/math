#!/usr/bin/env python3
"""(1) Ladder completeness: every equal-time merge of a residual source with its child h_3 (universal state after the
two-run) is a ladder collision: the last distinct odd values z_u, z_v satisfy z_v = 4^i z_u + (4^i-1)/3 (child ladder)
or z_u = 4^i z_v + (4^i-1)/3 (source ladder), i >= 1, final letters differing by 2i, merge at post-run Terras depth sum(u)+1,
and |v| = |u| + 3, F_v(1) = 4^i F_u(1) + (4^i-1)/3 (or symmetric).
(2) Bi-anchored drift: if u - c = 3^k (v - c') with c, c' cycle points of words w, w', the runs proceed in lockstep and the
debt drifts by |w|/sum(w) - |w'|/sum(w') per Terras step (checked on letters a, a' in 1..4 and on words (1,2), (1,3))."""
import random
from fractions import Fraction as Fr
def v2(x): return (x & -x).bit_length() - 1
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
def U(x):
    y = 3*x + 1; return y >> v2(y)
def F(word, x):
    for a in word: x = (3*x + 1) / Fr(2**a)
    return x
rnd = random.Random(2026)
CH = 0
# ---- (1) ladder completeness on actual sources ----
stats = {}
for trial in range(3000):
    K = rnd.randint(5, 40); J = rnd.randint(3, 12); r = rnd.choice((1, 2)); j = 2*J + r
    M = 1 << (j - 1)
    t0 = pow(3, -(K-1), M)
    while True:
        t = t0 + M * rnd.getrandbits(200)
        if t & 1 and v2(3**(K-1)*t - 1) == j - 1: break
    x = 2*3**(K-1)*t - 1; y = (x + 1)//27 - 1
    # advance both by the two-run (2J Terras steps): universal state
    u, v, k = x, y, 3
    for s in range(2*J):
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
    assert k == 3 and u - 1 == 27*(v - 1)
    X, Y = u, v
    # follow the chain up to 3000 steps; record the merge
    hist_u = [u]; hist_v = [v]; merged = None
    for s in range(1, 3001):
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
        hist_u.append(u); hist_v.append(v)
        if u == v and k == 0: merged = s; break
    if merged is None: continue
    # U-words of X and Y up to the merge
    def oddlist(h):
        return [h[i] for i in range(len(h)) if h[i] & 1]
    ou, ov = oddlist(hist_u[:merged]), oddlist(hist_v[:merged])
    zu, zv = ou[-1], ov[-1]
    assert zu != zv and U(zu) == U(zv), "last odd values are distinct U-preimages"
    if zv > zu:
        i = 0; z = zu
        while z < zv: z = 4*z + 1; i += 1
        assert z == zv; kind = 'child'
    else:
        i = 0; z = zv
        while z < zu: z = 4*z + 1; i += 1
        assert z == zu; kind = 'source'
    # words (heads without final letter)
    wu = [v2(3*a + 1) for a in ou[:-1]]; wv = [v2(3*a + 1) for a in ov[:-1]]
    assert len(wv) == len(wu) + 3
    cu, cv = v2(3*zu + 1), v2(3*zv + 1)
    if kind == 'child':
        assert cv - cu == 2*i and sum(wu) == sum(wv) + 2*i
        assert F(wv, Fr(1)) == 4**i * F(wu, Fr(1)) + Fr(4**i - 1, 3)
    else:
        assert cu - cv == 2*i and sum(wv) == sum(wu) + 2*i
        assert F(wu, Fr(1)) == 4**i * F(wv, Fr(1)) + Fr(4**i - 1, 3)
    # merge time = sum of source head + 1 for child ladders (source reaches (3zu+1)/2 first)
    if kind == 'child': assert merged == sum(wu) + 1
    stats[(kind, i)] = stats.get((kind, i), 0) + 1
    CH += 1
print("(1) ladder completeness: merges observed by (kind, index):", dict(sorted(stats.items())), f"; {CH} merges, all ladders")

# ---- (2) bi-anchored drift ----
def cycle_point(w):
    m = len(w); S = sum(w); B = 0; Q = 1
    for a in w:
        B = 3*B + Q; Q *= 2**a
    return Fr(B, 2**S - 3**m)
def two_adic(q, bits):  # q rational with odd denominator -> integer mod 2^bits
    M = 1 << bits
    return (q.numerator * pow(q.denominator, -1, M)) % M
tests = 0
for w in [(1,), (2,), (3,), (4,), (1, 2), (1, 3)]:
    for wp in [(1,), (2,), (3,), (4,), (1, 2), (1, 3)]:
        c, cp = cycle_point(w), cycle_point(wp)
        P, Pp = sum(w), sum(wp)
        from math import lcm
        Lam = lcm(P, Pp)
        rho = Fr(len(w), P) - Fr(len(wp), Pp)
        for _ in range(20):
            k = rnd.randint(0, 6); nu = 6*Lam + rnd.randint(1, 5)
            bits = nu + 60
            Mb = 1 << bits
            cpm = two_adic(cp, bits); cm = two_adic(c, bits)
            s_ = rnd.getrandbits(40) | 1
            v = (cpm + (1 << nu) * s_) % Mb
            u = (cm + pow(3, k) * (v - cpm)) % Mb
            # lift to positive integers with the right classes
            v += Mb * rnd.getrandbits(10); u += Mb * rnd.getrandbits(10)
            kk = k; uu, vv = u, v
            for step in range(Lam * 5):
                kk += (uu & 1) - (vv & 1); uu, vv = T(uu), T(vv)
            assert kk - k == rho * Lam * 5, (w, wp, kk - k, rho * Lam * 5)
            # relation persists modulo remaining precision: u' - c = 3^kk (v' - c')  (cycle points rotate back after full periods)
            rem = bits - Lam * 5 - 2
            assert (uu - two_adic(c, rem) - pow(3, kk) * (vv - two_adic(cp, rem))) % (1 << rem) == 0 if kk >= 0 else True
            tests += 1
print(f"(2) bi-anchored drift: {tests} exact tests over 36 (word, word') pairs; debt change = (|w|/sum w - |w'|/sum w') * time")
print("ALL CHECKS PASSED")
