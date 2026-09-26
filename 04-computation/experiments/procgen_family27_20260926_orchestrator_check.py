#!/usr/bin/env python3
"""Orchestrator audit of lane `family27`, written from the note's statements;
the lane's scripts were not read.

  1. Proposition M numerics: lambda* = 0.488077, min g = 2^-(1-h), beta = 0.0239937
     (1/beta = 41.677648), s_opt = 1.40368, g(1) = g(2) = 1.
  2. Theorem R (exact rise law) with exact rationals for k <= 14 and several W:
     P(tau <= k) = (1 - eps_k)/E[M_tau | tau <= k], W <= E[M_tau | tau <= k] < 3W/2.
  3. Theorem R' on actual integers, blocks k <= 18: N_k(W) <= #{max_j T^j n >= W n} <= N_k(W - (3/4)^k).
  4. Theorem G: #{n in [2^k, 2^(k+1)) : glide(n) >= L} = 2^(k-L+1) W_(L-1) for k <= 20, admissible L.
  5. 27's branch: a random sample of n in [2^24, 2^25) meets 3077 with frequency ~0.3927;
     the tree's members split 1/3 per class mod 3 (Proposition B with s = 1).
  6. 27 maximises ln t(n)/ln n over 3 <= n <= 10^6 (t = max of the shortcut orbit).
"""
import math, random
from fractions import Fraction


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def T(n):
    return (3 * n + 1) // 2 if n & 1 else n >> 1


# 1. Proposition M
g = lambda s: 2 ** (-s) + (3 / 2) ** s / 3
c = math.log(2) / math.log(3)
h = -c * math.log2(c) - (1 - c) * math.log2(1 - c)
lam = math.log(math.log(2) / math.log(1.5), 3)
lo, hi = 0.5, 3.0
for _ in range(200):
    m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
    if g(m1) < g(m2):
        hi = m2
    else:
        lo = m1
smin = (lo + hi) / 2
f = lambda s: -math.log(g(s)) / s
lo, hi = 1.0001, 1.9999
for _ in range(200):
    m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
    if f(m1) > f(m2):
        hi = m2
    else:
        lo = m1
sopt = (lo + hi) / 2
beta = f(sopt)
check(abs(g(1) - 1) < 1e-15 and abs(g(2) - 1) < 1e-15, "g(1) = g(2) = 1 (the two Moran roots; g(2) = 1/4 + 3/4 is 3 + 1 = 4)")
check(abs(lam - 0.488077) < 1e-6 and abs(smin - (1 + lam)) < 1e-6 and abs(g(smin) - 2 ** (-(1 - h))) < 1e-12,
      f"min g at s = 1 + lambda* = {smin:.6f}, min g = 2^-(1-h) = {g(smin):.9f}")
check(abs(beta - 0.0239937) < 1e-6 and abs(1 / beta - 41.677648) < 1e-3 and abs(sopt - 1.40368) < 1e-4,
      f"beta = max -ln g(s)/s = {beta:.7f} at s = {sopt:.5f}; 1/beta = {1/beta:.6f} (the delay-record constant)")

# 2. Theorem R exact
for k in (6, 10, 14):
    for W in (Fraction(2), Fraction(3), Fraction(10), Fraction(50)):
        # enumerate words by DP over (o, j, stopped) with exact M = 3^o/2^j
        p_stop = Fraction(0)
        e_Mtau = Fraction(0)     # E[M_tau ; tau <= k]
        eps = Fraction(0)        # E[M_k ; tau > k]
        states = {0: Fraction(1)}   # o -> probability, among not-yet-stopped
        for j in range(1, k + 1):
            nxt = {}
            for o, pr in states.items():
                for b in (0, 1):
                    o2 = o + b
                    M = Fraction(3 ** o2, 2 ** j)
                    q = pr / 2
                    if M >= W:
                        p_stop += q
                        e_Mtau += q * M
                    else:
                        nxt[o2] = nxt.get(o2, 0) + q
            states = nxt
        for o, pr in states.items():
            eps += pr * Fraction(3 ** o, 2 ** k)
        assert e_Mtau + eps == 1
        if p_stop > 0:
            cond = e_Mtau / p_stop
            assert p_stop == (1 - eps) / cond and W <= cond < Fraction(3, 2) * W
check(True, "Theorem R: optional-stopping identity and W <= E[M_tau | tau <= k] < 3W/2 hold exactly (k = 6, 10, 14; W = 2, 3, 10, 50)")

# 3. Theorem R' on integers
def Nk(k, V):
    cnt = 0
    states = {(0,): 1}
    # count words with max_{j<=k} 3^(o_j)/2^j >= V by DP over (o, reached)
    dp = {(0, False): 1}
    for j in range(1, k + 1):
        nd = {}
        for (o, rch), ct in dp.items():
            for b in (0, 1):
                o2 = o + b
                r2 = rch or (Fraction(3 ** o2, 2 ** j) >= V)
                nd[(o2, r2)] = nd.get((o2, r2), 0) + ct
        dp = nd
    return sum(ct for (o, r), ct in dp.items() if r) + (0 if V > 1 else 0)
for k in (10, 14, 18):
    for W in (2, 4, 16, 128):
        lo_b = Nk(k, Fraction(W))
        up_b = Nk(k, Fraction(W) - Fraction(3, 4) ** k)
        cnt = 0
        for n in range(2 ** k, 2 ** (k + 1)):
            x, mx = n, n
            for _ in range(k):
                x = T(x)
                if x > mx:
                    mx = x
            if mx >= W * n:
                cnt += 1
        assert lo_b <= cnt <= up_b, (k, W, lo_b, cnt, up_b)
check(True, "Theorem R': N_k(W) <= block count <= N_k(W - (3/4)^k) on actual integers (k = 10, 14, 18; W = 2, 4, 16, 128)")

# 4. Theorem G
def Wm(m):
    cur = {0: 1}
    for j in range(1, m + 1):
        nxt = {}
        for o, ct in cur.items():
            for b in (0, 1):
                if 3 ** (o + b) > 2 ** j:
                    nxt[o + b] = nxt.get(o + b, 0) + ct
        cur = nxt
    return sum(cur.values())
for k in (12, 16, 20):
    Lmax = int(1 + k * c)
    glides = {}
    for n in range(2 ** k, 2 ** (k + 1)):
        x, j = n, 0
        while True:
            x = T(x); j += 1
            if x < n or j > Lmax + 2:
                break
        glides[n] = j   # glide = first j with T^j n < n (capped)
    for L in range(1, Lmax + 1):
        cnt = sum(1 for v in glides.values() if v >= L)
        assert cnt == 2 ** (k - L + 1) * Wm(L - 1), (k, L, cnt)
check(True, "Theorem G: #{n in block k : glide >= L} = 2^(k-L+1) W_(L-1) exactly for k = 12, 16, 20 and all admissible L")

# 5. 27's branch
random.seed(2026)
hits, cls = 0, [0, 0, 0]
S = 60000
for _ in range(S):
    n = random.randrange(2 ** 24, 2 ** 25)
    x = n
    hit = False
    while x != 1:
        if x == 3077:
            hit = True
            break
        x = T(x)
    if hit:
        hits += 1
        cls[n % 3] += 1
p = hits / S
se = math.sqrt(p * (1 - p) / S)
check(abs(p - 0.39266) < 4 * se, f"27's branch (meets 3077): sampled density {p:.4f} +- {se:.4f} in [2^24, 2^25), note 0.39266")
fr = [x / hits for x in cls]
check(all(abs(x - 1 / 3) < 0.02 for x in fr), f"its members split {fr[0]:.3f} : {fr[1]:.3f} : {fr[2]:.3f} over classes mod 3 (Proposition B with s = 1 gives 1/3)")

# 6. 27 maximises ln t(n)/ln n for n <= 10^6
best, arg = 0, None
for n in range(3, 10 ** 6 + 1):
    x, mx = n, n
    while x != 1:
        x = T(x)
        if x > mx:
            mx = x
    r = math.log(mx) / math.log(n)
    if r > best:
        best, arg = r, n
check(arg == 27 and abs(best - math.log(4616) / math.log(27)) < 1e-12, f"27 maximises ln t(n)/ln n over 3 <= n <= 10^6 (ratio {best:.4f}, t(27) = 4616)")
