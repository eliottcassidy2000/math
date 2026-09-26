#!/usr/bin/env python3
"""Orchestrator audit of lane `floor` (min-max cycle density of sign strategies),
written from the note's statements; the lane's scripts were not read.

  1. rho*(q,k) = min over sign strategies of rho_max, by brute force with exact
     Karp, for q = 5, 7 and k = 2..4 (all strategies): equals 1/2.
  2. Theorem N: rho_max(sigma) >= log_(q+1) 2 for random strategies
     (q = 5, 7, 9, 11, 13; k = 3..7).
  3. Corollary 5 (q = 5): the potential Phi(u) = u^2 (Phi(3) = 8) satisfies
     Phi(v) <= 8 Phi(u) on odd steps and 4 Phi(v) <= Phi(u) on even steps of the
     negative-integer graph U_k (u -> u/2; u -> rho((5u +- 1)/2)), for k <= 16.
  4. Proposition F: the level-2 max-halving strategy has stationary law
     (1/3, 1/6, 1/3, 1/6) (odd frequency 1/3) and rho_max = 1/2, for q = 3..23.
"""
import math, random, itertools
from fractions import Fraction


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def targets(k, q, sig):
    N, H = 1 << k, 1 << (k - 1)
    return [((s // 2) % H) if s % 2 == 0 else (((q * s + sig[s // 2]) // 2) % H) for s in range(N)], H


def karp(k, q, sig):
    N = 1 << k
    t, H = targets(k, q, sig)
    NEG = -10 ** 9
    D = [[0] * N] + [[NEG] * N for _ in range(N)]
    for j in range(1, N + 1):
        Dj, Dp = D[j], D[j - 1]
        for u in range(N):
            if Dp[u] > NEG:
                val = Dp[u] + (u & 1)
                for v in (t[u], t[u] + H):
                    if val > Dj[v]:
                        Dj[v] = val
    best = None
    for v in range(N):
        if D[N][v] <= NEG:
            continue
        worst = min(Fraction(D[N][v] - D[j][v], N - j) for j in range(N) if D[j][v] > NEG)
        if best is None or worst > best:
            best = worst
    return best


# 1. brute force rho* for small k
for q in (5, 7):
    for k in (2, 3, 4):
        H = 1 << (k - 1)
        vals = [karp(k, q, list(sig)) for sig in itertools.product((1, -1), repeat=H)]
        assert min(vals) == Fraction(1, 2), (q, k, min(vals))
check(True, "rho*(q,k) = 1/2 for q = 5, 7 and k = 2, 3, 4 (exhaustive, exact Karp)")

# 2. Theorem N on random strategies
random.seed(26)
cnt = 0
for q in (5, 7, 9, 11, 13):
    bound = math.log(2) / math.log(q + 1)
    for k in range(3, 8):
        H = 1 << (k - 1)
        for _ in range(12 if k <= 6 else 4):
            sig = [random.choice((1, -1)) for _ in range(H)]
            r = karp(k, q, sig)
            assert float(r) >= bound - 1e-12, (q, k, r, bound)
            cnt += 1
check(True, f"Theorem N: rho_max >= log_(q+1) 2 on {cnt} random strategies (q = 5..13, k = 3..7)")

# 3. Corollary 5 potential inequalities on U_k (q = 5)
def Phi(u):
    return 8 if u == 3 else u * u
for k in range(2, 17):
    H = 1 << (k - 1)
    for u in range(1, H + 1):
        if u % 2 == 0:
            v = u // 2
            assert 4 * Phi(v) <= Phi(u), (k, u)
        else:
            for sgn in (1, -1):
                m = (5 * u + sgn) // 2
                v = ((m - 1) % H) + 1          # representative in [1, H]
                assert Phi(v) <= 8 * Phi(u), (k, u, sgn, v)
check(True, "Corollary 5: Phi(v) <= 8 Phi(u) (odd) and 4 Phi(v) <= Phi(u) (even) on U_k for all k <= 16, so every cycle has 5a >= 2p")
# the sporadic 5x+1 cycle attains 2/5
cyc, y = [], 1
while True:
    cyc.append(y); y = (5 * y + 1) // 2 if y % 2 else y // 2
    if y == 1:
        break
check(cyc == [1, 3, 8, 4, 2] and Fraction(sum(c % 2 for c in cyc), len(cyc)) == Fraction(2, 5),
      "the 5x+1 cycle (1,3,8,4,2) has density exactly 2/5")

# 4. Proposition F
for q in range(3, 24, 2):
    sig = []
    for r in (1, 3):
        sig.append(1 if (q * r + 1) % 4 == 0 else -1)   # flip r iff 4 does not divide q r + 1
    t, H = targets(2, q, sig)
    # uniform-lift chain on Z/4; stationary law by exact linear algebra (power iteration on fractions)
    P = [[Fraction(0)] * 4 for _ in range(4)]
    for s in range(4):
        for tg in (t[s], t[s] + H):
            P[s][tg] += Fraction(1, 2)
    pi = [Fraction(1, 3), Fraction(1, 6), Fraction(1, 3), Fraction(1, 6)]
    nxt = [sum(pi[s] * P[s][j] for s in range(4)) for j in range(4)]
    assert nxt == pi, (q, nxt)
    assert pi[1] + pi[3] == Fraction(1, 3) and karp(2, q, sig) == Fraction(1, 2)
check(True, "Prop F: the max-halving strategy has stationary law (1/3,1/6,1/3,1/6), odd frequency 1/3, rho_max 1/2 (q = 3..23)")
