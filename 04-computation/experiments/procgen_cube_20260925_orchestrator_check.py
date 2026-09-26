#!/usr/bin/env python3
"""Orchestrator's independent check of THM-4474 (strategy cube), written without reading the lane's code.
  S1  class (i) counts at levels 2, 3, 4 via the Theorem A criterion (max cycle odd-density < log_3 2, Karp's
      algorithm on the parity graph): 1, 1, 16 (totals, lifts included).
  S2  every class-(i) strategy at levels 2-4 sends every n < 3000 to {1,2}.
  S3  Proposition F: every class-(i) strategy has sigma(1) = + and sigma(-1) = -.
  S4  Theorem E: max density is invariant under nu (sigma'(m) = -sigma(-m mod 2^k)) at levels 2-4.
  S5  Collatz and 3n-1 have max cycle density 1 at levels 2-6 (the all-odd loops at -1, resp. 1).
"""
import math, itertools
c = math.log(2) / math.log(3)
def maxdensity(k, sig):
    M = 2**k; n = M; NEG = -10**9
    adj = []
    for s in range(M):
        t = (s // 2) % (M // 2) if s % 2 == 0 else ((3*s + sig[s // 2]) // 2) % (M // 2)
        adj.append((t, t + M // 2))
    w = [s % 2 for s in range(M)]
    D = [[0]*n] + [[NEG]*n for _ in range(n)]
    for j in range(1, n + 1):
        Dj, Dp = D[j], D[j-1]
        for u in range(n):
            if Dp[u] > NEG:
                val = Dp[u] + w[u]
                for v in adj[u]:
                    if val > Dj[v]: Dj[v] = val
    best = -1e9
    for v in range(n):
        if D[n][v] <= NEG: continue
        worst = min((D[n][v] - D[j][v]) / (n - j) for j in range(n) if D[j][v] > NEG)
        best = max(best, worst)
    return best
def converges(k, sig, N=3000, steps=5000):
    M = 2**k
    T = lambda n: n // 2 if n % 2 == 0 else (3*n + sig[(n % M) // 2]) // 2
    for n0 in range(1, N):
        n = n0; seen = set()
        for _ in range(steps):
            if n in (1, 2): break
            if n in seen: return False
            seen.add(n); n = T(n)
        else:
            return False
    return True
def nu(k, sig):
    M = 2**k
    return [-sig[((-(2*i + 1)) % M) // 2] for i in range(len(sig))]
expect = {2: 1, 3: 1, 4: 16}
for k in (2, 3, 4):
    cls = []
    for bits in itertools.product((1, -1), repeat=2**(k-1)):
        sig = list(bits); rho = maxdensity(k, sig)
        assert abs(rho - maxdensity(k, nu(k, sig))) < 1e-12          # S4
        if rho < c - 1e-12: cls.append(sig)
    assert len(cls) == expect[k], (k, len(cls))                      # S1
    assert all(converges(k, s) for s in cls)                         # S2
    assert all(s[0] == 1 and s[-1] == -1 for s in cls)               # S3: sigma(1)=+, sigma(2^k-1)=-
    print(f"S1-S4 level {k}: class (i) = {len(cls)} strategies; all converge (n<3000); all have sigma(1)=+, sigma(-1)=-; nu-invariant densities: ok")
for k in range(2, 7):
    assert maxdensity(k, [1]*(2**(k-1))) == 1.0 and maxdensity(k, [-1]*(2**(k-1))) == 1.0
print("S5  Collatz and 3n-1 have max cycle density 1 at levels 2-6: ok")
print("ALL CHECKS PASSED")
