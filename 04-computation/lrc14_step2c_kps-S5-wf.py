"""Heavy step: F(k) ceiling + cap-14 exhaustive minima k=10..13. Standalone for speed."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb
import sys, time

def mu_exact(E):
    E = sorted(set(int(e) for e in E))
    k = len(E)
    if k <= 1:
        return F(1) if k == 1 else F(0)
    TH = F(2, 7)
    bp = set([F(0), F(1)])
    diffs = set()
    for i in range(k):
        for j in range(i + 1, k):
            diffs.add(E[j] - E[i])
    for d in diffs:
        for m in range(0, d + 1):
            bp.add(F(m, d))
    obp = sorted(b for b in bp if F(0) <= b <= F(1))
    refined = set(obp)
    for a, b in zip(obp, obp[1:]):
        mid = (a + b) / 2
        floors = {e: (e * mid).__floor__() for e in E}
        order = sorted(E, key=lambda e: e * mid - floors[e])
        for t in range(k):
            if t == k - 1:
                ef, el = order[0], order[-1]
                slope = ef - el
                const = F(1) - floors[ef] + floors[el]
            else:
                eh, elo = order[t + 1], order[t]
                slope = eh - elo
                const = -floors[eh] + floors[elo]
            if slope != 0:
                xb = (TH - const) / slope
                if a < xb < b:
                    refined.add(xb)
    refined = sorted(refined)
    tot = F(0)
    for a, b in zip(refined, refined[1:]):
        mid = (a + b) / 2
        pts = sorted(set((e * mid) % 1 for e in E))
        if len(pts) == 1:
            mg = F(1)
        else:
            gaps = [pts[t + 1] - pts[t] for t in range(len(pts) - 1)]
            gaps.append(pts[0] + 1 - pts[-1])
            mg = max(gaps)
        if mg > TH:
            tot += (b - a)
    return tot

def F_ceiling(k):
    s = F(0); j = 1
    while True:
        base = F(1) - F(2*j, 7)
        if base <= 0:
            break
        s += (-1)**(j+1) * comb(k, j) * base**(k-1)
        j += 1
    return s

def exhaustive_min(k, cap):
    best = None; bestE = None
    for combo in combinations(range(1, cap+1), k-1):
        E = (0,) + combo
        g = 0
        for e in E:
            g = gcd(g, e)
        if g > 1:
            continue
        m = mu_exact(E)
        if best is None or m < best:
            best = m; bestE = E
    return best, bestE

if __name__ == "__main__":
    print("F(k) ceiling:", flush=True)
    claimedF = {12:F(574246018,1977326743),13:F(3132376013,13841287201)}
    for k in [12,13]:
        f = F_ceiling(k)
        print(f"  F({k})={f}={float(f):.6f} claim {claimedF[k]} {'OK' if f==claimedF[k] else 'MISMATCH'}", flush=True)
    claimed = {10:F(468,2695),11:F(409,2548),12:F(5367,35035),13:F(5367,35035)}
    for k in [10,11,12,13]:
        t0 = time.time()
        best, bestE = exhaustive_min(k, 14)
        c = claimed[k]
        print(f"  k={k} cap14 min={best}={float(best):.6f} at {bestE} claim {c} {'OK' if best==c else 'MISMATCH'} [{time.time()-t0:.0f}s]", flush=True)
