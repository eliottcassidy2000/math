# Basins of the signed descent D on odd integers, and the density function F.
#   D(h) = (h+1)/2 (h = 1 mod 4), (3h+3)/4 (h = 3 mod 8), (3h-1)/4 (h = 7 mod 8).   n = (h-1)/2.
import numpy as np, math, sys
X = 2 * 10**7
h = np.arange(1, X + 1, 2, dtype=np.int64)            # index i <-> h = 2i+1
def Dv(h):
    return np.where(h % 4 == 1, (h + 1) // 2, np.where(h % 8 == 3, (3 * h + 3) // 4, (3 * h - 1) // 4))
par = (Dv(h) - 1) // 2                                # index of D(h)
M = len(h)
# basin sizes: count[i] = # of h <= X whose orbit passes through 2i+1  (process from the top down)
cnt = np.ones(M, dtype=np.int64)
# D(h) < h for h > 3, so a single descending pass suffices; vectorise by chunks from the top
B = 1 << 16
for hi in range(M, 2, -B):
    lo = max(2, hi - B)
    # inside a chunk, parents are below lo only when h > 2*chunk top ... do it exactly with a loop over sub-blocks
    idx = np.arange(lo, hi)
    # parents of idx are < idx; those inside [lo, hi) need sequential handling: resolve by descending order in pure numpy via repeated passes
    p = par[lo:hi]
    inside = p >= lo
    if inside.any():
        order = np.argsort(-idx)                       # descending
        for k in order:
            cnt[p[k]] += cnt[lo + k]
    else:
        np.add.at(cnt, p, cnt[lo:hi])
print("X =", X)
def rho(t): return cnt[(t - 1) // 2] / M
print("basin densities: 21: %.4f  27: %.4f  41: %.4f  7: %.6f  13: %.6f" % (rho(21), rho(27), rho(41), rho(7), rho(13)))
# F(n) = t * rho(t) for branching t (t = 3, 5 mod 6), t <= 2001
ts = [t for t in range(5, 2001, 2) if t % 6 != 1]
F = {t: t * rho(t) for t in ts}
vals = np.array([F[t] for t in ts])
print("F = t*rho(t) over branching t in [5, 2000]: mean %.3f  min %.3f (t=%d)  max %.3f (t=%d)" % (vals.mean(), vals.min(), ts[int(vals.argmin())], vals.max(), ts[int(vals.argmax())]))
print("heuristic mean  (visits per dyadic window)/(ln 2 * 1/3) = %.3f" % ((1 / (0.5 + 0.5 * math.log2(4 / 3))) * 3 / math.log(2)))
# functional equation F(n) ~ F(2n)/2 + (3/4) F(rise)  in t-coordinates: rho(t) = rho(2t-1) + rho(rise(t)) exactly up to the cut-off
errs = []
for t in range(5, 400, 2):
    if t % 6 == 1: continue
    r = (4 * t - 3) // 3 if t % 6 == 3 else (4 * t + 1) // 3
    errs.append(abs(cnt[(t - 1) // 2] - 1 - cnt[t - 1] - cnt[(r - 1) // 2]))
print("exact additivity  basin(t) = {t} + basin(2t-1) + basin(rise t): max defect", max(errs))
# 3-adic dependence: F by n mod 9 class
for mod in (3, 9, 27):
    by = {}
    for t in ts:
        n = (t - 1) // 2; by.setdefault(n % mod, []).append(F[t])
    print("mean F by n mod %d:" % mod, {k: round(float(np.mean(v)), 2) for k, v in sorted(by.items())})
# landing multiplicity: number of orbit points in the dyadic window [Y, 2Y) for random h
import random
rng = random.Random(1); Y = 1 << 20; dist = {}
for _ in range(200000):
    x = rng.randrange(1 << 59, 1 << 60) | 1; c = 0
    while x >= Y:
        if x < 2 * Y: c += 1
        x = (x + 1) // 2 if x % 4 == 1 else ((3 * x + 3) // 4 if x % 8 == 3 else (3 * x - 1) // 4)
    dist[c] = dist.get(c, 0) + 1
tot = sum(dist.values())
print("orbit points in a dyadic window [Y, 2Y): distribution", {k: round(v / tot, 4) for k, v in sorted(dist.items())}, " mean %.4f (predicted 1/(1/2 + log2(4/3)/2) = %.4f)" % (sum(k * v for k, v in dist.items()) / tot, 1 / (0.5 + 0.5 * math.log2(4 / 3))))
# level sizes of the basin of 21 and of the whole tree above 3
def levels(root, K):
    lev = [root]; out = []
    for k in range(K):
        nxt = []
        for d in lev:
            nxt.append(2 * d - 1)
            if d % 6 == 3: nxt.append((4 * d - 3) // 3)
            elif d % 6 == 5: nxt.append((4 * d + 1) // 3)
        lev = nxt; out.append(len(lev))
    return out
L = levels(21, 30)
print("level sizes above 21:", L[:16], "... ratio of last two: %.4f" % (L[-1] / L[-2]))
L = levels(27, 30); print("level sizes above 27:", L[:12], " ratio %.4f" % (L[-1] / L[-2]))
