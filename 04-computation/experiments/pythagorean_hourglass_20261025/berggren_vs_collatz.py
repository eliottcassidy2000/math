# The Pythagorean "Collatz problem" that is a theorem: the Berggren parent map on primitive triples strictly lowers the
# hypotenuse c (a Lorentzian potential), so every primitive triple descends to (3,4,5).  Compare its depth statistics
# with Collatz total stopping times: both grow like a constant times log of the size, with Gaussian fluctuations
# (for Berggren a CLT for the even-continued-fraction cocycle of Romik's dynamical system is classical-type; for
# Collatz it is the conjectural random-walk model).  mac-mini-2026-10-25-pyth.
import math, statistics
from math import gcd, isqrt

def parent(a, b, c):
    """inverse Berggren step: the unique one of A^-1, B^-1, C^-1 giving a positive primitive triple"""
    for (x, y, w) in ((a - 2 * b + 2 * c, 2 * a - b + 2 * c, 2 * a - 2 * b + 3 * c),     # A^-1 (signs as needed)
                      (a + 2 * b - 2 * c, -2 * a - b + 2 * c, -2 * a - 2 * b + 3 * c),
                      (-a - 2 * b + 2 * c, -2 * a - b + 2 * c, -2 * a - 2 * b + 3 * c),
                      (a - 2 * b + 2 * c, 2 * a - b + 2 * c, -(2 * a - 2 * b + 3 * c)),
                      (-a + 2 * b + 2 * c, -2 * a + b + 2 * c, -2 * a + 2 * b + 3 * c),
                      (a + 2 * b + 2 * c, 2 * a + b + 2 * c, 2 * a + 2 * b + 3 * c)):
        pass
    # standard inverse matrices (Hall): the parent of (a,b,c) is |M| (a,b,c) for M one of the three
    cands = [(a - 2 * b + 2 * c, 2 * a - b + 2 * c, 2 * a - 2 * b + 3 * c),
             (a + 2 * b - 2 * c, 2 * a + b - 2 * c, 2 * a + 2 * b - 3 * c),
             (-a + 2 * b - 2 * c, -2 * a + b - 2 * c, -2 * a + 2 * b - 3 * c)]
    # the parent has strictly smaller hypotenuse: c' = |3c - 2(+-a +-b)|; choose the reduction (|a|,|b|,|c|) minimising c
    best = None
    for sa in (1, -1):
        for sb in (1, -1):
            x = abs(a * 1 + 0)  # placeholder (overwritten below)
            aa, bb = sa * a, sb * b
            np_ = (abs(aa - 2 * bb + 2 * c), abs(2 * aa - bb + 2 * c), abs(2 * aa - 2 * bb + 3 * c))
            if best is None or np_[2] < best[2]: best = np_
    return best

def descend(a, b, c):
    steps = 0
    while (a, b, c) not in ((3, 4, 5), (4, 3, 5)):
        na, nb, nc = parent(a, b, c)
        assert nc < c, (a, b, c, na, nb, nc)
        assert na * na + nb * nb == nc * nc
        a, b, c = na, nb, nc; steps += 1
        if steps > 10 ** 6: raise RuntimeError
    return steps

# all primitive triples with c <= Cmax via Euclid parameters
Cmax = 200000
depth = {}
for m in range(2, isqrt(Cmax) + 1):
    for n in range(1, m):
        if (m - n) % 2 == 1 and gcd(m, n) == 1:
            c = m * m + n * n
            if c > Cmax: break
            a, b = m * m - n * n, 2 * m * n
            depth[(a, b, c)] = descend(a, b, c)
print(f"Berggren descent: {len(depth)} primitive triples with c <= {Cmax}: every one reaches (3,4,5) with c strictly decreasing")
# depth vs log c in windows
for lo in (10 ** 3, 10 ** 4, 10 ** 5):
    ds = [d for (a, b, c), d in depth.items() if lo <= c < 2 * lo]
    if ds:
        mu = statistics.mean(ds); sd = statistics.pstdev(ds)
        print(f"  c in [{lo}, {2 * lo}): {len(ds)} triples, depth mean {mu:.2f} = {mu / math.log(1.5 * lo):.3f} log c, sd {sd:.2f} (sd^2/log c = {sd * sd / math.log(1.5 * lo):.3f}), max {max(ds)}")
# Collatz total stopping times for comparison (shortcut map T(n) = n/2 or (3n+1)/2)
def tst(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2; s += 1
    return s
for lo in (10 ** 3, 10 ** 4, 10 ** 5):
    ts = [tst(n) for n in range(lo, 2 * lo)]
    mu = statistics.mean(ts); sd = statistics.pstdev(ts)
    print(f"  Collatz n in [{lo}, {2 * lo}): stopping time mean {mu:.2f} = {mu / math.log(1.5 * lo):.3f} log n, sd {sd:.2f} (sd^2/log n = {sd * sd / math.log(1.5 * lo):.3f}), max {max(ts)}")
print("  heuristic Collatz constant: mean T-steps ~ log n / log(2/sqrt 3) =", f"{1 / math.log(2 / math.sqrt(3)):.3f} log n")

# Depth = a sum of continued-fraction partial quotients?  For Euclid parameters (m, n), m > n, the Berggren children
# of (m, n) are (2m - n, m), (2m + n, m), (m + 2n, n), so the parent map is a subtractive Euclidean algorithm.
# Compare the Berggren depth with S(m/n) = sum of the partial quotients of m/n, and the window means with the
# Yao-Knuth law for sums of partial quotients, mean ~ (6/pi^2) (log M)^2 with M ~ sqrt(c).
def cf_sum(m, n):
    s = 0
    while n:
        q, r = divmod(m, n); s += q; m, n = n, r
    return s
ratios = []; rows = []
for m in range(2, isqrt(Cmax) + 1):
    for n in range(1, m):
        if (m - n) % 2 == 1 and gcd(m, n) == 1:
            c = m * m + n * n
            if c > Cmax: break
            a, b = m * m - n * n, 2 * m * n
            rows.append((c, depth[(a, b, c)], cf_sum(m, n)))
import numpy as np
R = np.array(rows, dtype=float)
coef = np.polyfit(R[:, 2], R[:, 1], 1)
print(f"\nBerggren depth vs sum of partial quotients S(m/n): depth ~ {coef[0]:.4f} S + {coef[1]:.3f}; corr {np.corrcoef(R[:, 1], R[:, 2])[0, 1]:.5f}; "
      f"max |depth - (S - 2)| = {np.max(np.abs(R[:, 1] - (R[:, 2] - 2))):.0f}")
for lo in (10 ** 3, 10 ** 4, 10 ** 5):
    sel = (R[:, 0] >= lo) & (R[:, 0] < 2 * lo)
    L = math.log(math.sqrt(1.5 * lo))
    print(f"  c in [{lo}, {2 * lo}): mean depth {R[sel, 1].mean():.2f}; (6/pi^2)(log sqrt c)^2 = {6 / math.pi ** 2 * L * L:.2f}; mean S {R[sel, 2].mean():.2f}")
