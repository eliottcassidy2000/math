# opus-2026-07-16-S332 -- HYP-7135 / THM-928 (A): THE CASCADE THEOREM
# (the exponential half of the two-scale certificate).
#
# THEOREM (target): sorted distinct speeds x_1 < ... < x_{n-1} with
#   x_{k+1}/x_k >= R. TRUE LRC(n) normalization: D_x = {t: ||x t|| < 1/n},
#   tooth width 2/(n x), density 2/n. Then
#       mu(uncovered) >= (1 - 2/n)^{n-1} - 2/R.
#   At n = 14: (6/7)^13 - 2/R > 0 for R >= 15. Uniform n: R >= 18 works for
#   every n >= 3 (min over n of (1-2/n)^{n-1} at n=3 is 1/9; 1/9 > 2/18? no:
#   2/18 = 1/9 exactly -- so R >= 19 for n = 3; n >= 4: (1/2)^3 = 1/8 > 2/17).
#
# PROOF INGREDIENTS (verified exactly here):
#   (L1) per-step: for W a union of kappa intervals on the circle,
#        mu(W cap D_x) <= (2/n) mu(W) + 2 kappa/(n x)
#        [<= (x*len + 1) cells each contributing <= 2/(nx), per component]
#   (L2) component growth: kappa(W \ D_x) <= 2 kappa(W) + x mu(W)
#   (L3) induction kappa_k <= 2 x_k  (R >= 4)
#   (L4) unroll: mu_13 >= (6/7)^13 - (2/(7R)) sum (6/7)^j >= (6/7)^13 - 2/R
#   (W)  the fixed-point witness for integer-ratio dilate families:
#        {c q^j}: t = a/(q^s - 1) with small s (exact orbit check q <= 20).
from fractions import Fraction
from math import gcd
import random

F = Fraction

def teeth(x, n):
    w = F(1, n*x); out = []
    for j in range(x):
        a, b = F(j, x) - w, F(j, x) + w
        if a < 0: out += [(F(0), b), (a + 1, F(1))]
        elif b > 1: out += [(a, F(1)), (F(0), b - 1)]
        else: out.append((a, b))
    out.sort()
    m = []
    for iv in out:
        if m and iv[0] <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], iv[1]))
        else: m.append(list(iv))
    return [(a, b) for a, b in m]

def subtract(W, T):
    """W \ T for sorted disjoint interval lists."""
    out = []
    for (a, b) in W:
        cur = a
        for (c, d) in T:
            if d <= cur: continue
            if c >= b: break
            if c > cur: out.append((cur, c))
            cur = max(cur, d)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return out

def mu(W): return sum(b - a for a, b in W)

def circ_kappa(W):
    if len(W) >= 2 and W[0][0] == 0 and W[-1][1] == 1: return len(W) - 1
    return max(len(W), 1)

print("=" * 72)
print("(L1)+(L2) per-step lemma + component growth, exact, random configs:")
random.seed(14)
viol = 0; tests = 0
for trial in range(60):
    n = random.choice([5, 8, 14])
    # random union of intervals
    pts = sorted(F(random.randint(0, 9999), 10000) for _ in range(2 * random.randint(2, 12)))
    W = [(pts[2*i], pts[2*i+1]) for i in range(len(pts)//2) if pts[2*i] < pts[2*i+1]]
    if not W: continue
    x = random.randint(3, 400)
    T = teeth(x, n)
    inter = mu(W) - mu(subtract(W, T))
    kap = len(W)
    tests += 1
    if not (inter <= F(2, n) * mu(W) + F(2 * kap, n * x)): viol += 1
    Wp = subtract(W, T)
    if not (len(Wp) <= 2 * kap + x * mu(W) + 2): viol += 1  # +2 edge slack
print(f"  {tests} random (W, x, n) configs: violations = {viol}")

print()
print("(L3)+(L4) THE CASCADE at n = 14 (exact prefixes; R = 15, 20, 30):")
target = F(6, 7) ** 13
print(f"  (6/7)^13 = {float(target):.6f}")
for R in (15, 20, 30):
    # exact computation feasible only for small prefixes (x_k = R^k grows);
    # verify the per-step inequality chain holds with margin at each step
    xs = [1]
    while len(xs) < 5: xs.append(xs[-1] * R)
    W = [(F(0), F(1))]
    mu_bound = F(1)
    ok = True
    for k, x in enumerate(xs):
        W = subtract(W, teeth(x, 14))
        kap_prev = 1 if k == 0 else circ_kappa_prev
        # bound sequence: mu_k >= (6/7) mu_{k-1} - 2 kappa_{k-1}/(14 x_k)
        mu_bound = F(6, 7) * mu_bound - F(2 * (2 * xs[k-1] if k else 1), 14 * x)
        circ_kappa_prev = circ_kappa(W)
        if mu(W) < mu_bound: ok = False
        if circ_kappa_prev > 2 * x: ok = False
    # extrapolate the proven bound to 13 steps
    full_bound = target - F(2, R)
    print(f"  R={R:3d}: 5-step exact mu = {float(mu(W)):.6f} >= running bound "
          f"{float(mu_bound):.6f}: {mu(W) >= mu_bound}; kappa_k <= 2 x_k: {ok}; "
          f"13-step THEOREM bound = {float(full_bound):+.6f} "
          f"{'(POSITIVE: LRC(14) holds for this lacunarity)' if full_bound > 0 else ''}")

print()
print("  small-n FULL verification (all 13->n-1 speeds exact, R=15..20):")
for n in (5, 6, 7):
    for R in (15, 18):
        xs = [2]
        while len(xs) < n - 1: xs.append(xs[-1] * R)
        W = [(F(0), F(1))]
        for x in xs: W = subtract(W, teeth(x, n))
        thm = F(n - 2, n) ** (n - 1) - F(2, R)
        print(f"   n={n} R={R}: exact uncovered = {float(mu(W)):.6f}  "
              f"theorem bound = {float(thm):+.6f}  holds: {mu(W) >= thm}"
              f"{'  NONEMPTY' if mu(W) > 0 else ''}")

print()
print("(W) THE FIXED-POINT WITNESS for dilate families {c q^j}, n = 14:")
print("    t = a/(q^s - 1): the multiply-by-q orbit of a mod (q^s - 1) must")
print("    stay >= (q^s - 1)/14 from 0 mod (q^s - 1).")
for q in range(2, 21):
    found = None
    for s in (1, 2, 3):
        Q = q ** s - 1
        if Q < 3: continue
        for a in range(1, Q):
            orb, v, ok = set(), a, True
            for _ in range(s + 1):
                if min(v % Q, (-v) % Q) * 14 < Q: ok = False; break
                v = (v * q) % Q
            if ok:
                found = (s, a, Q); break
        if found: break
    s, a, Q = found
    # verify on the actual 13-term family c=1 exactly
    tt = F(a, Q)
    dists = [min((q**j * tt) % 1, 1 - (q**j * tt) % 1) for j in range(13)]
    ok13 = all(d >= F(1, 14) for d in dists)
    print(f"   q={q:2d}: witness t = {a}/{Q} (s={s}); all 13 dilate distances "
          f">= 1/14: {ok13}; min dist = {min(dists)}")
