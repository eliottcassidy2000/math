# opus-2026-07-15-S313 -- HYP-6920 stream 2: Lemma A (skew-Koksma rate) + exact
# AP floors + positivity mechanisms + exceptional scan.
# (A) EXACT integral_E a_dvec(t) dt by breakpoint sweep (piecewise linear,
#     trapezoid at exact Fractions -- replaces the withdrawn quadrature).
# (B) Lemma A empirics: |uncovered(N) - integral| * N over ranges of N.
# (C) B1 (q<=6 pigeonhole window) / B2 (q=13 adjacency dilation) mechanisms:
#     per-prefix scan; exceptional (P,R) list.
from fractions import Fraction
from math import gcd
import itertools

W = Fraction(1, 13)   # half-width of danger arc

def safe_set(P):
    ivs = [(Fraction(0), Fraction(1))]
    for q in P:
        bands = [(Fraction(13*k+1, 13*q), Fraction(13*(k+1)-1, 13*q)) for k in range(q)]
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def mu(ivs): return sum(b - a for (a, b) in ivs)

def union_measure_at(ds, t):
    # mu(union_i (W - d_i t)) at exact t: arcs centered -d_i t, half-width 1/13
    arcs = []
    for d in ds:
        c = (-d * t) % 1
        lo, hi = c - W, c + W
        if lo < 0: arcs.extend([(lo % 1, Fraction(1)), (Fraction(0), hi)])
        elif hi > 1: arcs.extend([(lo, Fraction(1)), (Fraction(0), hi % 1)])
        else: arcs.append((lo, hi))
    arcs.sort()
    tot, cur = Fraction(0), Fraction(0)
    for (lo, hi) in arcs:
        if hi <= cur: continue
        tot += hi - max(lo, cur)
        cur = hi
    return tot

def exact_integral_a(E, ds):
    # breakpoints: endpoints of arcs collide/wrap: (d_i - d_j) t = k +- 2/13, k
    # and d_i t = k +- 1/13 (wrap events). Collect within each E-component.
    bps = set()
    for (a, b) in E: bps.update([a, b])
    diffs = set()
    for di, dj in itertools.combinations(ds, 2): diffs.add(abs(di - dj))
    for d in ds: diffs.add(d)
    for d in diffs:
        if d == 0: continue
        # t = (k + r)/d for r in {0, +-2/13, +-1/13}
        for k in range(0, d + 1):
            for r in (Fraction(0), Fraction(2,13), Fraction(-2,13),
                      Fraction(1,13), Fraction(-1,13)):
                t = Fraction(k + 0, 1)  # placeholder
                t = (Fraction(k) + r) / d
                if 0 < t < 1: bps.add(t)
    bps = sorted(bps)
    total = Fraction(0)
    for (lo, hi) in E:
        cells = [lo] + [t for t in bps if lo < t < hi] + [hi]
        for c0, c1 in zip(cells, cells[1:]):
            # a(t) = 1 - union: linear on the cell; trapezoid exact
            a0 = 1 - union_measure_at(ds, c0 + (c1-c0)/1000)   # avoid corner
            a1 = 1 - union_measure_at(ds, c1 - (c1-c0)/1000)
            # use interior sample endpoints (the corner value may differ by 0-measure)
            total += (a0 + a1) / 2 * (c1 - c0)
    return total

def comb_teeth_in(x, a, b):
    import math
    w = Fraction(1, 13*x)
    out = []
    for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
        lo, hi = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def subtract_comb(ivs, x):
    out = []
    for (a, b) in ivs:
        cur = a
        for (lo, hi) in sorted(comb_teeth_in(x, a, b)):
            if lo > cur: out.append((cur, min(lo, b)))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return [iv for iv in out if iv[0] < iv[1]]

# ---- (A)+(B) for prefix {1..5}, consecutive pattern
P = [1, 2, 3, 4, 5]
E = safe_set(P)
ds = [0, 1, 2, 3, 4, 5, 6]
Ia = exact_integral_a(E, ds)
print(f"(A) EXACT integral_E a_(0..6) = {Ia} ~ {float(Ia):.6f}  "
      f"(S312 quadrature said ~0.0558)")

print("(B) |uncovered(N) - integral| * N   (residues: x_i = N + d_i):")
for N in (100, 200, 400, 800, 1600):
    U = list(E)
    for d in ds: U = subtract_comb(U, N + d)
    err = mu(U) - Ia
    print(f"    N={N:5d}: uncovered = {float(mu(U)):.6f}  err*N = {float(err*N):+.3f}")

# ---- (C) mechanisms scan over all 792 prefixes
def b1_windows(P):
    outs = []
    for q in range(2, 7):
        for p in range(1, q):
            if gcd(p, q) != 1: continue
            t = Fraction(p, q)
            ok = all(min((qq*t) % 1, 1 - (qq*t) % 1) >= Fraction(2,13) for qq in P)
            if ok: outs.append((p, q))
    return outs

def b2_ps(P, R):
    # allowed p: q'p not in {0, +-1} mod 13 for q' in P ; adjacency: exists
    # c,c' in S1: p(c - c') = +-1 mod 13, S1 = {r - r0 mod 13}
    r0 = R[0]
    S1 = {(r - r0) % 13 for r in R}
    D = {(x - y) % 13 for x in S1 for y in S1 if x != y}
    ok = []
    for p in range(1, 13):
        if any((qq * p) % 13 in (0, 1, 12) for qq in P): continue
        if any((p * d) % 13 in (1, 12) for d in D): ok.append(p)
    return ok

exceptions = []
for Pset in itertools.combinations(range(1, 13), 5):
    P5 = list(Pset)
    R7 = [r for r in range(1, 13) if r not in Pset]
    if b1_windows(P5): continue
    if b2_ps(P5, R7): continue
    exceptions.append((P5, R7))
print(f"\n(C) prefixes with NEITHER mechanism: {len(exceptions)} of 792")
for (P5, R7) in exceptions[:10]: print(f"    P={P5} R={R7}")
