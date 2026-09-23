#!/usr/bin/env python3
"""collatz_procgen_20260922_mirror_duality.py -- forward/backward numerator duality test (Q1 mirror lane).

Inputs (dimension-lane dumps, see collatz_procgen_20260922_exceptional_dimension.md sec. 4.5, 4.7):
  FWD_CENSUS  hostile_z13_M61_J25.txt   the 343 certified hostile -p/3^j, j <= 25, 1 <= p/3^j < 3/2
  FWD_BAD     fwd_bad_m61.txt, fwd_bad_m63.txt   exceptional classes mod 2^61 / 2^63 (forward)
  BWD_BAD     bwd_bad_r41.txt               exceptional classes mod 3^42 (backward)
Steps:
 (1) forward numerator set: the 343 census points plus the j = 26, 27 points of Bad_63 with 1 <= |x| < 3/2
     (39 points, certified hostile in the dimension lane), and the forward candidates with 3/2 < |x| <= 1.6
     (j <= 27) for reference;
 (2) backward census: every positive dyadic p/2^e (e <= 45) with value in [0.3, 2] lying in a class of
     Bad_41 (the dimension lane used e <= 40 and found 98 points);
 (3) every backward census point is classified by the backward prover (finite criterion below 3/2,
     slow-descent search above 3/2);
 (4) shared numerators and their status on both sides; the generation-1 form N = 3^j + 2^(e-1) and the
     straddle identity (|x| - 3/2)(y - 3/2) = -(rho - 1)^2/(2 rho), rho = 3^j/2^e;
 (5) the generation-1 backward points 1/2 + 3^j/2^e at every clock with 1 < 3^j/2^e < 1.06, e <= 200:
     slow descents found or not (their forward partners -1 - 2^(e-1)/3^j are hostile by Theorem F with
     alpha(e-1) = a0(e-1), which the loop DP gives for 10 <= e-1 <= 4000).
Usage: python3 ..._mirror_duality.py DIMDIR   (DIMDIR contains fwd30/ and bwd17/ as in scratch/procgen_dim)
"""
import sys, os, math
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("bp", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                               "collatz_procgen_20260922_mirror_bwd_prover.py"))
bp = importlib.util.module_from_spec(spec); spec.loader.exec_module(bp)

def load_classes(fn): return [int(l.split()[0]) for l in open(fn) if l.strip()]

def fwd_points(bad, M, jmin, jmax, lo, hi):
    """-p/3^j (jmin<=j<=jmax) with lo <= p/3^j < hi lying in a class of bad (mod 2^M)"""
    mod = 1 << M; out = []
    for j in range(jmin, jmax + 1):
        t = 3 ** j
        for c in bad:
            p = (-c * t) % mod                         # x = -p/3^j = c mod 2^M  <=>  p = -c 3^j
            if lo * t <= p < hi * t: out.append(F(-p, t))
    return sorted(set(out))

def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "scratch/procgen_dim"
    cen = [F(l.strip()) for l in open(f"{d}/fwd30/hostile_z13_M61_J25.txt") if l.strip()]
    bad63 = load_classes(f"{d}/fwd30/fwd_bad_m63.txt"); bad61 = load_classes(f"{d}/fwd30/fwd_bad_m61.txt")
    ext = fwd_points(bad63, 63, 26, 27, 1, F(3, 2))
    above = fwd_points(bad61, 61, 0, 25, F(3, 2), F(8, 5)) + fwd_points(bad63, 63, 26, 27, F(3, 2), F(8, 5))
    fwd = sorted(set(cen) | set(ext))
    above = sorted(set(above))
    print(f"(1) forward hostile census: {len(cen)} points (j<=25) + {len(set(ext) - set(cen))} new with j=26,27"
          f" = {len(fwd)}; forward candidates with 3/2 < |x| <= 1.6 (j <= 27): {len(above)}")
    fwdN = {x.numerator * -1: x for x in fwd}; aboveN = {x.numerator * -1: x for x in above}
    # (2) backward census
    bb = load_classes(f"{d}/bwd17/bwd_bad_r41.txt"); MODB = 3 ** 42; bpts = set()
    for e in range(0, 46):
        q = 2 ** e
        for cc in bb:
            p = (cc * q) % MODB
            if q * 3 // 10 <= p <= 2 * q and (p % 2 == 1 or e == 0): bpts.add(F(p, q))
    bpts = sorted(bpts)
    n40 = sum(1 for y in bpts if y.denominator <= 2 ** 40)
    print(f"(2) backward census (Bad_41, value in [0.3,2]): {len(bpts)} points with e <= 45 ({n40} with e <= 40);"
          f" {sum(1 for y in bpts if y > F(3,2))} above 3/2")
    # (3) classify
    status = {}
    for y in bpts:
        if y < F(3, 2):
            r = bp.hostile_finite(y)
            status[y] = r[0] if r[0] != 'descends' else f"descends{r[1][:6]}"
        else:
            r = bp.slow_descent(y, K1MAX=400)
            if r[0] == 'none': status[y] = f"above 3/2, no greedy descent (k1<=400)"
            else:
                k1, s, K, ks = r
                assert bp.check_path(y, ks)
                status[y] = f"above 3/2, DESCENDS (k1={k1}, 2^{K}<3^{s})"
    from collections import Counter
    cnt = Counter(v.split(' (')[0].split('[')[0] for v in status.values())
    print("(3) backward census statuses:", dict(cnt))
    for y in bpts:
        if not status[y].startswith('hostile'):
            print(f"      {str(y):>30} = {float(y):.6f}: {status[y]}")
    # (4) shared numerators
    print("(4) numerators shared between the forward census (hostile, j<=27) or forward candidates above 3/2,"
          " and the backward census:")
    shared = []
    for y in bpts:
        N = y.numerator
        side = 'fwd hostile' if N in fwdN else ('fwd candidate >3/2' if N in aboveN else None)
        if side:
            x = fwdN.get(N, aboveN.get(N)); j = int(round(math.log(x.denominator, 3))) if x.denominator > 1 else 0
            e = y.denominator.bit_length() - 1
            gen1 = (N == 3 ** j + 2 ** (e - 1)) if e >= 1 else False
            shared.append((N, x, y, side, status[y], gen1, j, e))
    for N, x, y, side, st, gen1, j, e in shared:
        rho = F(3 ** j, 2 ** e) if e >= 0 else None
        strad = ""
        if gen1:
            lhs = (abs(x) - F(3, 2)) * (y - F(3, 2)); rhs = -(rho - 1) ** 2 / (2 * rho)
            strad = f"; gen-1 N=3^{j}+2^{e-1}, rho=3^{j}/2^{e}={float(rho):.5f}, straddle identity {lhs == rhs}"
        print(f"      N={N}: forward {x} (|x|={float(abs(x)):.5f}, {side}); backward {y} (= {float(y):.5f}): {st}{strad}")
    # also: every census point p/2^e with p odd -- does -p/3^j lie in a forward exceptional class at all?
    # (5) generation-1 backward points at the lower-convergent clocks
    print("(5) generation-1 backward points 1/2 + 3^j/2^e with 1 < 3^j/2^e < 1.06, e <= 200, search k1 <= 2500 (forward partner"
          " -1 - 2^(e-1)/3^j is hostile by Theorem F, alpha(e-1) = a0(e-1)):")
    for e in range(2, 201):
        j = math.ceil(e * math.log(2) / math.log(3))
        while 3 ** j <= 2 ** e: j += 1
        rho = F(3 ** j, 2 ** e)
        if not (1 < rho < F(106, 100)): continue
        y = F(1, 2) + rho
        r = bp.slow_descent(y, K1MAX=2500)
        if r[0] == "none": st = f"no greedy descent for k1 <= 2500 (closest final R {r[1]:.4f})"
        else:
            k1, s, K, ks = r; ok = bp.check_path(y, ks)
            st = f"DESCENDS: k1={k1}, then greedy, 2^{K} < 3^{s} (R={float(F(2**K,3**s)):.6f}), replay {ok}"
        print(f"      (e,j)=({e},{j}) rho={float(rho):.6f} y={float(y):.6f}: {st}")

def fwd_slow_descent(x, AEXTRA=400, STEPS=20000):
    """forward route for x = -N/3^j with |x| > 3/2: M^a (a >= j, the value becomes an integer), then the
    Collatz map on negative integers (halve even, 3v+1 odd) until a cycle is entered; a descent is a
    halving at which 3^(a_total) < 2^b.  Returns (a, a_total, b) of the first descent found, else None."""
    j = 0; d = x.denominator
    while d % 3 == 0: d //= 3; j += 1
    for a in range(j, j + AEXTRA):
        v = x * 3 ** a + F(3 ** a - 1, 2)
        assert v.denominator == 1
        v = v.numerator; at = a; b = 0
        for _ in range(STEPS):
            if v in (-1, -5, -7, -17): break
            if v % 2 == 0:
                v //= 2; b += 1
                if 3 ** at < 2 ** b: return (a, at, b)
            else:
                v = 3 * v + 1; at += 1
    return None

def upper_clocks(EMAX=45):
    print(f"(6) the reverse straddle: clocks with 0.94 < rho = 3^j/2^e < 1, e <= {EMAX}: backward point 1/2 + rho < 3/2"
          " (backward prover), forward partner -1 - 2^(e-1)/3^j above 3/2 (route M^a + negative Collatz, a <= j+400):")
    for e in range(2, EMAX + 1):
        j = math.floor(e * math.log(2) / math.log(3))
        while 3 ** (j + 1) < 2 ** e: j += 1
        rho = F(3 ** j, 2 ** e)
        if not (F(94, 100) < rho < 1): continue
        y = F(1, 2) + rho; x = -(1 + F(2 ** (e - 1), 3 ** j))
        rb = bp.hostile_finite(y)
        rf = fwd_slow_descent(x)
        fst = f"DESCENDS (a={rf[0]}, 3^{rf[1]} < 2^{rf[2]})" if rf else "no route descent found"
        print(f"      (e,j)=({e},{j}) rho={float(rho):.6f}: backward {y} = {float(y):.6f}: {rb[0]}"
              f" ({rb[1] if rb[0] != 'descends' else rb[1][:8]}); forward {x} (|x|={float(-x):.6f}): {fst}")

def forward_above(d, AEXTRA=400):
    """(7) every forward candidate -p/3^j (j <= 27) with 3/2 < |x| <= 1.6 in Bad_61 / Bad_63: route search."""
    bad63 = load_classes(f"{d}/fwd30/fwd_bad_m63.txt"); bad61 = load_classes(f"{d}/fwd30/fwd_bad_m61.txt")
    above = sorted(set(fwd_points(bad61, 61, 0, 25, F(3, 2), F(8, 5)) + fwd_points(bad63, 63, 26, 27, F(3, 2), F(8, 5))))
    found = []; none = []
    for x in above:
        r = fwd_slow_descent(x, AEXTRA)
        (found if r else none).append((x, r))
    print(f"(7) forward candidates with 3/2 < |x| <= 1.6 (j <= 27): {len(above)}; explicit route descents "
          f"(M^a, a <= j+{AEXTRA}, then negative Collatz): {len(found)}; none found: {len(none)}")
    for x, r in found:
        print(f"      {str(x):>28} (|x|={float(-x):.6f}): DESCENDS a={r[0]}, 3^{r[1]} < 2^{r[2]}")
    for x, r in none:
        print(f"      {str(x):>28} (|x|={float(-x):.6f}): no route descent found")

if __name__ == '__main__':
    main()
    upper_clocks(70)
    forward_above(sys.argv[1] if len(sys.argv) > 1 else "scratch/procgen_dim")
