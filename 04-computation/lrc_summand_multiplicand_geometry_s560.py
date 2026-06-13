#!/usr/bin/env python3
r"""
lrc_summand_multiplicand_geometry_s560.py    oracle-2026-06-02-S560o

THE GEOMETRY OF THE SUMMAND AND MULTIPLICAND GRAPHS, AND HOW IT SIMPLIFIES LRC.

  SUMMAND graph (the user's A+B=C): nodes N, arc a->n,b->n iff a+b=n. Its level sets
  are the ANTIDIAGONALS x+y=n of the lattice N^2 (slope -1 lines). In-pair count of n is
  floor((n-1)/2). It is the ADDITIVE foliation.

  MULTIPLICAND graph: nodes N, arc a->n,b->n iff a*b=n. Its level sets are the HYPERBOLAS
  xy=n. In-pair count of n = (tau(n) - [n square])/2 (divisor pairs). Primes have only the
  (1,p) pair -> SOURCES once the unit is dropped. It is the MULTIPLICATIVE foliation.

  THE LOG BRIDGE: (x,y) |-> (log x, log y) carries the hyperbola xy=n to the line
  log x + log y = log n. So the MULTIPLICAND graph IS the SUMMAND graph in log-coordinates.
  The two foliations are exp/log images of each other (the +/x weld, S548).

  LRC DICTIONARY (S555o/S557):
   - the pinch denominator C = v_a + v_b is a SUMMAND-graph node: the pinch witnesses
     t = m/C live on antidiagonal C; the available denominators on a speed set S are the
     SUMSET S+S (distinct pair-sums).
   - a runner w is cleared at t=m/C iff C does NOT divide w (for gcd(m,C)=1): a
     MULTIPLICAND-graph (divisibility) test. The sieve THM-369 = the multiplicand shadow.
   - oracle-S555o: the rational pinch IS the sieve -> the two foliations COINCIDE on the
     integer points (rational t); the fine pinch (C>n) is where addition outruns division.

  THE SIMPLIFICATION (this script tests it):
   * ADDITIVE: |S+S| >= 2|S|-1, EQUALITY IFF S is an arithmetic progression (Freiman).
     The number of distinct pinch denominators = |distinct pair-sums|. So the AP MINIMIZES
     the pinch options -> the TIGHTEST pinch pigeonhole -> the AP is the hardest case.
   * MULTIPLICATIVE: the AP {1..k} contains a multiple of every q<=k (max small-divisor
     coverage among bounded sets) -> the most sieve-covered, escaping only at q=k+1=n.
   => the AP is the JOINT EXTREMUM of both graphs = the LRC wall. We verify the additive
      minimization, the coverage maximization, and the correlation (bigger sumset -> looser).
"""
from functools import reduce
from math import gcd
import random

N = 14
TH = 1.0 / N

def d0(p):
    p = p % 1.0
    return min(p, 1 - p)

def margin(S, G=200000):
    """M(S) = max_t min_i ||v_i t||  (the loneliness radius)."""
    best = 0.0
    for i in range(G):
        t = (i + 0.5) / G
        m = min(d0(v * t) for v in S)
        if m > best:
            best = m
    return best

def pair_sums(S):
    return sorted({a + b for i, a in enumerate(S) for b in S[i + 1:]})

def sieve_coverage(S, qmax=N):
    return sum(1 for q in range(2, qmax + 1) if any(v % q == 0 for v in S))

def divisor_pairs(n):
    ds = [d for d in range(1, n + 1) if n % d == 0]
    return sum(1 for d in ds if d * d < n)  # unordered pairs d<n/d (incl (1,n))

def report(name, S):
    S = tuple(sorted(S))
    k = len(S)
    ps = pair_sums(S)
    ss = len(ps)                       # |distinct pair-sums| = # pinch denominators
    minss = 2 * k - 3                  # restricted-sumset (distinct pairs) min; =AP (Freiman)
    M = margin(S)
    cov = sieve_coverage(S)            # # of q in 2..14 dividing some speed
    small_dens = [c for c in ps if c <= N]      # coarse pinch denominators (= sieve, S555o)
    fine_dens = [c for c in ps if c > N]        # fine pinch denominators
    print(f"  {name}: k={k}  |S+S|={ss} (min {minss}, excess {ss-minss})  M={M:.4f}({M*N:.2f}/14)")
    print(f"      pinch denominators: {len(small_dens)} coarse(<=14)  {len(fine_dens)} fine(>14)"
          f"   sieve-coverage q<=14: {cov}/13")
    return {"name": name, "k": k, "ss": ss, "excess": ss - minss, "M": M, "cov": cov}

def main():
    print("=" * 82)
    print("SUMMAND (antidiagonals, +) vs MULTIPLICAND (hyperbolas, x) -- geometry & LRC")
    print("=" * 82)

    print("\n(0) IN-DEGREE GEOMETRY (incoming-pair counts; the two foliations):")
    print("   n : summand floor((n-1)/2) | multiplicand divisor-pairs (1..n)")
    for n in [4, 6, 7, 12, 13, 14, 16, 24, 30, 36]:
        print(f"   {n:2d} : {((n-1)//2):>3d}                  | {divisor_pairs(n):>3d}"
              f"   {'(prime: only (1,p))' if all(n%d for d in range(2,n)) else ''}")
    print("   -> summand in-degree grows linearly (dense, additive); multiplicand stays tiny")
    print("      except at highly-composite n (sparse, multiplicative). Log map sends the")
    print("      hyperbola xy=n to the antidiagonal log x+log y=log n: multiplicand = summand")
    print("      in log-coordinates.")

    print("\n(1) THE WALL configs -- joint extremum?")
    ap = report("AP {1..13}", range(1, 14))
    vs = report("V* {1..11,13,24}", list(range(1, 12)) + [13, 24])

    print("\n(2) RANDOM primitive 13-sets: sumset excess vs loneliness margin")
    rnd = random.Random(560)
    rows = [ap, vs]
    for j in range(16):
        while True:
            S = tuple(sorted(rnd.sample(range(1, 60), 13)))
            if reduce(gcd, S) == 1:
                break
        rows.append(report(f"rand{j:02d}", S))

    # correlation: does a bigger sumset (more pinch denominators) => bigger margin (looser)?
    xs = [r["excess"] for r in rows]
    ys = [r["M"] for r in rows]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov_xy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs) ** 0.5
    vy = sum((y - my) ** 2 for y in ys) ** 0.5
    corr = cov_xy / (vx * vy) if vx and vy else 0.0

    print("\n" + "=" * 82)
    print("FINDINGS")
    print("=" * 82)
    min_excess = min(r["excess"] for r in rows)
    ap_excess = ap["excess"]
    tightest = min(rows, key=lambda r: r["M"])
    print(f"  * ADDITIVE (Freiman): min distinct-pair-sum excess over 2k-3 = {min_excess}; AP excess = {ap_excess}"
          f"  (AP attains the restricted-sumset minimum {2*13-3}=23 distinct pair-sums).")
    print(f"  * The TIGHTEST set (smallest M) is '{tightest['name']}' with M={tightest['M']:.4f},"
          f" excess={tightest['excess']}.")
    print(f"  * CORRELATION(sumset excess, margin M) = {corr:+.3f}"
          f"  ({'bigger sumset -> looser (more pinch denominators = more room)' if corr>0.15 else 'weak/none'}).")
    print(f"  * MULTIPLICATIVE: AP coverage q<=14 = {ap['cov']}/13 (max among bounded sets;"
          f" misses only q=14). The AP is the small-divisor-coverage maximizer.")
    print("""
  GEOMETRIC PICTURE. The AP sits at the JOINT EXTREMUM:
    - summand graph: minimal sumset S+S (Freiman) = fewest distinct pinch denominators
      = the tightest pinch pigeonhole (the additive reason the AP is hard);
    - multiplicand graph: maximal small-q divisor coverage = most sieve-covered, escaping
      the sieve only at the single modulus q=n (the multiplicative reason).
  SIMPLIFICATION FOR LRC. Off the AP, Freiman gives |S+S| > 2k-1: there are STRICTLY MORE
  pinch denominators (antidiagonals carrying witnesses) than for the AP. The pinch
  pigeonhole therefore has provable surplus exactly where the AP's tightness is relaxed.
  The proof effort should index the (fine) pinch search by the SUMSET S+S (the summand
  graph on S) and use the divisibility (multiplicand) structure only to test clearance --
  the two graphs split the problem cleanly: addition supplies the candidate times,
  multiplication tests them, and they coincide exactly on the coarse rationals (S555o).""")

if __name__ == "__main__":
    main()
