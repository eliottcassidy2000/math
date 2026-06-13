#!/usr/bin/env python3
"""
unit_distance_product_cap_s1.py
monad-explorer-2026-06-07-S1  (deep-research)

THESIS.  Average degree is ADDITIVE under the Cartesian/Minkowski product of
unit-distance graphs.  Consequently the *product* (Erdos-Minkowski) family
CANNOT beat 3N below a hard threshold, and in particular cannot beat it
anywhere inside the crossover window N* in [25,28] (THM-431).  Therefore N*
is necessarily a NON-PRODUCT (irreducible / Moser-lattice) phenomenon, and the
clean tie at N = 27 = 3^3 is exactly the Cartesian cube K3 [] K3 [] K3.

Key elementary facts (all exact integer arithmetic here):
  For graphs G, H the Cartesian product G[]H (= generic-angle Minkowski sum,
  the standard Erdos unit-distance product construction) has
      n(G[]H) = n(G) n(H)
      e(G[]H) = e(G) n(H) + n(G) e(H)
  hence
      avgdeg(G[]H) = 2 e/n = 2e(G)/n(G) + 2e(H)/n(H) = avgdeg(G) + avgdeg(H).
  Since A[]B is a unit-distance graph whenever A, B are (extra coincidences
  only RAISE the count), this yields the rigorous LOWER bound
      u(n(G) n(H)) >= e(G) n(H) + n(G) e(H).

We compute, over ALL factorizations N = prod n_i (n_i >= 2), the best product
construction edge count bestE[N], using the PROVEN small-n optimum u(n) (AMP
arXiv:2412.11914, exact for n<=21; proven lower bounds 22..30 for comparison),
and compare to 3N.
"""

# Proven Erdos unit-distance maxima u(n), n=0..21 (AMP arXiv:2412.11914, Thm 1; OEIS A186705)
U_EXACT = {0:0,1:0,2:1,3:3,4:5,5:7,6:9,7:12,8:14,9:18,10:20,11:23,
           12:27,13:30,14:33,15:37,16:41,17:43,18:46,19:50,20:54,21:57}
# Proven realizable lower bounds (Schade/Engel "Moser lattice") for 22..30 (AMP Table 1)
U_LOWER = dict(U_EXACT)
U_LOWER.update({22:60,23:64,24:72,25:75,26:78,27:81,28:85,29:89,30:93})
# NB: 24,25,26,27 lower bounds equal/exceed 3n at some n -> Moser ties/beats; these are
# the NON-product values.  We test the PRODUCT family against them.

NMAX = 42

def divisors(n):
    return [d for d in range(2, n) if n % d == 0]

# ---- product-construction DP --------------------------------------------------
# bestE[N] = max edge count of a unit-distance graph on N points obtainable as an
# iterated Cartesian product of OPTIMAL atomic factors u(n_i).
# Recurrence: bestE[N] = max( u(N) [atomic, if known],
#                             max_{d|N,2<=d<N} u(d)*(N/d) + bestE[N/d]*d )
bestE = {1:0}
best_factorization = {1:()}
for N in range(2, NMAX+1):
    cand = []
    # atomic
    if N in U_LOWER:
        cand.append((U_LOWER[N], (N,)))         # use the FULL (Moser-inclusive) optimum atomically
    # also the product-only atomic (no non-product help): handled by factorization below
    for d in divisors(N):
        b = N // d
        if b in bestE and d in U_LOWER:
            e = U_LOWER[d]*b + bestE[b]*d
            cand.append((e, (d,)+best_factorization[b]))
    if not cand:                                  # prime N outside the known table
        cand.append((-1, (N,)))
    e_best, fac = max(cand)
    bestE[N] = e_best
    best_factorization[N] = fac

# Separate: PRODUCT-ONLY cap = best you can do using products of factors each of size <= 21
# (i.e. NOT allowing the unknown/Moser atomic value for composite N itself, only as factors).
# This is the honest "is the crossover reachable by products of small optima?" question.
bestEprod = {1:0}
facprod = {1:()}
for N in range(2, NMAX+1):
    cand = []
    if N <= 21:
        cand.append((U_EXACT[N], (N,)))          # atomic, proven exact
    for d in divisors(N):
        b = N//d
        if b in bestEprod and d <= 21:
            cand.append((U_EXACT[d]*b + bestEprod[b]*d, (d,)+facprod[b]))
    if not cand:                                  # prime N > 21: NO product construction exists
        cand.append((-1, (N,)))
    e_best, fac = max(cand)
    bestEprod[N] = e_best
    facprod[N] = fac

# ---- report -------------------------------------------------------------------
def fr(n,d):
    from math import gcd
    g=gcd(n,d); return f"{n//g}/{d//g}"

print("="*78)
print("AVERAGE DEGREE IS ADDITIVE UNDER [] :  avgdeg(G[]H) = avgdeg(G)+avgdeg(H)")
print("Atomic avg degrees a(n) = 2u(n)/n  (proven u(n), AMP):")
for n in range(3,17):
    print(f"   a({n:2d}) = 2*{U_EXACT[n]}/{n} = {fr(2*U_EXACT[n],n):>6}  = {2*U_EXACT[n]/n:.4f}")
print()
print("="*78)
print("PRODUCT-ONLY construction cap vs 3N  (factors are proven small optima, n<=21):")
print(f"{'N':>3} {'3N':>4} {'prod_cap':>9} {'def=cap-3N':>11}  factorization (best)")
first_prod_beat = None
prod_ties = []
for N in range(3, NMAX+1):
    cap = bestEprod[N]
    if cap < 0:
        print(f"{N:>3} {3*N:>4} {'?(prime>21)':>9}")
        continue
    d = cap - 3*N
    fac = facprod[N]
    facs = "[]".join(f"u{x}" for x in sorted(fac, reverse=True))
    tag = ""
    if d > 0 and first_prod_beat is None:
        first_prod_beat = N; tag = "  <== FIRST PRODUCT BEAT"
    if d == 0:
        prod_ties.append(N); tag = "  (TIE)"
    print(f"{N:>3} {3*N:>4} {cap:>9} {d:>+11}  {facs}{tag}")
print()
print(f"First N where a PRODUCT construction strictly beats 3N: {first_prod_beat}")
print(f"Product TIES (cap == 3N) at N = {prod_ties}")
print()
print("="*78)
print("THE N=27 = 3^3 TIE, explicitly:")
for fac,name in [((3,3,3),"K3 [] K3 [] K3  (Cartesian CUBE of the triangle)"),
                 ((9,3),"G9 [] K3        (optimal 9-pt graph x triangle)")]:
    # edges via additivity:  avgdeg = sum a(n_i); edges = N/2 * avgdeg
    N=1; ad_num=0
    # compute edges exactly by folding
    e=0; n=1
    for x in fac:
        e = U_EXACT[x]*n + e*x
        n = n*x
    ad = f"{fr(2*e,n)}"
    print(f"   {name:<48}  N={n}, edges={e}={'3*27' if e==81 else e}, avgdeg={ad}")
print()
print("Both hit avgdeg = 6 EXACTLY (=3+3 or 2+2+2) -> u(27) >= 81 by products,")
print("and NO product can exceed it (additivity caps avgdeg at 2+2+2=6 for 3^3).")
print()
print("="*78)
print("CONSEQUENCE for N* in [25,28] (THM-431):")
for N in [25,26,27,28]:
    print(f"   N={N}: product cap = {bestEprod[N]} {'<=' if bestEprod[N]<=3*N else '>'} 3N={3*N}"
          f"   | proven u(N) lower bound (Moser, NON-product) = {U_LOWER[N]}")
print()
print("Since the product cap <= 3N for every N <= 31 (tie at most), ANY graph with")
print("u(N) > 3N at N <= 31 -- in particular at the crossover N* in [25,28] -- is")
print("NECESSARILY NON-PRODUCT (irreducible).  N* is a Moser-lattice phenomenon;")
print("the Erdos product family alone first beats 3N only at N =", first_prod_beat, ".")
print()
print("TANGENCY + TAIL:")
print("  P(N)-3N grazes 0 (tangent from below) at N=27,30 (the ties), dips to -1 at")
print("  N=28 between them, then crosses positive at 32 and stays positive:")
tail_ok = all((bestEprod[N] < 0) or (bestEprod[N] > 3*N) for N in range(32, NMAX+1))
print("  every composite N in [32,%d] strictly beats 3N: %s" % (NMAX, tail_ok))
print("  -> the crossover N* in [25,28] sits exactly in the grazing region (tie@27,")
print("     dip@28) where the product family touches but cannot cross 3N.")
