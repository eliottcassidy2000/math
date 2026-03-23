#!/usr/bin/env python3
"""
meta_graph_n78_s192.py -- opus-2026-03-22-S192

CHIPPING AWAY AT G_7 AND G_8: ANALYTICAL RESULTS

We can't enumerate all 2^21 or 2^28 tournaments,
but we CAN compute many properties analytically:

1. Burnside sums (V, SC, T, T_anti) — EXACT to any n
2. Fiber fraction — EXACT
3. Edge estimates — from T/2 and V*m*(1-f)/2
4. Degree predictions — from spectral gap
5. Score class structure — from the partition lattice
6. Specific tournament families (transitive, Paley, circulant)

ALSO: properties that DON'T need enumeration:
- Width = C(n-2, floor((n-2)/2))
- Diameter = n-2 (conjectured)
- Spectral gap = 2/n (conjectured)
- Sources = 1 (transitive class)
- Down edges = 0 (DAG property)

Author: opus-2026-03-22-S192
"""

import sys
from math import comb, factorial, gcd, sqrt, pi, log
from collections import Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def gen_partitions(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0: yield []; return
    for first in range(min(n, max_part), 0, -1):
        for rest in gen_partitions(n - first, first):
            yield [first] + rest

def cycle_type_count(n, ct):
    c = Counter(ct)
    r = factorial(n)
    for l, k in c.items(): r //= (l**k) * factorial(k)
    return r

def fix_tournaments(ct):
    for c in ct:
        if c % 2 == 0: return 0
    exp = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)):
            exp += gcd(ct[i], ct[j])
    return 2**exp

def fix_anti(ct, n):
    f = ct.count(1)
    if f >= 2: return 0
    perm = [0]*n; pos = 0
    for c in ct:
        for i in range(c): perm[pos+i] = pos + (i+1)%c
        pos += c
    visited = set(); free = 0
    for i in range(n):
        for j in range(n):
            if i==j or (i,j) in visited: continue
            orbit = []; a,b = i,j
            while (a,b) not in visited:
                visited.add((a,b)); orbit.append((a,b))
                a,b = perm[a],perm[b]
            L = len(orbit); orbit_set = set(orbit)
            rev = (j,i)
            if rev in orbit_set:
                k = orbit.index(rev)
                if k%2==1: free += 1
                else: return 0
            else:
                if rev in visited: pass
                else:
                    if L%2==0: free += 1
                    else: return 0
    return 2**free

def fix_arcs(ct):
    f = ct.count(1)
    t = ct.count(2)
    return comb(f,2) + t

def fiber_fraction(n):
    k = n - 2
    return Fraction(comb(2*k, k), 4**k)

# ============================================================
# PART 1: ANALYTICAL PROPERTIES OF G_n FOR n=3..13
# ============================================================

banner("PART 1: COMPLETE ANALYTICAL TABLE n=3..13")

print("  All computed from cycle-index algebra (no tournament enumeration).\n")

E_known = {3:1, 4:5, 5:30, 6:290, 7:4086}

print("  %-3s | %-10s | %-8s | %-12s | %-8s | %-8s | %-12s | %-8s | %-8s" %
      ("n", "V", "SC", "T_orb", "E_exact", "E~T/2", "V*m*(1-f)/2", "Width", "f(n)"))
print("  " + "-" * 100)

for n in range(3, 14):
    nf = factorial(n); m = comb(n,2)
    V_s=0; SC_s=0; T_s=0
    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        mult = cycle_type_count(n, ct)
        V_s += mult * fix_tournaments(ct)
        SC_s += mult * fix_anti(ct, n)
        T_s += mult * fix_tournaments(ct) * fix_arcs(ct)
    V = V_s // nf; SC = SC_s // nf; T = T_s // nf

    f = fiber_fraction(n)
    approx = float(Fraction(V*m, 2) * (1-f))
    width = comb(n-2, (n-2)//2)
    E_exact = E_known.get(n, "")
    E_T2 = T // 2

    print("  %-3d | %-10d | %-8d | %-12d | %-8s | %-8d | %-12.0f | %-8d | %-8s" %
          (n, V, SC, T, E_exact, E_T2, approx, width, str(f)[:8]))

# ============================================================
# PART 2: n=7 DEEP DIVE
# ============================================================

banner("PART 2: G_7 DEEP DIVE")

n = 7; m = comb(n,2); nf = factorial(n)
V_s=0; SC_s=0; T_s=0; TA_s=0
for ct_list in gen_partitions(n):
    ct = list(ct_list)
    mult = cycle_type_count(n, ct)
    ft = fix_tournaments(ct)
    fa = fix_anti(ct, n)
    farcs = fix_arcs(ct)
    V_s += mult*ft; SC_s += mult*fa
    T_s += mult*ft*farcs; TA_s += mult*fa*farcs

V=V_s//nf; SC=SC_s//nf; T=T_s//nf; TA=TA_s//nf

print("  V(G_7) = %d = A000568(7)" % V)
print("  SC(7) = %d" % SC)
print("  NSC(7) = %d, NSC pairs = %d" % (V-SC, (V-SC)//2))
print("  V_merged(7) = %d" % ((V+SC)//2))
print("  T_orb(7) = %d" % T)
print("  T_anti(7) = %d" % TA)
print("  T_merged(7) = %d" % ((T+TA)//2))
print("  m = C(7,2) = %d arcs" % m)
print("  f(7) = %s = %.6f" % (fiber_fraction(n), float(fiber_fraction(n))))
print()

print("  KNOWN FROM COMPUTATION (opus S212, kind-pasteur S20bz):")
print("  E(G_7) = 4086")
print("  avg_deg = 2*4086/456 = %.4f" % (2*4086/456))
print("  avg_deg/m = %.4f" % (2*4086/(456*21)))
print("  T/(2E) = %d/%d = %.4f" % (T, 2*4086, T/(2*4086)))
print()

# Predictions
print("  PREDICTIONS FOR G_7:")
print("  Width = C(5, 2) = %d" % comb(5, 2))
print("  Diameter = n-2 = %d (CONJECTURE)" % (n-2))
print("  Spectral gap = 2/n = %.4f (CONJECTURE)" % (2.0/n))
print("  Sources = 1 (transitive class)")
print("  Sinks = 1 (regular class, CONJECTURE)")
print()

# Score sequences at n=7
# A000571(7) = 22 distinct score sequences
print("  Score sequences: A000571(7) = 22")
print("  Classes per score (avg): %d/22 = %.1f" % (V, V/22))
print("  Max splitting score: the PoS (1,3,3,3,3,3,6)-like class")
print()

# Specific families at n=7
print("  SPECIFIC TOURNAMENTS AT n=7:")
print("  Transitive: |Aut|=1, H=1, 1 tiling, degree=m-neutral=%d-??" % m)
print("  Paley T_7: |Aut|=%d, H=189, %d tilings" % (7*6//2, 189//(7*6//2)))
print("  Circulant C_7^{1,2,3}: |Aut|>=7, H=175, %d tilings" % (175//7))

# ============================================================
# PART 3: n=8 ANALYTICAL PROFILE
# ============================================================

banner("PART 3: G_8 ANALYTICAL PROFILE")

n = 8; m = comb(n,2); nf = factorial(n)
V_s=0; SC_s=0; T_s=0; TA_s=0
for ct_list in gen_partitions(n):
    ct = list(ct_list)
    mult = cycle_type_count(n, ct)
    ft = fix_tournaments(ct)
    fa = fix_anti(ct, n)
    farcs = fix_arcs(ct)
    V_s += mult*ft; SC_s += mult*fa
    T_s += mult*ft*farcs; TA_s += mult*fa*farcs

V=V_s//nf; SC=SC_s//nf; T=T_s//nf; TA=TA_s//nf

print("  V(G_8) = %d = A000568(8)" % V)
print("  SC(8) = %d" % SC)
print("  NSC(8) = %d, NSC pairs = %d" % (V-SC, (V-SC)//2))
print("  V_merged(8) = %d" % ((V+SC)//2))
print("  T_orb(8) = %d" % T)
print("  T_anti(8) = %d" % TA)
print("  T_merged(8) = %d" % ((T+TA)//2))
print("  m = C(8,2) = %d arcs" % m)
print("  f(8) = %s = %.6f" % (fiber_fraction(n), float(fiber_fraction(n))))
print()

print("  EDGE PREDICTIONS FOR G_8:")
# Method 1: T/2
E_T2 = T // 2
print("  Method 1 (T/2): %d" % E_T2)
# Method 2: V*m*(1-f)/2
f8 = fiber_fraction(8)
E_approx = float(Fraction(V*m, 2) * (1-f8))
print("  Method 2 (V*m*(1-f)/2): %.0f" % E_approx)
# Method 3: T/(2 * ratio) where ratio ~ 1.07 at n=7
E_corrected = T / (2 * 1.05)  # ratio at n=8 estimated
print("  Method 3 (T/2 / 1.05): %.0f" % E_corrected)
print("  PREDICTED RANGE: %.0f - %.0f" % (E_approx, E_T2))
print()

print("  OTHER G_8 PROPERTIES:")
print("  Width = C(6, 3) = %d" % comb(6, 3))
print("  Diameter = n-2 = %d (CONJECTURE)" % (n-2))
print("  Spectral gap = 2/n = %.4f (CONJECTURE)" % (2.0/n))
print("  Score sequences: A000571(8) = 59")
print("  Classes per score (avg): %d/59 = %.1f" % (V, V/59))
print("  SC fraction: %.4f%%" % (100*SC/V))
print()

print("  EVEN-n SPECIFICS (n=8):")
print("  Blueself EXISTS at even n (THM-023)")
print("  The blue skeleton is NOT bipartite at even n")
print("  The transfer matrix tr(M) = 0 at even n")
print("  H_max at n=8: conjectured from regular tournaments")

# ============================================================
# PART 4: n=9 AND BEYOND
# ============================================================

banner("PART 4: n=9+ ANALYTICAL PROFILE")

for n in [9, 10, 11]:
    m = comb(n,2); nf = factorial(n)
    V_s=0; SC_s=0; T_s=0
    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        mult = cycle_type_count(n, ct)
        V_s += mult * fix_tournaments(ct)
        SC_s += mult * fix_anti(ct, n)
        T_s += mult * fix_tournaments(ct) * fix_arcs(ct)
    V = V_s//nf; SC = SC_s//nf; T = T_s//nf
    f = fiber_fraction(n)
    width = comb(n-2, (n-2)//2)

    print("  n=%d:" % n)
    print("    V = %d" % V)
    print("    SC = %d (%.4f%% of V)" % (SC, 100*SC/V))
    print("    T_orb = %d" % T)
    print("    E ~ T/2 = %d" % (T//2))
    print("    m = %d, f = %.6f" % (m, float(f)))
    print("    Width = %d" % width)
    print("    avg_deg ~ m*(1-1/(n-2)^2.4) = %.1f" % (m * (1 - 1/(n-2)**2.4)))
    print()

# ============================================================
# PART 5: GROWTH RATES AND RATIOS
# ============================================================

banner("PART 5: GROWTH RATES")

print("  V(n)/V(n-1) = descent ratio:\n")
prev_V = None
for n in range(3, 14):
    nf = factorial(n)
    V_s = 0
    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        V_s += cycle_type_count(n, ct) * fix_tournaments(ct)
    V = V_s // nf
    if prev_V:
        ratio = V / prev_V
        print("  V(%d)/V(%d) = %d/%d = %.4f" % (n, n-1, V, prev_V, ratio))
    prev_V = V

print("\n  T(n)/T(n-1) = arc-orbit growth ratio:\n")
prev_T = None
for n in range(3, 14):
    nf = factorial(n)
    T_s = 0
    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        T_s += cycle_type_count(n, ct) * fix_tournaments(ct) * fix_arcs(ct)
    T = T_s // nf
    if prev_T:
        ratio = T / prev_T
        print("  T(%d)/T(%d) = %.4f" % (n, n-1, ratio))
    prev_T = T

print("\n  T(n)/(V(n)*m(n)) = avg orbits per class per arc:\n")
for n in range(3, 14):
    nf = factorial(n); m = comb(n,2)
    V_s = 0; T_s = 0
    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        mult = cycle_type_count(n, ct)
        ft = fix_tournaments(ct)
        V_s += mult * ft
        T_s += mult * ft * fix_arcs(ct)
    V = V_s // nf; T = T_s // nf
    ratio = Fraction(T, V * m)
    print("  n=%2d: T/(V*m) = %s = %.6f" % (n, ratio, float(ratio)))

# ============================================================
# PART 6: THE |Aut| DISTRIBUTION AT n=7,8
# ============================================================

banner("PART 6: AUTOMORPHISM DISTRIBUTION PREDICTIONS")

print("""
  From Erdos-Renyi: P(|Aut|>1) ~ n * 2^{-(n-1)}.

  At n=7: P(|Aut|>1) ~ 7 * 2^{-6} = 7/64 = 10.9%
  At n=8: P(|Aut|>1) ~ 8 * 2^{-7} = 8/128 = 6.25%
  At n=9: P(|Aut|>1) ~ 9 * 2^{-8} = 9/256 = 3.52%

  Fraction of iso CLASSES with |Aut|>1 is DIFFERENT from
  fraction of TOURNAMENTS with |Aut|>1. Classes with |Aut|>1
  have FEWER labeled tournaments (size = n!/|Aut| < n!).

  From our data:
  n=5: 5/12 = 42% of classes have |Aut|>1
  n=6: should be ~30% of classes
  n=7: should be ~20% of classes
  n=8: should be ~15% of classes

  The |Aut| values that appear:
  n=7: |Aut| in {1, 3, 5, 7, 21, ...}
    7 divides n, so cyclic Z/7Z is possible
    21 = 7*3 for the Paley tournament (|Aut| = p(p-1)/2)
  n=8: |Aut| in {1, 2, 3, 4, 6, 8, ...}
    Even values appear at even n (from transposition automorphisms)
    8 for cyclic Z/8Z
""")

# ============================================================
# PART 7: THE EDGE FORMULA CONVERGENCE
# ============================================================

banner("PART 7: T/(2E) CONVERGENCE ANALYSIS")

E_data = {3:1, 4:5, 5:30, 6:290, 7:4086}

print("  The ratio T/(2E) converges to 1:\n")
for n in sorted(E_data.keys()):
    nf = factorial(n); m = comb(n,2)
    T_s = 0
    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        T_s += cycle_type_count(n, ct) * fix_tournaments(ct) * fix_arcs(ct)
    T = T_s // nf
    E = E_data[n]
    ratio = T / (2*E)
    deg = T - 2*E  # degeneracy
    print("  n=%d: T=%d, 2E=%d, T/(2E)=%.4f, degeneracy=%d, deg/T=%.4f" %
          (n, T, 2*E, ratio, deg, deg/T))

print("""
  The degeneracy fraction deg/T:
    n=3: 0.5000
    n=4: 0.3750
    n=5: 0.3182
    n=6: 0.1761
    n=7: 0.0828

  This decreases roughly as 1/n:
""")

for n in sorted(E_data.keys()):
    nf = factorial(n); m = comb(n,2)
    T_s = 0
    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        T_s += cycle_type_count(n, ct) * fix_tournaments(ct) * fix_arcs(ct)
    T = T_s // nf
    E = E_data[n]
    deg_frac = (T - 2*E) / T
    print("  n=%d: deg/T = %.4f, 1/n = %.4f, deg/T * n = %.4f" %
          (n, deg_frac, 1/n, deg_frac * n))

print("""
  deg/T * n sequence: 1.5, 1.5, 1.59, 1.06, 0.58
  NOT simply 1/n. But the product DECREASES, suggesting
  deg/T ~ C/n^alpha with alpha > 1.

  If deg/T ~ C/n^1.5:
  n=5: C * 5^{-1.5} = 0.318 -> C = 3.56
  n=7: C * 7^{-1.5} = 3.56/18.5 = 0.192 (actual: 0.0828)
  Does NOT fit. The decay is FASTER than n^{-1.5}.

  PREDICTION for n=8:
  If deg/T = 0.08 * (7/8)^2 ~ 0.06:
  T(8) = 188288
  2*E(8) = T * (1 - 0.06) = 188288 * 0.94 = 176991
  E(8) ~ 88,500

  If deg/T = 0.05 (faster decay):
  E(8) ~ 188288 * 0.95 / 2 = 89,437
""")

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: WHAT WE KNOW ABOUT G_7 AND G_8")

print("""
  G_7 (COMPREHENSIVE):
    V = 456, E = 4086, SC = 88, NSC = 368
    m = 21, avg_deg = 17.92, avg_deg/m = 0.853
    T_orb = 8912, T/(2E) = 1.091
    Width = 10, Diameter = 5 (conjectured)
    Spectral gap = 2/7 (conjectured)
    f(7) = 63/256 = 0.246
    Blue fraction ~ 77%
    SC fraction = 19.3% (declining fast)
    |Aut| values: {1, 3, 5, 7, 21, ...}
    Score sequences: 22

  G_8 (PREDICTED):
    V = 6880, SC = 176, NSC = 6704
    m = 28, T_orb = 188288
    E ~ 88,000 (range: 75,000 - 95,000)
    avg_deg ~ 25.6 (from E ~ 88000)
    avg_deg/m ~ 0.91
    Width = 20
    Diameter = 6 (conjectured)
    Spectral gap = 1/4 (conjectured)
    f(8) = 231/1024 = 0.226
    SC fraction = 2.56% (even n: very few SC classes)
    Blueself EXISTS (even n)
    Score sequences: 59
    Classes per score: ~117 (much more splitting)

  ASYMPTOTIC (large n):
    V ~ 2^{C(n,2)} / n!
    SC / V -> 0
    T / (V * m) -> 1
    E ~ T/2 ~ V * m / 2
    avg_deg / m -> 1
    G_n becomes asymptotically complete
    f(n) -> 0 (fibers thin)
    Width ~ 2^{n-2} / sqrt(pi*(n-2)/2)
""")

print("Done. opus-2026-03-22-S192")
