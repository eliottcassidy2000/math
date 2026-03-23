#!/usr/bin/env python3
"""
three_burnsides_s190.py -- opus-2026-03-22-S190

THE THREE BURNSIDES OF TOURNAMENT THEORY

Three different Burnside-type formulas control G_n:

1. CLASSICAL BURNSIDE: V = (1/n!) sum Fix(sigma)
   Fix(sigma) = 2^p(sigma) where p = pair orbits (with parity constraint)
   COUNTS: iso classes = vertices of G_n
   GIVES: A000568(n)

2. ANTI-AUTOMORPHISM BURNSIDE: SC = (1/n!) sum Fix_anti(sigma)
   Fix_anti(sigma) from directed arc orbit analysis
   COUNTS: self-complementary iso classes
   GIVES: 2, 2, 8, 12, 88, 176, 2752, ...

3. EDGE BURNSIDE: T = (1/n!) sum Fix(sigma) * Fix_arcs(sigma)
   T = total arc-orbits (weighted Burnside)
   CONSTRAINS: edge count via E ~ T/2 asymptotically
   GIVES: 4, 16, 88, 704, 8912, 188288, ...

Feng's DUAL BURNSIDE PROCESS (arXiv:2510.25202) unifies these:
  Classical = K = BA kernel (state-to-group)
  Dual = Q = AB kernel (group-to-state)
  Both share nonzero eigenvalues

This script computes all three to n=13 and derives the
combined implications for the edge count.

Author: opus-2026-03-22-S190
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

# Burnside 1: Classical
def fix_tournaments(ct):
    for c in ct:
        if c % 2 == 0: return 0
    exp = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)):
            exp += gcd(ct[i], ct[j])
    return 2**exp

# Burnside 2: Anti-automorphism
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

# Burnside 3: Arc-weighted (for edge count)
def fix_arcs(ct):
    f = ct.count(1)
    t = ct.count(2)
    return comb(f,2) + t

# ============================================================
# COMPUTE ALL THREE TO n=13
# ============================================================

banner("THE THREE BURNSIDE SUMS")

print("  %-3s | %-12s | %-12s | %-12s | %-12s | %-10s | %-10s" %
      ("n", "V (cls)", "SC (anti)", "T (arc-wt)", "T_anti_wt", "V_merged", "T_merged"))
print("  " + "-" * 85)

all_data = {}
for n in range(3, 14):
    nf = factorial(n); m = comb(n,2)
    V_s=0; SC_s=0; T_s=0; TA_s=0

    for ct_list in gen_partitions(n):
        ct = list(ct_list)
        mult = cycle_type_count(n, ct)
        ft = fix_tournaments(ct)
        fa = fix_anti(ct, n)
        farcs = fix_arcs(ct)

        V_s += mult * ft
        SC_s += mult * fa
        T_s += mult * ft * farcs
        TA_s += mult * fa * farcs

    V = V_s // nf; SC = SC_s // nf
    T = T_s // nf; TA = TA_s // nf
    Vm = (V+SC)//2; Tm = (T+TA)//2

    all_data[n] = {'V':V, 'SC':SC, 'T':T, 'TA':TA, 'Vm':Vm, 'Tm':Tm, 'm':m}

    print("  %-3d | %-12d | %-12d | %-12d | %-12d | %-10d | %-10d" %
          (n, V, SC, T, TA, Vm, Tm))

# ============================================================
# THE EDGE COUNT FROM T
# ============================================================

banner("EDGE COUNT ESTIMATES FROM T")

E_known = {3:1, 4:5, 5:30, 6:290, 7:4086}

print("  %-3s | %-8s | %-8s | %-8s | %-8s | %-10s | %-8s" %
      ("n", "T", "T/2", "E_actual", "T/(2E)", "V*m*(1-f)/2", "E_merged"))
print("  " + "-" * 70)

for n in range(3, 9):
    d = all_data[n]; m = d['m']
    T = d['T']; V = d['V']
    E = E_known.get(n, None)

    # Fiber fraction
    k = n-2; f_num=1; f_den=1
    for i in range(k): f_num*=(2*i+1); f_den*=(2*i+2)
    f = Fraction(f_num, f_den)
    approx = float(Fraction(V*m,2)*(1-f))

    # Merged estimate
    Tm = d['Tm']

    if E:
        ratio = T/(2*E)
        print("  %-3d | %-8d | %-8d | %-8d | %-8.4f | %-10.1f | %-8d" %
              (n, T, T//2, E, ratio, approx, Tm))
    else:
        print("  %-3d | %-8d | %-8d | %-8s | %-8s | %-10.1f | %-8d" %
              (n, T, T//2, "?", "?", approx, Tm))

# ============================================================
# THE KEY RELATIONSHIPS
# ============================================================

banner("KEY RELATIONSHIPS BETWEEN THE THREE BURNSIDES")

print("  RELATIONSHIP 1: T_anti / SC -> floor(n/2)\n")
for n in range(3, 14):
    d = all_data[n]
    if d['SC'] > 0:
        ratio = d['TA'] / d['SC']
        print("    n=%2d: T_anti/SC = %d/%d = %.4f, floor(n/2) = %d" %
              (n, d['TA'], d['SC'], ratio, n//2))

print("\n  RELATIONSHIP 2: V_merged / V -> 1/2\n")
for n in range(3, 14):
    d = all_data[n]
    ratio = d['Vm'] / d['V']
    print("    n=%2d: V_merged/V = %.6f" % (n, ratio))

print("\n  RELATIONSHIP 3: SC / V -> 0\n")
for n in range(3, 14):
    d = all_data[n]
    ratio = d['SC'] / d['V']
    print("    n=%2d: SC/V = %.6f" % (n, ratio))

print("\n  RELATIONSHIP 4: T / (V * m) = avg arc-orbits per class per arc\n")
for n in range(3, 10):
    d = all_data[n]
    ratio = Fraction(d['T'], d['V'] * d['m'])
    print("    n=%2d: T/(V*m) = %s = %.6f" % (n, ratio, float(ratio)))

# ============================================================
# THE ORBIFOLD EULER CHARACTERISTIC
# ============================================================

banner("ORBIFOLD EULER CHARACTERISTIC (from kind-pasteur S20ca)")

print("  chi_orb = sum_C 1/|Aut(C)| = V / <|Aut|>\n")
print("  From Burnside: sum 1/|Aut| = 2^m / n! (orbit-stabilizer).\n")

for n in range(3, 10):
    m = comb(n,2)
    chi_orb = Fraction(2**m, factorial(n))
    print("  n=%d: chi_orb = 2^%d/%d! = %s = %.6f" %
          (n, m, n, chi_orb, float(chi_orb)))

# ============================================================
# THE TRIANGLE UNIFICATION (kind-pasteur S20cm)
# ============================================================

banner("THE TRIANGLE UNIFICATION")

print("""
  From kind-pasteur S20cm: EVERYTHING IS THE TRIANGLE.

  THE STAIRCASE delta_{n-2} IS a right isosceles triangle.
  Its three sides encode the three Burnside formulas:

  HYPOTENUSE (anti-diagonal, i+j = const):
    -> Hook lengths 2(k-d)-1 (odd numbers)
    -> Fiber fraction f(n) = (1/2)_k / k! = [x^k](1-x)^{-1/2}
    -> Controls the TOPOLOGY (two-sheeted cover, pi)
    -> BURNSIDE 2 (anti-automorphism) lives here
       (SC classes = tilings on the branch locus)

  VERTICAL LEG (column j = const, same target vertex):
    -> Score contributions
    -> OCR = projection onto score space
    -> Controls the HIERARCHY (transitive at top)
    -> BURNSIDE 1 (classical) lives here
       (vertex count = orbits of tiling space)

  HORIZONTAL LEG (row i = const, same source vertex):
    -> Arc range = "skip distance"
    -> Controls the DYNAMICS (arc reversal = Glauber)
    -> BURNSIDE 3 (edge/arc-weighted) lives here
       (edge count ~ arc-orbit count)

  THE FOUR CONSTANTS:
    sqrt(2) = hypotenuse/leg = the "pseudo-doubling" ratio
    pi = from fiber fractions (Wallis product)
    e = from Stirling/Burnside asymptotics
    gamma = Euler-Mascheroni from harmonic corrections

  FENG'S DUAL BURNSIDE (arXiv:2510.25202):
    Classical K = BA acts on group elements (cycle types)
    Dual Q = AB acts on states (tournaments)
    Both share eigenvalues -> the THREE Burnsides are spectral
    projections of the SAME underlying Markov chain.
""")

# ============================================================
# THE EDGE FORMULA: FINAL STATUS
# ============================================================

banner("THE EDGE FORMULA: DEFINITIVE STATUS")

print("""
  |E(G_n)| = 1, 5, 30, 290, 4086

  THREE BURNSIDE FORMULAS give THREE estimates:

  1. T/2 (from BURNSIDE 3):
     E ~ T/2 where T = sum Fix(sigma)*Fix_arcs(sigma) / n!
     T/(2E) -> 1 as n -> inf. Currently 1.09 at n=7.

  2. V*m*(1-f)/2 (from BURNSIDE 1 + fiber fraction):
     E ~ A000568(n) * C(n,2) * (1 - (1/2)_{n-2}/(n-2)!) / 2
     95%+ for n >= 6. Overshoots then undershoots.

  3. Tm/Em (from BURNSIDE 2 + merged graph):
     E_merged via T_merged = (T + T_anti) / 2
     E_orig = E_merged + collapsed + twin

  THE EXACT FORMULA:
  E(G_n) = (1/2) * sum_C degree(C)
  where degree(C) = #{distinct classes among cross-orbits of Aut(C)}
  This is COMPUTABLE but requires knowing Aut(C) for each class.

  THE ASYMPTOTIC FORMULA:
  E(G_n) ~ A000568(n) * C(n,2) / 2 as n -> inf
  Equivalently: avg_deg -> m = C(n,2)
  The graph becomes ASYMPTOTICALLY COMPLETE.

  THE REMAINING GAP:
  A CLOSED-FORM cycle-index expression that gives E(G_n)
  exactly for all n, without computing individual classes.
  This would require a Burnside formula that accounts for
  the DEGENERACY (multiple arc-orbits hitting the same class).

  CURRENT BEST ESTIMATE for n=8:
  T(8) = 188288, so E(8) ~ 188288/2 * (1/1.07) ~ 88,000
  V*m*(1-f)/2 = 74,592
  Range: 75,000 - 90,000
""")

# ============================================================
# THE FIVE EXACT SEQUENCES
# ============================================================

banner("FIVE EXACT SEQUENCES FROM THE THREE BURNSIDES")

print("  All computable to any n via cycle-index algebra:\n")

print("  1. V(n) = A000568 (tournament iso classes):")
v_seq = [all_data[n]['V'] for n in range(3,14)]
print("     %s\n" % v_seq)

print("  2. SC(n) (self-complementary classes):")
sc_seq = [all_data[n]['SC'] for n in range(3,14)]
print("     %s\n" % sc_seq)

print("  3. T(n) (total arc-orbits):")
t_seq = [all_data[n]['T'] for n in range(3,14)]
print("     %s\n" % t_seq)

print("  4. V_merged(n) = (V + SC) / 2:")
vm_seq = [all_data[n]['Vm'] for n in range(3,14)]
print("     %s\n" % vm_seq)

print("  5. T_merged(n) = (T + T_anti) / 2:")
tm_seq = [all_data[n]['Tm'] for n in range(3,14)]
print("     %s\n" % tm_seq)

print("  Additional derived sequences:")
print("  NSC(n) = V - SC: %s" % [all_data[n]['V']-all_data[n]['SC'] for n in range(3,14)])
print("  T_anti(n): %s" % [all_data[n]['TA'] for n in range(3,14)])

print("\nDone. opus-2026-03-22-S190")
