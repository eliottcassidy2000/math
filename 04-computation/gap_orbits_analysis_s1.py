"""
Analyze the gap_orbits sequence: 2, 5, 20, 86, 490, 3703, 47889
Looking for patterns, generating functions, or Burnside-type formulas.

kind-pasteur-2026-03-25-S1
"""
import sys
from math import factorial, comb, gcd
from fractions import Fraction

sys.stdout.reconfigure(line_buffering=True)

# The sequence (n=3..9)
gap = {3: 2, 4: 5, 5: 20, 6: 86, 7: 490, 8: 3703, 9: 47889}

# Related sequences
A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536}  # tournament iso classes
E_Gn = {3: 1, 4: 5, 5: 30, 6: 290, 7: 4086, 8: 91161, 9: 3380751}  # metagraph edges
m_tiles = {n: comb(n-1, 2) for n in range(3, 10)}
T_tilings = {n: 2**m_tiles[n] for n in range(3, 10)}  # total tilings

# D(n) from OPEN-Q-039 (n=3..10)
D = {3: 2, 4: 6, 5: 16, 6: 60, 7: 328, 8: 3160, 9: 54928, 10: 1722992}

print("=" * 60)
print("GAP ORBITS ANALYSIS")
print("=" * 60)

print("\n1. Basic data:")
print(f"{'n':>3} {'gap':>8} {'V(Gn)':>8} {'E(Gn)':>10} {'m':>4} {'2^m':>12} {'D(n)':>8}")
for n in range(3, 10):
    print(f"{n:>3} {gap[n]:>8} {A000568[n]:>8} {E_Gn[n]:>10} {m_tiles[n]:>4} {T_tilings[n]:>12} {D.get(n, '?'):>8}")

print("\n2. Ratios gap(n+1)/gap(n):")
for n in range(3, 9):
    r = gap[n+1] / gap[n]
    print(f"  gap({n+1})/gap({n}) = {r:.4f}")

print("\n3. gap(n) / (n-2)!:")
for n in range(3, 10):
    r = Fraction(gap[n], factorial(n-2))
    print(f"  gap({n}) / {n-2}! = {float(r):.6f} = {r}")

print("\n4. gap(n) / A000568(n):")
for n in range(3, 10):
    r = Fraction(gap[n], A000568[n])
    print(f"  gap({n}) / V({n}) = {float(r):.4f} = {r}")

print("\n5. gap(n) / D(n) (D = even-cycle Burnside formula):")
for n in range(3, 10):
    if n in D:
        r = Fraction(gap[n], D[n])
        print(f"  gap({n}) / D({n}) = {float(r):.4f} = {r}")

print("\n6. Check recurrence a(n) = c1*a(n-1) + c2*a(n-2):")
# Try to solve: a(5) = c1*a(4) + c2*a(3) and a(6) = c1*a(5) + c2*a(4)
# 20 = 5*c1 + 2*c2
# 86 = 20*c1 + 5*c2
# From first: c2 = (20 - 5*c1)/2
# Sub into second: 86 = 20*c1 + 5*(20 - 5*c1)/2 = 20*c1 + 50 - 12.5*c1 = 7.5*c1 + 50
# c1 = 36/7.5 = 4.8
# c2 = (20 - 24)/2 = -2
c1 = Fraction(36, 75) * 10  # = 48/10 = 24/5
c2 = Fraction(20 - 5 * c1, 2)
print(f"  Trying c1 = {c1}, c2 = {c2}")
for n in range(5, 10):
    pred = c1 * gap[n-1] + c2 * gap[n-2]
    print(f"  a({n}) = {c1}*{gap[n-1]} + {c2}*{gap[n-2]} = {pred}, actual = {gap[n]}, {'OK' if pred == gap[n] else 'FAIL'}")

print("\n7. Try a(n) = f(n)*a(n-1) + g(n)*a(n-2):")
# For each n >= 5, solve: a(n) = f(n)*a(n-1) + g(n)*a(n-2)
# We have 5 equations and need f(n), g(n) as functions of n
for n in range(5, 10):
    # One equation, two unknowns. Try g(n) = -1:
    f_val = Fraction(gap[n] + gap[n-2], gap[n-1])
    print(f"  n={n}: if g={-1}, f = (a(n)+a(n-2))/a(n-1) = {float(f_val):.4f}")

print("\n8. Try a(n) = (n-1)*a(n-1) - something:")
for n in range(4, 10):
    residual = gap[n] - (n-1)*gap[n-1]
    print(f"  a({n}) - {n-1}*a({n-1}) = {residual}")

print("\n9. Try a(n) = n*a(n-1) - something:")
for n in range(4, 10):
    residual = gap[n] - n*gap[n-1]
    print(f"  a({n}) - {n}*a({n-1}) = {residual}")

print("\n10. Check: is edge_orbits = T_n/2 + (n-2)! correct?")
for n in range(3, 10):
    edge_orb = T_tilings[n] // 2 + factorial(n-2)
    gap_computed = edge_orb - E_Gn[n]
    print(f"  n={n}: edge_orbits = {T_tilings[n]//2}+{factorial(n-2)} = {edge_orb}, gap = {edge_orb}-{E_Gn[n]} = {gap_computed}, expected gap = {gap[n]}, {'OK' if gap_computed == gap[n] else 'FAIL'}")

print("\n11. Factorizations:")
for n in range(3, 10):
    g = gap[n]
    factors = []
    temp = g
    for p in range(2, temp + 1):
        while temp % p == 0:
            factors.append(p)
            temp //= p
        if temp == 1:
            break
    print(f"  gap({n}) = {g} = {'*'.join(map(str, factors))}")

print("\n12. EGF coefficients a(n)/n!:")
for n in range(3, 10):
    r = Fraction(gap[n], factorial(n))
    print(f"  gap({n})/{n}! = {float(r):.8f}")

print("\n13. Check a(n) = sum over cycle types in S_n of something:")
# Burnside: if gap_orbits = (1/n!) * sum_{sigma in S_n} Fix_gap(sigma),
# then n! * gap(n) should be an integer
for n in range(3, 10):
    print(f"  {n}! * gap({n}) = {factorial(n) * gap[n]}")
