"""
Search for nullity formula from known data.

P_7 (m=3, p=7):
  d:  0  1  2  3  4   5   6
  n:  0  0  0  0  9  33  42

P_11 (m=5, p=11):
  d:  0  1  2  3   4    5     6
  n:  0  0  0  0  25  250  1620

Also: Omega values:
P_7:  [1, 3, 6, 9, 9, 6, 3]
P_11: [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]

Betti: β_d = Omega_d - rank(∂_d on Omega) - rank(∂_{d+1} on Omega)
For the TOTAL complex (summing over eigenspaces).

opus-2026-03-13-S71b
"""
from math import comb

# Known data
data = {
    7: {  # p=7, m=3
        'A': [1, 3, 9, 21, 39, 45, 27],
        'junk': [0, 0, 3, 12, 39, 72, 66],
        'rank': [0, 0, 3, 12, 30, 39, 24],
        'nullity': [0, 0, 0, 0, 9, 33, 42],
        'Omega': [1, 3, 6, 9, 9, 6, 3],
    },
    11: {  # p=11, m=5
        'A': [1, 5, 25, 110, 430, 1430, 3970],
        'junk': [0, 0, 5, 40, 250, 1220, 4890],
        'rank': [0, 0, 5, 40, 225, 970, 3270],
        'nullity': [0, 0, 0, 0, 25, 250, 1620],
        'Omega': [1, 5, 20, 70, 205, 460, 700],
    }
}

# Additional known: P_19 d=4 nullity = 81, m=9
# Additional known: P_23 d=4 nullity = 121, m=11

# Look for patterns in nullity
for p, d_data in data.items():
    m = (p - 1) // 2
    print(f"\nP_{p} (m={m}):")
    for d in range(len(d_data['nullity'])):
        n = d_data['nullity'][d]
        if n == 0:
            continue
        # Try various formulas
        candidates = []
        # Check m^k for various k
        for k in range(1, 6):
            if n == m**k:
                candidates.append(f"m^{k}")
        # Check C(d-1,k) * m^j
        for k in range(0, d):
            for j in range(1, 5):
                if comb(d-1, k) * m**j == n:
                    candidates.append(f"C({d-1},{k})*m^{j}")
        # Check C(m,k) patterns
        for k in range(1, m+2):
            if comb(m, k) == n:
                candidates.append(f"C({m},{k})")
            for j in range(1, 5):
                if comb(m, k) * m**j == n:
                    candidates.append(f"C({m},{k})*m^{j}")
        # Check product formulas
        for a in range(1, 20):
            for b in range(1, 20):
                if a * m**2 + b * m == n:
                    candidates.append(f"{a}*m² + {b}*m")
                if a * m**2 - b * m == n:
                    candidates.append(f"{a}*m² - {b}*m")
                if a * m**3 + b * m**2 == n:
                    candidates.append(f"{a}*m³ + {b}*m²")
                if a * m**3 - b * m**2 == n:
                    candidates.append(f"{a}*m³ - {b}*m²")
        # Check rank formula: rank = A - Omega, and if rank factors nicely
        r = d_data['rank'][d]
        a = d_data['A'][d]
        omega = d_data['Omega'][d]
        print(f"  d={d}: nullity={n}, rank={r}, A={a}, Omega={omega}")
        if candidates:
            print(f"    Matches: {', '.join(candidates[:10])}")

# Cross-prime comparison at fixed d
print("\n\nCross-prime comparison:")
for d in [4, 5, 6]:
    print(f"\n  d={d}:")
    for p in [7, 11]:
        m = (p - 1) // 2
        if d < len(data[p]['nullity']):
            n = data[p]['nullity'][d]
            print(f"    P_{p} (m={m}): nullity={n}")

# Try: is nullity(d,m) a polynomial in m?
# d=4: n(3)=9, n(5)=25, n(9)=81, n(11)=121 → n = m² ✓
# d=5: n(3)=33, n(5)=250
#   Try n = a*m⁴ + b*m³ + c*m² + d*m
#   33 = 81a + 27b + 9c + 3d
#   250 = 625a + 125b + 25c + 5d
# Subtract: 217 = 544a + 98b + 16c + 2d
# Need more data points. But let me try some specific forms:
print("\n\nTrying polynomial forms for d=5 nullity:")
m_vals = [3, 5]
n_vals = [33, 250]

# Try n = a*m³ + b*m²
# 27a + 9b = 33, 125a + 25b = 250
# From first: 3a + b = 33/9 = 11/3 → NOT integer. So not this form.

# Try n = a*m⁴ + b*m²
# 81a + 9b = 33, 625a + 25b = 250
# From first: 9a + b = 33/9 → not integer.

# Try n = m * f(m)
# 33/3 = 11, 250/5 = 50
# f(3) = 11, f(5) = 50
# Try f = a*m² + b*m + c:
# 9a + 3b + c = 11
# 25a + 5b + c = 50
# Subtract: 16a + 2b = 39. Not integer coefficients.

# Try f = a*m³ + b*m:
# 27a + 3b = 11 → 9a + b = 11/3. Not integer.

# Try n = m² * g(m):
# g(3) = 33/9 = 11/3. Not integer!

# So it's not a simple polynomial in m. Let me try with p.
# d=5: n(7)=33, n(11)=250
# n = a*p² + b*p:  49a + 7b = 33, 121a + 11b = 250
# From first: 7a + b = 33/7 → not integer.

# Try n = (something involving both m and d)
# What about: sum over k of nullity contribution from each same-sign pair?
# At d=5: same-sign pairs (1,3) and (2,4).
# If (1,3) coupling gives nullity = m² (as at d=4)
# and (2,4) coupling gives additional nullity
# and maybe (1,3)+(2,4)+(1,2,3,4) 4-way gives more...

# Let me look at the OMEGA pattern more carefully
print("\n\nOmega values:")
for p in [7, 11]:
    m = (p - 1) // 2
    print(f"P_{p}: Omega = {data[p]['Omega']}")
    # Check: Omega_d = C(p-1, d) - something?
    print(f"  C(p-1,d) = {[comb(p-1, d) for d in range(len(data[p]['Omega']))]}")
    # Ratio Omega_d / C(p-1,d):
    ratios = [data[p]['Omega'][d] / comb(p-1, d) if comb(p-1, d) > 0 else 0
              for d in range(len(data[p]['Omega']))]
    print(f"  Omega/C(p-1,d) = {[f'{r:.4f}' for r in ratios]}")

# Check: Omega_d = A_d * product
print("\n\nOmega/A ratios:")
for p in [7, 11]:
    m = (p - 1) // 2
    ratios = [data[p]['Omega'][d] / data[p]['A'][d] if data[p]['A'][d] > 0 else 0
              for d in range(len(data[p]['Omega']))]
    print(f"P_{p}: {[f'{r:.4f}' for r in ratios]}")
    # Check if these are (m-d+1)/m or similar
    for d in range(len(ratios)):
        r = ratios[d]
        # Check if r = 1 - d/p
        expected = 1 - d/p if d <= len(ratios) else 0
        print(f"  d={d}: Omega/A = {r:.4f}, 1-d/p = {expected:.4f}, "
              f"1-(d-1)/m = {1-(d-1)/m:.4f}")
