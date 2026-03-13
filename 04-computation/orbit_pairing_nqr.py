#!/usr/bin/env python3
"""
orbit_pairing_nqr.py — Why do H-values come in pairs of orbits?

At p=17, every H value has exactly 2 QR orbits of 8.
HYPOTHESIS: The two orbits are related by NQR multiplication.

For b ∈ NQR: P_b maps QR orbits to QR orbits, and since
b maps T_Paley → T_Paley^op, we have H(P_b σ) = H(σ).
So P_b pairs orbits at the same H level.

But P_b is NOT in the QR group, so it maps one orbit to a DIFFERENT one.

This should give exactly the pairing we observe.

Author: opus-2026-03-12-S63
"""

import numpy as np
from itertools import product

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

def make_signed_permutation(a, p):
    m = (p-1)//2
    P = np.zeros((m, m))
    for k in range(1, m+1):
        ak = (a * k) % p
        if ak <= m:
            P[ak-1, k-1] = 1
        else:
            P[(p-ak)-1, k-1] = -1
    return P

# p=17
p = 17
m = 8
qr = sorted(a for a in range(1, p) if legendre(a, p) == 1)
nqr = sorted(a for a in range(1, p) if legendre(a, p) == -1)

print(f"p={p}, m={m}")
print(f"QR = {qr}")
print(f"NQR = {nqr}")

# Build QR signed permutations
P_qr = {a: make_signed_permutation(a, p) for a in qr}
P_nqr = {b: make_signed_permutation(b, p) for b in nqr}

# The two orbits at H_max from the previous analysis
orbit_0 = [
    (-1, -1, -1, -1, -1, -1, -1, -1),
    (-1, -1, -1, -1, 1, 1, 1, 1),
    (-1, -1, 1, 1, -1, -1, 1, 1),
    (-1, 1, -1, 1, -1, 1, -1, 1),
    (1, -1, 1, -1, 1, -1, 1, -1),
    (1, 1, -1, -1, 1, 1, -1, -1),
    (1, 1, 1, 1, -1, -1, -1, -1),
    (1, 1, 1, 1, 1, 1, 1, 1),
]

orbit_28 = [
    (-1, -1, 1, 1, 1, -1, -1, -1),
    (-1, 1, -1, -1, 1, -1, -1, 1),
    (-1, 1, -1, 1, -1, -1, 1, -1),
    (-1, 1, 1, -1, -1, 1, -1, -1),
    (1, -1, -1, 1, 1, -1, 1, 1),
    (1, -1, 1, -1, 1, 1, -1, 1),
    (1, -1, 1, 1, -1, 1, 1, -1),
    (1, 1, -1, -1, -1, 1, 1, 1),
]

# Test: does NQR element b=3 map orbit_0 to orbit_28?
print(f"\n--- NQR ORBIT PAIRING ---")
for b in nqr[:3]:
    Pb = P_nqr[b]
    print(f"\nP_{b} (NQR element {b}):")
    mapped = []
    for sigma in orbit_0:
        new_s = tuple(int(x) for x in Pb @ np.array(sigma, dtype=float))
        in_28 = new_s in orbit_28
        mapped.append(new_s)
        print(f"  {sigma} → {new_s}  in orbit_28: {in_28}")

    # Check if ALL of orbit_0 maps to orbit_28
    all_in_28 = all(s in orbit_28 for s in mapped)
    print(f"  All orbit_0 → orbit_28: {all_in_28}")

# Also check: P_b maps orbit_28 to orbit_0?
print(f"\nReverse: P_3 maps orbit_28 → orbit_0?")
Pb = P_nqr[3]
mapped_back = []
for sigma in orbit_28:
    new_s = tuple(int(x) for x in Pb @ np.array(sigma, dtype=float))
    in_0 = new_s in orbit_0
    mapped_back.append(new_s)
    print(f"  {sigma} → {new_s}  in orbit_0: {in_0}")

print(f"  All orbit_28 → orbit_0: {all(s in orbit_0 for s in mapped_back)}")

# Now let's understand WHY this pairing preserves H.
print(f"\n--- WHY NQR PAIRING PRESERVES H ---")
print(f"""
For b ∈ NQR, the map x → bx on Z_p is an ANTI-automorphism of T_Paley.
On orientations: P_b σ gives a tournament T(P_b σ) = T(σ)^op.
Since H(T) = H(T^op), we have H(P_b σ) = H(σ).

But P_b is NOT in the QR group {{P_a : a ∈ QR}}.
So P_b maps each QR orbit to a DIFFERENT QR orbit.
All QR orbits in the same H-level are paired by ANY NQR element.

KEY: The FULL multiplicative group (Z/pZ)* = QR ∪ NQR acts on orientations.
QR orbits have |QR| = m elements each.
(Z/pZ)* orbits have |QR| + |NQR| = p-1 elements, but orientations
are in {{±1}}^m (only 2^m elements), so the orbit sizes are smaller.

The pairing: each (Z/pZ)* orbit = exactly 2 QR orbits (related by NQR).
Since (Z/pZ)* includes BOTH QR and NQR, and both preserve H (one by
automorphism, the other by anti-automorphism + H=H(T^op)).
""")

# For p=13, verify the same pairing
print(f"\n--- VERIFICATION AT p=13 ---")
p13 = 13
m13 = 6
nqr13 = sorted(a for a in range(1, p13) if legendre(a, p13) == -1)
P13_nqr = {b: make_signed_permutation(b, p13) for b in nqr13}

orbit_13_0 = [
    (-1, -1, -1, -1, -1, -1),
    (-1, -1, 1, 1, -1, -1),
    (-1, 1, 1, -1, 1, 1),
    (1, -1, -1, 1, -1, -1),
    (1, 1, -1, -1, 1, 1),
    (1, 1, 1, 1, 1, 1),
]

orbit_13_1 = [
    (-1, -1, -1, 1, 1, 1),
    (-1, 1, -1, 1, -1, 1),
    (-1, 1, -1, 1, 1, -1),
    (1, -1, 1, -1, -1, 1),
    (1, -1, 1, -1, 1, -1),
    (1, 1, 1, -1, -1, -1),
]

b13 = nqr13[0]
Pb13 = P13_nqr[b13]
print(f"NQR element b={b13} at p=13:")
all_maps = True
for sigma in orbit_13_0:
    new_s = tuple(int(x) for x in Pb13 @ np.array(sigma, dtype=float))
    if new_s not in orbit_13_1:
        all_maps = False
        print(f"  FAIL: {sigma} → {new_s} NOT in orbit_1")
if all_maps:
    print(f"  All orbit_0 → orbit_1: True ✓")

# The structure theorem
print(f"\n{'='*70}")
print("STRUCTURE THEOREM FOR p ≡ 1 mod 4")
print("="*70)
print(f"""
THEOREM: For p ≡ 1 mod 4, the orientation cube {{±1}}^m decomposes into:
  - (2^m / (2m)) QR orbits, each of size m (Interval lies in one such orbit)
  - Wait, 2^m / m is not always an integer. Let's check...

  p=13: 2^6 / 6 = 64/6 = 10.67... but we found 12 orbits of size 6 and
    a few of size 2. So orbit sizes aren't uniform.

Actually from the data:
  p=13: 12 QR orbits. 10 orbits of size 6, plus 2 orbits of size 2?
    12*6 = 72 ≠ 64, so NOT all size 6. Let me recheck.
    From the output: each H value has 12 orientations (mostly). 12*5 + 4 = 64.
    So: 5 H-values with 12 orientations + 1 H-value with 4 = 64. ✓
    Orbits: 5 × 2 orbits of 6 + 1 × 2 orbits of 2 = 10 + 2 = 12 orbits. ✓

  p=17: 32 QR orbits. 16 H-values × 2 orbits = 32. Each orbit size 8.
    32*8 = 256 = 2^8. ✓

  So orbit sizes CAN be different (6 or 2 at p=13, but uniform 8 at p=17).

  The uniform size at p=17 (all orbits size 8 = m) happens because:
  - The QR action on chord types is transitive
  - No chord type is fixed by any non-identity QR element
  - (Equivalently: no a ∈ QR with a ≠ 1 has ak ≡ ±k mod p for any chord k)

  At p=13: some orbits are smaller because the QR action has fixed points
  on certain orientations. Specifically, (1, -1, 1, 1, -1, -1) = σ_P is
  fixed by all QR elements → orbit size = 1 (but paired with -σ_P → 2).

CORRECT FORMULATION:
  The orientation cube decomposes into QR orbits.
  Each orbit has size dividing m.
  NQR multiplication pairs orbits: each (Z/pZ)* orbit = 2 QR orbits.
  H is constant on (Z/pZ)* orbits (since both QR and NQR preserve H).
  So each H-level is a union of (Z/pZ)* orbits.

  For p ≡ 3 mod 4, this is DIFFERENT:
  - NQR includes -1
  - P_{-1} = -I maps σ → -σ
  - σ and -σ may be in DIFFERENT QR orbits (since -1 ∉ QR)
  - But they give the same H (since H(σ) = H(-σ))
  - So the pairing is just complement symmetry
  - The Paley orientation σ_P is in a QR orbit of size 1
    (since it's fixed by all P_a) → exceptional behavior!
""")

# Count which orientations at p=13 have orbit size < m
print("p=13: Checking orbit sizes...")
all_sigmas = list(product([1, -1], repeat=m13))
qr13 = sorted(a for a in range(1, p13) if legendre(a, p13) == 1)
P13_qr = {a: make_signed_permutation(a, p13) for a in qr13}

visited = set()
for sigma in all_sigmas:
    if sigma in visited:
        continue
    orbit = set()
    queue = [sigma]
    while queue:
        s = queue.pop()
        if s in orbit:
            continue
        orbit.add(s)
        for a in qr13:
            Pa = P13_qr[a]
            new_s = tuple(int(x) for x in Pa @ np.array(s, dtype=float))
            if new_s not in orbit:
                queue.append(new_s)
    visited.update(orbit)
    if len(orbit) < m13:
        print(f"  Small orbit (size {len(orbit)}): {sorted(orbit)}")

print("\nDONE.")
