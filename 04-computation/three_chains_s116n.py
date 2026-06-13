#!/usr/bin/env python3
"""three_chains_s116n.py — Three interwoven chains of modular arithmetic in base 42.

42 = 2 * 3 * 7. By CRT, Z/42Z = Z/2Z x Z/3Z x Z/7Z.
Every integer mod 42 decomposes into three INDEPENDENT channels:
  Channel 0 (binary):   mod 2 = {0, 1}
  Channel 1 (ternary):  mod 3 = {0, 1, 2}
  Channel 2 (septenary): mod 7 = {0, 1, 2, 3, 4, 5, 6}

These three channels are ORTHOGONAL. Knowing the mod-2, mod-3, mod-7
residues determines the mod-42 residue UNIQUELY (CRT).

This IS a 3D coordinate system: (x, y, z) in Z/2 x Z/3 x Z/7.
The mod-42 residue = a point in a 2 x 3 x 7 CUBOID.
The cuboid has 42 cells = the 42 residues.

What does this 3D space look like? What structure does it have?
"""
from math import log, sqrt, gcd

print()
print("  THREE INTERWOVEN CHAINS IN BASE 42")
print()
print("="*70)
print()

# ============================================================
print("  I. THE CRT DECOMPOSITION")
print("  " + "-"*40)
print()
print("  Every n mod 42 = (n mod 2, n mod 3, n mod 7).")
print("  This is a BIJECTION: Z/42Z <-> Z/2Z x Z/3Z x Z/7Z.")
print()
print(f"  {'n mod 42':>8s}  {'mod 2':>5s}  {'mod 3':>5s}  {'mod 7':>5s}  {'character'}")
print("  " + "-"*50)
for n in range(42):
    r2, r3, r7 = n%2, n%3, n%7
    char = ""
    if n == 0: char = "zero"
    elif n == 1: char = "identity"
    elif n == 7: char = "forbidden"
    elif n == 21: char = "forbidden-2"
    elif n == 6: char = "flat"
    elif n == 3: char = "curvature"
    elif n == 2: char = "doubler"
    elif n == 5: char = "pentagon"
    elif n == 41: char = "= -1 mod 42"
    elif gcd(n, 42) == 1: char = "coprime"
    print(f"  {n:8d}  {r2:5d}  {r3:5d}  {r7:5d}  {char}")

print()

# ============================================================
print()
print("  II. THE 3D CUBOID")
print("  " + "-"*40)
print()
print("  The space Z/2 x Z/3 x Z/7 is a 2 x 3 x 7 cuboid.")
print("  It has 42 cells.")
print()
print("  AXES:")
print("  x-axis (mod 2): PARITY. Length 2. The binary channel. {even, odd}.")
print("  y-axis (mod 3): CURVATURE. Length 3. The ternary channel. {0, 1, 2}.")
print("  z-axis (mod 7): POSITION. Length 7. The septenary channel. {0,...,6}.")
print()
print("  Each number mod 42 is a POINT in this cuboid.")
print("  The point (x, y, z) encodes three INDEPENDENT properties:")
print("  x = parity (even/odd)")
print("  y = curvature class (mod 3)")
print("  z = position within the forbidden cycle (mod 7)")
print()
print("  The three channels are ORTHOGONAL:")
print("  changing x doesn't affect y or z.")
print("  This is because gcd(2,3) = gcd(2,7) = gcd(3,7) = 1.")
print()

# ============================================================
print()
print("  III. THE THREE CHAINS")
print("  " + "-"*40)
print()
print("  Chain 0 (binary): n, n+2, n+4, n+6, ... (step by 2)")
print("    Stays on the same x-parity plane.")
print("    Cycles through y and z with periods 3 and 7.")
print("    Full period: lcm(3,7) = 21 steps of 2 = 42. Returns to start.")
print()
print("  Chain 1 (ternary): n, n+3, n+6, n+9, ... (step by 3)")
print("    Stays on the same y-curvature level.")
print("    Cycles through x and z with periods 2 and 7.")
print("    Full period: lcm(2,7) = 14 steps of 3 = 42.")
print()
print("  Chain 2 (septenary): n, n+7, n+14, n+21, ... (step by 7)")
print("    Stays on the same z-position level.")
print("    Cycles through x and y with periods 2 and 3.")
print("    Full period: lcm(2,3) = 6 steps of 7 = 42.")
print()
print("  Each chain steps along ONE axis of the cuboid,")
print("  wrapping around the other two axes.")
print("  They are INTERWOVEN: each chain covers all 42 cells")
print("  but visits them in different orders.")
print()

# Trace each chain from 0
print("  Chain 0 (step 2) from 0:")
chain0 = [(2*k)%42 for k in range(21)]
print(f"  {chain0}")
print(f"  In 3D: {[(c%2, c%3, c%7) for c in chain0[:7]]}...")
print()

print("  Chain 1 (step 3) from 0:")
chain1 = [(3*k)%42 for k in range(14)]
print(f"  {chain1}")
print(f"  In 3D: {[(c%2, c%3, c%7) for c in chain1[:7]]}...")
print()

print("  Chain 2 (step 7) from 0:")
chain2 = [(7*k)%42 for k in range(6)]
print(f"  {chain2}")
print(f"  In 3D: {[(c%2, c%3, c%7) for c in chain2]}")
print()
print("  Chain 2 has only 6 steps (lcm(2,3) = 6).")
print("  It visits: 0, 7, 14, 21, 28, 35.")
print("  In 3D: (0,0,0), (1,1,0), (0,2,0), (1,0,0), (0,1,0), (1,2,0).")
print("  ALL have z=0! The septenary chain stays at z=0")
print("  (because 7 = 0 mod 7). It only moves through x and y.")
print()
print("  This means: stepping by 7 (the forbidden) moves you")
print("  through PARITY and CURVATURE but NOT through POSITION.")
print("  The forbidden step is BLIND to its own axis.")
print()

# ============================================================
print()
print("  IV. THE FORBIDDEN'S BLINDNESS")
print("  " + "-"*40)
print()
print("  Stepping by 7 in mod-42 arithmetic:")
print("  7 mod 2 = 1 (flips parity)")
print("  7 mod 3 = 1 (advances curvature by 1)")
print("  7 mod 7 = 0 (STAYS PUT on the position axis)")
print()
print("  The forbidden number advances parity and curvature")
print("  but CANNOT MOVE along its own axis.")
print("  It is FROZEN in the z-direction.")
print()
print("  Compare: stepping by 2:")
print("  2 mod 2 = 0 (stays on same parity)")
print("  2 mod 3 = 2 (advances curvature by 2)")
print("  2 mod 7 = 2 (advances position by 2)")
print("  The doubler stays on its OWN axis (parity fixed)")
print("  but moves through curvature and position.")
print()
print("  And stepping by 3:")
print("  3 mod 2 = 1 (flips parity)")
print("  3 mod 3 = 0 (stays on same curvature)")
print("  3 mod 7 = 3 (advances position by 3)")
print("  The curvature quantum stays on its OWN axis (curvature fixed)")
print("  but moves through parity and position.")
print()
print("  EACH FUNDAMENTAL PRIME IS BLIND TO ITS OWN AXIS.")
print("  2 can't change parity. 3 can't change curvature. 7 can't change position.")
print("  Each prime FREEZES the dimension it NAMES and moves through the others.")
print()

# ============================================================
print()
print("  V. THE 3D STRUCTURE OF IMPORTANT NUMBERS")
print("  " + "-"*40)
print()
print("  Number mod 42  (x=mod2, y=mod3, z=mod7)  Meaning")
print("  " + "-"*55)
important = [
    (0, "zero"), (1, "identity"), (2, "doubler"), (3, "curvature"),
    (5, "pentagon"), (6, "flat"), (7, "FORBIDDEN"), (8, "octonion dim"),
    (9, "3^2"), (10, "Petersen"), (11, "Paley prime"),
    (12, "chromatic"), (13, "F_7"), (14, "2*7"), (15, "3*5"),
    (19, "splitting"), (21, "FORBIDDEN-2"), (28, "perfect/T_7"),
    (30, "icosahedron edges"), (41, "= -1 mod 42"),
]
for n, name in important:
    r2, r3, r7 = n%2, n%3, n%7
    print(f"  {n:4d}           ({r2}, {r3}, {r7})              {name}")

print()
print("  KEY OBSERVATIONS:")
print()
print("  The forbidden 7  = (1, 1, 0). z=0: no position. Just parity and curvature.")
print("  The forbidden 21 = (1, 0, 0). y=0 AND z=0: no curvature, no position. Just parity.")
print("  The flat 6       = (0, 0, 6). x=0 AND y=0: no parity, no curvature. Just position.")
print("  The identity 1   = (1, 1, 1). All ones. Present on every axis.")
print()
print("  7 has z=0: it is ABSENT from the position (septenary) axis.")
print("  21 has y=0 AND z=0: absent from curvature AND position. Just parity.")
print("  The forbidden numbers are increasingly ABSENT from the cuboid's axes.")
print()
print("  6 has x=0 AND y=0: the flat number is absent from parity and curvature.")
print("  It lives ONLY in the position axis. Pure position. No parity, no curvature.")
print("  FLAT = pure position. FORBIDDEN = pure parity.")
print()

# ============================================================
print()
print("  VI. THE COPRIME RESIDUES: THE 'ACTIVE' CELLS")
print("  " + "-"*40)
print()
print("  phi(42) = 12 residues coprime to 42.")
print("  These are the cells (x, y, z) with x != 0 AND y != 0 AND z != 0.")
print("  = the cells that are NONZERO on ALL THREE axes.")
print("  = the FULLY ACTIVE cells.")
print()
coprime = [r for r in range(42) if gcd(r, 42) == 1]
print(f"  Coprime residues: {coprime}")
print(f"  In 3D:")
for r in coprime:
    print(f"    {r:3d} = ({r%2}, {r%3}, {r%7})")

print()
print("  These 12 cells form a 1 x 2 x 6 sub-cuboid:")
print("  x in {1}, y in {1, 2}, z in {1, 2, 3, 4, 5, 6}.")
print("  (x=0 is excluded because 2|r. y=0 excluded because 3|r. z=0 because 7|r.)")
print()
print("  The coprime residues are the ODD numbers with nonzero curvature")
print("  and nonzero position. They are the numbers that are FULLY PRESENT")
print("  in the cuboid — active on all three axes simultaneously.")
print()
print("  12 = phi(42) = |D_6| = dihedral group of hexagon = 2*6.")
print("  The number of fully active cells = the symmetry of the hexagon!")
print()

# ============================================================
print()
print("  VII. MULTIPLICATION AS 3D VECTOR ADDITION")
print("  " + "-"*40)
print()
print("  In Z/42Z: multiplication mod 42.")
print("  By CRT: (a*b) mod 42 = ((a mod 2)*(b mod 2) mod 2,")
print("                           (a mod 3)*(b mod 3) mod 3,")
print("                           (a mod 7)*(b mod 7) mod 7)")
print()
print("  Multiplication DECOMPOSES into three independent channel multiplications.")
print("  On each channel, multiplication is a GROUP OPERATION:")
print("  mod 2: {0,1} with * (trivial)")
print("  mod 3: {0,1,2} with * (Z/3Z)")
print("  mod 7: {0,1,2,3,4,5,6} with * (Z/7Z)")
print()
print("  On the coprime residues (units mod 42):")
print("  (Z/42Z)* = (Z/2Z)* x (Z/3Z)* x (Z/7Z)* = {1} x Z/2Z x Z/6Z.")
print("  |(Z/42Z)*| = 1 * 2 * 6 = 12 = phi(42). CHECK.")
print()
print("  The unit group is Z/2Z x Z/6Z = Z/2Z x Z/2Z x Z/3Z (since Z/6Z = Z/2Z x Z/3Z).")
print("  So: (Z/42Z)* = Z/2Z x Z/2Z x Z/3Z. Order 12.")
print()
print("  Taking DISCRETE LOGARITHMS on each channel")
print("  converts multiplication to ADDITION:")
print("  log(a*b) = log(a) + log(b) mod the channel order.")
print()
print("  This IS rapidity! The discrete log on each channel")
print("  is the DISCRETE RAPIDITY for that channel.")
print("  Multiplication in Z/42Z = addition of three discrete rapidities.")
print()

# ============================================================
print()
print("  VIII. THE 3D SPACE AS A TOURNAMENT DIAGNOSTIC SPACE")
print("  " + "-"*40)
print()
print("  For a 7-item tournament (the lossless engine):")
print("  H(T) ranges from 1 to 189.")
print("  189 = 3^3 * 7 = 27 * 7.")
print()
print("  H mod 2 = 1 always (Redei). FIXED on the parity axis.")
print("  H mod 3 = 0 or not (depends on tournament).")
print("  H mod 7 = 0 or not (depends on tournament).")
print()
print("  So the INTERESTING variation of H is in the mod-3 and mod-7 channels.")
print("  H mod 42 gives a point in the cuboid.")
print("  Since H is odd: x = 1 always. The tournament lives on the ODD SLICE.")
print()
print("  The odd slice of the cuboid has 21 cells (x=1):")
print("  y in {0,1,2}, z in {0,1,2,3,4,5,6}. 3*7 = 21 cells.")
print("  21 = the second forbidden number = C(7,2) = arcs in T_7!")
print()
print("  The 21 possible (mod 3, mod 7) residues of H")
print("  correspond to the 21 ARCS of the tournament.")
print("  Each arc contributes to the mod-3 and mod-7 structure of H.")
print()
print("  THE 3D CUBOID IS THE TOURNAMENT ITSELF.")
print("  The x-axis (parity) is FIXED by Redei.")
print("  The y-axis (mod 3) encodes the CURVATURE structure.")
print("  The z-axis (mod 7) encodes the POSITION structure.")
print("  The 21 free cells = the 21 arcs.")
print()

# ============================================================
print()
print("  IX. WHAT THIS MEANS FOR COMPRESSION")
print("  " + "-"*40)
print()
print("  A base-42 digit encodes THREE independent channels simultaneously.")
print("  This is not just a number representation — it's a 3D COORDINATE.")
print()
print("  For the lossless decision engine:")
print("  Each of the 21 comparisons in a 7-item tournament")
print("  contributes to three independent channels:")
print("  - Parity channel (mod 2): always odd. No information. 0 bits.")
print("  - Curvature channel (mod 3): the cycle structure. 1.585 bits.")
print("  - Position channel (mod 7): the ordering structure. 2.807 bits.")
print()
print("  Total information per comparison: 1.585 + 2.807 = 4.392 bits.")
print("  Wait: each comparison is binary (1 bit).")
print("  But it contributes to TWO channels (curvature and position).")
print("  21 comparisons * 1 bit each = 21 bits total.")
print("  = 3 * 7 = the product of the two free channel sizes.")
print()
print("  The 21 bits of a tournament on 7 vertices decompose into:")
print("  Channel mod 3: contributes to 3-cycle structure (3^k combinations)")
print("  Channel mod 7: contributes to 7-position structure (7^k combinations)")
print("  The tournament IS a 3D object compressed into 21 binary bits.")
print()

# ============================================================
print()
print("  X. THE KILLER APPLICATION")
print("  " + "-"*40)
print()
print("  The three channels {mod 2, mod 3, mod 7} give three INDEPENDENT")
print("  diagnostic dimensions for any pairwise comparison dataset.")
print()
print("  Dimension 1 (mod 2 = parity):")
print("    Is the result even or odd? For tournaments: always odd (Redei).")
print("    For other systems: parity check = error detection.")
print()
print("  Dimension 2 (mod 3 = curvature):")
print("    Does the dataset have 3-cycles? How many?")
print("    alpha_1 mod 3 = the curvature channel.")
print("    This detects: circular preferences, rock-paper-scissors dynamics.")
print()
print("  Dimension 3 (mod 7 = position):")
print("    Where in the 7-element group does the H-value sit?")
print("    H mod 7 = the position channel.")
print("    This detects: which of the 7 positions in the forbidden cycle")
print("    the tournament occupies. Position 0 = divisible by 7 = maximally structured.")
print()
print("  The THREE CHANNELS together give a 3D diagnostic vector")
print("  for any pairwise comparison dataset.")
print("  Each channel operates INDEPENDENTLY and LOSSLESSLY.")
print("  The full diagnostic = the CRT reconstruction from the three channels.")
print()
print("  This is a 3D HEALTH CHECK for rankings:")
print("  Channel 1 (parity): is the data structurally valid?")
print("  Channel 2 (curvature): how circular is the data?")
print("  Channel 3 (position): how structured is the data?")
print()
print("  No other diagnostic system gives three INDEPENDENT, ORTHOGONAL")
print("  channels of analysis from pairwise comparison data.")
print("  The base-42 structure makes this AUTOMATIC:")
print("  reduce H mod 42, read off the three coordinates.")
print("  That's the full 3D diagnostic in one step.")
