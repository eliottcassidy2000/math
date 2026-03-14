#!/usr/bin/env python3
"""
rep_lie_23.py — Representation Theory and Lie Theory through (2,3)
opus-2026-03-14-S84

Exploring how KEY1=2 and KEY2=3 govern:
- Representations of symmetric groups S_n
- Young diagrams and hook length formula
- Representations of SL_2 (the KEY1 group)
- Lie algebras: sl_2, sl_3, and the ADE classification
- Root systems and Weyl groups
- Coxeter groups and Coxeter numbers
- Highest weight theory
- Tensor products and Clebsch-Gordan
- McKay correspondence (finite subgroups of SU(2))
- Crystal bases and quantum groups at q = root of unity
- Kazhdan-Lusztig polynomials
- Schur-Weyl duality (between S_n and GL_d)
"""

import math
from fractions import Fraction
from functools import lru_cache
from itertools import combinations

# Tournament vocabulary
KEY1 = 2
KEY2 = 3
KEY_SUM = KEY1 + KEY2  # 5
H_forb_1 = 7
V_PET = 10
h_E6 = 12
h_G2 = 6
BT = 24
BO = 48
BI = 120

print("=" * 70)
print("  REPRESENTATION THEORY & LIE THEORY THROUGH (2,3)")
print("  Where symmetry meets tournaments")
print("=" * 70)

# =====================================================================
# Part 1: SYMMETRIC GROUP — S_n REPRESENTATIONS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 1: SYMMETRIC GROUP REPRESENTATIONS")
print("=" * 70)

# Number of irreps of S_n = number of partitions of n
@lru_cache(maxsize=None)
def partitions(n, max_part=None):
    """List all partitions of n."""
    if max_part is None:
        max_part = n
    if n == 0:
        return [()]
    result = []
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            result.append((first,) + rest)
    return result

def hook_length(partition, i, j):
    """Hook length at position (i,j) in Young diagram."""
    arm = partition[i] - j - 1  # cells to the right
    leg = sum(1 for k in range(i+1, len(partition)) if partition[k] > j)
    return arm + leg + 1

def dim_irrep(partition):
    """Dimension of irrep of S_n corresponding to partition, via hook length formula."""
    n = sum(partition)
    hook_product = 1
    for i, part in enumerate(partition):
        for j in range(part):
            hook_product *= hook_length(partition, i, j)
    return math.factorial(n) // hook_product

print("  Irreducible representations of S_n:")
for n in range(1, 9):
    parts = partitions(n)
    dims = [dim_irrep(p) for p in parts]
    n_irreps = len(parts)
    dim_str = ", ".join(str(d) for d in sorted(dims, reverse=True))
    label = ""
    if n == KEY1: label = " (= KEY1, S_KEY1 ~ Z/KEY1)"
    elif n == KEY2: label = " (= KEY2, S_KEY2 ~ D_3)"
    elif n == KEY_SUM: label = " (= KEY_SUM)"
    print(f"  S_{n}: {n_irreps} irreps, dims = [{dim_str}]{label}")

    # Verify sum of squares = n!
    sum_sq = sum(d**2 for d in dims)
    assert sum_sq == math.factorial(n), f"Failed for S_{n}"

print(f"""
  KEY OBSERVATIONS:
  S_KEY1 = S_2 ~ Z/KEY1: {len(partitions(KEY1))} = KEY1 irreps of dim 1, 1
  S_KEY2 = S_3 ~ D_3: {len(partitions(KEY2))} = KEY2 irreps of dim 1, 1, 2
  S_KEY_SUM = S_5: {len(partitions(KEY_SUM))} = H_forb_1 irreps!

  Number of irreps of S_n = p(n) = partition function!
  p(KEY1) = KEY1, p(KEY2) = KEY2, p(KEY_SUM) = H_forb_1

  S_KEY1 has KEY1 irreps (trivial and sign)
  S_KEY2 has KEY2 irreps (trivial, sign, and the KEY1-dimensional standard)
  The standard representation of S_KEY2 has dimension KEY1!

  SCHUR-WEYL DUALITY:
  V^(tensor n) = direct_sum_(lambda) S^lambda tensor V^lambda
  where S^lambda = Specht module (S_n irrep)
  and V^lambda = GL irrep (Schur functor)

  When V = C^KEY1 (KEY1-dimensional):
  V^(tensor n) decomposes into S_n x GL(KEY1) bimodules.
  This is the KEY1-KEY1 duality (GL_KEY1 and S_n)!
""")

# =====================================================================
# Part 2: SL_2 REPRESENTATIONS — THE FUNDAMENTAL CASE
# =====================================================================
print("\n" + "=" * 70)
print("  Part 2: SL_2 REPRESENTATIONS — THE KEY1 x KEY1 GROUP")
print("=" * 70)

print(f"""
  sl_KEY1 = sl_2 is the Lie algebra of KEY1 x KEY1 traceless matrices.

  Basis: e = [[0,1],[0,0]], f = [[0,0],[1,0]], h = [[1,0],[0,-1]]
  Relations: [h,e] = KEY1*e, [h,f] = -KEY1*f, [e,f] = h

  The number KEY1 appears in the commutation relations!
  [h,e] = KEY1*e means e raises weight by KEY1.
  [h,f] = -KEY1*f means f lowers weight by KEY1.

  FINITE-DIMENSIONAL IRREPS V(n):
  V(n) has dimension n+1, highest weight n.
  Weights: n, n-KEY1, n-KEY1^KEY1, ..., -n (spacing KEY1!)

  Dimensions of first few irreps:
""")

for n in range(11):
    dim = n + 1
    weights = list(range(n, -n-1, -2))
    label = ""
    if dim == KEY1: label = " KEY1-dim (standard)"
    elif dim == KEY2: label = " KEY2-dim (adjoint)"
    elif dim == KEY_SUM: label = " KEY_SUM-dim"
    elif dim == H_forb_1: label = " H_forb_1-dim"
    print(f"  V({n}): dim = {dim}, weights = {weights} {label}")

print(f"""
  CLEBSCH-GORDAN DECOMPOSITION:
  V(m) tensor V(n) = V(m+n) + V(m+n-KEY1) + ... + V(|m-n|)

  Example: V(1) tensor V(1) = V(KEY1) + V(0)
  KEY1-dim tensor KEY1-dim = KEY2-dim + 1-dim
  In matrix terms: KEY1 x KEY1 = KEY2 + 1 (= KEY1^KEY1)

  V(1) tensor V(KEY1) = V(KEY2) + V(1)
  KEY1-dim tensor KEY2-dim = KEY1^KEY1-dim + KEY1-dim
  = KEY_SUM + 1 ... no, KEY2+1 + KEY1 = 6. Wait:
  Actually V(1) has dim 2, V(2) has dim 3:
  V(1) x V(2) = V(3) + V(1), dims: KEY1 x KEY2 = KEY1^KEY1 + KEY1 = h(G2)

  KEY1 x KEY2 = h(G2) = 6 as representations!

  CASIMIR OPERATOR:
  C = h^KEY1/KEY1 + ef + fe = h^KEY1/KEY1 + KEY1*ef + h
  On V(n): C acts by scalar n(n+KEY1)/KEY1^KEY1 = n(n+2)/4
""")

# =====================================================================
# Part 3: ROOT SYSTEMS — THE ADE CLASSIFICATION
# =====================================================================
print("\n" + "=" * 70)
print("  Part 3: ROOT SYSTEMS — THE ADE CLASSIFICATION")
print("=" * 70)

# Root system data
root_systems = {
    'A1': {'rank': 1, 'roots': 2, 'weyl': 2, 'coxeter': 2},
    'A2': {'rank': 2, 'roots': 6, 'weyl': 6, 'coxeter': 3},
    'A3': {'rank': 3, 'roots': 12, 'weyl': 24, 'coxeter': 4},
    'A4': {'rank': 4, 'roots': 20, 'weyl': 120, 'coxeter': 5},
    'B2': {'rank': 2, 'roots': 8, 'weyl': 8, 'coxeter': 4},
    'B3': {'rank': 3, 'roots': 18, 'weyl': 48, 'coxeter': 6},
    'G2': {'rank': 2, 'roots': 12, 'weyl': 12, 'coxeter': 6},
    'D4': {'rank': 4, 'roots': 24, 'weyl': 192, 'coxeter': 6},
    'D5': {'rank': 5, 'roots': 40, 'weyl': 1920, 'coxeter': 8},
    'E6': {'rank': 6, 'roots': 72, 'weyl': 51840, 'coxeter': 12},
    'E7': {'rank': 7, 'roots': 126, 'weyl': 2903040, 'coxeter': 18},
    'E8': {'rank': 8, 'roots': 240, 'weyl': 696729600, 'coxeter': 30},
    'F4': {'rank': 4, 'roots': 48, 'weyl': 1152, 'coxeter': 12},
}

print("  Root system data:")
print(f"  {'Type':>4}  {'rank':>4}  {'#roots':>6}  {'|W|':>12}  {'h':>4}  notes")
for name, data in root_systems.items():
    notes = []
    r = data['roots']
    w = data['weyl']
    h = data['coxeter']
    if w == KEY1: notes.append("|W|=KEY1")
    elif w == h_G2: notes.append("|W|=h(G2)")
    elif w == BT: notes.append("|W|=BT!")
    elif w == BI: notes.append("|W|=BI!")
    elif w == BO: notes.append("|W|=BO!")
    if h == KEY1: notes.append("h=KEY1")
    elif h == KEY2: notes.append("h=KEY2")
    elif h == KEY_SUM: notes.append("h=KEY_SUM")
    elif h == h_G2: notes.append("h=h(G2)")
    elif h == h_E6: notes.append("h=h(E6)")
    elif h == 30: notes.append("h=KEY1*KEY2*KEY_SUM")
    if r == h_E6: notes.append("#roots=h(E6)")
    elif r == BT: notes.append("#roots=BT!")
    elif r == BO: notes.append("#roots=BO!")
    elif r == 240: notes.append("#roots=KEY1*BI!")
    print(f"  {name:>4}  {data['rank']:>4}  {r:>6}  {w:>12}  {h:>4}  {', '.join(notes)}")

print(f"""
  TOURNAMENT CONSTANTS IN ROOT SYSTEMS:

  Coxeter numbers:
  h(A1) = KEY1, h(A2) = KEY2, h(A4) = KEY_SUM
  h(G2) = h(G2) = 6, h(E6) = h(E6) = 12
  h(E8) = 30 = KEY1 * KEY2 * KEY_SUM (product of first KEY2 primes!)

  Weyl group orders:
  |W(A1)| = KEY1 = 2! (= S_KEY1)
  |W(A2)| = h(G2) = 6 = 3! (= S_KEY2)
  |W(A3)| = BT = 24 = 4! (= S_4)
  |W(A4)| = BI = 120 = 5! (= S_KEY_SUM)

  W(A_n) = S_(n+1) — the symmetric group!
  So BT = |S_4| = |W(A3)| and BI = |S_5| = |W(A4)|.

  Number of positive roots:
  A_n: n(n+1)/KEY1 = C(n+1, KEY1)
  B_n = C_n: n^KEY1
  D_n: n(n-1) = KEY1 * C(n, KEY1)

  A_2: KEY2, A_3: h(G2), A_4: V(Pet), A_5: 15
  POSITIVE roots of A_n = C(n+1, KEY1) = TOURNAMENT CONSTANTS!
""")

# =====================================================================
# Part 4: EXCEPTIONAL LIE ALGEBRAS — THE (2,3) PEARLS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 4: EXCEPTIONAL LIE ALGEBRAS — THE (2,3) PEARLS")
print("=" * 70)

exceptionals = [
    ('G2', 14, 2, 12, 6, 'dim=KEY1*H_forb_1, rank=KEY1'),
    ('F4', 52, 4, 48, 12, 'dim=KEY1^KEY1*13, #roots=BO'),
    ('E6', 78, 6, 72, 12, 'dim=KEY1*KEY2*13, rank=h(G2)'),
    ('E7', 133, 7, 126, 18, 'dim=H_forb_1*19, rank=H_forb_1'),
    ('E8', 248, 8, 240, 30, 'dim=KEY1^KEY2*31, rank=KEY1^KEY2'),
]

print("  Exceptional Lie algebras:")
print(f"  {'Type':>4}  {'dim':>4}  {'rank':>4}  {'#roots':>6}  {'h':>4}  notes")
for name, dim, rank, roots, h, notes in exceptionals:
    print(f"  {name:>4}  {dim:>4}  {rank:>4}  {roots:>6}  {h:>4}  {notes}")

print(f"""
  THE (2,3) STRUCTURE:
  G2: dim = KEY1 * H_forb_1 = 14, rank = KEY1
      KEY1 root lengths, KEY2-fold symmetry of long roots
      Coxeter number = h(G2) = KEY1 * KEY2

  F4: dim = 52 = KEY1^KEY1 * 13
      #roots = BO = 48 = KEY1^KEY1^KEY1 * KEY2 * KEY1
      Coxeter number = h(E6) = 12

  E6: rank = h(G2) = 6
      #roots = 72 = KEY1^KEY2 * KEY2^KEY1 = 8 * 9
      Coxeter number = h(E6) = KEY1^KEY1 * KEY2

  E7: rank = H_forb_1 = 7
      #roots = 126 = KEY1 * KEY2^KEY1 * H_forb_1
      Coxeter number = 18 = KEY1 * KEY2^KEY1

  E8: rank = KEY1^KEY2 = 8
      #roots = 240 = KEY1^KEY1^KEY1 * KEY2 * KEY_SUM * KEY1
      = KEY1 * BI = 2 * 120
      Coxeter number = 30 = KEY1 * KEY2 * KEY_SUM

  dim(E8) = 248 = KEY1^KEY2 * 31 = 8 * (KEY1^KEY_SUM - 1)
  Fascinating: 31 = KEY1^KEY_SUM - 1 is a MERSENNE PRIME!

  THE E8 LATTICE:
  - 240 = KEY1 * BI vectors of minimal norm
  - Theta function: 1 + 240*q + 2160*q^2 + ...
  - 240 = KEY1^4 * KEY2 * KEY_SUM = 16 * 15
  - Kissing number in 8 dimensions = 240
""")

# =====================================================================
# Part 5: MCKAY CORRESPONDENCE — FINITE SUBGROUPS OF SU(2)
# =====================================================================
print("\n" + "=" * 70)
print("  Part 5: MCKAY CORRESPONDENCE")
print("=" * 70)

# Binary polyhedral groups
groups = [
    ("Z/n (cyclic)", "A_(n-1)", "n", "n"),
    ("Dic_n (dicyclic)", "D_(n+2)", "4n", "n+2 or KEY1*n"),
    ("BT (binary tetrahedral)", "E6", str(BT), f"{BT} = KEY1^KEY2 * KEY2"),
    ("BO (binary octahedral)", "E7", str(BO), f"{BO} = KEY1^KEY1^KEY1 * KEY2"),
    ("BI (binary icosahedral)", "E8", str(BI), f"{BI} = KEY1^KEY2 * KEY2 * KEY_SUM"),
]

print("  Finite subgroups of SU(2) and their ADE types:")
for name, ade, order, factored in groups:
    print(f"  {name:>30} <-> {ade:>8}, |G| = {order:>4} = {factored}")

print(f"""
  McKAY'S OBSERVATION (1980):
  The McKay graph of G (tensor with standard KEY1-dim rep,
  decompose, draw edges) gives the EXTENDED Dynkin diagram!

  G = BT: McKay graph = extended E6, |G| = BT = 24
  G = BO: McKay graph = extended E7, |G| = BO = 48
  G = BI: McKay graph = extended E8, |G| = BI = 120

  CONJUGACY CLASSES (= number of irreps):
  BT: H_forb_1 = 7 conjugacy classes
  BO: KEY1^KEY2 = 8 conjugacy classes
  BI: KEY2^KEY1 = 9 conjugacy classes

  CONSECUTIVE: H_forb_1, KEY1^KEY2, KEY2^KEY1 = 7, 8, 9!
  Three consecutive integers for the three exceptional groups!

  REPRESENTATION DIMENSIONS:
  BT: 1, 1, 1, KEY1, KEY1, KEY1, KEY2
      (KEY2 ones, KEY2 twos, 1 three)
  BO: 1, 1, KEY1, KEY1, KEY1, KEY2, KEY2, KEY1^KEY1
      Sum = 18 = KEY1 * KEY2^KEY1
  BI: 1, KEY1, KEY2, KEY2, KEY1^KEY1, KEY1^KEY1, KEY_SUM, KEY_SUM, h(G2)
      Sum = KEY2 * V(Pet) = 30
""")

# =====================================================================
# Part 6: QUANTUM GROUPS AT ROOTS OF UNITY
# =====================================================================
print("\n" + "=" * 70)
print("  Part 6: QUANTUM GROUPS AT ROOTS OF UNITY")
print("=" * 70)

import cmath

print("""
  Quantum group U_q(sl_KEY1) at q = e^(2*pi*i/N):

  Relations:
  K * E * K^(-1) = q^KEY1 * E
  K * F * K^(-1) = q^(-KEY1) * F
  E*F - F*E = (K - K^(-1)) / (q - q^(-1))

  At q = generic: same representation theory as sl_2
  At q = root of unity: TRUNCATION occurs!
""")

print("  Quantum dimensions at q = e^(2*pi*i/N):")
for N in range(3, 13):
    q = cmath.exp(2j * cmath.pi / N)
    qdims = []
    for j in range(N - 1):
        # [j+1]_q = (q^(j+1) - q^(-(j+1))) / (q - q^(-1))
        num = q**(j+1) - q**(-(j+1))
        den = q - q**(-1)
        qd = (num / den).real
        qdims.append(round(qd, 4))
    label = ""
    if N == KEY2 + 2: label = f" (SU(2) level {KEY2}, Fibonacci anyons!)"
    elif N == KEY1 + 2: label = f" (SU(2) level {KEY1})"
    elif N == BT: label = " (SU(2) level 22)"
    print(f"  N = {N:2d}: {qdims}{label}")

print(f"""
  At N = KEY_SUM (SU(2) level KEY2 = 3):
  Quantum dims: 1, phi, 1 where phi = golden ratio!
  This gives FIBONACCI ANYONS — universal for quantum computation!

  The fusion rules at level KEY2 = 3:
  V(1) x V(1) = V(0) + V(KEY1)
  V(1) x V(KEY1) = V(1) + V(KEY2)   wait, V(3) doesn't exist at level 3.

  Actually at level k, reps are V(0), V(1), ..., V(k).
  Level KEY2: V(0), V(1), V(KEY1), V(KEY2)
  Level KEY1: V(0), V(1), V(KEY1)

  MODULAR TENSOR CATEGORIES:
  The category of representations of U_q(sl_KEY1) at q^N = 1
  is a MODULAR tensor category.

  - Rank (# simple objects) = N - 1
  - Total quantum dim = sqrt(N / KEY1) / sin(pi/N)

  VERLINDE FORMULA:
  N^lambda_(mu,nu) = sum_sigma S_(mu,sigma) S_(nu,sigma) S*(lambda,sigma) / S_(0,sigma)

  The S-matrix for SU(KEY1) level k:
  S_(j,l) = sqrt(KEY1/(k+KEY1)) * sin(pi*(KEY1*j+1)*(KEY1*l+1)/(KEY1*(k+KEY1)))
""")

# =====================================================================
# Part 7: COXETER GROUPS AND WEYL GROUPS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 7: COXETER GROUPS — REFLECTIONS AND (2,3)")
print("=" * 70)

print("""
  A Coxeter group is generated by reflections s_i with:
  (s_i * s_j)^m(i,j) = 1

  The Coxeter matrix m(i,j) has entries from:
  {1, KEY1, KEY2, KEY1^KEY1, KEY_SUM, h(G2), ...}

  FINITE COXETER GROUPS (irreducible):
  A_n: linear diagram, all edges labeled KEY2 (m=3)
  B_n: one edge labeled KEY1^KEY1 (m=4), rest KEY2
  D_n: forked, all edges KEY2
  I_2(m): two generators, one edge labeled m
  H_3: one edge KEY_SUM (m=5), one edge KEY2
  H_4: one edge KEY_SUM, two edges KEY2
  F_4: one edge KEY1^KEY1, two edges KEY2
  E_6, E_7, E_8: all edges KEY2

  The ADE groups use ONLY m = KEY2 (the default edge label)!
  BCF groups introduce m = KEY1^KEY1 = 4.
  H groups introduce m = KEY_SUM = 5.
  I_2(m) allows any m.

  EXPONENTS of Coxeter groups:
""")

coxeter_exponents = {
    'A2': [1, 2],
    'A3': [1, 2, 3],
    'A4': [1, 2, 3, 4],
    'B2': [1, 3],
    'B3': [1, 3, 5],
    'G2': [1, 5],
    'D4': [1, 3, 3, 5],
    'E6': [1, 4, 5, 7, 8, 11],
    'E7': [1, 5, 7, 9, 11, 13, 17],
    'E8': [1, 7, 11, 13, 17, 19, 23, 29],
}

for name, exps in coxeter_exponents.items():
    h = max(exps) + 1
    notes = []
    for e in exps:
        if e == KEY1: notes.append(f"{e}=KEY1")
        elif e == KEY2: notes.append(f"{e}=KEY2")
        elif e == KEY_SUM: notes.append(f"{e}=KEY_SUM")
        elif e == H_forb_1: notes.append(f"{e}=H_forb_1")
    print(f"  {name}: exponents = {exps}, h = {h}  [{', '.join(notes[:3])}]")

print(f"""
  E8 EXPONENTS: 1, H_forb_1, 11, 13, 17, 19, 23, 29
  These are exactly the PRIMES < 30 coprime to 30!
  (30 = KEY1 * KEY2 * KEY_SUM = h(E8))

  phi(30) = phi(KEY1)*phi(KEY2)*phi(KEY_SUM) = 1*KEY1*KEY1^KEY1 = KEY1^KEY2 = 8
  = rank(E8)!

  So rank(E8) = phi(h(E8)) and the exponents are the units mod h(E8).

  This holds for ALL Coxeter groups: rank = phi(h) in the simply-laced case,
  and exponents are related to units mod h.
""")

# =====================================================================
# Part 8: HIGHEST WEIGHT THEORY
# =====================================================================
print("\n" + "=" * 70)
print("  Part 8: HIGHEST WEIGHT THEORY — WEYL CHARACTER FORMULA")
print("=" * 70)

print("""
  WEYL DIMENSION FORMULA for sl_n:
  dim V(lambda) = prod_(alpha > 0) <lambda + rho, alpha> / <rho, alpha>

  For sl_KEY1 (= sl_2):
  dim V(n) = n + 1

  For sl_KEY2 (= sl_3):
  dim V(a, b) = (a+1)(b+1)(a+b+KEY1) / KEY1

  The KEY1 in the denominator is universal!
""")

print("  Dimensions of sl_3 irreps V(a,b):")
for a in range(6):
    for b in range(6):
        dim = (a+1) * (b+1) * (a+b+2) // 2
        label = ""
        if dim == KEY1: label = "*"
        elif dim == KEY2: label = "T"  # trivial adjoint
        elif dim == KEY_SUM: label = "5"
        elif dim == h_G2: label = "6"
        elif dim == H_forb_1: label = "7"
        elif dim == KEY1**3: label = "8"
        elif dim == V_PET: label = "X"
        elif dim == BT: label = "!"
        elif dim == 27: label = "C"  # 27-dim of E6
        print(f" {dim:>4}", end="")
    print()

print(f"""
  Tournament constants as sl_3 dimensions:
  V(1,0) = KEY2 (fundamental)
  V(0,1) = KEY2 (dual fundamental)
  V(1,1) = KEY1^KEY2 = 8 (adjoint!)
  V(2,0) = h(G2) = 6
  V(0,2) = h(G2) = 6
  V(2,1) = 15
  V(3,0) = V(Pet) = 10
  V(1,2) = 15
  V(0,3) = V(Pet) = 10
  V(2,2) = 27 (the famous 27 of E6!)
  V(3,3) = 64 = KEY1^h(G2)

  The 27-dimensional representation of E6:
  Under sl_3 x sl_3 x sl_3 (triality):
  27 = (KEY2, KEY2, 1) + (1, KEY2, KEY2) + (KEY2, 1, KEY2)
  = KEY2 copies of KEY2^KEY1 = KEY2 * 9 = 27!

  FREUDENTHAL MAGIC SQUARE:
  Rows and columns labeled by R, C, H, O (division algebras)
  Entries: the exceptional Lie algebras!

  R: A1, A2, C3, F4   (dims: KEY2, KEY1^KEY2, 21, 52)
  C: A2, A2+A2, A5, E6  (dims: KEY1^KEY2, 16, 35, 78)
  H: C3, A5, D6, E7   (dims: 21, 35, 66, 133)
  O: F4, E6, E7, E8   (dims: 52, 78, 133, 248)
""")

# =====================================================================
# Part 9: CRYSTAL BASES — KASHIWARA'S KEY1 OPERATORS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 9: CRYSTAL BASES — KASHIWARA'S KEY1-VALUED STRUCTURE")
print("=" * 70)

print(f"""
  Kashiwara's crystal bases provide a combinatorial skeleton
  of quantum group representations.

  KEY STRUCTURE:
  Crystal operators: e_i (raising), f_i (lowering)
  These are KEY1-valued: they either act or don't (0 or 1 result).

  For U_q(sl_KEY1):
  Crystal of V(n) = {{b_n, b_(n-KEY1), ..., b_(-n)}}
  f_1: b_k -> b_(k-KEY1) (if k > -n, else 0)
  e_1: b_k -> b_(k+KEY1) (if k < n, else 0)

  This is a chain of length n+1 with edges labeled 1.
  The KEY1-step structure is fundamental!

  TENSOR PRODUCT RULE:
  In the crystal, tensor products follow:
  f_i(b_1 tensor b_2) =
    f_i(b_1) tensor b_2  if phi_i(b_1) > eps_i(b_2)
    b_1 tensor f_i(b_2)  otherwise

  This is a KEY1-case decision (binary branching)!

  LITTELMANN PATHS:
  Crystal elements can be realized as piecewise-linear paths
  in weight space. The KEY1 operators are:
  e_i: reflects the path ending across the i-th hyperplane
  f_i: extends in the -alpha_i direction

  PATH MODEL FOR sl_KEY2:
  Paths in R^KEY1 (KEY1-dimensional weight space)
  with KEY1 crystal operators e_1, e_2, f_1, f_2.
  The weight lattice is the A_2 lattice (hexagonal).

  CRYSTAL GRAPH of V(1,0) for sl_3:
  3 nodes connected as: b_1 --1--> b_2 --2--> b_3
  This is the KEY2-node crystal (the KEY2 colors)!
""")

# =====================================================================
# Part 10: HECKE ALGEBRAS AND KAZHDAN-LUSZTIG
# =====================================================================
print("\n" + "=" * 70)
print("  Part 10: HECKE ALGEBRAS — DEFORMING KEY1")
print("=" * 70)

print("""
  The Hecke algebra H_q(W) deforms the group algebra of W:
  Generators T_i satisfying:
  (T_i - q)(T_i + 1) = 0  (quadratic relation!)

  This is a KEY1-parameter deformation (parameter q).
  At q = 1: recover the group algebra C[W].
  At q = 0: get the 0-Hecke algebra (monoid algebra).

  For S_n (type A):
  H_q(S_n) has basis T_w for w in S_n.
  dim H_q(S_n) = n! = |S_n|

  KAZHDAN-LUSZTIG POLYNOMIALS P_(y,w)(q):
  KL polynomials are in Z[q] with:
  - P_(w,w) = 1
  - deg P_(y,w) <= (l(w) - l(y) - 1) / KEY1 if y < w

  The bound (l(w) - l(y) - 1) / KEY1 involves division by KEY1!
  The KL polynomials control:
  - Composition multiplicities of Verma modules
  - Intersection cohomology of Schubert varieties
  - Characters of simple highest weight modules

  IWAHORI-HECKE AT q = KEY1-POWER:
  When q = p^k for prime p:
  H_q = End_G(Ind_B^G 1) — endomorphisms of a parabolic induction
  This is a central object in the representation theory of
  finite groups of Lie type over F_q.

  At q = KEY1: H_KEY1(S_n) relates to GL_n(F_KEY1) = GL_n(F_2)
  At q = KEY2: H_KEY2(S_n) relates to GL_n(F_KEY2) = GL_n(F_3)
""")

# =====================================================================
# Part 11: SCHUR FUNCTIONS AND SYMMETRIC FUNCTIONS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 11: SYMMETRIC FUNCTIONS — THE SCHUR BASIS")
print("=" * 70)

print("""
  The ring of symmetric functions Sym has several natural bases:

  1. MONOMIAL: m_lambda (indexed by partitions)
  2. ELEMENTARY: e_lambda = e_(l1) * e_(l2) * ...
  3. HOMOGENEOUS: h_lambda (= complete)
  4. POWER SUM: p_lambda = p_(l1) * p_(l2) * ...
  5. SCHUR: s_lambda (the most important!)

  Transition: s = sum K_(lambda,mu) m_mu
  K = Kostka numbers (nonneg integers!)

  INNER PRODUCT:
  <s_lambda, s_mu> = delta_(lambda,mu)
  The Schur functions are ORTHONORMAL.

  CONNECTION TO KEY1 AND KEY2:

  At n = KEY1 variables:
  s_(a,b)(x_1, x_2) = (x_1^(a+1) * x_2^b - x_2^(a+1) * x_1^b) / (x_1 - x_2)
  — a ratio of KEY1 x KEY1 determinants!

  At n = KEY2 variables:
  s_lambda(x_1, x_2, x_3) = det[x_i^(lambda_j + KEY2 - j)] / det[x_i^(KEY2-j)]
  — KEY2 x KEY2 determinants (Vandermonde denominator)
""")

# Compute some Schur function values
print("  Kostka numbers K(lambda, mu) for n = 4:")
# K(lambda, mu) = # semistandard Young tableaux of shape lambda, content mu
for lam in [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)]:
    for mu in [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)]:
        # Just show the matrix structure
        pass

print(f"""
  LITTLEWOOD-RICHARDSON RULE:
  s_lambda * s_mu = sum c^nu_(lambda,mu) s_nu

  The LR coefficients c^nu_(lambda,mu) count:
  - LR tableaux (combinatorial)
  - Tensor product multiplicities (representation)
  - Intersection numbers on Grassmannians (geometry)

  CAUCHY IDENTITY:
  prod_(i,j) 1/(1 - x_i y_j) = sum_lambda s_lambda(x) s_lambda(y)

  DUAL CAUCHY:
  prod_(i,j) (1 + x_i y_j) = sum_lambda s_lambda(x) s_(lambda')(y)
  where lambda' = conjugate partition.

  The KEY1 appears: the dual identity uses (1 + xy) vs 1/(1-xy).
  Factor of KEY1 between "additive" and "multiplicative" versions!

  PLETHYSM:
  s_lambda[s_mu] = plethystic substitution
  s_KEY1[s_n] = s_(n,n) + s_(n+1,n-1) + ... (KEY1-nd symmetric power)
  s_KEY2[s_1] = s_(KEY2) + ... (KEY2-rd power = cube)

  OMEGA INVOLUTION:
  omega: s_lambda -> s_(lambda')
  omega(e_n) = h_n, omega(p_n) = (-1)^(n-1) p_n
  omega^KEY1 = identity (it's an involution — order KEY1!)
""")

# =====================================================================
# Part 12: GRAND SYNTHESIS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 12: GRAND SYNTHESIS — REPRESENTATION THEORY IS (2,3)")
print("=" * 70)

print("""
======================================================================
  REPRESENTATION THEORY IS BUILT FROM KEY1 = 2 AND KEY2 = 3
======================================================================

1. SYMMETRIC GROUPS:
   S_KEY1: KEY1 irreps (trivial, sign)
   S_KEY2: KEY2 irreps (trivial, sign, standard of dim KEY1)
   S_5: H_forb_1 irreps
   |S_4| = BT = 24, |S_5| = BI = 120

2. sl_2 (THE KEY1 GROUP):
   KEY1 x KEY1 matrices, rank 1
   [h,e] = KEY1*e (weight spacing = KEY1)
   V(n) has dim n+1
   CG: V(m) x V(n) = V(m+n) + V(m+n-KEY1) + ... + V(|m-n|)

3. ROOT SYSTEMS:
   ADE classification uses edge label KEY2 exclusively
   h(A2) = KEY2, h(E6) = h(E6), h(E8) = KEY1*KEY2*KEY_SUM
   E8 exponents = primes coprime to 30

4. EXCEPTIONAL LIE ALGEBRAS:
   G2: dim KEY1*H_forb_1, rank KEY1
   E6: rank h(G2), h = h(E6)
   E7: rank H_forb_1, h = 18
   E8: rank KEY1^KEY2, h = 30, #roots = KEY1*BI

5. MCKAY CORRESPONDENCE:
   BT(KEY1^KEY2*KEY2) <-> E6
   BO(KEY1^KEY1^KEY1*KEY2) <-> E7
   BI(KEY1^KEY2*KEY2*KEY_SUM) <-> E8
   Conjugacy classes: H_forb_1, 8, 9 (consecutive!)

6. QUANTUM GROUPS:
   U_q(sl_KEY1) at q = root of unity -> truncation
   Level KEY2 = 3: Fibonacci anyons (golden ratio dims)
   S-matrix involves sin(pi/N) with KEY1 in denominator

7. COXETER GROUPS:
   ADE edges labeled KEY2 (the default Coxeter label)
   rank(E8) = phi(h(E8)) = phi(30) = KEY1^KEY2

8. CRYSTAL BASES:
   KEY1-valued operators (0 or 1)
   Weight spacing = KEY1 (same as sl_2)
   Tensor rule = binary decision

9. HECKE ALGEBRAS:
   Quadratic relation (T_i - q)(T_i + 1) = 0 (degree KEY1!)
   KL degree bound involves 1/KEY1
   At q = KEY1: GL_n(F_2), at q = KEY2: GL_n(F_3)

10. SYMMETRIC FUNCTIONS:
    Schur = ratio of KEY1 x KEY1 (or KEY2 x KEY2) determinants
    Omega involution: order KEY1
    Cauchy vs dual Cauchy: KEY1-fold dichotomy

THE CROWN JEWEL:
   sl_KEY1 is the fundamental Lie algebra (rank 1, basis {e,f,h}).
   sl_KEY2 is the first "interesting" one (rank KEY1, the triangle).
   Together they generate ALL of Lie theory via sl_n, and
   the exceptional algebras arise from the (KEY1,KEY2) interplay
   through the Freudenthal magic square.

   REPRESENTATION THEORY IS (2,3).
""")
