"""
Decompose H(T_7) = 189 Hamiltonian paths into irreducible representations of D_14.

T_7 = Paley tournament on Z_7, connection set QR_7 = {1, 2, 4}.
D_14 = <r, s | r^7 = s^2 = 1, srs = r^{-1}> acts on T_7.
  r: v -> v+1 (mod 7)  [automorphism]
  s: v -> -v (mod 7)   [anti-automorphism, so reverses + negates]

D_14 irreps:
  rho_0: trivial (1-dim)
  rho_1: sign repr, r->1, s->-1 (1-dim)
  rho_k (k=1,2,3): 2-dim, r -> rotation by 2*pi*k/7, s -> reflection

Author: opus session script
"""

import numpy as np
from itertools import permutations
from collections import defaultdict

# === Step 1: Define T_7 and enumerate all 189 HPs ===

N = 7
QR = {1, 2, 4}  # quadratic residues mod 7

def has_edge(u, v):
    """T_7 has edge u->v iff (v-u) mod 7 in QR."""
    return (v - u) % N in QR

# Adjacency matrix
adj = [[has_edge(i, j) for j in range(N)] for i in range(N)]

def enumerate_hps():
    """Enumerate all Hamiltonian paths of T_7 by brute force."""
    hps = []
    for perm in permutations(range(N)):
        valid = True
        for i in range(N - 1):
            if not adj[perm[i]][perm[i+1]]:
                valid = False
                break
        if valid:
            hps.append(perm)
    return hps

print("=" * 70)
print("D_14 Irreducible Decomposition of H(T_7)")
print("=" * 70)

print("\nStep 1: Enumerating Hamiltonian paths of T_7 (Paley tournament on Z_7)...")
print(f"  QR_7 = {QR}")
print(f"  Edge rule: u -> v iff (v - u) mod 7 in QR_7")

all_hps = enumerate_hps()
print(f"  Found {len(all_hps)} Hamiltonian paths.")
assert len(all_hps) == 189, f"Expected 189, got {len(all_hps)}"

# Index HPs for fast lookup
hp_set = set(all_hps)
hp_to_idx = {hp: i for i, hp in enumerate(all_hps)}

# === Step 2: Group action on HPs ===

def rotate_hp(hp, k):
    """Apply r^k: v -> v + k (mod 7) to an HP."""
    return tuple((v + k) % N for v in hp)

def reflect_hp(hp):
    """Apply s: reverse the path and negate each vertex.
    Given P = (v_0, ..., v_6), s(P) = ((-v_6) mod 7, ..., (-v_0) mod 7).
    This works because s is an anti-automorphism (reverses edges) and T_7 = T_7^op."""
    return tuple((-v) % N for v in reversed(hp))

# Verify that rotation and reflection map HPs to HPs
print("\nStep 2: Verifying D_14 action on HPs...")
for hp in all_hps:
    for k in range(N):
        rhp = rotate_hp(hp, k)
        assert rhp in hp_set, f"Rotation of {hp} by {k} gives {rhp}, not an HP!"
    shp = reflect_hp(hp)
    assert shp in hp_set, f"Reflection of {hp} gives {shp}, not an HP!"
print("  All rotations and reflections map HPs to HPs. Confirmed.")

# === Step 3: Compute orbits under Z_7 rotations ===

print("\nStep 3: Computing Z_7 orbits...")
visited = set()
orbits = []
for hp in all_hps:
    if hp in visited:
        continue
    orbit = []
    for k in range(N):
        rhp = rotate_hp(hp, k)
        orbit.append(rhp)
        visited.add(rhp)
    orbits.append(orbit)

print(f"  Found {len(orbits)} orbits of size 7 each.")
print(f"  Total HPs in orbits: {sum(len(o) for o in orbits)}")
assert len(orbits) == 27
assert all(len(o) == 7 for o in orbits)

# === Step 4: Reflection action on orbits ===

print("\nStep 4: Analyzing reflection action on orbits...")

# Map each HP to its orbit index
hp_to_orbit = {}
for oi, orbit in enumerate(orbits):
    for hp in orbit:
        hp_to_orbit[hp] = oi

# For each orbit, find where reflection sends it
orbit_self_paired = []  # orbits mapped to themselves by s
orbit_pairs = []        # pairs of orbits swapped by s
paired = set()

for oi, orbit in enumerate(orbits):
    if oi in paired:
        continue
    rep = orbit[0]
    srep = reflect_hp(rep)
    target_orbit = hp_to_orbit[srep]
    if target_orbit == oi:
        orbit_self_paired.append(oi)
    else:
        orbit_pairs.append((oi, target_orbit))
        paired.add(oi)
        paired.add(target_orbit)

print(f"  Self-paired orbits (s maps orbit to itself): {len(orbit_self_paired)}")
print(f"  Swapped orbit pairs (s swaps two orbits): {len(orbit_pairs)}")
print(f"  Check: {len(orbit_self_paired)} + 2 * {len(orbit_pairs)} = {len(orbit_self_paired) + 2 * len(orbit_pairs)} orbits")

# For self-paired orbits, find the shift: s sends orbit[0] to orbit[j] for some j
print("\n  Self-paired orbit details:")
for oi in orbit_self_paired:
    orbit = orbits[oi]
    rep = orbit[0]
    srep = reflect_hp(rep)
    # Find which rotation of rep gives srep
    for k in range(N):
        if rotate_hp(rep, k) == srep:
            print(f"    Orbit {oi}: rep={rep}, s(rep) = r^{k}(rep)")
            break

# === Step 5: Character computation ===

print("\n" + "=" * 70)
print("Step 5: Computing the character of the permutation representation")
print("=" * 70)

# The representation is on C^189 (permutation representation).
# Character chi(g) = number of HPs fixed by g.

def apply_group_element(hp, k, has_s):
    """Apply r^k (if has_s=False) or s*r^k (if has_s=True) to an HP."""
    result = rotate_hp(hp, k)
    if has_s:
        result = reflect_hp(result)
    return result

# Compute character for each conjugacy class of D_14
# Conjugacy classes of D_14 (n=7 odd):
#   {e}, {r, r^6}, {r^2, r^5}, {r^3, r^4}, {s, sr, sr^2, ..., sr^6}
# So 5 conjugacy classes.

conj_classes = [
    ("e", [(0, False)]),
    ("r,r^6", [(1, False), (6, False)]),
    ("r^2,r^5", [(2, False), (5, False)]),
    ("r^3,r^4", [(3, False), (4, False)]),
    ("s,sr,...,sr^6", [(k, True) for k in range(N)]),
]

print("\nConjugacy classes of D_14:")
characters = {}
for name, elements in conj_classes:
    # Count fixed points for one representative
    k, has_s = elements[0]
    fixed = sum(1 for hp in all_hps if apply_group_element(hp, k, has_s) == hp)
    characters[name] = fixed
    print(f"  {name:20s}: |class| = {len(elements):2d}, chi = {fixed}")

print(f"\n  Verification: sum of multiplicities * irrep_dim should equal 189")
print(f"  Burnside orbit count: (1/|D_14|) * sum |class|*chi = ", end="")
total = sum(len(elts) * characters[name] for name, elts in conj_classes)
print(f"(1/14) * {total} = {total / 14:.0f} D_14-orbits")

# === Step 6: Decompose into irreps ===

print("\n" + "=" * 70)
print("Step 6: Decomposition into irreducible representations")
print("=" * 70)

# Irreducible characters of D_14 (order 14, n=7 odd)
# Classes: e, {r,r^6}, {r^2,r^5}, {r^3,r^4}, {s,...,sr^6}
# Class sizes: 1, 2, 2, 2, 7

class_sizes = [1, 2, 2, 2, 7]
class_names = ["e", "r,r^6", "r^2,r^5", "r^3,r^4", "s,sr,...,sr^6"]

# Irrep characters (rows = irreps, cols = conjugacy classes)
# rho_0 (trivial): 1, 1, 1, 1, 1
# rho_1 (sign): 1, 1, 1, 1, -1
# rho_k (k=1,2,3, 2-dim): 2, 2*cos(2*pi*k/7), 2*cos(4*pi*k/7), 2*cos(6*pi*k/7), 0

omega = 2 * np.pi / 7
irrep_chars = {
    "rho_0 (trivial, dim 1)": [1, 1, 1, 1, 1],
    "rho_1 (sign, dim 1)":    [1, 1, 1, 1, -1],
}
for k in range(1, 4):
    chars = [2,
             2 * np.cos(k * omega),
             2 * np.cos(2 * k * omega),
             2 * np.cos(3 * k * omega),
             0]
    irrep_chars[f"rho_{k+1} (dim 2, k={k})"] = chars

# Build chi vector from our representation
chi_vec = [characters[name] for name, _ in conj_classes]
print(f"\nPermutation character: {chi_vec}")

print("\nIrrep character table of D_14:")
print(f"  {'Irrep':30s} | {'e':>4s} | {'r,r6':>8s} | {'r2,r5':>8s} | {'r3,r4':>8s} | {'s..sr6':>8s}")
print(f"  {'-'*30}-+-{'-'*4}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}")
for name, chars in irrep_chars.items():
    cstrs = [f"{c:8.4f}" for c in chars]
    print(f"  {name:30s} | {cstrs[0]} | {cstrs[1]} | {cstrs[2]} | {cstrs[3]} | {cstrs[4]}")

print("\nMultiplicity computation: m_i = (1/|G|) * sum_C |C| * chi_V(C) * chi_i(C)*")
print()

G_order = 14
multiplicities = {}
total_dim = 0

for irr_name, irr_chars in irrep_chars.items():
    m = 0
    detail_parts = []
    for j in range(5):
        contrib = class_sizes[j] * chi_vec[j] * irr_chars[j]
        m += contrib
        detail_parts.append(f"{class_sizes[j]}*{chi_vec[j]}*{irr_chars[j]:.4f}")
    m /= G_order
    multiplicities[irr_name] = m
    dim = 1 if "dim 1" in irr_name else 2
    total_dim += round(m) * dim
    print(f"  {irr_name}:")
    print(f"    m = (1/{G_order}) * ({' + '.join(detail_parts)})")
    print(f"    m = (1/{G_order}) * {m * G_order:.4f} = {m:.4f}")
    print(f"    Multiplicity: {round(m)}")
    print()

print(f"  Dimension check: ", end="")
parts = []
for irr_name, m in multiplicities.items():
    dim = 1 if "dim 1" in irr_name else 2
    parts.append(f"{round(m)}*{dim}")
print(" + ".join(parts) + f" = {total_dim}")
assert total_dim == 189, f"Dimension mismatch: {total_dim} != 189"
print(f"  Matches H(T_7) = 189. Confirmed!")

# === Step 7: Interpretation ===

print("\n" + "=" * 70)
print("Step 7: Interpretation")
print("=" * 70)

m_trivial = round(multiplicities["rho_0 (trivial, dim 1)"])
m_sign = round(multiplicities["rho_1 (sign, dim 1)"])
print(f"\n  Trivial irrep multiplicity: {m_trivial}")
print(f"    -> There are {m_trivial} D_14-invariant linear combinations of HPs.")
print(f"    -> These are the HP combinations invariant under BOTH rotation AND reflection.")
print(f"\n  Sign irrep multiplicity: {m_sign}")
print(f"    -> There are {m_sign} combinations that are rotation-invariant but flip sign under reflection.")

for k in range(1, 4):
    name = f"rho_{k+1} (dim 2, k={k})"
    m = round(multiplicities[name])
    print(f"\n  {name}: multiplicity {m}")
    print(f"    -> {m} copies of the 2-dim representation with rotation eigenvalue e^{{2*pi*i*{k}/7}}")

# === Step 8: Additional analysis — orbit structure ===

print("\n" + "=" * 70)
print("Step 8: Orbit structure details")
print("=" * 70)

print(f"\n  Total Z_7 orbits: {len(orbits)}")
print(f"  Self-paired under s: {len(orbit_self_paired)}")
print(f"  Swapped pairs: {len(orbit_pairs)}")

# For self-paired orbits, determine the fixed-point structure under s
print("\n  Self-paired orbits and their s-action:")
for oi in orbit_self_paired:
    orbit = orbits[oi]
    # Count fixed points of s within this orbit
    fixed_in_orbit = sum(1 for hp in orbit if reflect_hp(hp) == hp)
    # Find the shift parameter
    rep = orbit[0]
    srep = reflect_hp(rep)
    for k in range(N):
        if rotate_hp(rep, k) == srep:
            shift = k
            break
    print(f"    Orbit {oi}: shift={shift}, fixed by s: {fixed_in_orbit}, rep={orbit[0]}")

# Check: total fixed by s
total_fixed_s = characters["s,sr,...,sr^6"]
print(f"\n  Total HPs fixed by s (identity element of reflection class): {total_fixed_s}")

# Show a few example orbits
print("\n  Example orbits:")
for i in range(min(5, len(orbits))):
    print(f"    Orbit {i}: {orbits[i][0]} -> ... (size {len(orbits[i])})")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"\n  H(T_7) = 189 decomposes under D_14 as:")
parts_summary = []
for irr_name, m in multiplicities.items():
    mi = round(m)
    short = irr_name.split("(")[0].strip()
    parts_summary.append(f"{mi}*{short}")
print(f"    {' + '.join(parts_summary)}")
print(f"\n  Or in terms of dimensions:")
dim_parts = []
for irr_name, m in multiplicities.items():
    mi = round(m)
    dim = 1 if "dim 1" in irr_name else 2
    dim_parts.append(f"{mi}x{dim}")
print(f"    {' + '.join(dim_parts)} = {total_dim}")
print()
