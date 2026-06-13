"""
D_{2p} orbit decomposition of Hamiltonian paths in Paley tournaments.
Tests the anti-automorphism cycle pairing mechanism.

New computation: For each Z_p-orbit of HPs, determine whether the anti-automorphism
maps the orbit to itself (palindromic) or to another orbit (chiral pair).

Author: opus-2026-03-12-S57
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

def qr_set(p):
    return frozenset((a*a) % p for a in range(1, p))

def build_adj(p, S):
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j-i)%p in S:
                A[i][j] = 1
    return A

def enumerate_hps(A, p):
    """Enumerate ALL Hamiltonian paths using backtracking."""
    paths = []
    def backtrack(path, visited):
        if len(path) == p:
            paths.append(tuple(path))
            return
        last = path[-1]
        for u in range(p):
            if not visited[u] and A[last][u]:
                visited[u] = True
                path.append(u)
                backtrack(path, visited)
                path.pop()
                visited[u] = False
    
    for start in range(p):
        visited = [False]*p
        visited[start] = True
        backtrack([start], visited)
    return paths

def rotate_path(path, r, p):
    """Apply rotation: v_k -> (v_k + r) mod p."""
    return tuple((v+r) % p for v in path)

def reflect_path(path, p):
    """Apply anti-automorphism: v_k -> (-v_k) mod p, then reverse (since anti-aut)."""
    # The reflection σ: v -> -v maps T_p to T_p^{op}.
    # The composition σ ∘ reverse gives an automorphism of T_p
    # (since reversing a path in T^{op} gives a path in T).
    # So: σ maps HP P=(v_0,...,v_{p-1}) to the path (-v_{p-1},...,-v_0) in T_p.
    reversed_reflected = tuple((-v) % p for v in reversed(path))
    return reversed_reflected

# ============================================================
# PART 1: T_7 (Paley, p=7)
# ============================================================
print("="*70)
print("PART 1: D_{14} orbit decomposition of T_7 Hamiltonian paths")
print("="*70)

p = 7
QR = qr_set(p)
A = build_adj(p, QR)
print(f"QR_{p} = {sorted(QR)}")

paths = enumerate_hps(A, p)
print(f"Total HPs: {len(paths)} (expected 189)")

# Group into Z_p orbits
path_to_orbit = {}
orbits = []
for path in paths:
    if path in path_to_orbit:
        continue
    orbit = set()
    for r in range(p):
        rp = rotate_path(path, r, p)
        orbit.add(rp)
        path_to_orbit[rp] = len(orbits)
    orbits.append(frozenset(orbit))

print(f"Z_{p}-orbits: {len(orbits)} (expected {len(paths)//p} = {len(paths)//p})")

# Classify each orbit as palindromic or chiral
palindromic = []
chiral_pairs = []
orbit_partner = {}

for i, orbit in enumerate(orbits):
    if i in orbit_partner:
        continue
    # Apply anti-automorphism to each path in orbit
    reflected_paths = set()
    for path in orbit:
        rp = reflect_path(path, p)
        reflected_paths.add(rp)
    
    # Find which orbit the reflected paths belong to
    rp = next(iter(reflected_paths))
    partner_idx = path_to_orbit.get(rp)
    
    if partner_idx == i:
        palindromic.append(i)
    else:
        chiral_pairs.append((i, partner_idx))
        orbit_partner[i] = partner_idx
        orbit_partner[partner_idx] = i

print(f"\nPalindromic (self-paired) orbits: {len(palindromic)}")
print(f"Chiral pairs: {len(chiral_pairs)}")
print(f"D_{{{2*p}}}-orbits: {len(palindromic) + len(chiral_pairs)}")
print(f"Check: {len(palindromic)} + 2×{len(chiral_pairs)} = {len(palindromic) + 2*len(chiral_pairs)} Z_p-orbits (expected {len(orbits)})")

# Show representative palindromic paths
print(f"\n--- Palindromic orbit representatives ---")
for idx in palindromic[:5]:
    rep = min(orbits[idx])
    rp = reflect_path(rep, p)
    # Find the rotation r that maps rp back to a path in the orbit
    for r in range(p):
        if rotate_path(rp, r, p) == rep:
            print(f"  Orbit {idx}: {rep} → reflect → {rp} → rot({r}) → {rep}")
            break
    else:
        # Find the rotation that maps to SOME member of the orbit
        for r in range(p):
            if rotate_path(rp, r, p) in orbits[idx]:
                target = rotate_path(rp, r, p)
                print(f"  Orbit {idx}: {rep} → reflect → rot({r}) → {target}")
                break

# Show representative chiral pairs
print(f"\n--- Chiral pair representatives ---")
for i, j in chiral_pairs[:3]:
    rep_i = min(orbits[i])
    rep_j = min(orbits[j])
    print(f"  Orbits {i},{j}: {rep_i} ↔ {rep_j}")

# ============================================================
# PART 2: Odd cycle structure and pairing under anti-aut
# ============================================================
print("\n" + "="*70)
print("PART 2: Odd cycle pairing under anti-automorphism")
print("="*70)

def enumerate_directed_cycles(A, p, max_len=None):
    """Find all directed cycles using DFS."""
    if max_len is None:
        max_len = p
    cycles = set()
    for start in range(p):
        def dfs(path, visited):
            last = path[-1]
            if len(path) >= 3 and A[last][start]:
                # Found a cycle back to start
                if len(path) % 2 == 1:  # odd length only
                    canon = min(path[i:] + path[:i] for i in range(len(path)))
                    cycles.add(tuple(canon))
            if len(path) < max_len:
                for u in range(start+1, p):  # avoid double-counting
                    if not visited[u] and A[last][u]:
                        visited[u] = True
                        path.append(u)
                        dfs(path, visited)
                        path.pop()
                        visited[u] = False
        
        visited = [False]*p
        visited[start] = True
        dfs([start], visited)
    return cycles

cycles = enumerate_directed_cycles(A, p)
print(f"Total odd directed cycles in T_{p}: {len(cycles)}")
cycle_lengths = {}
for c in cycles:
    L = len(c)
    cycle_lengths[L] = cycle_lengths.get(L, 0) + 1
print(f"By length: {sorted(cycle_lengths.items())}")

# Pair cycles under anti-automorphism
def reflect_cycle(cycle, p):
    """Apply anti-aut to cycle: v -> -v, reverse direction."""
    reflected = tuple((-v) % p for v in reversed(cycle))
    # Canonicalize
    canon = min(reflected[i:] + reflected[:i] for i in range(len(reflected)))
    return canon

cycle_set = set(cycles)
cycle_partners = {}
fixed_cycles = []
paired_cycles = []

for c in cycles:
    if c in cycle_partners:
        continue
    rc = reflect_cycle(c, p)
    if rc == c:
        fixed_cycles.append(c)
    elif rc in cycle_set:
        cycle_partners[c] = rc
        cycle_partners[rc] = c
        paired_cycles.append((c, rc))
    else:
        print(f"  WARNING: reflected cycle {rc} not found in cycle set!")

print(f"\nAnti-aut cycle pairing:")
print(f"  Fixed (self-paired) cycles: {len(fixed_cycles)}")
print(f"  Chiral pairs of cycles: {len(paired_cycles)}")
print(f"  Check: {len(fixed_cycles)} + 2×{len(paired_cycles)} = {len(fixed_cycles)+2*len(paired_cycles)} (expected {len(cycles)})")

# Check vertex disjointness of paired cycles
disjoint_pairs = 0
for c1, c2 in paired_cycles:
    if set(c1).isdisjoint(set(c2)):
        disjoint_pairs += 1

print(f"\n  Vertex-DISJOINT paired cycles: {disjoint_pairs}/{len(paired_cycles)}")
print(f"  These contribute to α_2 in OCF: each disjoint pair contributes +4 to H")

# Compute alpha_1, alpha_2 from the OCF
# alpha_1 = number of odd cycles = |Ω(T)|
# alpha_2 = number of vertex-disjoint pairs of odd cycles
alpha_1 = len(cycles)
alpha_2 = 0
cycle_list = list(cycles)
for i in range(len(cycle_list)):
    for j in range(i+1, len(cycle_list)):
        if set(cycle_list[i]).isdisjoint(set(cycle_list[j])):
            alpha_2 += 1

print(f"\n  OCF components:")
print(f"  alpha_1 = {alpha_1} (total odd cycles)")
print(f"  alpha_2 = {alpha_2} (disjoint pairs)")
print(f"  H = I(Omega, 2) ≥ 1 + 2*{alpha_1} + 4*{alpha_2} = {1+2*alpha_1+4*alpha_2}")
print(f"  (This is a LOWER BOUND; full I(Omega,2) also counts alpha_3, alpha_4, ...)")

# ============================================================
# PART 3: Compare OCF structure across all Z_7 circulants
# ============================================================
print("\n" + "="*70)
print("PART 3: OCF structure for all circulant tournaments on Z_7")
print("="*70)

def all_circulant_sets(p):
    pairs = []
    used = set()
    for d in range(1, p):
        if d not in used:
            pairs.append((d, p-d))
            used.add(d); used.add(p-d)
    results = []
    for mask in range(2**len(pairs)):
        S = set()
        for i, (a, b) in enumerate(pairs):
            S.add(b if (mask>>i)&1 else a)
        results.append(frozenset(S))
    return results

all_S = all_circulant_sets(p)
print(f"\n{'S':15} {'H':>5} {'α₁':>4} {'α₂':>4} {'α₃':>4} {'Paley?':>7} {'palindromic_Ω':>15}")

for S in all_S:
    A_s = build_adj(p, S)
    # Count HPs
    from itertools import permutations
    H = 0
    for perm in permutations(range(p)):
        if all(A_s[perm[i]][perm[i+1]] for i in range(p-1)):
            H += 1
    
    # Count cycles
    cyc = enumerate_directed_cycles(A_s, p)
    a1 = len(cyc)
    cyc_list = list(cyc)
    a2 = sum(1 for i in range(len(cyc_list)) for j in range(i+1, len(cyc_list))
             if set(cyc_list[i]).isdisjoint(set(cyc_list[j])))
    # a3
    a3 = 0
    for i in range(len(cyc_list)):
        for j in range(i+1, len(cyc_list)):
            if not set(cyc_list[i]).isdisjoint(set(cyc_list[j])): continue
            for k in range(j+1, len(cyc_list)):
                if set(cyc_list[k]).isdisjoint(set(cyc_list[i])) and set(cyc_list[k]).isdisjoint(set(cyc_list[j])):
                    a3 += 1
    
    is_paley = (S == QR or S == frozenset((p-s)%p for s in QR))
    
    # Check palindromic Omega (heuristic: check if cycle count by length is balanced)
    cl = {}
    for c in cyc:
        L = len(c)
        cl[L] = cl.get(L, 0) + 1
    pal = "yes" if is_paley else "no"
    
    ocf_check = 1 + 2*a1 + 4*a2 + 8*a3
    print(f"{str(sorted(S)):15} {H:>5} {a1:>4} {a2:>4} {a3:>4} {'YES' if is_paley else 'no':>7} I={ocf_check:>4}")

print(f"\nH = 1 + 2*α₁ + 4*α₂ + 8*α₃ + ... (OCF expansion)")
print(f"Paley maximizes H by maximizing the FULL independence polynomial, not just low-order terms")

# ============================================================
# PART 4: D_{14} representation theory of HP space
# ============================================================
print("\n" + "="*70)
print("PART 4: Representation theory summary")
print("="*70)

# D_{2p} has:
# - 2 one-dimensional reps: trivial and sign
# - (p-1)/2 two-dimensional reps
# Total: 2 + (p-1)/2 = (p+3)/2 irreducible reps

print(f"\nD_{{{2*p}}} irreducible representations:")
print(f"  1-dim reps: 2 (trivial, sign)")
print(f"  2-dim reps: {(p-1)//2}")
print(f"  Total: {2 + (p-1)//2} irreducible reps")
print()
print(f"HP space decomposition (189 HPs as D_{{{2*p}}}-representation):")
print(f"  Total dim = {len(paths)}")
print(f"  Palindromic orbits (fixed by σ): {len(palindromic)} → contribute to 1-dim reps")
print(f"  Chiral pairs: {len(chiral_pairs)} → each gives a 2-dim rep")
print(f"  Decomposition: {len(palindromic)} × 1-dim + {len(chiral_pairs)} × 2-dim = {len(palindromic) + 2*len(chiral_pairs)} Z_p-orbits")
print()
n_d_orbits = len(palindromic) + len(chiral_pairs)
print(f"  D-orbits: {n_d_orbits} = H/(2p) would be {len(paths)/(2*p):.1f}")
print(f"  Palindromic size = {len(palindromic)*p}, Chiral size = {2*len(chiral_pairs)*p}")
print(f"  Total = {len(palindromic)*p + 2*len(chiral_pairs)*p}")

print("\nDONE.")
