"""
Doubling/Blowup study for tournaments.
Computes H values, SC status, recursive patterns for both Lex and SC blowups.
"""

import itertools
from collections import defaultdict

# ============================================================
# Core tournament utilities
# ============================================================

def make_tournament(n, arcs):
    """arcs: set of (u,v) meaning u->v. Returns adjacency matrix."""
    A = [[0]*n for _ in range(n)]
    for u,v in arcs:
        A[u][v] = 1
    return A

def adj_to_arcs(A):
    n = len(A)
    return frozenset((i,j) for i in range(n) for j in range(n) if A[i][j]==1)

def arcs_to_matrix(n, arcs):
    A = [[0]*n for _ in range(n)]
    for u,v in arcs:
        A[u][v] = 1
    return A

def score_seq(A):
    n = len(A)
    return tuple(sorted(sum(A[i]) for i in range(n)))

def complement_arcs(n, arcs):
    """Reverse all arcs (T^op)."""
    return frozenset((v,u) for u,v in arcs)

def apply_perm(n, arcs, perm):
    """Relabel vertices: vertex i -> perm[i]."""
    return frozenset((perm[u], perm[v]) for u,v in arcs)

def canonical_form(n, arcs):
    """Canonical form: min over all relabelings (brute force for small n)."""
    best = None
    for perm in itertools.permutations(range(n)):
        relabeled = apply_perm(n, arcs, perm)
        if best is None or relabeled < best:
            best = relabeled
    return best

def is_isomorphic(n, arcs1, arcs2):
    return canonical_form(n, arcs1) == canonical_form(n, arcs2)

def is_SC(n, arcs):
    """Check if tournament is self-complementary (isomorphic to its reverse)."""
    rev = complement_arcs(n, arcs)
    return is_isomorphic(n, arcs, rev)

# ============================================================
# Hamiltonian path counting
# ============================================================

def count_hamiltonian_paths(n, arcs):
    """Count directed Hamiltonian paths."""
    arc_set = set(arcs)
    count = 0
    for perm in itertools.permutations(range(n)):
        if all((perm[i], perm[i+1]) in arc_set for i in range(n-1)):
            count += 1
    return count

def count_hamiltonian_paths_fast(n, arcs):
    """DP over subsets for Hamiltonian path counting."""
    arc_set = set(arcs)
    # dp[mask][v] = number of paths using exactly the vertices in mask, ending at v
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1<<w):
                    continue
                if (v,w) in arc_set:
                    dp[mask|(1<<w)][w] += dp[mask][v]
    full = (1<<n)-1
    return sum(dp[full][v] for v in range(n))

# ============================================================
# Blowup constructions
# ============================================================

def lex_blowup(n, arcs):
    """
    T[K2]: vertex v -> {2v, 2v+1}, internal arc 2v->2v+1.
    u->v in T gives all 4 arcs: ui->vj for i,j in {0,1}.
    Returns (2n, new_arcs).
    """
    N = 2*n
    new_arcs = set()
    # Internal arcs
    for v in range(n):
        new_arcs.add((2*v, 2*v+1))
    # Inter-pair arcs
    for u,v in arcs:
        for i in range(2):
            for j in range(2):
                new_arcs.add((2*u+i, 2*v+j))
    return N, frozenset(new_arcs)

def sc_blowup(n, arcs):
    """
    T_SC: vertex v -> {2v, 2v+1}, internal arc 2v->2v+1.
    u->v in T:
      Lane arcs: u0->v0, u1->v1 (same subscript, follow T)
      Cross arcs: v0->u1, v1->u0 (mixed subscript, follow T^op)
    Returns (2n, new_arcs).
    """
    N = 2*n
    new_arcs = set()
    # Internal arcs
    for v in range(n):
        new_arcs.add((2*v, 2*v+1))
    # Inter-pair arcs
    for u,v in arcs:
        # Lane: same subscript, follow T
        new_arcs.add((2*u,   2*v))    # u0->v0
        new_arcs.add((2*u+1, 2*v+1))  # u1->v1
        # Cross: T^op direction
        new_arcs.add((2*v,   2*u+1))  # v0->u1
        new_arcs.add((2*v+1, 2*u))    # v1->u0
    return N, frozenset(new_arcs)

# ============================================================
# Enumerate all non-isomorphic tournaments for small n
# ============================================================

def all_tournaments(n):
    """Generate one representative per isomorphism class."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()
    reps = []
    for choices in itertools.product([0,1], repeat=len(pairs)):
        arcs = set()
        for k,(i,j) in enumerate(pairs):
            if choices[k]==0:
                arcs.add((i,j))
            else:
                arcs.add((j,i))
        arcs = frozenset(arcs)
        can = canonical_form(n, arcs)
        if can not in seen:
            seen.add(can)
            reps.append(arcs)
    return reps

# ============================================================
# Main computations
# ============================================================

print("="*70)
print("PART 1: H values for all iso classes at n=3,4,5")
print("="*70)

results = {}
for n in [3, 4, 5]:
    print(f"\nn={n}:")
    print(f"{'Class':<6} {'Score seq':<20} {'H(T)':<8} {'H(Lex)':<8} {'H(SC)':<8} {'T is SC':<8} {'Lex is SC':<10} {'SC is SC':<8}")
    print("-"*80)

    classes = all_tournaments(n)
    results[n] = []

    for idx, arcs in enumerate(classes):
        ss = score_seq(arcs_to_matrix(n, list(arcs) + [(j,i) for i,j in arcs if (i,j) not in arcs]))
        # Recompute score seq properly
        A = arcs_to_matrix(n, arcs)
        ss = score_seq(A)

        H_T = count_hamiltonian_paths_fast(n, arcs)

        Nl, lex_arcs = lex_blowup(n, arcs)
        H_lex = count_hamiltonian_paths_fast(Nl, lex_arcs)

        Ns, sc_arcs = sc_blowup(n, arcs)
        H_sc = count_hamiltonian_paths_fast(Ns, sc_arcs)

        t_sc = is_SC(n, arcs)
        lex_sc = is_SC(Nl, lex_arcs)
        sc_sc = is_SC(Ns, sc_arcs)

        print(f"T{idx:<5} {str(ss):<20} {H_T:<8} {H_lex:<8} {H_sc:<8} {str(t_sc):<8} {str(lex_sc):<10} {str(sc_sc):<8}")

        results[n].append({
            'idx': idx,
            'arcs': arcs,
            'ss': ss,
            'H_T': H_T,
            'H_lex': H_lex,
            'H_sc': H_sc,
            't_sc': t_sc,
            'lex_sc': lex_sc,
            'sc_sc': sc_sc,
        })

print("\n" + "="*70)
print("PART 2: Ratio analysis - looking for H(T_SC) = f(H(T))")
print("="*70)
for n in [3, 4, 5]:
    print(f"\nn={n}:")
    for r in results[n]:
        ratio = r['H_sc'] / r['H_T'] if r['H_T'] > 0 else float('inf')
        print(f"  T{r['idx']}: H(T)={r['H_T']}, H(SC)={r['H_sc']}, ratio={ratio:.4f}, H(SC)/H(T)^2={r['H_sc']/r['H_T']**2:.4f}")

print("\n" + "="*70)
print("PART 3: Circulant tournaments")
print("="*70)

def circulant_tournament(n, S):
    """Circulant tournament C_n(S): i->j iff (j-i) mod n in S."""
    arcs = set()
    for i in range(n):
        for s in S:
            arcs.add((i, (i+s)%n))
    return frozenset(arcs)

def is_circulant(n, arcs):
    """Check if tournament is a circulant (admits a cyclic automorphism)."""
    # Check if rotation by 1 is an automorphism
    rotated = frozenset(((u+1)%n, (v+1)%n) for u,v in arcs)
    return rotated == arcs

# n=3: C3 = cyclic 3-cycle
c3 = circulant_tournament(3, {1})
print(f"\nC3 = circulant(3, {{1}}): arcs={c3}")
print(f"  is circulant: {is_circulant(3, c3)}")
H_c3 = count_hamiltonian_paths_fast(3, c3)
print(f"  H(C3) = {H_c3}")

N_c3sc, c3_sc_arcs = sc_blowup(3, c3)
print(f"\nC3_SC (n=6):")
print(f"  arcs = {sorted(c3_sc_arcs)}")
print(f"  is circulant (n=6): {is_circulant(6, c3_sc_arcs)}")
H_c3sc = count_hamiltonian_paths_fast(6, c3_sc_arcs)
print(f"  H(C3_SC) = {H_c3sc}")
print(f"  is SC tournament: {is_SC(6, c3_sc_arcs)}")
ss_c3sc = score_seq(arcs_to_matrix(6, c3_sc_arcs))
print(f"  score seq: {ss_c3sc}")

# n=5: C5(1,2)
c5 = circulant_tournament(5, {1,2})
print(f"\nC5(1,2): is circulant: {is_circulant(5, c5)}, H={count_hamiltonian_paths_fast(5,c5)}")
N_c5sc, c5_sc_arcs = sc_blowup(5, c5)
print(f"C5(1,2)_SC (n=10): is circulant: {is_circulant(10, c5_sc_arcs)}, H={count_hamiltonian_paths_fast(10, c5_sc_arcs)}")

# n=7: C7(1,2,4) - quadratic residues
c7 = circulant_tournament(7, {1,2,4})
print(f"\nC7(1,2,4) [Paley(7)]: is circulant: {is_circulant(7, c7)}, H={count_hamiltonian_paths_fast(7,c7)}")
N_c7sc, c7_sc_arcs = sc_blowup(7, c7)
print(f"C7(1,2,4)_SC (n=14): is circulant: {is_circulant(14, c7_sc_arcs)}, H={count_hamiltonian_paths_fast(14, c7_sc_arcs)}")

print("\n" + "="*70)
print("PART 4: Tower computation C3 -> C3_SC -> (C3_SC)_SC")
print("="*70)
print(f"H(C3) = {H_c3}")
print(f"H(C3_SC) = {H_c3sc}")

# (C3_SC)_SC
N2, c3_sc2_arcs = sc_blowup(N_c3sc, c3_sc_arcs)
H_c3sc2 = count_hamiltonian_paths_fast(N2, c3_sc2_arcs)
print(f"H((C3_SC)_SC) [n=12] = {H_c3sc2}")
print(f"Ratios: {H_c3sc}/{H_c3} = {H_c3sc/H_c3:.4f}, {H_c3sc2}/{H_c3sc} = {H_c3sc2/H_c3sc:.4f}")
print(f"H values: {H_c3}, {H_c3sc}, {H_c3sc2}")
print(f"If geometric: ratio = {H_c3sc/H_c3:.2f}, next predicted = {H_c3sc * (H_c3sc/H_c3):.0f}")

print("\n" + "="*70)
print("PART 5: Coarse path structure in SC blowup")
print("="*70)
# In T_SC (n=6 from C3), a Hamiltonian path visits 6 vertices.
# Pairs are {0,1},{2,3},{4,5}. A path visits each pair with v_even before v_odd (internal arc).
# What coarse patterns (permutations of pairs) appear?

def get_hp_coarse_patterns(n, sc_arcs):
    """For SC blowup of n-vertex tournament (giving 2n vertices),
    enumerate all Hamiltonian paths and record:
    - coarse permutation of pairs
    - whether within each pair the order is forward (0 before 1) or backward
    """
    N = 2*n
    arc_set = set(sc_arcs)

    coarse_patterns = defaultdict(int)
    order_patterns = defaultdict(int)

    for perm in itertools.permutations(range(N)):
        if all((perm[i], perm[i+1]) in arc_set for i in range(N-1)):
            # This is a Hamiltonian path
            # Extract pair visits
            pair_order = []
            for v in perm:
                pair = v // 2
                sub = v % 2  # 0 or 1
                pair_order.append((pair, sub))

            # Coarse: order of pairs
            pairs_visited = []
            for v in perm:
                p = v // 2
                if not pairs_visited or pairs_visited[-1] != p:
                    pairs_visited.append(p)
            coarse = tuple(pairs_visited)

            # Within each pair: is it (0,1) or (1,0) or split?
            # For each pair, when does 0 appear vs 1?
            pair_sub = {}
            for pos, v in enumerate(perm):
                p = v // 2
                s = v % 2
                if p not in pair_sub:
                    pair_sub[p] = [pos, s]
                else:
                    pair_sub[p].append(pos)
                    pair_sub[p].append(s)

            coarse_patterns[coarse] += 1

    return coarse_patterns

print("\nCoarse patterns for C3_SC (n=6):")
cp = get_hp_coarse_patterns(3, c3_sc_arcs)
total = sum(cp.values())
print(f"Total Hamiltonian paths: {total}")
print(f"Distinct coarse patterns: {len(cp)}")
for pat, cnt in sorted(cp.items()):
    print(f"  {pat}: {cnt} paths")

print("\nCoarse patterns for C5(1,2)_SC (n=10):")
cp5 = get_hp_coarse_patterns(5, c5_sc_arcs)
total5 = sum(cp5.values())
print(f"Total Hamiltonian paths: {total5}")
print(f"Distinct coarse patterns: {len(cp5)}")
# Group by count
from collections import Counter
cnt_dist = Counter(cp5.values())
print(f"Distribution of counts per coarse pattern: {dict(cnt_dist)}")
# Check which coarse patterns correspond to HP in C5
c5_hps = set()
c5_arc_set = set(c5)
for perm in itertools.permutations(range(5)):
    if all((perm[i],perm[i+1]) in c5_arc_set for i in range(4)):
        c5_hps.add(perm)
print(f"Hamiltonian paths in C5: {len(c5_hps)} = {sorted(c5_hps)}")
print(f"Coarse patterns that are HPs in C5: {sum(1 for pat in cp5 if pat in c5_hps)}")
print(f"Coarse patterns that are NOT HPs in C5: {sum(1 for pat in cp5 if pat not in c5_hps)}")

print("\n" + "="*70)
print("PART 6: Check SC blowup of SC tournament")
print("="*70)
for n in [3, 4, 5]:
    sc_tournaments = [r for r in results[n] if r['t_sc']]
    print(f"\nn={n}: {len(sc_tournaments)} SC tournament(s)")
    for r in sc_tournaments:
        print(f"  T{r['idx']} (ss={r['ss']}): H(T)={r['H_T']}, H(T_SC)={r['H_sc']}, T_SC is SC: {r['sc_sc']}")

print("\n" + "="*70)
print("PART 7: Detailed formula search")
print("="*70)
print("\nLooking at H(T_SC) vs H(T)^2 and H(T)*n!/(n)")
print(f"{'n':<4} {'T':<6} {'H(T)':<8} {'H(SC)':<10} {'H(SC)/H(T)^2':<15} {'H(SC)/H(T)':<12}")
for n in [3, 4, 5]:
    for r in results[n]:
        ratio1 = r['H_sc'] / r['H_T']**2 if r['H_T'] > 0 else 0
        ratio2 = r['H_sc'] / r['H_T'] if r['H_T'] > 0 else 0
        print(f"{n:<4} T{r['idx']:<5} {r['H_T']:<8} {r['H_sc']:<10} {ratio1:<15.4f} {ratio2:<12.4f}")

print("\n" + "="*70)
print("PART 8: Paley(7) analysis")
print("="*70)
# Paley(7) = C7({1,2,4})
print(f"Paley(7) = C7({{1,2,4}}): H={count_hamiltonian_paths_fast(7, c7)}")
print(f"Paley(7) is SC: {is_SC(7, c7)}")
H_p7sc = count_hamiltonian_paths_fast(14, c7_sc_arcs)
print(f"Paley(7)_SC (n=14): H={H_p7sc}")
print(f"Paley(7)_SC is SC: {is_SC(14, c7_sc_arcs)}")
print(f"Paley(7)_SC is circulant: {is_circulant(14, c7_sc_arcs)}")

# Check if Paley(13) = circulant(13, {1,3,4,9,10,12}) has similar structure
# Paley(13) connection set = QR mod 13
qr13 = {x**2 % 13 for x in range(1,13)} - {0}
qr13_half = {s for s in qr13 if s < 13-s}  # take canonical half
# Actually for Paley tournament: i->j iff (j-i) mod p is a QR
# QR mod 13: 1,3,4,9,10,12
print(f"\nQR mod 13: {sorted(qr13)}")
# The connection set for the tournament has size (p-1)/2
qr13_conn = {s for s in range(1,13) if pow(s, (13-1)//2, 13) == 1}
print(f"Paley(13) connection set: {sorted(qr13_conn)}")

print("\nDone!")
