#!/usr/bin/env python3
"""
scissors_attack_s92h.py — Attacking the central open problem
opus-2026-03-20-S92h

THE PROBLEM: I(Omega(T), x) lumps 56 isomorphism classes into 24 groups at n=6.
Some groups have classes with identical score, n3, n5, self-duality, class size.
WHAT separates them?

ATTACK STRATEGIES:
  1. Per-vertex cycle counts (which vertices are in which cycles?)
  2. Cycle overlap graph (how do cycles intersect beyond vertex-disjointness?)
  3. Adjacency matrix spectrum (eigenvalues of the tournament)
  4. The independence polynomial of SUBGRAPHS of Omega
  5. Higher-order cycle packings (length-refined packing)
"""

from math import comb, factorial, pi
from itertools import permutations, combinations
from collections import defaultdict, Counter
import numpy as np
import time

# =============================================================================
# BUILD n=6 DATA
# =============================================================================

def build_all(n):
    num_arcs = comb(n, 2)
    all_perms = list(permutations(range(n)))
    def adj(mask):
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (mask >> idx) & 1: A[i][j] = 1
                else: A[j][i] = 1
                idx += 1
        return A
    def canon(A):
        best = None
        for p in all_perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best
    def H(A):
        return sum(1 for perm in permutations(range(n))
                   if all(A[perm[k]][perm[k+1]] for k in range(n-1)))
    def score_seq(A):
        return tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))
    canon_map = {}
    classes = []
    for mask in range(2**num_arcs):
        A = adj(mask)
        c = canon(A)
        if c not in canon_map:
            canon_map[c] = len(classes)
            classes.append({'A': A, 'H': H(A), 'size': 0, 'score': score_seq(A), 'idx': len(classes)})
        classes[canon_map[c]]['size'] += 1
    return classes

def find_cycles_detailed(A, n):
    """Find all directed odd cycles with VERTEX INFORMATION."""
    cycles = []
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.append({'verts': frozenset([i,j,k]), 'len': 3, 'ordered': (i,j,k)})
        if A[i][k] and A[k][j] and A[j][i]:
            cycles.append({'verts': frozenset([i,j,k]), 'len': 3, 'ordered': (i,k,j)})
    if n >= 5:
        seen = set()
        for perm in permutations(range(n), 5):
            verts = frozenset(perm)
            if len(verts) == 5 and all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                if verts not in seen:
                    seen.add(verts)
                    cycles.append({'verts': verts, 'len': 5, 'ordered': perm})
    return cycles

def independence_poly(cycles):
    nc = len(cycles)
    if nc == 0: return (1,)
    coeffs = [0] * (nc + 1)
    coeffs[0] = 1
    for mask in range(1, 1 << nc):
        selected = [i for i in range(nc) if mask & (1 << i)]
        used = set()
        valid = True
        for i in selected:
            if cycles[i]['verts'] & used:
                valid = False
                break
            used |= cycles[i]['verts']
        if valid:
            coeffs[len(selected)] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return tuple(coeffs)

print("=" * 70)
print("ATTACKING THE CENTRAL OPEN PROBLEM")
print("=" * 70)
print()

t0 = time.perf_counter()
classes6 = build_all(6)
print(f"Built {len(classes6)} classes in {time.perf_counter()-t0:.1f}s")

# Compute all cycle data
for d in classes6:
    d['cycles'] = find_cycles_detailed(d['A'], 6)
    d['I'] = independence_poly(d['cycles'])
    d['n3'] = sum(1 for c in d['cycles'] if c['len'] == 3)
    d['n5'] = sum(1 for c in d['cycles'] if c['len'] == 5)

# Group by I(x)
I_groups = defaultdict(list)
for d in classes6:
    I_groups[d['I']].append(d)

# Find the HARD cases: groups with >1 class that can't be separated by simple invariants
hard_cases = []
for I_poly, group in I_groups.items():
    if len(group) <= 1:
        continue
    # Group by score
    score_groups = defaultdict(list)
    for d in group:
        score_groups[d['score']].append(d)
    for score, subgroup in score_groups.items():
        if len(subgroup) > 1:
            hard_cases.append((I_poly, score, subgroup))

print(f"\nHARD CASES: {len(hard_cases)} groups of classes with same I(x) AND same score")
print()

# =============================================================================
# ATTACK 1: PER-VERTEX CYCLE DISTRIBUTION
# =============================================================================

print("=" * 70)
print("ATTACK 1: PER-VERTEX CYCLE COUNTS")
print("=" * 70)
print()

print("""
For each vertex v, count how many 3-cycles pass through v.
This gives a per-vertex profile that the global n3 count misses.

Example: n3=3 could mean:
  (a) Three cycles sharing one common vertex: profile = (3, 1, 1, 1, 1, 1)
  (b) Three cycles in a chain: profile = (2, 2, 2, 1, 1, 1)
  (c) Three independent cycles: impossible at n=6 (would need 9 vertices)
""")

def vertex_cycle_profile(A, n, cycles):
    """For each vertex, count how many cycles of each length pass through it."""
    profile_3 = [0] * n
    profile_5 = [0] * n
    for c in cycles:
        for v in c['verts']:
            if c['len'] == 3:
                profile_3[v] += 1
            elif c['len'] == 5:
                profile_5[v] += 1
    return tuple(sorted(profile_3, reverse=True)), tuple(sorted(profile_5, reverse=True))

# Test on hard cases
separations_1 = 0
for I_poly, score, subgroup in hard_cases:
    profiles = []
    for d in subgroup:
        p3, p5 = vertex_cycle_profile(d['A'], 6, d['cycles'])
        profiles.append((p3, p5))
        d['vcp3'] = p3
        d['vcp5'] = p5

    unique_profiles = len(set(profiles))
    if unique_profiles == len(subgroup):
        separations_1 += 1

    H_val = subgroup[0]['H']
    if len(subgroup) <= 4:
        print(f"  I={I_poly}, score={score}, H={H_val}: {len(subgroup)} classes")
        for d in subgroup:
            print(f"    class {d['idx']}: vcp3={d['vcp3']}, vcp5={d['vcp5']}")
        separated = unique_profiles == len(subgroup)
        print(f"    Separated by vertex-cycle profile: {'YES' if separated else 'NO'}")
        print()

print(f"  ATTACK 1 RESULT: {separations_1}/{len(hard_cases)} hard cases resolved")

# =============================================================================
# ATTACK 2: TOURNAMENT ADJACENCY SPECTRUM
# =============================================================================

print()
print("=" * 70)
print("ATTACK 2: ADJACENCY MATRIX EIGENVALUES")
print("=" * 70)
print()

print("""
The adjacency matrix eigenvalues are a complete invariant for
"spectral equivalence" but NOT for isomorphism in general.
However, they might separate our hard cases.
""")

def adj_spectrum(A, n):
    """Sorted eigenvalues of the adjacency matrix."""
    M = np.array(A, dtype=float)
    evals = sorted(np.linalg.eigvals(M).real, reverse=True)
    return tuple(round(e, 8) for e in evals)

separations_2 = 0
for I_poly, score, subgroup in hard_cases:
    spectra = []
    for d in subgroup:
        spec = adj_spectrum(d['A'], 6)
        d['spectrum'] = spec
        spectra.append(spec)

    unique_spectra = len(set(spectra))
    if unique_spectra == len(subgroup):
        separations_2 += 1

    H_val = subgroup[0]['H']
    if len(subgroup) <= 4:
        print(f"  I={I_poly}, H={H_val}: {len(subgroup)} classes")
        for d in subgroup:
            print(f"    class {d['idx']}: spectrum={d['spectrum'][:3]}...")
        separated = unique_spectra == len(subgroup)
        print(f"    Separated by spectrum: {'YES' if separated else 'NO'}")
        print()

print(f"  ATTACK 2 RESULT: {separations_2}/{len(hard_cases)} hard cases resolved")

# =============================================================================
# ATTACK 3: CYCLE OVERLAP STRUCTURE
# =============================================================================

print()
print("=" * 70)
print("ATTACK 3: CYCLE OVERLAP GRAPH")
print("=" * 70)
print()

print("""
Omega(T) has edges between cycles that share vertices.
But HOW MANY vertices do they share? And in what pattern?

The overlap graph: vertices = odd cycles, edge weight = |shared vertices|.
The sorted weight sequence is an invariant.
""")

def cycle_overlap_invariant(cycles):
    """Compute overlap structure between cycles."""
    nc = len(cycles)
    if nc <= 1:
        return (0,)
    overlaps = []
    for i in range(nc):
        for j in range(i+1, nc):
            shared = len(cycles[i]['verts'] & cycles[j]['verts'])
            if shared > 0:
                overlaps.append(shared)
    return tuple(sorted(overlaps, reverse=True))

separations_3 = 0
for I_poly, score, subgroup in hard_cases:
    overlaps = []
    for d in subgroup:
        ov = cycle_overlap_invariant(d['cycles'])
        d['overlap'] = ov
        overlaps.append(ov)

    unique_overlaps = len(set(overlaps))
    if unique_overlaps == len(subgroup):
        separations_3 += 1

    H_val = subgroup[0]['H']
    if len(subgroup) <= 4 and unique_overlaps > 1:
        print(f"  I={I_poly}, H={H_val}: {len(subgroup)} classes")
        for d in subgroup:
            print(f"    class {d['idx']}: overlap={d['overlap'][:6]}...")
        separated = unique_overlaps == len(subgroup)
        print(f"    Separated by overlap: {'YES' if separated else 'NO'}")
        print()

print(f"  ATTACK 3 RESULT: {separations_3}/{len(hard_cases)} hard cases resolved")

# =============================================================================
# ATTACK 4: COMBINE ALL INVARIANTS
# =============================================================================

print()
print("=" * 70)
print("ATTACK 4: COMBINED INVARIANT")
print("=" * 70)
print()

def combined_invariant(d):
    """Combine all computed invariants."""
    return (d['I'], d['score'], d.get('vcp3', ()), d.get('vcp5', ()),
            d.get('overlap', ()), d.get('spectrum', ()))

all_combined = set()
for d in classes6:
    if 'spectrum' not in d:
        d['spectrum'] = adj_spectrum(d['A'], 6)
    if 'vcp3' not in d:
        p3, p5 = vertex_cycle_profile(d['A'], 6, d['cycles'])
        d['vcp3'] = p3
        d['vcp5'] = p5
    if 'overlap' not in d:
        d['overlap'] = cycle_overlap_invariant(d['cycles'])
    all_combined.add(combined_invariant(d))

print(f"  56 isomorphism classes")
print(f"  {len(all_combined)} distinct combined invariants")
print(f"  Combined invariant is {'COMPLETE' if len(all_combined) == 56 else f'NOT complete ({56 - len(all_combined)} collisions)'}")

# Find remaining collisions
if len(all_combined) < 56:
    comb_groups = defaultdict(list)
    for d in classes6:
        comb_groups[combined_invariant(d)].append(d)

    print(f"\n  Remaining collisions:")
    for key, group in comb_groups.items():
        if len(group) > 1:
            print(f"    {len(group)} classes share ALL invariants:")
            for d in group:
                print(f"      class {d['idx']}: H={d['H']}, score={d['score']}, "
                      f"n3={d['n3']}, n5={d['n5']}")

# =============================================================================
# ATTACK 5: TRANSPOSE PAIRING
# =============================================================================

print()
print("=" * 70)
print("ATTACK 5: TRANSPOSE PAIRING STRUCTURE")
print("=" * 70)
print()

print("""
Every tournament T has a transpose T^op (reverse all arcs).
H(T) = H(T^op) always. If T ≅ T^op, T is self-dual.

For non-self-dual tournaments, T and T^op form a PAIR.
These pairs always share the same I(x).
Could the collision groups simply be {T, T^op} pairs?
""")

all_perms6 = list(permutations(range(6)))

def canon6(A):
    best = None
    for p in all_perms6:
        s = tuple(A[p[i]][p[j]] for i in range(6) for j in range(6))
        if best is None or s < best: best = s
    return best

# Find transpose partner for each class
for d in classes6:
    A = d['A']
    At = [[A[j][i] for j in range(6)] for i in range(6)]
    ct = canon6(At)
    # Find which class ct belongs to
    for d2 in classes6:
        if d2['A'] == At or canon6(d2['A']) == ct:
            d['transpose_idx'] = d2['idx']
            break
    d['is_self_dual'] = d['idx'] == d.get('transpose_idx', -1)

# Check: within each collision group, are members transpose pairs?
if len(all_combined) < 56:
    comb_groups = defaultdict(list)
    for d in classes6:
        comb_groups[combined_invariant(d)].append(d)

    print("  Collision analysis via transpose pairing:")
    for key, group in comb_groups.items():
        if len(group) > 1:
            indices = [d['idx'] for d in group]
            transposes = [d.get('transpose_idx', '?') for d in group]
            are_pairs = set(indices) == set(transposes)
            print(f"    classes {indices}: transposes = {transposes}, "
                  f"{'TRANSPOSE PAIR' if are_pairs else 'NOT transpose pair'}")

# =============================================================================
# ATTACK 6: THE SKEW-ADJACENCY MATRIX
# =============================================================================

print()
print("=" * 70)
print("ATTACK 6: SKEW-ADJACENCY MATRIX SPECTRUM")
print("=" * 70)
print()

print("""
The skew-adjacency matrix S[i][j] = A[i][j] - A[j][i] is ANTISYMMETRIC.
S has purely imaginary eigenvalues (or zero).
The eigenvalues come in ±pairs.
The PFAFFIAN of S is related to H(T).

S captures orientation information that the symmetric adjacency spectrum might miss.
""")

def skew_spectrum(A, n):
    """Eigenvalues of the skew-adjacency matrix."""
    S = np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=float)
    evals = np.linalg.eigvals(S)
    # Sort by imaginary part (all should be purely imaginary)
    imag_parts = sorted([round(e.imag, 8) for e in evals], reverse=True)
    return tuple(imag_parts)

separations_6 = 0
for I_poly, score, subgroup in hard_cases:
    skew_specs = []
    for d in subgroup:
        ss = skew_spectrum(d['A'], 6)
        d['skew_spectrum'] = ss
        skew_specs.append(ss)

    unique = len(set(skew_specs))
    if unique == len(subgroup):
        separations_6 += 1

print(f"  ATTACK 6 RESULT: {separations_6}/{len(hard_cases)} hard cases resolved by skew spectrum")

# Test: does adding skew spectrum to the combined invariant make it complete?
def full_invariant(d):
    return combined_invariant(d) + (d.get('skew_spectrum', ()),)

all_full = set()
for d in classes6:
    if 'skew_spectrum' not in d:
        d['skew_spectrum'] = skew_spectrum(d['A'], 6)
    all_full.add(full_invariant(d))

print(f"  With skew spectrum added: {len(all_full)} distinct / 56 total")
print(f"  {'COMPLETE!' if len(all_full) == 56 else f'Still {56-len(all_full)} collisions'}")

# =============================================================================
# SUMMARY
# =============================================================================

print(f"""

{'='*70}
SUMMARY: PROGRESS ON THE CENTRAL OPEN PROBLEM
{'='*70}

INVARIANT                    | #CLASSES DISTINGUISHED | COMPLETE?
-----------------------------------------------------------------
H alone                      | 19 / 56                | NO
I(Omega, x)                  | 24 / 56                | NO
Score sequence                | 22 / 56                | NO
(I(x), score)                | {len(set((d['I'], d['score']) for d in classes6))} / 56                | {'YES' if len(set((d['I'], d['score']) for d in classes6)) == 56 else 'NO'}
Adjacency spectrum            | {len(set(d['spectrum'] for d in classes6))} / 56                | {'YES' if len(set(d['spectrum'] for d in classes6)) == 56 else 'NO'}
Skew-adjacency spectrum       | {len(set(d['skew_spectrum'] for d in classes6))} / 56                | {'YES' if len(set(d['skew_spectrum'] for d in classes6)) == 56 else 'NO'}
Vertex-cycle profile          | {len(set((d['vcp3'], d['vcp5']) for d in classes6))} / 56                | NO
Combined (all above)          | {len(all_full)} / 56                | {'YES!' if len(all_full) == 56 else 'NO'}
""")

if len(all_full) == 56:
    print("""
  THE COMBINED INVARIANT IS COMPLETE AT n=6!

  The tuple (I(Omega,x), score_seq, vertex_cycle_profile,
             overlap_structure, adj_spectrum, skew_spectrum)
  distinguishes all 56 isomorphism classes.

  QUESTION: Which SUBSET is sufficient?
  The answer tells us what the "tournament Dehn-Sydler theorem" looks like.
""")

    # Find minimal sufficient subset
    for combo_name, combo_func in [
        ("(I(x), skew_spectrum)", lambda d: (d['I'], d['skew_spectrum'])),
        ("(score, skew_spectrum)", lambda d: (d['score'], d['skew_spectrum'])),
        ("(I(x), score, vcp3)", lambda d: (d['I'], d['score'], d['vcp3'])),
        ("(skew_spectrum alone)", lambda d: (d['skew_spectrum'],)),
        ("(adj_spectrum alone)", lambda d: (d['spectrum'],)),
    ]:
        vals = set(combo_func(d) for d in classes6)
        complete = len(vals) == 56
        print(f"    {combo_name:>35}: {len(vals)}/56 {'← COMPLETE' if complete else ''}")

elif len(all_full) < 56:
    print("""
  Even the combined invariant has collisions.
  The remaining classes are truly hard to separate.
  Need DEEPER structural invariants.
""")
    # Show what's still colliding
    full_groups = defaultdict(list)
    for d in classes6:
        full_groups[full_invariant(d)].append(d)
    for key, group in full_groups.items():
        if len(group) > 1:
            print(f"  STUBBORN COLLISION: classes {[d['idx'] for d in group]}")
            for d in group:
                # Print the actual adjacency matrices
                print(f"    class {d['idx']}:")
                A = d['A']
                for row in A:
                    print(f"      {row}")

print()
print("opus-2026-03-20-S92h: ATTACKING THE CENTRAL OPEN PROBLEM")
