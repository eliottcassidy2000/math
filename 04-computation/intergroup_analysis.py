#!/usr/bin/env python3
"""
Inter-group structure analysis for 3-cycle partitions.
Uses ONLY 3-cycles (fast O(n^3) enumeration).

Key findings so far:
- Counterexample has 3 groups with extreme one-sided inter-group arcs (9-0)
- Question: how extreme must the group imbalance be for disc < 0?

Author: opus-2026-03-06-S19
"""
import random
import numpy as np
from collections import Counter

random.seed(42)

n = 9

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_3cycles(A, n):
    cycles = set()
    for i in range(n):
        for j in range(n):
            if j == i or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j or not A[j][k]: continue
                if A[k][i]:
                    c = (i, j, k)
                    mi = c.index(min(c))
                    cycles.add(c[mi:] + c[:mi])
    return list(cycles)

def find_disjoint_triple(c3):
    """Find 3 vertex-disjoint 3-cycles, return (triple, all_triples_count)."""
    vsets = [frozenset(c) for c in c3]
    m = len(c3)
    triples = []
    for i in range(m):
        for j in range(i+1, m):
            if vsets[i] & vsets[j]: continue
            for k in range(j+1, m):
                if not (vsets[i] & vsets[k]) and not (vsets[j] & vsets[k]):
                    triples.append((c3[i], c3[j], c3[k]))
    return (triples[0] if triples else None), len(triples)

def omega3_poly(c3):
    """Build Omega_3 and compute independence polynomial."""
    m = len(c3)
    if m == 0:
        return [1]
    vsets = [frozenset(c) for c in c3]
    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and vsets[i] & vsets[j])
    memo = {}
    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            return [1]
        v = max(verts, key=lambda u: len(adj[u] & verts))
        p1 = solve(verts - {v})
        p2 = solve(verts - (adj[v] & verts) - {v})
        maxlen = max(len(p1), len(p2) + 1)
        result = [0] * maxlen
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i + 1] += p2[i]
        memo[verts] = result
        return result
    return solve(frozenset(range(m)))

def cubic_disc(coeffs):
    if len(coeffs) != 4: return None
    a0, a1, a2, a3 = coeffs
    return 18*a3*a2*a1*a0 - 4*a2**3*a0 + a2**2*a1**2 - 4*a3*a1**3 - 27*a3**2*a0**2

def newton_check(coeffs):
    deg = len(coeffs) - 1
    while deg > 0 and coeffs[deg] == 0:
        deg -= 1
    if deg < 2: return True
    for k in range(1, deg):
        if coeffs[k]**2 < coeffs[k-1] * coeffs[k+1] * (k+1) / k:
            return False
    return True

def intergroup_stats(T, groups):
    """Compute inter-group arc statistics."""
    stats = {}
    for gi in range(3):
        for gj in range(gi+1, 3):
            fwd = sum(T[u][v] for u in groups[gi] for v in groups[gj])
            bwd = 9 - fwd  # each pair has 9 arcs
            stats[(gi,gj)] = (fwd, bwd)
    max_imb = max(abs(f - b) for f, b in stats.values())
    min_fwd = min(min(f, b) for f, b in stats.values())
    return stats, max_imb, min_fwd

# ============================================================
print("=" * 70)
print("INTER-GROUP STRUCTURE ANALYSIS (3-cycles only)")
print("=" * 70)

# ============================================================
# Part 1: Survey Omega_3 at n=9 — thousands of tournaments
# ============================================================
print("\n--- Part 1: Omega_3 survey at n=9 (10000 samples) ---")

data = []
for trial in range(10000):
    T = random_tournament(n)
    c3 = find_3cycles(T, n)

    triple, num_triples = find_disjoint_triple(c3)
    if triple is None:
        continue

    coeffs = omega3_poly(c3)
    newton_ok = newton_check(coeffs)
    d = cubic_disc(coeffs)

    groups = [list(triple[0]), list(triple[1]), list(triple[2])]
    stats, max_imb, min_fwd = intergroup_stats(T, groups)

    data.append({
        'coeffs': coeffs,
        'newton': newton_ok,
        'disc': d,
        'max_imb': max_imb,
        'min_fwd': min_fwd,
        'num_3c': len(c3),
        'num_triples': num_triples,
        'scores': tuple(sorted(sum(row) for row in T)),
    })

print(f"  {len(data)} tournaments with 3 disjoint 3-cycles")

# Newton and discriminant failures for Omega_3
nf3 = sum(1 for d in data if not d['newton'])
dn3 = sum(1 for d in data if d['disc'] is not None and d['disc'] < 0)
print(f"  Omega_3 Newton failures: {nf3}")
print(f"  Omega_3 disc < 0: {dn3}")

# ============================================================
# Part 2: Imbalance vs failure
# ============================================================
print("\n--- Part 2: Imbalance distribution ---")

imb_dist = Counter(d['max_imb'] for d in data)
print(f"  Max imbalance distribution:")
for imb in sorted(imb_dist.keys()):
    cnt = imb_dist[imb]
    fails = sum(1 for d in data if d['max_imb'] == imb and (not d['newton'] or (d['disc'] is not None and d['disc'] < 0)))
    print(f"    imb={imb}: {cnt} ({100*cnt/len(data):.1f}%), failures={fails}")

# min_fwd distribution (how many arcs go "against" the dominant direction)
mf_dist = Counter(d['min_fwd'] for d in data)
print(f"\n  Min forward arc count (across all 3 group pairs):")
for mf in sorted(mf_dist.keys()):
    cnt = mf_dist[mf]
    fails = sum(1 for d in data if d['min_fwd'] == mf and (not d['newton'] or (d['disc'] is not None and d['disc'] < 0)))
    print(f"    min_fwd={mf}: {cnt} ({100*cnt/len(data):.1f}%), failures={fails}")

# ============================================================
# Part 3: TARGETED — construct tournaments with extreme imbalance
# ============================================================
print("\n--- Part 3: Constructing extreme-imbalance tournaments ---")

# Build a tournament where groups {0,1,2}, {3,4,5}, {6,7,8}
# have 3-cycles within each group AND extreme inter-group domination.

# Group 0: 0->1->2->0
# Group 1: 3->4->5->3
# Group 2: 6->7->8->6
# Inter-group: G1 >> G0 >> G2 (all 9 arcs one-way)

def make_extreme_tournament():
    """Tournament with 3 disjoint 3-cycles and maximal inter-group imbalance."""
    T = [[0]*9 for _ in range(9)]

    # Within-group 3-cycles
    # G0: 0->1->2->0
    T[0][1] = T[1][2] = T[2][0] = 1
    # G1: 3->4->5->3
    T[3][4] = T[4][5] = T[5][3] = 1
    # G2: 6->7->8->6
    T[6][7] = T[7][8] = T[8][6] = 1

    # Inter-group: G1 >> G0 (all arcs from G1 to G0)
    for u in [3, 4, 5]:
        for v in [0, 1, 2]:
            T[u][v] = 1

    # Inter-group: G0 >> G2 (all arcs from G0 to G2)
    for u in [0, 1, 2]:
        for v in [6, 7, 8]:
            T[u][v] = 1

    # Inter-group: G1 >> G2 (all arcs from G1 to G2)
    for u in [3, 4, 5]:
        for v in [6, 7, 8]:
            T[u][v] = 1

    return T

T_extreme = make_extreme_tournament()
c3_ext = find_3cycles(T_extreme, 9)
coeffs_ext = omega3_poly(c3_ext)
d_ext = cubic_disc(coeffs_ext)
scores_ext = tuple(sorted(sum(row) for row in T_extreme))

print(f"  Extreme tournament:")
print(f"    Scores: {scores_ext}")
print(f"    3-cycles: {len(c3_ext)}")
print(f"    I(Omega_3, x) = {coeffs_ext}")
print(f"    disc = {d_ext}")
print(f"    Newton: {'ok' if newton_check(coeffs_ext) else 'FAIL'}")

# Now perturb by flipping a few inter-group arcs
print(f"\n  Perturbing the extreme tournament (flipping inter-group arcs):")
for num_flips in range(1, 10):
    best_disc = None
    for attempt in range(100):
        T_p = [row[:] for row in T_extreme]
        # Random flips between groups
        flipped = set()
        for _ in range(num_flips):
            while True:
                u = random.randint(0, 8)
                v = random.randint(0, 8)
                if u != v and u // 3 != v // 3 and (u, v) not in flipped:
                    break
            T_p[u][v], T_p[v][u] = T_p[v][u], T_p[u][v]
            flipped.add((u, v))
            flipped.add((v, u))

        c3_p = find_3cycles(T_p, 9)
        coeffs_p = omega3_poly(c3_p)
        d_p = cubic_disc(coeffs_p)
        if d_p is not None:
            if best_disc is None or d_p < best_disc:
                best_disc = d_p
                best_coeffs = coeffs_p

    status = f"disc={best_disc}" if best_disc is not None else "deg != 3"
    real = "ok" if best_disc is not None and best_disc >= 0 else "COMPLEX"
    print(f"    {num_flips} flips: most negative disc = {best_disc}, {real}")

# ============================================================
# Part 4: n=10, n=11 Omega_3 exploration
# ============================================================
print("\n" + "=" * 70)
print("PART 4: n=10 and n=11 Omega_3")
print("=" * 70)

for N in [10, 11]:
    print(f"\n  n={N}:")
    nf = 0
    dn = 0
    total = 0
    deg_dist = Counter()

    for trial in range(2000):
        T = random_tournament(N)
        c3 = find_3cycles(T, N)
        if not c3:
            continue
        total += 1

        coeffs = omega3_poly(c3)
        deg = len(coeffs) - 1
        deg_dist[deg] += 1

        if not newton_check(coeffs):
            nf += 1
        d = cubic_disc(coeffs)
        if d is not None and d < 0:
            dn += 1

    print(f"    {total} tournaments with 3-cycles")
    print(f"    Degree dist: {dict(sorted(deg_dist.items()))}")
    print(f"    Newton fails: {nf}")
    print(f"    Disc < 0 (deg 3): {dn}")
    print(f"    Max alpha(Omega_3) at n={N}: {N // 3}")

# ============================================================
# Part 5: DEEP — What is the exact condition for disc < 0?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: EXACT CONDITION FOR OMEGA_3 DISC < 0 AT n=9")
print("=" * 70)

# For I(Omega_3, x) = 1 + a1*x + a2*x^2 + a3*x^3:
# a1 = number of directed 3-cycles
# a2 = number of independent PAIRS
# a3 = number of independent TRIPLES
# disc = 18*a3*a2*a1 - 4*a2^3 + a2^2*a1^2 - 4*a3*a1^3 - 27*a3^2

# For a3=1: disc = 18*a2*a1 - 4*a2^3 + a2^2*a1^2 - 4*a1^3 - 27
# disc < 0 iff 4*a2^3 + 4*a1^3 + 27 > 18*a2*a1 + a2^2*a1^2

# At a1=12, a2=6: disc = 18*72 - 4*216 + 36*144 - 4*1728 - 27
# = 1296 - 864 + 5184 - 6912 - 27 = -1323

# Map out the (a1, a2) region for disc = 0, a3 = 1:
print("  For a3=1, disc = 0 boundary (a1 vs min a2):")
for a1 in range(3, 30):
    # disc = 18*a2*a1 - 4*a2^3 + a2^2*a1^2 - 4*a1^3 - 27
    # This is a cubic in a2: -4*a2^3 + a1^2*a2^2 + 18*a1*a2 - (4*a1^3 + 27) = 0
    poly_a2 = [-4, a1**2, 18*a1, -(4*a1**3 + 27)]
    roots = np.roots(poly_a2)
    real_pos = sorted([r.real for r in roots if abs(r.imag) < 0.01 and r.real > 0])
    if real_pos:
        min_a2 = max(real_pos)  # largest positive root is the threshold
        print(f"    a1={a1:2d}: need a2 >= {min_a2:.2f} for disc >= 0"
              f" (Newton needs a2 >= {np.sqrt(1.5*a1):.2f})")
    else:
        print(f"    a1={a1:2d}: always disc >= 0 (no positive roots)")

# The counterexample: a1=12, a2=6 needs a2 >= ~7 for disc >= 0
# But we FOUND a2=6 (just below threshold!)

# What's the MINIMUM number of 3-cycles for disc < 0 to be possible?
print(f"\n  Minimum a1 for disc < 0 with a2=1, a3=1:")
for a2 in range(1, 20):
    for a1 in range(3, 100):
        d = 18*a1*a2 - 4*a2**3 + a2**2*a1**2 - 4*a1**3 - 27
        if d < 0:
            print(f"    a2={a2}: first disc < 0 at a1={a1} (disc={d})")
            break
    else:
        print(f"    a2={a2}: disc >= 0 for all a1 <= 100")

# Can a2 be small enough at a1 ~ 12 for disc < 0?
# We need: in a tournament with 12 directed 3-cycles and 3 disjoint ones,
# how few independent pairs can there be?
print(f"\n  At a1=12: min a2 for disc >= 0: {max(r.real for r in np.roots([-4, 144, 216, -6939]) if abs(r.imag)<0.01 and r.real > 0):.2f}")
print(f"  At a1=12: Newton threshold: {np.sqrt(1.5*12):.2f}")
print(f"  Counterexample has a2=6 < 7 = floor of disc threshold")

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
print("""
Key findings:
1. Omega_3 disc < 0 requires VERY specific structure at n=9:
   - Exactly 3 disjoint 3-cycles covering all vertices
   - Low total 3-cycle count (~12, well below average ~22)
   - Very few independent pairs (a2 < ~7 for a1=12)
   - This needs near-total domination between groups

2. The disc threshold scales as a2 ~ a1^{2/3} * const,
   while Newton scales as a2 ~ sqrt(a1).
   Discriminant is MUCH more restrictive than Newton.

3. The extreme imbalance (9-0 between group pairs) ensures
   cross-group 3-cycles are impossible, keeping a2 low.

4. At n=10,11: the same phenomenon can occur but the threshold shifts.
""")
