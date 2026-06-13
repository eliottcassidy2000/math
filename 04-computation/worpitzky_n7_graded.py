#!/usr/bin/env python3
"""
worpitzky_n7_graded.py — Test graded Worpitzky-OCF pattern at n=7.

At n=6 we found:
  c_j = C(6,j) + delta_j where
  delta_5 = delta_4 = 0                    [universal]
  delta_3 = 8*t3                           [3-cycles]
  delta_2 = 12*t3                          [3-cycles]
  delta_1 = 8*t3 + 4*t5 + 8*alpha_2       [3-cycles, 5-cycles, pairs]
  delta_0 = 2*t3 + 2*t5 + 4*alpha_2 = H-1 [OCF]

QUESTION: At n=7, what is the pattern?
Expected universal levels: d=0, d=1 (delta = 0)
Expected t3-only levels: d=2, d=3
Expected t3+t5+alpha_2 levels: d=4?
New invariant at d=5?: t7 (7-cycles), alpha_3 (disjoint 3-cycle triples)?

Since 2^21 tournaments is too many, we sample.

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from fractions import Fraction
from collections import defaultdict
import random
import numpy as np

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def worpitzky_a(F, n, m):
    return sum(F[k] * comb(m + n - 1 - k, n - 1) for k in range(n))

def exact_worpitzky_coeffs(F, n):
    """Get exact Worpitzky polynomial coefficients using Fraction arithmetic.
    Returns coeffs[j] = coefficient of m^j (low to high)."""
    N = n
    pts = list(range(N))
    vals = [Fraction(worpitzky_a(F, n, m)) for m in pts]
    mat = [[Fraction(m**j) for j in range(N)] + [vals[i]] for i, m in enumerate(pts)]
    for col in range(N):
        for row in range(col, N):
            if mat[row][col] != 0:
                mat[col], mat[row] = mat[row], mat[col]
                break
        pivot = mat[col][col]
        for row in range(col+1, N):
            factor = mat[row][col] / pivot
            for k in range(N+1):
                mat[row][k] -= factor * mat[col][k]
    coeffs = [Fraction(0)] * N
    for row in range(N-1, -1, -1):
        val = mat[row][N]
        for k in range(row+1, N):
            val -= mat[row][k] * coeffs[k]
        coeffs[row] = val / mat[row][row]
    return coeffs

def find_3cycles(adj, n):
    cycles = []
    for i, j, k in combinations(range(n), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            cycles.append((i, j, k))
        elif adj[i][k] and adj[k][j] and adj[j][i]:
            cycles.append((i, k, j))
    return cycles

def count_kcycles(adj, n, k):
    count = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%k]] for i in range(k)):
                count += 1
    return count // k

def count_disjoint_3cycle_pairs(three_cycles):
    count = 0
    for i in range(len(three_cycles)):
        for j in range(i+1, len(three_cycles)):
            if set(three_cycles[i]).isdisjoint(set(three_cycles[j])):
                count += 1
    return count

def count_disjoint_3cycle_triples(three_cycles):
    """Count vertex-disjoint triples of 3-cycles (need n>=9, so 0 at n=7)."""
    count = 0
    for i in range(len(three_cycles)):
        si = set(three_cycles[i])
        for j in range(i+1, len(three_cycles)):
            sj = set(three_cycles[j])
            if not si.isdisjoint(sj):
                continue
            for k in range(j+1, len(three_cycles)):
                sk = set(three_cycles[k])
                if si.isdisjoint(sk) and sj.isdisjoint(sk):
                    count += 1
    return count

def count_disjoint_35_pairs(three_cycles, adj, n):
    """Count pairs: one 3-cycle + one 5-cycle, vertex-disjoint."""
    count = 0
    # Enumerate 5-cycles
    for combo5 in combinations(range(n), 5):
        s5 = set(combo5)
        for perm5 in permutations(combo5):
            if all(adj[perm5[i]][perm5[(i+1)%5]] for i in range(5)):
                # Found a 5-cycle; check against each 3-cycle
                for c3 in three_cycles:
                    if s5.isdisjoint(set(c3)):
                        count += 1
    # Each 5-cycle counted 5 times (rotations)
    return count // 5

# ============================================================
# MAIN: Sample n=7 tournaments
# ============================================================
n = 7
m_edges = n*(n-1)//2  # 21
random.seed(42)

seen_F = {}
all_data = []
NUM_SAMPLES = 500000

print(f"Sampling {NUM_SAMPLES} random tournaments at n={n}...")
print(f"(Also including special tournaments: transitive, rotational)")

# Include special tournaments
special_bits = [0]  # transitive
# Paley tournament for n=7 (QRs mod 7: {1,2,4})
# i->j iff (j-i) mod 7 in {1,2,4}
paley_adj = [[0]*7 for _ in range(7)]
qr = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in qr:
            paley_adj[i][j] = 1
# Convert to bits
paley_bits = 0
idx = 0
for i in range(7):
    for j in range(i+1, 7):
        if paley_adj[i][j]:
            paley_bits |= (1 << idx)
        idx += 1
special_bits.append(paley_bits)

trial = 0
for bits in special_bits:
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F[key] = True

    t3_cycles = find_3cycles(adj, n)
    t3 = len(t3_cycles)
    t5 = count_kcycles(adj, n, 5)
    t7 = count_kcycles(adj, n, 7)
    alpha_2 = count_disjoint_3cycle_pairs(t3_cycles)
    # At n=7, can have at most 2 disjoint 3-cycles (need 6 vertices)
    # so alpha_3 (triples) = 0 always

    coeffs = exact_worpitzky_coeffs(F, n)
    delta = [coeffs[j] - Fraction(comb(n, j)) for j in range(n)]

    all_data.append({
        't3': t3, 't5': t5, 't7': t7, 'alpha_2': alpha_2,
        'coeffs': coeffs, 'delta': delta, 'F': F, 'H': F[n-1]
    })

for trial in range(NUM_SAMPLES):
    bits = random.getrandbits(m_edges)
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F[key] = True

    t3_cycles = find_3cycles(adj, n)
    t3 = len(t3_cycles)
    t5 = count_kcycles(adj, n, 5)
    t7 = count_kcycles(adj, n, 7)
    alpha_2 = count_disjoint_3cycle_pairs(t3_cycles)

    coeffs = exact_worpitzky_coeffs(F, n)
    delta = [coeffs[j] - Fraction(comb(n, j)) for j in range(n)]

    all_data.append({
        't3': t3, 't5': t5, 't7': t7, 'alpha_2': alpha_2,
        'coeffs': coeffs, 'delta': delta, 'F': F, 'H': F[n-1]
    })

    if len(all_data) % 20 == 0:
        print(f"  Found {len(all_data)} distinct F-vectors so far (trial {trial})...")

print(f"\nTotal distinct F-vectors found: {len(all_data)}")

# ============================================================
# VERIFY UNIVERSAL LEVELS
# ============================================================
print("\n" + "=" * 70)
print("CHECKING WHICH LEVELS ARE UNIVERSAL")
print("=" * 70)

for j in range(n-1, -1, -1):
    vals = set(d['delta'][j] for d in all_data)
    if len(vals) == 1 and vals.pop() == 0:
        print(f"  delta_{j} (coeff of m^{j}): UNIVERSAL = 0")
    else:
        print(f"  delta_{j} (coeff of m^{j}): {len(vals)} distinct values, range [{min(vals)}, {max(vals)}]")

# ============================================================
# CHECK t3 DEPENDENCE
# ============================================================
print("\n" + "=" * 70)
print("CHECKING t3 DEPENDENCE")
print("=" * 70)

for j in range(n-1, -1, -1):
    vals = set(d['delta'][j] for d in all_data)
    if len(vals) <= 1:
        continue
    # Check if t3 alone determines this
    t3_to_delta = defaultdict(set)
    for d in all_data:
        t3_to_delta[d['t3']].add(d['delta'][j])
    determined = all(len(v) == 1 for v in t3_to_delta.values())
    if determined:
        print(f"  delta_{j}: DETERMINED by t3 alone")
        # Find linear formula
        items = [(t3, list(vals)[0]) for t3, vals in sorted(t3_to_delta.items())]
        if len(items) >= 2:
            slope = (items[1][1] - items[0][1]) / (items[1][0] - items[0][0])
            intercept = items[0][1] - slope * items[0][0]
            ok = all(abs(v - (intercept + slope * t)) < Fraction(1, 1000) for t, v in items)
            if ok:
                print(f"    Formula: delta_{j} = {intercept} + {slope}*t3")
            else:
                print(f"    Not linear in t3")
    else:
        n_ambiguous = sum(1 for v in t3_to_delta.values() if len(v) > 1)
        print(f"  delta_{j}: NOT determined by t3 alone ({n_ambiguous} ambiguous t3 values)")

# ============================================================
# LINEAR REGRESSION: delta_j = a + b*t3 + c*t5 + d*alpha_2 + e*t7
# ============================================================
print("\n" + "=" * 70)
print("LINEAR REGRESSION FOR EACH DELTA LEVEL")
print("=" * 70)

for j in range(n-1, -1, -1):
    vals = set(d['delta'][j] for d in all_data)
    if len(vals) <= 1:
        continue

    y = np.array([float(d['delta'][j]) for d in all_data])

    # Try: a + b*t3
    X1 = np.array([[1, d['t3']] for d in all_data], dtype=float)
    c1, res1, _, _ = np.linalg.lstsq(X1, y, rcond=None)
    err1 = max(abs(X1 @ c1 - y))

    # Try: a + b*t3 + c*t5 + d*alpha_2
    X2 = np.array([[1, d['t3'], d['t5'], d['alpha_2']] for d in all_data], dtype=float)
    c2, res2, _, _ = np.linalg.lstsq(X2, y, rcond=None)
    err2 = max(abs(X2 @ c2 - y))

    # Try: a + b*t3 + c*t5 + d*alpha_2 + e*t7
    X3 = np.array([[1, d['t3'], d['t5'], d['alpha_2'], d['t7']] for d in all_data], dtype=float)
    c3, res3, _, _ = np.linalg.lstsq(X3, y, rcond=None)
    err3 = max(abs(X3 @ c3 - y))

    # Try: a + b*t3 + c*t5 + d*alpha_2 + e*t7 + f*alpha_35
    # alpha_35 = disjoint (3-cycle, 5-cycle) pairs
    # This is expensive to compute, skip for now

    print(f"\n  delta_{j}:")
    print(f"    t3 only:           coeffs=[{', '.join(f'{x:.4f}' for x in c1)}], max_err={err1:.6f}")
    print(f"    t3+t5+alpha_2:     coeffs=[{', '.join(f'{x:.4f}' for x in c2)}], max_err={err2:.6f}")
    print(f"    t3+t5+alpha_2+t7:  coeffs=[{', '.join(f'{x:.4f}' for x in c3)}], max_err={err3:.6f}")

# ============================================================
# VERIFY c0 = H(T)
# ============================================================
print("\n" + "=" * 70)
print("VERIFYING c0 = H(T)")
print("=" * 70)

all_match = True
for d in all_data:
    c0 = d['coeffs'][0]
    H = Fraction(d['H'])
    if c0 != H:
        print(f"  FAIL: c0={c0}, H={H}")
        all_match = False
        break
if all_match:
    print(f"  ALL {len(all_data)} F-classes: c0 = H(T)")

# ============================================================
# CHECK: delta_0 = H-1 = OCF decomposition
# ============================================================
print("\n" + "=" * 70)
print("OCF DECOMPOSITION OF delta_0 = H-1")
print("=" * 70)

# At n=6: H-1 = 2*t3 + 2*t5 + 4*alpha_2
# At n=7: expect H-1 = 2*t3 + 2*t5 + 2*t7 + 4*alpha_2 + ???
# (each odd-L cycle contributes 2, each disjoint pair contributes 4?)

for d in all_data[:5]:
    H_minus_1 = d['H'] - 1
    formula = 2*d['t3'] + 2*d['t5'] + 2*d['t7'] + 4*d['alpha_2']
    print(f"  t3={d['t3']}, t5={d['t5']}, t7={d['t7']}, alpha_2={d['alpha_2']}, "
          f"H-1={H_minus_1}, 2(t3+t5+t7)+4*alpha_2={formula}, "
          f"match={'YES' if H_minus_1 == formula else 'NO (diff=' + str(H_minus_1-formula) + ')'}")

# Full check
n_match = sum(1 for d in all_data if d['H']-1 == 2*d['t3']+2*d['t5']+2*d['t7']+4*d['alpha_2'])
print(f"\n  Simple formula 2(t3+t5+t7)+4*alpha_2 = H-1: {n_match}/{len(all_data)} match")

if n_match < len(all_data):
    # Need more invariants. Try adding disjoint (3,5) pairs
    print("\n  Computing disjoint (3-cycle, 5-cycle) pairs...")
    for d in all_data:
        if d['H']-1 != 2*d['t3']+2*d['t5']+2*d['t7']+4*d['alpha_2']:
            diff = d['H']-1 - (2*d['t3']+2*d['t5']+2*d['t7']+4*d['alpha_2'])
            print(f"    Mismatch: t3={d['t3']}, t5={d['t5']}, t7={d['t7']}, alpha_2={d['alpha_2']}, diff={diff}")

# ============================================================
# FULL TABLE (sorted by t3)
# ============================================================
print("\n" + "=" * 70)
print(f"FULL TABLE: {len(all_data)} F-VECTOR CLASSES AT n={n}")
print("=" * 70)
print(f"{'#':>3} {'t3':>3} {'t5':>4} {'t7':>3} {'a2':>3} {'H':>5} "
      f"{'d6':>4} {'d5':>4} {'d4':>6} {'d3':>6} {'d2':>8} {'d1':>8} {'d0':>6}")

for i, d in enumerate(sorted(all_data, key=lambda x: (x['t3'], x['t5'], x['t7']))):
    dd = [int(x) for x in d['delta']]
    print(f"{i+1:>3} {d['t3']:>3} {d['t5']:>4} {d['t7']:>3} {d['alpha_2']:>3} {d['H']:>5} "
          f"{dd[6]:>4} {dd[5]:>4} {dd[4]:>6} {dd[3]:>6} {dd[2]:>8} {dd[1]:>8} {dd[0]:>6}")
