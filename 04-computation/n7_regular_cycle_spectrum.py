#!/usr/bin/env python3
"""
n7_regular_cycle_spectrum.py -- kind-pasteur-2026-03-13-S61

FOCUSED: Extract the cycle spectrum (c3, c5, c7) for all 2640 regular
n=7 tournaments, and analyze the Fourier-Cycle bridge.

At n=7 regular:
  - All scores = 3, so Var_s = 0 and H_2 is maximal
  - H = H_0 + H_2 + H_4 + H_6 (max degree = 2*floor(6/2) = 6)
  - Three rigid H classes: 171, 175, 189

Key questions:
  1. Are c5 and c7 each constant within an H class?
  2. If so, can we separate H_4 and H_6 using the bridge:
     H_4 = f(c5), H_6 = g(c7)?
  3. What is the overlap weight structure for each class?

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_subset(A, verts):
    k = len(verts)
    if k < 3:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


n = 7
m = n * (n - 1) // 2  # 21
total = 1 << m
EH = math.factorial(n) / (2 ** (n - 1))
c2 = math.factorial(n - 2) / (2 ** (n - 2))
H2_reg = c2 * n * (n - 1) / 2  # H_2 for regular tournaments

print("=" * 70)
print(f"REGULAR n={n} TOURNAMENT CYCLE SPECTRUM")
print("=" * 70)
print(f"m = {m}, total = {total}")
print(f"E[H] = {EH:.4f}")
print(f"H_2(regular) = {H2_reg:.4f}")
print(f"H_0 + H_2 = {EH + H2_reg:.4f}")
print()

# Scan all 2^21 tournaments for regular ones
regular = []
scanned = 0
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if any(s != 3 for s in scores):
        continue

    H = count_ham_paths(A, n)

    # Count directed 3-cycles
    c3_dir = 0
    c3_vsets = set()
    for a, b, c in combinations(range(n), 3):
        fwd = A[a][b] * A[b][c] * A[c][a]
        bwd = A[a][c] * A[c][b] * A[b][a]
        c3_dir += fwd + bwd
        if fwd or bwd:
            c3_vsets.add(frozenset([a, b, c]))

    # Count directed 5-cycles (on each 5-subset)
    c5_dir = 0
    c5_vsets = set()
    for subset in combinations(range(n), 5):
        cnt = count_directed_ham_cycles_subset(A, list(subset))
        c5_dir += cnt
        if cnt > 0:
            c5_vsets.add(frozenset(subset))

    # Count directed 7-cycles (Hamiltonian cycles)
    c7_dir = count_directed_ham_cycles_subset(A, list(range(n)))

    # Overlap analysis for 3-cycle vertex sets
    c3_list = list(c3_vsets)
    disj_33 = 0
    ov1_33 = 0
    ov2_33 = 0
    for i in range(len(c3_list)):
        for j in range(i+1, len(c3_list)):
            ov = len(c3_list[i] & c3_list[j])
            if ov == 0:
                disj_33 += 1
            elif ov == 1:
                ov1_33 += 1
            else:
                ov2_33 += 1

    # Also check 3-5 overlap and 5-5 overlap
    c5_list = list(c5_vsets)
    disj_35 = 0
    for c3s in c3_list:
        for c5s in c5_list:
            if not (c3s & c5s):
                disj_35 += 1

    disj_55 = 0
    for i in range(len(c5_list)):
        for j in range(i+1, len(c5_list)):
            if not (c5_list[i] & c5_list[j]):
                disj_55 += 1

    # AA^T off-diagonal variance
    off_diag = []
    for i in range(n):
        for j in range(i+1, n):
            off_diag.append(sum(A[i][k] * A[j][k] for k in range(n)))
    aat_var = sum((x - sum(off_diag)/len(off_diag))**2 for x in off_diag) / len(off_diag)

    regular.append({
        'bits': bits, 'H': H,
        'c3_dir': c3_dir, 'c3_vsets': len(c3_vsets),
        'c5_dir': c5_dir, 'c5_vsets': len(c5_vsets),
        'c7_dir': c7_dir,
        'disj_33': disj_33, 'ov1_33': ov1_33, 'ov2_33': ov2_33,
        'disj_35': disj_35, 'disj_55': disj_55,
        'aat_var': round(aat_var, 6)
    })

print(f"Found {len(regular)} regular tournaments")
print()

# Group by H
by_H = defaultdict(list)
for d in regular:
    by_H[d['H']].append(d)

print("=" * 70)
print("CLASS SUMMARY")
print("=" * 70)
print(f"{'H':>5s} | {'cnt':>5s} | {'c3_dir':>7s} | {'c3_vs':>6s} | {'c5_dir':>7s} | {'c5_vs':>6s} | {'c7_dir':>7s} | {'disj33':>7s} | {'disj35':>7s} | {'disj55':>7s} | {'AA^T var':>10s}")
print(f"{'':->5s}-+-{'':->5s}-+-{'':->7s}-+-{'':->6s}-+-{'':->7s}-+-{'':->6s}-+-{'':->7s}-+-{'':->7s}-+-{'':->7s}-+-{'':->7s}-+-{'':->10s}")

for H in sorted(by_H.keys()):
    g = by_H[H]
    def show(key):
        vals = sorted(set(d[key] for d in g))
        return str(vals[0]) if len(vals) == 1 else str(vals)
    print(f"  {H:>3d} | {len(g):>5d} | {show('c3_dir'):>7s} | {show('c3_vsets'):>6s} | "
          f"{show('c5_dir'):>7s} | {show('c5_vsets'):>6s} | {show('c7_dir'):>7s} | "
          f"{show('disj_33'):>7s} | {show('disj_35'):>7s} | {show('disj_55'):>7s} | {show('aat_var'):>10s}")

# Now try to decompose H_4 and H_6
print(f"\n{'='*70}")
print("FOURIER-CYCLE BRIDGE: SEPARATING H_4 AND H_6")
print("=" * 70)

print(f"\nH_0 + H_2 = {EH + H2_reg:.4f}")
print(f"So H_4 + H_6 = H - {EH + H2_reg:.4f}")

# If c3, c5, c7 are each constant per H class, we have:
# H = f(c3, c5, c7)
# And since c3 is score-determined (always same for regular), c3 drops out of H_4 + H_6
# We need: H_4 + H_6 = a*c5_dir + b*c7_dir + constant

data_pts = []
for H in sorted(by_H.keys()):
    g = by_H[H]
    c5 = list(set(d['c5_dir'] for d in g))
    c7 = list(set(d['c7_dir'] for d in g))
    disj33 = list(set(d['disj_33'] for d in g))
    H46 = H - EH - H2_reg

    print(f"\n  H={H}: H_4+H_6 = {H46:.4f}")
    print(f"    c5_dir = {c5}, c7_dir = {c7}, disj_33 = {disj33}")

    if len(c5) == 1 and len(c7) == 1:
        data_pts.append((c5[0], c7[0], H46, disj33[0]))

# Solve for coefficients
if len(data_pts) == 3:
    print(f"\n  System of equations (3 unknowns: a, b, c):")
    print(f"  H_4+H_6 = a*c5_dir + b*c7_dir + c")

    c5s = [p[0] for p in data_pts]
    c7s = [p[1] for p in data_pts]
    H46s = [p[2] for p in data_pts]
    d33s = [p[3] for p in data_pts]

    for i in range(3):
        print(f"    {H46s[i]:.4f} = a*{c5s[i]} + b*{c7s[i]} + c")

    # Solve 3x3 system manually
    dc5_01 = c5s[1] - c5s[0]
    dc7_01 = c7s[1] - c7s[0]
    dH_01 = H46s[1] - H46s[0]

    dc5_02 = c5s[2] - c5s[0]
    dc7_02 = c7s[2] - c7s[0]
    dH_02 = H46s[2] - H46s[0]

    det = dc5_01 * dc7_02 - dc5_02 * dc7_01
    if abs(det) > 1e-10:
        a = (dH_01 * dc7_02 - dH_02 * dc7_01) / det
        b = (dc5_01 * dH_02 - dc5_02 * dH_01) / det
        c_const = H46s[0] - a * c5s[0] - b * c7s[0]

        print(f"\n  Solution:")
        print(f"    a (c5 coeff) = {a:.8f}")
        print(f"    b (c7 coeff) = {b:.8f}")
        print(f"    c (constant) = {c_const:.8f}")

        # Verify
        print(f"\n  Verification:")
        for i in range(3):
            pred = a * c5s[i] + b * c7s[i] + c_const
            print(f"    H_4+H_6 predicted={pred:.8f}, actual={H46s[i]:.8f}, "
                  f"diff={abs(pred-H46s[i]):.10f}")

        # Check if a and b are nice fractions
        # a could be related to (n-2)!/2^{n-2} or similar
        c2_val = math.factorial(n-2) / (2 ** (n-2))
        c4_val = math.factorial(n-4) / (2 ** (n-4)) if n >= 4 else 0
        print(f"\n  Checking if coefficients are recognizable:")
        print(f"    a = {a:.8f}")
        print(f"    b = {b:.8f}")
        print(f"    c2 = (n-2)!/2^(n-2) = {c2_val:.8f}")
        print(f"    c4 = (n-4)!/2^(n-4) = {c4_val:.8f}")
        print(f"    a/c4 = {a/c4_val:.8f}" if c4_val != 0 else "    c4 = 0")
        print(f"    a * 5 = {a*5:.8f}")
        print(f"    a * 10 = {a*10:.8f}")
        print(f"    b * 7 = {b*7:.8f}")
        print(f"    b * 14 = {b*14:.8f}")
        print(f"    b * 42 = {b*42:.8f}")

        # Also try the model: H_4+H_6 = a'*disj_33 + b'*c7_dir + c'
        print(f"\n  Alternative model: H_4+H_6 = a'*disj_33 + b'*c7_dir + c'")
        dd_01 = d33s[1] - d33s[0]
        dd_02 = d33s[2] - d33s[0]
        det2 = dd_01 * dc7_02 - dd_02 * dc7_01
        if abs(det2) > 1e-10:
            a2 = (dH_01 * dc7_02 - dH_02 * dc7_01) / det2
            b2 = (dd_01 * dH_02 - dd_02 * dH_01) / det2
            c2_const = H46s[0] - a2 * d33s[0] - b2 * c7s[0]
            print(f"    a' (disj_33) = {a2:.8f}")
            print(f"    b' (c7_dir)  = {b2:.8f}")
            print(f"    c' (constant) = {c2_const:.8f}")
    else:
        print(f"  System is degenerate (det={det:.6f})")
        # Try just 2-param models
        print(f"\n  Trying 2-param model: H_4+H_6 = a*c5_dir + c")
        if dc5_01 != 0:
            a_simple = dH_01 / dc5_01
            c_simple = H46s[0] - a_simple * c5s[0]
            print(f"    a = {a_simple:.8f}, c = {c_simple:.8f}")
            for i in range(3):
                pred = a_simple * c5s[i] + c_simple
                print(f"    H_4+H_6 predicted={pred:.8f}, actual={H46s[i]:.8f}, diff={abs(pred-H46s[i]):.6f}")

        print(f"\n  Trying 2-param model: H_4+H_6 = b*c7_dir + c")
        dc7_01_val = c7s[1] - c7s[0]
        if dc7_01_val != 0:
            b_simple = dH_01 / dc7_01_val
            c_simple = H46s[0] - b_simple * c7s[0]
            print(f"    b = {b_simple:.8f}, c = {c_simple:.8f}")
            for i in range(3):
                pred = b_simple * c7s[i] + c_simple
                print(f"    H_4+H_6 predicted={pred:.8f}, actual={H46s[i]:.8f}, diff={abs(pred-H46s[i]):.6f}")


# Also: OCF decomposition for each class
print(f"\n{'='*70}")
print("OCF ALPHA DECOMPOSITION")
print("=" * 70)

for H in sorted(by_H.keys()):
    d = by_H[H][0]  # Representative
    bits = d['bits']
    A = binary_to_tournament(bits, n)

    # Get ALL odd-cycle vertex sets
    all_cycle_vsets = []
    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            cnt = count_directed_ham_cycles_subset(A, list(subset))
            if cnt > 0:
                all_cycle_vsets.append(frozenset(subset))

    nc = len(all_cycle_vsets)

    # Build adjacency of Omega
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if all_cycle_vsets[i] & all_cycle_vsets[j]:
                adj[i][j] = 1
                adj[j][i] = 1

    # Count independent sets by backtracking
    nbr = [0] * nc
    for i in range(nc):
        for j in range(nc):
            if adj[i][j]:
                nbr[i] |= (1 << j)

    alpha = [0] * (nc + 1)
    def backtrack(v, mask, size):
        alpha[size] += 1
        for w in range(v + 1, nc):
            if not (mask & (1 << w)):
                new_mask = mask | nbr[w]
                backtrack(w, new_mask, size + 1)

    backtrack(-1, 0, 0)

    # I(Omega, 2)
    I_val = sum(alpha[j] * (2**j) for j in range(nc + 1))

    max_j = max(j for j in range(nc + 1) if alpha[j] > 0)

    print(f"\n  H={H} (alpha_1={nc}):")
    for j in range(max_j + 1):
        if alpha[j] > 0:
            print(f"    alpha_{j} = {alpha[j]}")
    print(f"    I(Omega, 2) = {I_val} (should be {H}): {'MATCH' if I_val == H else 'MISMATCH!'}")

    # Show cycle composition
    c3_in = sum(1 for s in all_cycle_vsets if len(s) == 3)
    c5_in = sum(1 for s in all_cycle_vsets if len(s) == 5)
    c7_in = sum(1 for s in all_cycle_vsets if len(s) == 7)
    print(f"    Omega vertices: {c3_in} from c3 + {c5_in} from c5 + {c7_in} from c7 = {nc}")


print("\n" + "=" * 70)
print("DONE.")
