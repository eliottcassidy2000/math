#!/usr/bin/env python3
"""
gk_k4_investigation_s112.py — Deep investigation of k >= 4
kind-pasteur-2026-03-15-S112

Key questions:
1. Verify weight formula 2^c/(n)_L at n=10 (first n with 4 separated pairs)
2. Exact structure of corrections Delta_k = TM_g_k - opus_g_k
3. WHY corrections cancel in the CV^2 sum
4. Closed form for corrections
"""

from fractions import Fraction
from math import factorial, comb
from itertools import permutations, combinations

def falling_fact(n, k):
    r = Fraction(1)
    for i in range(k):
        r *= (n - i)
    return r

# ============================================================
# PART 1: Verify weight formula at n=10 for 4 separated pairs
# ============================================================

def verify_weight_n10():
    """Verify E[prod Z_j] = 2^c/(n)_L at n=10 for selected large configs."""
    n = 10
    m = n - 1  # 9 positions (0..8)
    nfact = factorial(n)

    # Precompute Z
    print("Computing Z values for all n=10 permutations...")
    z_all = []
    count = 0
    for perm in permutations(range(n)):
        z = []
        for j in range(m):
            x = 1 if perm[j+1] == perm[j] + 1 else 0
            y = 1 if perm[j+1] == perm[j] - 1 else 0
            z.append(x - y)
        z_all.append(z)
        count += 1
        if count % 500000 == 0:
            print(f"  {count}/{nfact} permutations processed...")

    print(f"Done. {len(z_all)} permutations.")

    # Test specific configurations
    test_configs = [
        # 4 separated pairs
        ({0,1, 3,4, 6,7}, "3 sep pairs, gaps 1,1"),
        ({0,1, 3,4, 6,7, 8}, "invalid - odd size"),  # skip odd
        ({0,1, 4,5, 7,8}, "3 sep pairs, gaps 2,1"),
        # Need 4 pairs = 8 positions from 0..8 (9 available)
        ({0,1, 3,4, 6,7}, "3 pairs"),
        # For 4 pairs in 9 positions: minimum is (0,1),(2,3),(4,5),(6,7) - 8 positions
        # But (2,3) shares position 2 with... wait, pairs are (j,j+1).
        # 4 matching edges from path with 8 edges (positions 0-8):
        # e.g., edges (0,1),(2,3),(4,5),(6,7) → positions {0,1,2,3,4,5,6,7}, L=8
        # connected? 0-1-2-3-4-5-6-7 is all connected! c=1.
        # e.g., edges (0,1),(2,3),(4,5),(7,8) → positions {0,...,5,7,8}, L=8
        # components: {0,...,5} and {7,8}. c=2.
        # e.g., edges (0,1),(3,4),(5,6),(8,?) - no edge 8 exists
        # Actually edges are (j,j+1) for j=0..7 at n=10 (positions 0..8, 8 edges)
    ]

    # Generate all 4-edge matchings of the path with 8 edges
    edges = list(range(8))  # edge i connects positions i and i+1
    matchings_4 = []
    for combo in combinations(edges, 4):
        # Check matching: no two edges share endpoint
        valid = True
        for i in range(len(combo)):
            for j in range(i+1, len(combo)):
                if abs(combo[i] - combo[j]) <= 1:
                    valid = False
                    break
            if not valid:
                break
        if valid:
            matchings_4.append(combo)

    print(f"\n4-edge matchings of P_9 (8 edges): {len(matchings_4)} total")

    for matching in matchings_4:
        positions = set()
        for e in matching:
            positions.add(e)
            positions.add(e + 1)
        S = frozenset(positions)
        L = len(S)

        # Components
        pos_sorted = sorted(S)
        comps = []
        cur = [pos_sorted[0]]
        for p in pos_sorted[1:]:
            if p == cur[-1] + 1:
                cur.append(p)
            else:
                comps.append(cur)
                cur = [p]
        comps.append(cur)
        c = len(comps)

        # Compute E[prod Z_j]
        total = Fraction(0)
        for z in z_all:
            prod = 1
            for j in S:
                prod *= z[j]
            total += prod
        E = total / nfact

        pred = Fraction(2**c, falling_fact(n, L))
        match = "OK" if E == pred else f"FAIL (E={E}, pred={pred})"

        comp_sizes = tuple(len(cc) for cc in comps)
        print(f"  edges={matching}, comps={comp_sizes}, c={c}, L={L}: {match}")

    # Also verify a few 3-pair configs
    print("\n3-edge matchings (spot check):")
    for combo in [(0, 2, 4), (0, 3, 6), (0, 2, 7), (1, 4, 7)]:
        positions = set()
        for e in combo:
            positions.add(e)
            positions.add(e + 1)
        S = frozenset(positions)
        L = len(S)

        pos_sorted = sorted(S)
        comps = []
        cur = [pos_sorted[0]]
        for p in pos_sorted[1:]:
            if p == cur[-1] + 1:
                cur.append(p)
            else:
                comps.append(cur)
                cur = [p]
        comps.append(cur)
        c = len(comps)

        total = Fraction(0)
        for z in z_all:
            prod = 1
            for j in S:
                prod *= z[j]
            total += prod
        E = total / nfact

        pred = Fraction(2**c, falling_fact(n, L))
        match = "OK" if E == pred else f"FAIL"
        comp_sizes = tuple(len(cc) for cc in comps)
        print(f"  edges={combo}, comps={comp_sizes}, c={c}, L={L}: {match}")

verify_weight_n10()

# ============================================================
# PART 2: Exact correction structure Delta_k = TM - opus
# ============================================================

print("\n" + "="*60)
print("PART 2: CORRECTION STRUCTURE")
print("="*60)

def transfer_gk_values(k, m_max):
    """Compute transfer matrix g_k(m) for m=0..m_max."""
    results = {}
    for m in range(0, m_max + 1):
        n = m + 2*k
        num_edges = n - 2
        k_max = (n - 1) // 2
        if k > k_max:
            continue
        state = [[Fraction(1)] + [Fraction(0)] * k_max,
                 [Fraction(0)] * (k_max + 1),
                 [Fraction(0)] * (k_max + 1)]
        for step in range(num_edges):
            A, B, C = state
            nA = [A[i] + C[i] for i in range(k_max + 1)]
            nB = [Fraction(0)] + [2*A[i] + C[i] for i in range(k_max)]
            nC = list(B)
            state = [nA, nB, nC]
        total = [state[0][i] + state[1][i] + state[2][i] for i in range(k_max + 1)]
        if k < len(total) and total[k] != 0:
            results[m] = total[k] / 2
        else:
            results[m] = Fraction(0)
    return results

def opus_gk(k, m):
    C = {3:(2,0,1,0), 4:(10,-33,50,-24), 5:(388,-2040,3431,-1776),
         6:(69660,-380445,653748,-342960), 7:(19826270,-109486152,189674605,-100014720),
         8:(7309726742,-40641958545,70757788486,-37425556680),
         9:(3262687720240,-18232387983408,31858349908595,-16888649645424)}
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m * m)
    if k not in C: return None
    a,b,c,d = C[k]
    return Fraction(a*m**3 + b*m**2 + c*m + d, 3)

# Compute corrections for k=4..9
for k in range(4, 10):
    print(f"\nk={k}: Delta_k(m) = TM_g_k(m) - opus_g_k(m)")
    tm = transfer_gk_values(k, 15)
    deltas = []
    for m in range(0, 13):
        if m not in tm:
            continue
        op = opus_gk(k, m)
        if op is None:
            continue
        d = tm[m] - op
        deltas.append((m, d))
        print(f"  m={m}: TM={tm[m]}, opus={op}, Delta={d}")

    # Forward differences of Delta
    vals = [d for _, d in deltas if _ >= 2]
    ms = [m for m, _ in deltas if m >= 2]
    if len(vals) >= 6:
        d1 = [vals[i+1]-vals[i] for i in range(len(vals)-1)]
        d2 = [d1[i+1]-d1[i] for i in range(len(d1)-1)]
        d3 = [d2[i+1]-d2[i] for i in range(len(d2)-1)]
        d4 = [d3[i+1]-d3[i] for i in range(len(d3)-1)]
        d5 = [d4[i+1]-d4[i] for i in range(len(d4)-1)]
        print(f"  Forward diffs of Delta (from m=2):")
        print(f"    d1: {d1[:6]}")
        print(f"    d2: {d2[:6]}")
        print(f"    d3: {d3[:6]}")
        print(f"    d4: {d4[:6]}")
        print(f"    d5: {d5[:5]}")

    # Check if Delta = c * C(m-1, k) for some c
    print(f"  Ratio Delta/C(m-1,{k}):")
    for m, d in deltas:
        if m >= k+1:
            binom = comb(m-1, k)
            if binom > 0:
                ratio = d / binom
                print(f"    m={m}: Delta={d}, C(m-1,{k})={binom}, ratio={ratio}")

# ============================================================
# PART 3: Binomial basis decomposition of Delta_k
# ============================================================

print("\n" + "="*60)
print("PART 3: BINOMIAL BASIS OF CORRECTIONS")
print("="*60)

for k in range(4, 8):
    tm = transfer_gk_values(k, 15)
    deltas = {}
    for m in range(0, 13):
        if m not in tm:
            continue
        op = opus_gk(k, m)
        if op is not None:
            deltas[m] = tm[m] - op

    # Compute forward differences from m=0 to get binomial basis coefficients
    # Delta_k(m) = sum_{j>=0} alpha_j * C(m, j)
    # alpha_j = j-th forward difference of Delta_k at m=0
    vals = [deltas.get(m, Fraction(0)) for m in range(k + 5)]
    alphas = []
    cur = list(vals)
    for j in range(len(cur)):
        alphas.append(cur[0])
        cur = [cur[i+1] - cur[i] for i in range(len(cur) - 1)]
        if not cur:
            break

    print(f"\nk={k}: Delta_{k}(m) in binomial basis C(m,j):")
    for j, a in enumerate(alphas):
        if a != 0:
            print(f"  alpha_{j} = {a}")

    # So Delta_k(m) = sum alpha_j * C(m,j) with nonzero alpha only for j >= 4
    # This means Delta_k is a polynomial of degree >= 4 that vanishes at m=0,1,2,3

print("\n" + "="*60)
print("PART 4: PATTERN IN BINOMIAL COEFFICIENTS alpha_{k,j}")
print("="*60)

# Collect alpha_{k,j} for k=4..9, j=4..k+1
alpha_table = {}
for k in range(4, 10):
    tm = transfer_gk_values(k, k + 6)
    deltas = {}
    for m in range(0, k + 5):
        if m not in tm:
            continue
        op = opus_gk(k, m)
        if op is not None:
            deltas[m] = tm[m] - op

    vals = [deltas.get(m, Fraction(0)) for m in range(k + 5)]
    cur = list(vals)
    alphas = []
    for j in range(len(cur)):
        alphas.append(cur[0])
        cur = [cur[i+1] - cur[i] for i in range(len(cur) - 1)]
        if not cur:
            break
    alpha_table[k] = alphas

print("alpha_{k,j} table:")
print(f"{'k':>3} | ", end="")
for j in range(12):
    print(f"{'j='+str(j):>14}", end="")
print()
print("-" * 100)
for k in sorted(alpha_table.keys()):
    print(f"{k:>3} | ", end="")
    for j in range(min(12, len(alpha_table[k]))):
        print(f"{str(alpha_table[k][j]):>14}", end="")
    print()

# Check ratios alpha_{k,k}/alpha_{k,4}
print("\nRatio alpha_{k,k}/alpha_{k,4}:")
for k in sorted(alpha_table.keys()):
    if k < len(alpha_table[k]) and alpha_table[k][4] != 0:
        if k < len(alpha_table[k]) and alpha_table[k][k] != 0:
            ratio = alpha_table[k][k] / alpha_table[k][4]
            print(f"  k={k}: alpha_{k},{k}={alpha_table[k][k]}, alpha_{k},4={alpha_table[k][4]}, ratio={ratio}")

# Check if alpha_{k,4} follows a pattern
print("\nalpha_{k,4} sequence:")
for k in sorted(alpha_table.keys()):
    if len(alpha_table[k]) > 4:
        a4 = alpha_table[k][4]
        # Check against 8*C(k-1,3) or similar
        for formula_name, formula in [
            ("8*C(k-1,3)", 8*comb(k-1,3)),
            ("2^(k-1)", 2**(k-1)),
            ("2^k", 2**k),
            ("-8*C(k-1,3)", -8*comb(k-1,3)),
        ]:
            if a4 == formula:
                print(f"  k={k}: alpha_4 = {a4} = {formula_name}")
                break
        else:
            print(f"  k={k}: alpha_4 = {a4}")

print("\nDone!")
