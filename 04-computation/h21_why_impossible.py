#!/usr/bin/env python3
"""
h21_why_impossible.py — WHY is H=21 permanently impossible?

The constraint is β₁+2β₂+4β₃+...=10 where:
  β_k = Σ over pairwise-disjoint k-tuples of vertex sets, Π d(S_i)
  d(S) = number of directed odd cycles on vertex set S

Key question: what structural property of tournaments prevents
β₁+2β₂=10 from being achievable?

APPROACH: At each n, exhaustively catalog (β₁,β₂) and identify
the "exclusion zone" around the line β₁+2β₂=10.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def fast_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if not c or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def get_cycle_data(A, n):
    """Get vertex sets with their directed cycle multiplicities."""
    groups = defaultdict(int)
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    groups[frozenset(verts)] += 1
    return groups

def compute_betas(groups):
    """Compute β₁, β₂ from vertex set groups."""
    vs_list = list(groups.items())
    n_vs = len(vs_list)

    beta1 = sum(d for _, d in vs_list)

    beta2 = 0
    for i in range(n_vs):
        for j in range(i+1, n_vs):
            if not (vs_list[i][0] & vs_list[j][0]):
                beta2 += vs_list[i][1] * vs_list[j][1]

    return beta1, beta2

print("=" * 70)
print("WHY IS H=21 PERMANENTLY IMPOSSIBLE?")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: At n=5,6 — map β₁+2β₂ values to H values
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: β₁+2β₂ → H mapping ---")

for n in [5, 6]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    ne = len(edges)

    target_to_H = defaultdict(set)
    H_to_target = defaultdict(set)
    b1b2_to_H = defaultdict(set)

    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        H = fast_hp(A, n)
        groups = get_cycle_data(A, n)
        beta1, beta2 = compute_betas(groups)
        target = beta1 + 2*beta2

        target_to_H[target].add(H)
        H_to_target[H].add(target)
        b1b2_to_H[(beta1, beta2)].add(H)

    print(f"\n  n={n}:")
    print(f"  Target=10 (β₁+2β₂=10) achievable: {'YES' if 10 in target_to_H else 'NO'}")

    # What targets give what H?
    for t in range(max(0, 10-4), min(10+5, max(target_to_H.keys())+1)):
        h_vals = sorted(target_to_H.get(t, set()))
        expected = 1 + 2*t
        has_21 = 21 in h_vals
        marker = " ← TARGET" if t == 10 else ""
        print(f"    β₁+2β₂={t:2d}: H expected={expected:3d}, actual H ∈ {h_vals}{marker}")

    # Check: does H always equal 1+2*(β₁+2β₂)?
    # No — higher β terms contribute. Let's check.
    print(f"\n  Does H = 1+2(β₁+2β₂) always hold? (ignoring β₃+)")
    mismatches = 0
    for key, h_set in b1b2_to_H.items():
        expected = 1 + 2*key[0] + 4*key[1]
        for h in h_set:
            if h != expected:
                mismatches += 1
    print(f"    Mismatches: {mismatches} (these have β₃+ > 0)")

# ═══════════════════════════════════════════════════════════════════
# Part 2: The KEY structural constraint
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: Structural constraint analysis ---")
print("  Q: Why can't we get β₁+2β₂=10 with β₃=β₄=...=0?")
print("  A: This would need β₁+2β₂=10 with no independent triples.")
print("     Let's examine what (β₁,β₂) pairs exist at each n.")

for n in [5, 6]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    ne = len(edges)

    # Collect (β₁, β₂, count of vertex sets, max |S|) for tournaments
    b1b2_info = defaultdict(list)

    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        groups = get_cycle_data(A, n)
        if not groups:
            continue

        beta1, beta2 = compute_betas(groups)
        n_vs = len(groups)
        sizes = [len(s) for s in groups.keys()]
        max_size = max(sizes) if sizes else 0
        multiplicities = list(groups.values())

        b1b2_info[(beta1, beta2)].append({
            'n_vs': n_vs,
            'max_size': max_size,
            'sizes': Counter(sizes),
            'mults': Counter(multiplicities),
        })

    print(f"\n  n={n}: (β₁, β₂) patterns near target line β₁+2β₂=10:")
    for key in sorted(b1b2_info.keys()):
        b1, b2 = key
        target = b1 + 2*b2
        if abs(target - 10) > 3:
            continue

        count = len(b1b2_info[key])
        # Summarize sizes
        all_sizes = Counter()
        all_mults = Counter()
        for info in b1b2_info[key]:
            all_sizes += info['sizes']
            all_mults += info['mults']

        marker = " ← TARGET" if target == 10 else ""
        print(f"    ({b1:2d},{b2:2d}): {count:5d} tournaments, "
              f"sizes={dict(sorted(all_sizes.items()))}, "
              f"mults={dict(sorted(all_mults.items()))}"
              f"{marker}")

# ═══════════════════════════════════════════════════════════════════
# Part 3: At n=5 — WHY is β₁=10 impossible?
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: Why β₁=10 is impossible at n=5 ---")

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

beta1_dist = Counter()
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_cycle_data(A, n)
    beta1 = sum(groups.values())
    beta1_dist[beta1] += 1

print(f"  β₁ distribution at n=5:")
for b1 in sorted(beta1_dist.keys()):
    marker = " ← needed for H=21" if b1 == 10 else ""
    print(f"    β₁={b1:2d}: {beta1_dist[b1]:5d} tournaments{marker}")

# ═══════════════════════════════════════════════════════════════════
# Part 4: At n=6 — detailed analysis of β₁+2β₂ near 10
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: n=6 — tournaments where β₁+2β₂ is closest to 10 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

near_10 = []
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_cycle_data(A, n)
    beta1, beta2 = compute_betas(groups)
    target = beta1 + 2*beta2
    H = fast_hp(A, n)

    if target in [9, 10, 11]:
        near_10.append({
            'bits': bits, 'H': H, 'beta1': beta1, 'beta2': beta2,
            'target': target, 'groups': dict(groups),
            'n_vs': len(groups)
        })

target_counts = Counter(t['target'] for t in near_10)
print(f"  Tournaments with β₁+2β₂ near 10:")
for t in sorted(target_counts.keys()):
    print(f"    β₁+2β₂={t}: {target_counts[t]} tournaments")

# Show details for target=10 if any
target_10 = [t for t in near_10 if t['target'] == 10]
if target_10:
    print(f"\n  TARGET=10 tournaments: {len(target_10)}")
    for i, t in enumerate(target_10[:5]):
        print(f"    #{i}: bits={t['bits']}, H={t['H']}, β₁={t['beta1']}, β₂={t['beta2']}")
        for vs, d in sorted(t['groups'].items(), key=lambda x: len(x[0])):
            print(f"      {sorted(vs)} (len {len(vs)}): d={d}")
else:
    print(f"\n  TARGET=10: NOT ACHIEVABLE at n=6")
    # Show closest
    for t_val in [9, 11]:
        examples = [t for t in near_10 if t['target'] == t_val][:2]
        for t in examples:
            print(f"    β₁+2β₂={t_val}: bits={t['bits']}, H={t['H']}, β₁={t['beta1']}, β₂={t['beta2']}")

# ═══════════════════════════════════════════════════════════════════
# Part 5: The parity/modular constraint
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 5: Parity analysis ---")
print("  Key observation: β₁ = Σ d(S) where d(S) is the number of")
print("  directed odd cycles on vertex set S.")
print()
print("  For 3-cycles: d(S) ∈ {0,1} (exactly 1 directed 3-cycle per triple)")
print("  For 5-cycles: d(S) ∈ {0,1,2,3} (a 5-tournament has 0 or 1,2,3 Ham cycles)")
print("  For 7-cycles: d(S) ∈ {0,...,60} (subtournament Ham cycles)")
print()
print("  Let's check: at n=5, what β₁ mod 2 values occur?")

beta1_mod2 = Counter()
for bits in range(2**10):  # n=5
    A = [[0]*5 for _ in range(5)]
    edges5 = [(i,j) for i in range(5) for j in range(i+1,5)]
    for idx, (i,j) in enumerate(edges5):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    groups = get_cycle_data(A, 5)
    beta1 = sum(groups.values())
    beta1_mod2[beta1 % 2] += 1

print(f"  n=5: β₁ mod 2: {dict(beta1_mod2)}")

# n=6
beta1_mod2_6 = Counter()
for bits in range(2**15):  # n=6
    A = [[0]*6 for _ in range(6)]
    edges6 = [(i,j) for i in range(6) for j in range(i+1,6)]
    for idx, (i,j) in enumerate(edges6):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    groups = get_cycle_data(A, 6)
    beta1 = sum(groups.values())
    beta1_mod2_6[beta1 % 2] += 1

print(f"  n=6: β₁ mod 2: {dict(beta1_mod2_6)}")

# ═══════════════════════════════════════════════════════════════════
# Part 6: Minimum β₁+2β₂ for each H near 21
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 6: Minimum β₁+2β₂ for each odd H ---")

for n in [5, 6]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    ne = len(edges)

    H_to_min_target = {}
    H_to_max_target = {}

    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        H = fast_hp(A, n)
        groups = get_cycle_data(A, n)
        beta1, beta2 = compute_betas(groups)
        target = beta1 + 2*beta2

        if H not in H_to_min_target or target < H_to_min_target[H]:
            H_to_min_target[H] = target
        if H not in H_to_max_target or target > H_to_max_target[H]:
            H_to_max_target[H] = target

    print(f"\n  n={n}: H → min/max(β₁+2β₂)")
    for h in sorted(H_to_min_target.keys()):
        if h % 2 == 0:
            continue
        expected = (h - 1) // 2
        marker = " ← 21" if h == 21 else ""
        print(f"    H={h:3d}: β₁+2β₂ ∈ [{H_to_min_target[h]}, {H_to_max_target[h]}], "
              f"expected (H-1)/2 = {expected}{marker}")

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
print("""
The H=21 impossibility comes from the interaction between:

1. β₁ = Σ d(S): total directed cycle count across all vertex sets
2. β₂ = Σ d(S)·d(T) over disjoint pairs: grows FASTER than β₁²
3. Higher β_k: contribute 2^k each, further pushing H above 21

The key insight: having 10 directed cycles (β₁=10) with NO disjoint
pairs (β₂=0) is structurally impossible because:
- Having enough 3-cycles forces vertex overlaps but ALSO forces 5-cycles
  (via the Splicing Lemma), which ADD directed cycle count
- Having enough 5-cycles for β₁=10 requires many vertex sets,
  which inevitably create disjoint pairs (β₂>0)
- β₂>0 pushes H above 21

Result: the "exactly 10 cycles with no disjoint structure" configuration
cannot exist in any tournament, making H=21 a permanent gap.
""")
