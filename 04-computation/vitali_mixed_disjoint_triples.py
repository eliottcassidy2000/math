"""
vitali_mixed_disjoint_triples.py -- kind-pasteur-2026-03-13-S61
At n=10: disjoint triples of type (c3,c3,c5) become possible (3+3+5=11>10... no!)
Actually: need disjoint, so 3+3+5 = 11 > 10. NOT possible at n=10.
At n=11: (c3,c3,c5) needs 3+3+5 = 11 = n. FIRST TIME.

Wait, but (c3,c5) pairs already work at n=8 (3+5=8=n).
And (c3,c3) at n=6 (3+3=6=n).
The i_k channels are:
  i_1 = number of directed odd cycles = alpha_1
  i_2 = sum of products over disjoint PAIRS
  i_3 = sum of products over disjoint TRIPLES

i_3 includes MIXED triples like (c3,c3,c5), (c3,c5,c5), etc.
At n=10: can we have disjoint triples?
  (c3,c3,c3): 3+3+3=9 <= 10. YES.
  (c3,c3,c5): 3+3+5=11 > 10. NO.
  Everything else: even larger. NO.

So at n=10, delta_i3 counts only disjoint c3 triples.
And THM-171 guarantees delta_(c3,c3,c3) = 0.
Therefore delta_i3 = 0 at n=10.

At n=11: (c3,c3,c5) triples become possible. 3+3+5=11=n.
These are NOT purely c3, so THM-171 doesn't force delta=0.
The c5 multiplicities CAN change, so the product m_c3 * m_c3 * m_c5
can change. This would open the i3 channel!

But dc5 = 0 is still hypothesized... and the c3 identities are
lambda-determined. So the (c3,c3,c5) triple product changes
would come from c5 multiplicity changes on sets disjoint from
the c3 pair. This is exactly the i2-type mechanism but one level up.

Let me check: at n=10, does the i3 channel vanish?
(It should: only c3 triples possible, and THM-171 forces delta=0.)
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import time

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_cycles_by_set(A, n, length):
    result = {}
    for combo in combinations(range(n), length):
        count = 0
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for k in range(length):
                if A[path[k]][path[(k+1) % length]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
        if count > 0:
            result[frozenset(combo)] = count
    return result

def count_disjoint_directed_triples_by_type(cd, lengths, n):
    """Count disjoint directed triples broken down by type."""
    results = {}
    for i, l1 in enumerate(lengths):
        for j, l2 in enumerate(lengths):
            if l2 < l1:
                continue
            for k, l3 in enumerate(lengths):
                if l3 < l2:
                    continue
                if l1 + l2 + l3 > n:
                    continue

                d1 = cd.get(l1, {})
                d2 = cd.get(l2, {})
                d3 = cd.get(l3, {})
                vs1 = list(d1.keys())
                vs2 = list(d2.keys())
                vs3 = list(d3.keys())

                total = 0
                for a in range(len(vs1)):
                    for b in range(len(vs2)):
                        if l1 == l2 and b <= a:
                            continue
                        if vs1[a] & vs2[b]:
                            continue
                        for c in range(len(vs3)):
                            if l2 == l3 and c <= b:
                                continue
                            if l1 == l3 and l2 != l3 and c <= a:
                                continue
                            if (vs1[a] & vs3[c]) or (vs2[b] & vs3[c]):
                                continue
                            total += d1[vs1[a]] * d2[vs2[b]] * d3[vs3[c]]

                if total > 0 or (l1, l2, l3) in [(3,3,3)]:
                    results[(l1, l2, l3)] = total

    return results

n = 10
total_bits = n * (n - 1) // 2

print("=" * 70)
print(f"MIXED DISJOINT TRIPLES AT n={n}")
print(f"Possible triple types:")
for l1 in [3, 5, 7]:
    for l2 in [l1, l1+2, l1+4]:
        if l2 > 9:
            continue
        for l3 in [l2, l2+2, l2+4]:
            if l3 > 9:
                continue
            if l1 + l2 + l3 <= n:
                print(f"  ({l1},{l2},{l3}): needs {l1+l2+l3} vertices {'<=' if l1+l2+l3<=n else '>'} {n}")
print("=" * 70)

np.random.seed(2026)
data = []
t0 = time.time()

for trial in range(2000):
    bits = 0
    for chunk in range(0, total_bits, 30):
        chunk_size = min(30, total_bits - chunk)
        bits |= int(np.random.randint(0, 1 << chunk_size)) << chunk

    A = bits_to_adj(bits, n)
    lam = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(lam, lambda_graph(B, n)):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        if H_A == H_B:
            continue

        lengths = [3, 5, 7, 9]
        cd_A = {}
        cd_B = {}
        dc = {}
        for l in lengths:
            c_A = count_directed_cycles_by_set(A, n, l)
            c_B = count_directed_cycles_by_set(B, n, l)
            cd_A[l] = c_A
            cd_B[l] = c_B
            dc[l] = sum(c_B.values()) - sum(c_A.values())

        # Count i2
        def total_disjoint_pairs(cd):
            all_vs = {}
            for l, d in cd.items():
                for vs, m in d.items():
                    all_vs[vs] = m
            vlist = list(all_vs.keys())
            total = 0
            for a in range(len(vlist)):
                for b in range(a+1, len(vlist)):
                    if not (vlist[a] & vlist[b]):
                        total += all_vs[vlist[a]] * all_vs[vlist[b]]
            return total

        i2_A = total_disjoint_pairs(cd_A)
        i2_B = total_disjoint_pairs(cd_B)
        di2 = i2_B - i2_A

        # Count i3 by type
        i3_by_type_A = count_disjoint_directed_triples_by_type(cd_A, [3, 5, 7], n)
        i3_by_type_B = count_disjoint_directed_triples_by_type(cd_B, [3, 5, 7], n)

        all_types = set(i3_by_type_A.keys()) | set(i3_by_type_B.keys())
        di3_by_type = {}
        for t in all_types:
            di3_by_type[t] = i3_by_type_B.get(t, 0) - i3_by_type_A.get(t, 0)

        total_di3 = sum(di3_by_type.values())

        delta_H = H_B - H_A
        predicted = 2*(dc[7]+dc[9]) + 4*di2 + 8*total_di3

        data.append({
            'delta_H': delta_H,
            'dc': dc,
            'di2': di2,
            'di3_by_type': di3_by_type,
            'total_di3': total_di3,
            'predicted': predicted,
            'i3_A': i3_by_type_A,
            'i3_B': i3_by_type_B,
        })
        break

    if trial % 500 == 0 and trial > 0:
        elapsed = time.time() - t0
        print(f"  ... {trial}/2000, found {len(data)}, elapsed {elapsed:.0f}s")

elapsed = time.time() - t0
print(f"\nTotal: {len(data)} H-changing examples in {elapsed:.0f}s")

if data:
    # Formula check
    matches = sum(1 for d in data if d['predicted'] == d['delta_H'])
    print(f"\nFormula delta_H = 2*(dc7+dc9) + 4*di2 + 8*di3:")
    print(f"  Match rate: {matches}/{len(data)}")

    # Without di3:
    matches2 = sum(1 for d in data if 2*(d['dc'][7]+d['dc'][9]) + 4*d['di2'] == d['delta_H'])
    print(f"Without di3 term: {matches2}/{len(data)}")

    # i3 analysis
    print(f"\n  di3 total distribution: {dict(sorted(Counter(d['total_di3'] for d in data).items()))}")

    for t in sorted(set(t for d in data for t in d['di3_by_type'])):
        vals = [d['di3_by_type'].get(t, 0) for d in data]
        nonzero = sum(1 for v in vals if v != 0)
        print(f"  di3 type {t}: nonzero {nonzero}/{len(data)}, dist={dict(sorted(Counter(vals).items()))}")

    # Show i3 absolute values
    for t in sorted(set(t for d in data for t in d['i3_A'])):
        vals_A = [d['i3_A'].get(t, 0) for d in data]
        print(f"  i3_A type {t}: min={min(vals_A)}, max={max(vals_A)}, mean={np.mean(vals_A):.1f}")

print("\nDone.")
