#!/usr/bin/env python3
"""
CHARACTERIZE: Which tournaments have M = (H/n)*I?

Findings so far:
  - ALL circulant tournaments: YES
  - ALL VT tournaments at n=3,5,7: YES
  - Some non-circulant regular with H=175: YES
  - Regular with H=171: NO

Is M scalar <==> vertex-transitive?

kind-pasteur-2026-03-06-S25b (continuation)
"""

from itertools import permutations, combinations
import numpy as np

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    ea += prod
            bb2 = 0
            if len(R_set) == 1:
                bb2 = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    bb2 += prod
            val += sign * ea * bb2
        return val

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        M[a, a] = compute_M_entry(T, n, a, a)
        for b in range(a+1, n):
            M[a, b] = compute_M_entry(T, n, a, b)
            M[b, a] = M[a, b]
    return M

def check_automorphism(T, n, sigma):
    """Check if sigma is an automorphism of T."""
    for i in range(n):
        for j in range(n):
            if i != j:
                if T.get((i,j), 0) != T.get((sigma[i], sigma[j]), 0):
                    return False
    return True

def aut_group_size(T, n):
    """Count automorphisms of T (brute force)."""
    count = 0
    for perm in permutations(range(n)):
        if check_automorphism(T, n, perm):
            count += 1
    return count

def is_vertex_transitive(T, n):
    """Check if Aut(T) acts transitively on vertices."""
    # For each pair (0, v), check if some automorphism maps 0 to v
    for v in range(1, n):
        found = False
        for perm in permutations(range(n)):
            if perm[0] != v:
                continue
            if check_automorphism(T, n, perm):
                found = True
                break
        if not found:
            return False
    return True


# ============================================================
# n=5: Exhaustive classification
# ============================================================
print("=" * 70)
print("n=5: M scalar vs vertex-transitive (exhaustive)")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

# Group by isomorphism class (by canonical form)
seen = set()
results = []

for bits in range(1 << len(pairs)):
    b_list = [(bits >> k) & 1 for k in range(len(pairs))]
    T = {}
    idx = 0
    for (i,j) in pairs:
        if b_list[idx]:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
        idx += 1

    # Canonical: score sequence + H
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    H = count_H(T, n)
    M = compute_M(T, n)
    scalar = H // n if H % n == 0 else None
    is_scalar = scalar is not None and np.array_equal(M, scalar * np.eye(n, dtype=int))

    # Only process one representative per isomorphism class
    key = (scores, H, is_scalar)
    if key in seen:
        continue
    seen.add(key)

    vt = is_vertex_transitive(T, n)
    aut = aut_group_size(T, n)

    results.append({
        'scores': scores, 'H': H, 'is_scalar': is_scalar,
        'vt': vt, 'aut': aut, 'M': M
    })

print(f"\n  {len(results)} isomorphism classes:")
for r in sorted(results, key=lambda x: x['H']):
    vt_str = "VT" if r['vt'] else "not-VT"
    scalar_str = f"M={r['H']//r['scores'].count(r['scores'][0])}*I" if r['is_scalar'] else "M not scalar"
    print(f"    scores={r['scores']}: H={r['H']}, |Aut|={r['aut']}, {vt_str}, {scalar_str}")


# ============================================================
# n=7: Regular tournaments — VT vs scalar M
# ============================================================
print("\n" + "=" * 70)
print("n=7: Regular (3-regular) tournaments — VT vs scalar M")
print("=" * 70)

import random
random.seed(12345)

n = 7
pairs_7 = [(i,j) for i in range(n) for j in range(i+1, n)]

tested = 0
vt_scalar = 0
vt_not_scalar = 0
not_vt_scalar = 0
not_vt_not_scalar = 0
examples = {}

for trial in range(200000):
    T = {}
    for (i,j) in pairs_7:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    scores = tuple(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n))
    if not all(s == 3 for s in scores):
        continue

    H = count_H(T, n)

    # Quick scalar check: just compute M[0,0] and M[0,1]
    M00 = compute_M_entry(T, n, 0, 0)
    M01 = compute_M_entry(T, n, 0, 1)

    if H % n == 0 and M01 == 0:
        # Likely scalar, verify fully
        M = compute_M(T, n)
        is_scalar = np.array_equal(M, (H // n) * np.eye(n, dtype=int))
    else:
        is_scalar = False

    vt = is_vertex_transitive(T, n)

    if vt and is_scalar:
        vt_scalar += 1
    elif vt and not is_scalar:
        vt_not_scalar += 1
        if 'vt_not_scalar' not in examples:
            examples['vt_not_scalar'] = (H, T.copy())
            print(f"  COUNTEREXAMPLE: VT but M not scalar! H={H}")
    elif not vt and is_scalar:
        not_vt_scalar += 1
        if 'not_vt_scalar' not in examples:
            examples['not_vt_scalar'] = (H, T.copy())
            print(f"  COUNTEREXAMPLE: scalar M but not VT! H={H}")
    else:
        not_vt_not_scalar += 1

    tested += 1
    if tested >= 100:
        break

print(f"\n  Tested {tested} regular tournaments at n=7:")
print(f"    VT + scalar:      {vt_scalar}")
print(f"    VT + not scalar:  {vt_not_scalar}")
print(f"    not VT + scalar:  {not_vt_scalar}")
print(f"    not VT + not scalar: {not_vt_not_scalar}")

if vt_not_scalar == 0 and not_vt_scalar == 0:
    print("\n  CONJECTURE CONFIRMED: M scalar <==> VT (at n=7, 100 samples)")
elif vt_not_scalar == 0:
    print(f"\n  VT => scalar: TRUE (but scalar does NOT imply VT)")
elif not_vt_scalar == 0:
    print(f"\n  scalar => VT: TRUE (but VT does NOT imply scalar)")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
