#!/usr/bin/env python3
"""
MOMENT HIERARCHY at n=5: complete analysis.

At n=5, all 12 isomorphism classes can be enumerated.
sum_P f^j for each j: what invariant does it depend on?

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

def count_5_cycles(A, n):
    count = 0
    for p in permutations(range(n)):
        if all(A[p[i]][p[(i+1)%n]] for i in range(n)):
            count += 1
    return count // n

def ham_count(A, n):
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] for i in range(n-1)))

n = 5
print("=" * 70)
print(f"MOMENT HIERARCHY at n={n}")
print("=" * 70)

# Enumerate ALL tournaments up to isomorphism by checking all 2^10 = 1024
seen = set()
iso_classes = []
for bits in range(1 << 10):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles(A, n)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    H = ham_count(A, n)

    moments = [0] * (n)
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        fpow = 1
        for j in range(n):
            moments[j] += fpow
            fpow *= f

    key = (scores, t3, t5, H, tuple(moments))
    if key not in seen:
        seen.add(key)
        iso_classes.append({
            'scores': scores, 't3': t3, 't5': t5, 'H': H,
            'moments': moments
        })

print(f"\n{len(iso_classes)} distinct (score, t3, t5, H, moments) tuples")
print(f"\n{'scores':<16} t3  t5   H  | m0    m1     m2       m3        m4")
for d in sorted(iso_classes, key=lambda x: (x['t3'], x['t5'])):
    m = d['moments']
    print(f"{str(d['scores']):<16} {d['t3']:2d}  {d['t5']:2d}  {d['H']:2d}  | "
          f"{m[0]:3d}  {m[1]:4d}  {m[2]:5d}  {m[3]:7d}  {m[4]:8d}")

# Check what m2, m3, m4 depend on
print(f"\nDependence analysis:")
import numpy as np

for j in range(2, n):
    y = np.array([d['moments'][j] for d in iso_classes], dtype=float)

    # Linear in t3
    X_t3 = np.column_stack([
        [d['t3'] for d in iso_classes],
        [1 for _ in iso_classes],
    ])
    coeffs_t3, _, _, _ = np.linalg.lstsq(X_t3, y, rcond=None)
    err_t3 = max(abs(X_t3 @ coeffs_t3 - y))

    # Linear in (t3, t5)
    X_t35 = np.column_stack([
        [d['t3'] for d in iso_classes],
        [d['t5'] for d in iso_classes],
        [1 for _ in iso_classes],
    ])
    coeffs_t35, _, _, _ = np.linalg.lstsq(X_t35, y, rcond=None)
    err_t35 = max(abs(X_t35 @ coeffs_t35 - y))

    # Linear in (t3, t5, H)
    X_full = np.column_stack([
        [d['t3'] for d in iso_classes],
        [d['t5'] for d in iso_classes],
        [d['H'] for d in iso_classes],
        [1 for _ in iso_classes],
    ])
    coeffs_full, _, _, _ = np.linalg.lstsq(X_full, y, rcond=None)
    err_full = max(abs(X_full @ coeffs_full - y))

    print(f"\n  j={j}: sum_P f^{j}")
    print(f"    f(t3): err={err_t3:.4f}")
    if err_t3 < 0.01:
        print(f"      = {coeffs_t3[0]:.0f}*t3 + {coeffs_t3[1]:.0f}")
    else:
        print(f"    f(t3,t5): err={err_t35:.4f}")
        if err_t35 < 0.01:
            print(f"      = {coeffs_t35[0]:.0f}*t3 + {coeffs_t35[1]:.0f}*t5 + {coeffs_t35[2]:.0f}")
        else:
            print(f"    f(t3,t5,H): err={err_full:.4f}")
            if err_full < 0.01:
                print(f"      = {coeffs_full[0]:.0f}*t3 + {coeffs_full[1]:.0f}*t5 + {coeffs_full[2]:.0f}*H + {coeffs_full[3]:.0f}")

# KEY: at n=5, f^4 = f*(f-1)*(f-2)*(f-3) + lower terms
# The only way to get product of ALL 4 distinct positions is {0,1,2,3}
# which = directed Hamiltonian path indicator
# So the "pure 4th moment" contribution is sum_P T_0*T_1*T_2*T_3 = H
# with multinomial coefficient 4!/1^4 = 24
print(f"\n{'='*70}")
print("ALGEBRAIC STRUCTURE")
print(f"{'='*70}")
print(f"""
  At n=5: n-1=4 edges per permutation.
  f^4 = (T_0+T_1+T_2+T_3)^4

  The product T_0*T_1*T_2*T_3 has multinomial coefficient 4! = 24.
  sum_P T_0*T_1*T_2*T_3 = H (Ham path count).

  So: m4 = 24*H + (lower terms depending on t3)

  Verify: m4 - 24*H should depend only on t3.
""")

y_adj = np.array([d['moments'][4] - 24*d['H'] for d in iso_classes], dtype=float)
X_t3_only = np.column_stack([
    [d['t3'] for d in iso_classes],
    [1 for _ in iso_classes],
])
coeffs_adj, _, _, _ = np.linalg.lstsq(X_t3_only, y_adj, rcond=None)
err_adj = max(abs(X_t3_only @ coeffs_adj - y_adj))
print(f"  m4 - 24*H = {coeffs_adj[0]:.0f}*t3 + {coeffs_adj[1]:.0f} (err={err_adj:.4f})")
print(f"  => m4 = 24*H + {coeffs_adj[0]:.0f}*t3 + {coeffs_adj[1]:.0f}")

# Verify
for d in iso_classes[:5]:
    m4 = d['moments'][4]
    pred = 24*d['H'] + coeffs_adj[0]*d['t3'] + coeffs_adj[1]
    print(f"  t3={d['t3']}, t5={d['t5']}, H={d['H']}: m4={m4}, pred={pred:.0f}")

# The full picture at n=5
print(f"\nFULL MOMENT TABLE at n=5:")
print(f"  m0 = {factorial(n)} (universal)")
print(f"  m1 = {iso_classes[0]['moments'][1]} (universal)")

# m2 coefficient
m2_coeffs = np.linalg.lstsq(X_t3_only,
    np.array([d['moments'][2] for d in iso_classes], dtype=float), rcond=None)[0]
print(f"  m2 = {m2_coeffs[0]:.0f}*t3 + {m2_coeffs[1]:.0f}")

# m3 coefficient
m3_coeffs = np.linalg.lstsq(X_t3_only,
    np.array([d['moments'][3] for d in iso_classes], dtype=float), rcond=None)[0]
print(f"  m3 = {m3_coeffs[0]:.0f}*t3 + {m3_coeffs[1]:.0f}")

print(f"  m4 = 24*H + {coeffs_adj[0]:.0f}*t3 + {coeffs_adj[1]:.0f}")

# At n=5, H = 1 + 2*(t3+t5) by OCF (since alpha_2=0)
# So m4 = 24*(1+2*t3+2*t5) + a*t3 + b = 48*t3 + 48*t5 + (24+b) + a*t3
# = (48+a)*t3 + 48*t5 + (24+b)
print(f"\n  Since H = 1 + 2*(t3+t5) at n=5:")
a, b = coeffs_adj[0], coeffs_adj[1]
print(f"  m4 = 24*(1+2*t3+2*t5) + {a:.0f}*t3 + {b:.0f}")
print(f"     = {48+a:.0f}*t3 + 48*t5 + {24+b:.0f}")

# Compare with direct fit
m4_fit = np.linalg.lstsq(
    np.column_stack([
        [d['t3'] for d in iso_classes],
        [d['t5'] for d in iso_classes],
        [1 for _ in iso_classes],
    ]),
    np.array([d['moments'][4] for d in iso_classes], dtype=float),
    rcond=None
)[0]
print(f"  Direct fit: m4 = {m4_fit[0]:.0f}*t3 + {m4_fit[1]:.0f}*t5 + {m4_fit[2]:.0f}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
