#!/usr/bin/env python3
"""
worpitzky_coefficients.py — Worpitzky-type expansion of F(T,x).

F(T,x) / (1-x)^n = sum_{m>=0} a_m x^m  (formal power series)

where a_m = sum_{k=0}^{n-1} F_k * C(m+n-1-k, n-1)

Since C(m+n-1-k, n-1) is polynomial in m of degree n-1, a_m is polynomial in m.

DISCOVERY (S46): The polynomial a_m has UNIVERSAL top coefficients:
  - Leading coeff (m^{n-1}): always n (related to n!/n!)
  - Second coeff (m^{n-2}): always C(n,2)
  - Third coeff (m^{n-3}): C(n,2) + 6*t3  (depends on 3-cycle count!)

This means the Worpitzky polynomial encodes tournament invariants.

QUESTION: Does the full polynomial determine F(T,x) uniquely?
(Answer: yes — it IS F(T,x) in a different basis, so same info.)
But which invariants appear at each level?

Author: opus-2026-03-07-S46
"""
from itertools import permutations, combinations
from math import comb, factorial
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

def count_3cycles(adj, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

def worpitzky_a(F, n, m):
    """Compute a_m = sum_k F_k * C(m+n-1-k, n-1)."""
    return sum(F[k] * comb(m + n - 1 - k, n - 1) for k in range(n))

# ============================================================
# VERIFY WORPITZKY POLYNOMIAL STRUCTURE
# ============================================================
print("=" * 60)
print("WORPITZKY POLYNOMIAL a_m = P(m) for F(T,x)")
print("=" * 60)

for n in [4, 5, 6]:
    m_vals = n*(n-1)//2
    seen = set()
    results = []

    import random
    random.seed(42)
    num = min(1 << m_vals, 100000)

    for trial in range(num):
        if n <= 5:
            bits = trial
        else:
            bits = random.getrandbits(m_vals)

        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        t3 = count_3cycles(adj, n)

        # Compute a_m at enough points to determine polynomial
        # a_m is degree n-1 polynomial in m, need n points
        m_points = list(range(n + 2))
        a_vals = [worpitzky_a(F, n, m) for m in m_points]

        # Fit polynomial using numpy
        coeffs = np.polyfit(m_points, a_vals, n - 1)
        # coeffs[0] = leading, coeffs[-1] = constant

        # Normalize by 1/(n-1)! to get nice coefficients
        norm_coeffs = [c / factorial(n-1) for c in coeffs]

        results.append({
            'F': F, 't3': t3, 'coeffs': coeffs, 'norm': norm_coeffs,
            'a_vals': a_vals[:5]
        })

    print(f"\nn={n}: {len(results)} distinct F-vectors")

    # Analyze which coefficients are universal
    for deg_idx in range(n):
        vals = set(round(r['coeffs'][deg_idx], 4) for r in results)
        if len(vals) == 1:
            print(f"  coeff of m^{n-1-deg_idx}: UNIVERSAL = {vals.pop():.4f}")
        else:
            # Check if it depends on t3
            t3_to_val = {}
            for r in results:
                t3 = r['t3']
                val = round(r['coeffs'][deg_idx], 4)
                if t3 not in t3_to_val:
                    t3_to_val[t3] = set()
                t3_to_val[t3].add(val)

            # Check if t3 determines this coefficient
            unique_per_t3 = all(len(v) == 1 for v in t3_to_val.values())
            if unique_per_t3:
                print(f"  coeff of m^{n-1-deg_idx}: DETERMINED BY t3")
                for t3_val in sorted(t3_to_val.keys()):
                    c_val = t3_to_val[t3_val].pop()
                    print(f"    t3={t3_val}: coeff={c_val}")
                # Try linear fit: coeff = A + B*t3
                t3_list = sorted(t3_to_val.keys())
                if len(t3_list) >= 2:
                    c0 = list(t3_to_val[t3_list[0]])[0] if t3_to_val[t3_list[0]] else 0
                    c1 = list(t3_to_val[t3_list[1]])[0] if t3_to_val[t3_list[1]] else 0
                    if t3_list[1] != t3_list[0]:
                        slope = (c1 - c0) / (t3_list[1] - t3_list[0])
                        intercept = c0 - slope * t3_list[0]
                        # Verify fit for all t3 values
                        fit_ok = all(
                            abs(list(t3_to_val[t])[0] - (intercept + slope * t)) < 0.01
                            for t in t3_list if t3_to_val[t]
                        )
                        if fit_ok:
                            print(f"    Linear fit: coeff = {intercept:.4f} + {slope:.4f}*t3 (EXACT)")
                        else:
                            print(f"    Linear fit: coeff ~ {intercept:.4f} + {slope:.4f}*t3 (approximate)")
            else:
                print(f"  coeff of m^{n-1-deg_idx}: NOT determined by t3 alone ({len(vals)} distinct values)")
                # Show how many distinct values per t3
                for t3_val in sorted(t3_to_val.keys()):
                    n_distinct = len(t3_to_val[t3_val])
                    if n_distinct > 1:
                        print(f"    t3={t3_val}: {n_distinct} distinct values")

    # Show examples
    print("\n  Examples:")
    for r in results[:4]:
        print(f"    F={r['F']}, t3={r['t3']}")
        print(f"    a(0..4)={r['a_vals']}")
        print(f"    poly coeffs (high to low): {[f'{c:.2f}' for c in r['coeffs']]}")

# ============================================================
# DEEPER ANALYSIS: WHAT INVARIANT AT EACH LEVEL?
# ============================================================
print("\n" + "=" * 60)
print("WHICH INVARIANTS DETERMINE EACH COEFFICIENT?")
print("=" * 60)

# For n=5, compute additional tournament invariants
n = 5
m_vals = n*(n-1)//2
seen = set()
full_data = []

for bits in range(1 << m_vals):
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    t3 = count_3cycles(adj, n)

    # Score sequence (sorted out-degrees)
    scores = sorted([sum(adj[i][j] for j in range(n) if j != i) for i in range(n)])

    # Worpitzky polynomial coefficients
    m_points = list(range(n + 2))
    a_vals = [worpitzky_a(F, n, m) for m in m_points]
    coeffs = np.polyfit(m_points, a_vals, n - 1)

    full_data.append({
        'F': F, 't3': t3, 'scores': tuple(scores),
        'coeffs': coeffs, 'key': key
    })

print(f"\nn=5: {len(full_data)} F-classes")
print("\nFull coefficient table:")
print(f"{'F':<30} {'t3':>3} {'scores':<16} {'c4':>8} {'c3':>8} {'c2':>8} {'c1':>8} {'c0':>8}")
for d in sorted(full_data, key=lambda x: x['t3']):
    F_str = str(d['F'])
    c = d['coeffs']
    print(f"{F_str:<30} {d['t3']:>3} {str(list(d['scores'])):<16} {c[0]:>8.2f} {c[1]:>8.2f} {c[2]:>8.2f} {c[3]:>8.2f} {c[4]:>8.2f}")

# ============================================================
# WORPITZKY + DELETION-CONTRACTION CONNECTION
# ============================================================
print("\n" + "=" * 60)
print("WORPITZKY + DELETION-CONTRACTION")
print("=" * 60)

# If F_T(x) = F_{T\e}(x) + (x-1)*F(T/e, x), what does this mean
# for the Worpitzky coefficients?
#
# a_m(T) = sum_k F_T[k] * C(m+n-1-k, n-1)
# a_m(T\e) uses same formula but F_{T\e}[k] coefficients
# F(T/e, x) has degree n-2, so its Worpitzky uses C(m+n-2-k, n-2)
#
# The (x-1) multiplication in the F-world translates to...
# In Worpitzky: x * sum a_m x^m = sum a_{m-1} x^m
# So (x-1) * sum b_m x^m = sum (b_{m-1} - b_m) x^m
#
# Therefore: a_m(T) = a_m(T\e) + a_{m-1}(T/e) - a_m(T/e)

print("Verifying: a_m(T) = a_m(T\\e) + a_{m-1}(T/e) - a_m(T/e)")

n = 5
ok = 0
total = 0
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)
    key = tuple(F_T)

    # Test for first edge (0->? or ?->0)
    u, v = 0, 1
    if not adj[u][v]:
        u, v = v, u

    # Deletion: remove edge u->v
    adj_del = [row[:] for row in adj]
    adj_del[u][v] = 0

    # Contraction: merge u,v
    # w inherits IN from u, OUT from v
    remaining = [i for i in range(n) if i != u and i != v]
    nc = n - 1
    adj_con = [[0]*nc for _ in range(nc)]
    # vertex 0 = w (merged), vertices 1..nc-1 = remaining
    for ii, a in enumerate(remaining):
        for jj, b in enumerate(remaining):
            if a != b:
                adj_con[ii+1][jj+1] = adj[a][b]
        # edges involving w
        adj_con[0][ii+1] = adj[v][a]  # w->a = v->a (OUT from v)
        adj_con[ii+1][0] = adj[a][u]  # a->w = a->u (IN from u)

    F_del = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj_del[P[i]][P[i+1]])
        F_del[fwd] += 1

    F_con = compute_F(adj_con, nc)

    # Check Worpitzky identity for several m values
    for m in [0, 1, 2, 3, 5]:
        a_T = worpitzky_a(F_T, n, m)
        a_del = worpitzky_a(F_del, n, m)
        # For contraction, use n-1
        a_con_m = worpitzky_a(F_con, nc, m)
        a_con_m1 = worpitzky_a(F_con, nc, m-1) if m > 0 else 0

        lhs = a_T
        rhs = a_del + a_con_m1 - a_con_m
        total += 1
        if abs(lhs - rhs) < 0.01:
            ok += 1

    if bits >= 100:
        break

print(f"  Checked {ok}/{total} (first 100 tournaments)")

# ============================================================
# EULERIAN NUMBER CONNECTION
# ============================================================
print("\n" + "=" * 60)
print("COMPARISON: EULERIAN vs TOURNAMENT WORPITZKY")
print("=" * 60)

# The STANDARD Worpitzky identity: sum_{k=0}^{n-1} A(n,k+1) * C(m+n-1-k, n-1) = (m+1)^n - 1? No.
# Actually: x^n = sum_k A(n,k) C(x+n-1-k, n) ... wrong form.
# Standard: sum_k A(n,k) C(x+k, n) = x^n
# Or: (m+1)^n = sum_{k=0}^{n-1} A(n,k+1) * C(m+k, n) -- hm
#
# Actually the standard Worpitzky identity is:
# sum_{k=0}^{n-1} <n,k> * C(x+k, n-1) = x^n  (Eulerian numbers <n,k>)
# where <n,k> = sum_{j=0}^k (-1)^j C(n+1,j) (k+1-j)^n
#
# Let's just compute: what is sum_k A_n[k] * C(m+n-1-k, n-1) ?

for n in [4, 5]:
    # Eulerian numbers A_n[k] for k=0,...,n-1
    # A_n[k] = number of perms of [n] with exactly k descents
    from math import perm
    A_n = [0]*n
    for P in permutations(range(n)):
        desc = sum(1 for i in range(n-1) if P[i] > P[i+1])
        A_n[desc] += 1

    print(f"\nn={n}: Eulerian numbers = {A_n}")

    # Worpitzky of Eulerian = ?
    for m in range(6):
        val = sum(A_n[k] * comb(m + n - 1 - k, n - 1) for k in range(n))
        print(f"  a({m}) = {val}  (compare (m+1)^n = {(m+1)**n})")

# ============================================================
# NORMALIZED WORPITZKY: F(T) vs Eulerian
# ============================================================
print("\n" + "=" * 60)
print("NORMALIZED: a_m(T) / a_m(Eulerian)")
print("=" * 60)

n = 5
A_n = [0]*n
for P in permutations(range(n)):
    desc = sum(1 for i in range(n-1) if P[i] > P[i+1])
    A_n[desc] += 1

print(f"Eulerian A_5 = {A_n}")
print(f"Sum = {sum(A_n)} = 5! = {factorial(5)}")
# Worpitzky of Eulerian
eul_a = [sum(A_n[k] * comb(m + n - 1 - k, n - 1) for k in range(n)) for m in range(8)]
print(f"Eulerian a(m) = {eul_a}")

for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)

    if key == tuple(A_n):
        print(f"\n  *** Found tournament with F = Eulerian! bits={bits}")
        # Is there such a tournament?
        break
else:
    print("\n  No tournament at n=5 has F = Eulerian polynomial")
    print(f"  Eulerian: {A_n}")
    print(f"  Tournament F-vectors:")
    seen = set()
    for bits in range(1 << (n*(n-1)//2)):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key not in seen:
            seen.add(key)
            print(f"    {list(key)}")
