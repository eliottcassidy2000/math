#!/usr/bin/env python3
"""
SEESAW BOUND: Why β₁+β₃ ≤ 1?
================================
We proved: β₁+β₃ = ker(d₁) + Ω₃ - Ω₂ - im(d₄)
         = [C(n,2)-n+1] + Ω₃ - Ω₂ - im(d₄)

So β₁+β₃ ≤ 1 iff Ω₂ - Ω₃ + im(d₄) ≥ C(n,2)-n = E - n + 1 - 1 = E - n

where E = C(n,2) = number of edges.

Key: Ω₂ - Ω₃ = dim(Ω₂) - dim(Ω₃) relates to the Euler characteristic
of the chain complex at levels 2-3.

Let's also compute the ALTERNATING SUM (Euler char of Ω complex)
and see if there's a universal formula for Ω₂ - Ω₃.

We know: rank(d₂) + ker(d₂) = Ω₂
         rank(d₃) + ker(d₃) = Ω₃
         β₂ = 0 => ker(d₂) = rank(d₃) = im(d₃)

So: Ω₂ - Ω₃ = rank(d₂) + ker(d₂) - rank(d₃) - ker(d₃)
             = rank(d₂) + rank(d₃) - rank(d₃) - ker(d₃)
             = rank(d₂) - ker(d₃)

And: β₁+β₃ = ker(d₁) + Ω₃ - Ω₂ - im(d₄)
            = ker(d₁) - rank(d₂) + ker(d₃) - im(d₄)
            = β₁ + β₃  ✓  (this is just the definition!)

Hmm, so the identity is tautological given β₂=0. The question is whether
the RHS is always ≤ 1.

Let me try a DIFFERENT approach: track how dim(Ω₂) and dim(Ω₃) vary
with the tournament, and see if there's a universal bound.

opus-2026-03-15-S72d
"""
import numpy as np
from collections import Counter
import time

# Inline path homology (same as before)
def enum_paths(A, n, p):
    if p == 0: return [(v,) for v in range(n)]
    if p < 0: return []
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j]: adj[i].append(j)
    paths = []; stack = []
    for s in range(n):
        stack.append(([s], 1<<s))
        while stack:
            path, vis = stack.pop()
            if len(path)==p+1: paths.append(tuple(path)); continue
            for u in adj[path[-1]]:
                if not (vis & (1<<u)):
                    stack.append((path+[u], vis|(1<<u)))
    return paths

def bcoeffs(path):
    return [((-1)**i, path[:i]+path[i+1:]) for i in range(len(path))]

def omega_dim(A, n, p):
    """Compute just dim(Ω_p)."""
    al_p = enum_paths(A, n, p)
    al_pm = enum_paths(A, n, p-1)
    d = len(al_p)
    if d == 0: return 0
    if p == 0: return d
    aset = set(al_pm); naf = {}; c = 0
    for j, path in enumerate(al_p):
        for s, f in bcoeffs(path):
            if len(set(f)) == len(f) and f not in aset:
                if f not in naf: naf[f] = c; c += 1
    if c == 0: return d
    P = np.zeros((c, d))
    for j, path in enumerate(al_p):
        for s, f in bcoeffs(path):
            if f in naf: P[naf[f], j] += s
    S = np.linalg.svd(P, compute_uv=False)
    r = sum(s > 1e-10 for s in S)
    return d - r

def tourn(n, bits):
    A = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1<<idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                w = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
                if max(w) < 2: c3 += 1
    return c3

# === MAIN ===
print("=" * 70)
print("SEESAW BOUND ANALYSIS")
print("=" * 70)

# At n=6, check if there's a formula for Ω₂ and Ω₃
n = 6; E = n*(n-1)//2
print(f"\nn={n}: E={E}, C(n,2)-n={E-n}")

t0 = time.time()
data = []

for bits in range(1 << E):
    A = tourn(n, bits)
    o2 = omega_dim(A, n, 2)
    o3 = omega_dim(A, n, 3)
    c3 = count_3cycles(A, n)
    data.append((bits, o2, o3, c3))

    if (bits+1) % 10000 == 0:
        print(f"  {bits+1}/{1<<E}, t={time.time()-t0:.0f}s", flush=True)

print(f"  Done in {time.time()-t0:.1f}s")

# Analyze: is Ω₂ a function of c₃?
print("\n--- Is Ω₂ determined by c₃? ---")
o2_by_c3 = {}
for bits, o2, o3, c3 in data:
    if c3 not in o2_by_c3: o2_by_c3[c3] = Counter()
    o2_by_c3[c3][o2] += 1

for c3 in sorted(o2_by_c3.keys()):
    vals = o2_by_c3[c3]
    if len(vals) == 1:
        v = list(vals.keys())[0]
        print(f"  c₃={c3}: Ω₂={v} always ({sum(vals.values())} tournaments)")
    else:
        print(f"  c₃={c3}: Ω₂={dict(sorted(vals.items()))} (VARIES)")

# Check: Ω₂ = E + c₃ ?  or Ω₂ = some formula?
print("\n--- Ω₂ formula check ---")
for bits, o2, o3, c3 in data[:20]:
    print(f"  bits={bits}: Ω₂={o2}, c₃={c3}, Ω₂-E={o2-E}, Ω₂-c₃={o2-c3}")

# Check O2 = E + something?
o2_minus_E = Counter()
for bits, o2, o3, c3 in data:
    o2_minus_E[(o2 - E, c3)] += 1

# Print unique (Ω₂-E) values
print(f"\n  Ω₂ - E values: {sorted(set(o2-E for _, o2, _, _ in data))}")
print(f"  Ω₃ values: {sorted(set(o3 for _, _, o3, _ in data))}")

# Check: is Ω₂ - Ω₃ + something universal?
gap_vals = Counter()
for bits, o2, o3, c3 in data:
    gap_vals[(o2-o3, c3)] += 1

print(f"\n  (Ω₂-Ω₃) values: {sorted(set(o2-o3 for _, o2, o3, _ in data))}")

# Check Ω₂ = f(c₃) directly
print("\n--- Direct formula search ---")
for bits, o2, o3, c3 in data[:30]:
    # Ω₁ = E (all edges are allowed 1-paths in a tournament)
    # Ω₂ = ? Let's check Ω₂ = E + c₃ - something
    print(f"  c₃={c3}: Ω₂={o2}, Ω₃={o3}, Ω₂-c₃={o2-c3}, 2*c₃-Ω₂={2*c3-o2}")

# Check: does Ω₂ = c₃ + E - n?
print("\n--- Ω₂ = c₃ + E - n? ---")
all_match = all(o2 == c3 + E - n for _, o2, _, c3 in data)
print(f"  Ω₂ = c₃ + E - n = c₃ + {E-n}? {all_match}")

# Try: Ω₂ = c₃ + constant?
first_diff = data[0][1] - data[0][3]
all_match2 = all(o2 - c3 == first_diff for _, o2, _, c3 in data)
print(f"  Ω₂ - c₃ = constant? diff₀={first_diff}, all_match={all_match2}")

# Maybe Ω₂ depends on more than c₃. Let's check if it's determined by the score sequence.
print("\n--- Ω₂ by score sequence ---")
o2_by_score = {}
for bits, o2, o3, c3 in data:
    A = tourn(n, bits)
    sc = tuple(sorted([sum(A[i]) for i in range(n)]))
    if sc not in o2_by_score: o2_by_score[sc] = Counter()
    o2_by_score[sc][o2] += 1

for sc in sorted(o2_by_score.keys()):
    vals = o2_by_score[sc]
    if len(vals) > 1:
        print(f"  score={list(sc)}: Ω₂={dict(sorted(vals.items()))}")

print("\nDone.")
