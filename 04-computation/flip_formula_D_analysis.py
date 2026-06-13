#!/usr/bin/env python3
"""
flip_formula_D_analysis.py — Analyzing D(x) = G_uv(x) - G_vu(x)

PROVED (n=4,5):
  F(T,x) - F(T',x) = (x-1) * D(x)
  D(x) is anti-palindromic: D_k = -D_{n-2-k}
  G_uv + G_vu = 2*F(T/e,x)
  D(1) = 0

What determines D(x)?
- D measures the "asymmetry" of the arc e within T.
- For the transitive tournament (all arcs aligned), D is always 0
  (by symmetry, or because G_uv = G_vu = F(T/e)).

QUESTION: Is D determined by local structure (score of u,v, common neighbors)?

Author: opus-2026-03-07-S45
"""
from itertools import permutations, combinations
import math

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

def compute_D(adj, n, u, v):
    """Compute D(x) = G_uv(x) - G_vu(x) for arc u->v."""
    G_uv = [0]*(n-1)
    G_vu = [0]*(n-1)
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        for i in range(n-1):
            if P[i]==u and P[i+1]==v:
                G_uv[fwd-1] += 1
                break
            elif P[i]==v and P[i+1]==u:
                if fwd < n-1:
                    G_vu[fwd] += 1
                break
    return [G_uv[k] - G_vu[k] for k in range(n-1)]

# ============================================================
# D(x) FOR ALL ARCS — WHAT DETERMINES IT?
# ============================================================
print("=" * 60)
print("D(x) = G_uv(x) - G_vu(x) — structure analysis")
print("=" * 60)

n = 5
m = n*(n-1)//2

# Collect (D, local_info) for all arcs in all tournaments
data = []
seen_F = set()

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)
    key_F = tuple(F_T)

    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue

            D = compute_D(adj, n, u, v)

            # Local info about arc (u,v)
            score_u = sum(adj[u])
            score_v = sum(adj[v])
            # How many others x have u->x->v (common "out of u, in to v")?
            others = [x for x in range(n) if x != u and x != v]
            common_uv = sum(1 for x in others if adj[u][x] and adj[x][v])
            common_vu = sum(1 for x in others if adj[v][x] and adj[x][u])

            data.append({
                'D': tuple(D), 'bits': bits, 'u': u, 'v': v,
                'score_u': score_u, 'score_v': score_v,
                'common_uv': common_uv, 'common_vu': common_vu,
                'F': key_F
            })

# Group by D values
from collections import Counter, defaultdict
D_counter = Counter(d['D'] for d in data)
print(f"\nn=5: {len(data)} total (arc, tournament) pairs")
print(f"Distinct D values: {len(D_counter)}")
for D_val, count in D_counter.most_common(20):
    print(f"  D={list(D_val)}: {count} times")

# Does (score_u, score_v, common_uv, common_vu) determine D?
sig_to_D = defaultdict(set)
for d in data:
    sig = (d['score_u'], d['score_v'], d['common_uv'], d['common_vu'])
    sig_to_D[sig].add(d['D'])

print(f"\nDoes (score_u, score_v, common_uv, common_vu) determine D?")
multi = 0
for sig, Ds in sorted(sig_to_D.items()):
    if len(Ds) > 1:
        multi += 1
        if multi <= 5:
            print(f"  sig={sig}: {len(Ds)} distinct D values")
print(f"  Signatures with multiple D: {multi}/{len(sig_to_D)}")

# Does D depend on the FULL tournament or just on the arc type?
# D = 0 when? (symmetric arcs)
zero_D = [d for d in data if all(x == 0 for x in d['D'])]
print(f"\nD = 0 (symmetric arc): {len(zero_D)}/{len(data)} = {100*len(zero_D)/len(data):.1f}%")
# Characterize arcs with D=0
if zero_D:
    zero_scores = Counter((d['score_u'], d['score_v']) for d in zero_D)
    print(f"  Score pairs for D=0: {dict(zero_scores)}")

# ============================================================
# D(x) AND THE TRANSITIVE TOURNAMENT
# ============================================================
print("\n" + "=" * 60)
print("D(x) FOR TRANSITIVE TOURNAMENT (all arcs i->j for i>j)")
print("=" * 60)

for n in [4, 5]:
    adj = tournament_from_bits(n, 0)  # transitive: lower index beats higher
    print(f"\nn={n}, transitive:")
    for u in range(n):
        for v in range(n):
            if adj[u][v]:
                D = compute_D(adj, n, u, v)
                print(f"  arc {u}->{v}: D={D}")

# ============================================================
# D(x) AND SCORE SEQUENCE
# ============================================================
print("\n" + "=" * 60)
print("D(x) DETERMINED BY SCORE DIFFERENCE?")
print("=" * 60)

n = 5
# For each arc, score_diff = score_u - score_v
diff_to_D = defaultdict(set)
for d in data:
    sdiff = d['score_u'] - d['score_v']
    diff_to_D[sdiff].add(d['D'])

for sdiff in sorted(diff_to_D.keys()):
    Ds = diff_to_D[sdiff]
    print(f"  score_diff={sdiff:+d}: {len(Ds)} distinct D values")
    if len(Ds) <= 3:
        for D in sorted(Ds):
            print(f"    D={list(D)}")

# ============================================================
# D(x) COEFFICIENTS — LINEAR IN SCORE DIFFERENCE?
# ============================================================
print("\n" + "=" * 60)
print("D COEFFICIENTS vs SCORE DIFFERENCE")
print("=" * 60)

# Since D is anti-palindromic of degree n-2=3, D = [a, b, -b, -a]
# Only need a and b. Let's see if they depend on score difference.

n = 5
score_diff_to_ab = defaultdict(set)
for d in data:
    sdiff = d['score_u'] - d['score_v']
    D = d['D']
    a, b = D[0], D[1]
    score_diff_to_ab[sdiff].add((a, b))

for sdiff in sorted(score_diff_to_ab.keys()):
    abs_values = score_diff_to_ab[sdiff]
    print(f"  score_diff={sdiff:+d}: {len(abs_values)} distinct (a,b) pairs: {sorted(abs_values)[:10]}")

# ============================================================
# FLIP FORMULA: F(T) - F(T') and the "descent change"
# ============================================================
print("\n" + "=" * 60)
print("FLIP FORMULA IMPLICATIONS")
print("=" * 60)

# F(T,x) - F(T',x) = (x-1) * D(x)
# At x=1: H(T) - H(T') = 0 * D(1) = 0. So H(T) = H(T')?
# No! H can change under flip. But (x-1)*D(x) at x=1 requires D(1)=0
# for the formula to give 0. Let's check.

n = 5
m = n*(n-1)//2

for bits in [0, 10, 100, 500]:
    if bits >= (1 << m):
        continue
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)
    H_T = sum(F_T)

    for u in range(n):
        for v in range(n):
            if adj[u][v]:
                break
        else:
            continue
        break

    flip_adj = [row[:] for row in adj]
    flip_adj[u][v] = 0
    flip_adj[v][u] = 1
    F_flip = compute_F(flip_adj, n)
    H_flip = sum(F_flip)

    D = compute_D(adj, n, u, v)
    D_at_1 = sum(D)

    diff_at_1 = H_T - H_flip
    xm1_D_at_1 = 0  # (x-1)*D at x=1 = 0*D(1) -- but really it's a limit

    # Actually: F(T,1) - F(T',1) = (1-1)*D(1) = 0 only if D has no pole.
    # But (x-1)*D(x) at x=1 gives: lim = 0 * D(1) = 0 if D(1) finite.
    # We showed D(1) = 0 (sum of anti-palindromic). So (x-1)*D has zero at x=1.
    # More precisely: (x-1)*D(x) at x=1 = 0. And H(T)-H(T') = sum(F(T)-F(T')).
    # So H(T) = H(T') ALWAYS under a single arc flip? That can't be right...

    print(f"  bits={bits} arc {u}->{v}: H_T={H_T}, H_flip={H_flip}, diff={H_T-H_flip}")
    print(f"    D={D}, D(1)={D_at_1}")

# Wait, (x-1)*D(x) evaluated at x=1:
# sum_k [(x-1)*D(x)]_k * 1^k = sum_k coefficient_k = value at x=1
# coefficient_k of (x-1)*D(x) = D[k-1] - D[k]
# sum = D[-1] - D[n-2] + D[0] - D[1] + ... = ?
# Actually, sum of (x-1)*D(x) coefficients = D(1) * (1-1) + ... no.
# Polynomial product: if f(x) = (x-1)*D(x), then f(1) = (1-1)*D(1) = 0.
# So f(1) = sum of f's coefficients = 0.
# And f(1) = F_T(1) - F_flip(1) = H_T - H_flip.
# So H_T = H_flip for every arc flip!!

# THIS IS A CONSEQUENCE: single arc flip preserves H(T).
# But wait, is this true? Let me check...
