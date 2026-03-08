#!/usr/bin/env python3
"""
flip_reversal_proof.py — Algebraic proof of the arc-flip formula and
computational verification of the reversal involution.

THEOREM (Arc-Flip Factorization):
  F(T,x) - F(T',x) = (x-1) * D(x)
  where D(x) = G_uv(x) - G_vu(x) and T' is T with arc u->v flipped to v->u.

THEOREM (Anti-Palindromicity of D):
  D_k = -D_{n-2-k}  (i.e., D is anti-palindromic of degree n-2).

PROOF SKETCH:
  The reversal involution P -> P^rev sends fwd(P) to (n-1)-fwd(P).
  If P contains ...u,v... at some position, then P^rev contains ...v,u...
  So reversal maps type-(a) paths (containing ...u,v...) with fwd=f
  to type-(b) paths (containing ...v,u...) with fwd=(n-1-f).

  In shifted coordinates:
    G_uv[k] counts type-(a) paths with fwd = k+1
    These map to type-(b) paths with fwd = (n-1)-(k+1) = n-2-k
    So G_uv[k] = G_vu[n-2-k]

  Therefore D[k] = G_uv[k] - G_vu[k] = G_vu[n-2-k] - G_vu[k] = -D[n-2-k].

This script verifies:
  1. G_uv[k] = G_vu[n-2-k] for all tournaments at n=4,5
  2. D is anti-palindromic
  3. F(T)-F(T') = (x-1)*D(x)

Author: opus-2026-03-07
"""
from itertools import permutations

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

def compute_G_uv_vu(adj, n, u, v):
    """
    G_uv[k] = # perms P containing ...u,v... consecutively with fwd(P) = k+1
               (shifted down by 1 because u->v is a forward step).
    G_vu[k] = # perms P containing ...v,u... consecutively with fwd(P) = k.
    """
    G_uv = [0]*(n-1)
    G_vu = [0]*(n-1)
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        for i in range(n-1):
            if P[i] == u and P[i+1] == v:
                G_uv[fwd - 1] += 1
                break
            elif P[i] == v and P[i+1] == u:
                G_vu[fwd] += 1
                break
    return G_uv, G_vu


# ============================================================
# TEST 1: G_uv[k] = G_vu[n-2-k] (reversal involution identity)
# ============================================================
print("=" * 70)
print("TEST 1: Reversal involution identity G_uv[k] = G_vu[n-2-k]")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 0
    passes = 0

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue
                total += 1

                G_uv, G_vu = compute_G_uv_vu(adj, n, u, v)

                # Check: G_uv[k] = G_vu[n-2-k] for all k
                d = n - 2  # degree of D
                ok = all(G_uv[k] == G_vu[d - k] for k in range(d + 1))
                if ok:
                    passes += 1
                else:
                    print(f"  FAIL: n={n} bits={bits} arc {u}->{v}")
                    print(f"    G_uv = {G_uv}")
                    print(f"    G_vu = {G_vu}")
                    print(f"    G_vu reversed = {G_vu[::-1]}")

    print(f"  n={n}: {passes}/{total} pass")


# ============================================================
# TEST 2: D is anti-palindromic: D[k] = -D[n-2-k]
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Anti-palindromicity D[k] = -D[n-2-k]")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 0
    passes = 0

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue
                total += 1

                G_uv, G_vu = compute_G_uv_vu(adj, n, u, v)
                D = [G_uv[k] - G_vu[k] for k in range(n-1)]
                d = n - 2

                anti_pal = all(D[k] == -D[d - k] for k in range(d + 1))
                if anti_pal:
                    passes += 1
                else:
                    print(f"  FAIL: n={n} bits={bits} arc {u}->{v}")
                    print(f"    D = {D}")

    print(f"  n={n}: {passes}/{total} pass")


# ============================================================
# TEST 3: F(T) - F(T') = (x-1) * D(x)
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: F(T,x) - F(T',x) = (x-1) * D(x)")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 0
    passes = 0

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F_T = compute_F(adj, n)

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue
                total += 1

                # Flip arc u->v to v->u
                adj_flip = [row[:] for row in adj]
                adj_flip[u][v] = 0
                adj_flip[v][u] = 1
                F_flip = compute_F(adj_flip, n)

                diff = [F_T[k] - F_flip[k] for k in range(n)]

                G_uv, G_vu = compute_G_uv_vu(adj, n, u, v)
                D = [G_uv[k] - G_vu[k] for k in range(n-1)]

                # (x-1) * D(x): coefficient of x^k is D[k-1] - D[k]
                prod = [0]*n
                for k in range(n-1):
                    prod[k] -= D[k]      # -1 * x^k term
                    prod[k+1] += D[k]    # x * x^k term

                if diff == prod:
                    passes += 1
                else:
                    print(f"  FAIL: n={n} bits={bits} arc {u}->{v}")
                    print(f"    diff = {diff}")
                    print(f"    (x-1)*D = {prod}")

    print(f"  n={n}: {passes}/{total} pass")


# ============================================================
# TEST 4: Explicit reversal involution demonstration
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Explicit reversal involution (one example)")
print("=" * 70)

n = 4
adj = tournament_from_bits(n, 0)  # transitive tournament
u, v = 0, 1
if not adj[u][v]:
    # find an arc
    for u in range(n):
        for v in range(n):
            if adj[u][v]:
                break
        else:
            continue
        break

print(f"n={n}, arc {u}->{v}")
print(f"Adjacency matrix:")
for i in range(n):
    print(f"  {adj[i]}")

# Show the bijection explicitly
type_a = []  # paths containing ...u,v...
type_b = []  # paths containing ...v,u...

for P in permutations(range(n)):
    fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
    for i in range(n-1):
        if P[i] == u and P[i+1] == v:
            type_a.append((list(P), fwd))
            break
        elif P[i] == v and P[i+1] == u:
            type_b.append((list(P), fwd))
            break

print(f"\nType (a) paths (containing ...{u},{v}...):")
for P, f in sorted(type_a, key=lambda x: x[1]):
    P_rev = list(reversed(P))
    fwd_rev = sum(1 for i in range(n-1) if adj[P_rev[i]][P_rev[i+1]])
    print(f"  P={P}, fwd={f}, P^rev={P_rev}, fwd(P^rev)={fwd_rev}")

print(f"\nType (b) paths (containing ...{v},{u}...):")
for P, f in sorted(type_b, key=lambda x: x[1]):
    P_rev = list(reversed(P))
    fwd_rev = sum(1 for i in range(n-1) if adj[P_rev[i]][P_rev[i+1]])
    print(f"  P={P}, fwd={f}, P^rev={P_rev}, fwd(P^rev)={fwd_rev}")

G_uv, G_vu = compute_G_uv_vu(adj, n, u, v)
D = [G_uv[k] - G_vu[k] for k in range(n-1)]
print(f"\nG_uv = {G_uv}")
print(f"G_vu = {G_vu}")
print(f"G_vu reversed = {G_vu[::-1]}")
print(f"G_uv == G_vu reversed: {G_uv == G_vu[::-1]}")
print(f"D = {D}")
print(f"D is anti-palindromic: {all(D[k] == -D[n-2-k] for k in range(n-1))}")


# ============================================================
# COROLLARY: H(T) is invariant under single arc flip
# ============================================================
print("\n" + "=" * 70)
print("COROLLARY: H(T) = H(T') under single arc flip")
print("  (Since F(T,1) - F(T',1) = (1-1)*D(1) = 0)")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 0
    passes = 0

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        H_T = sum(compute_F(adj, n))

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue
                total += 1

                adj_flip = [row[:] for row in adj]
                adj_flip[u][v] = 0
                adj_flip[v][u] = 1
                H_flip = sum(compute_F(adj_flip, n))

                if H_T == H_flip:
                    passes += 1
                else:
                    print(f"  COUNTEREXAMPLE: n={n} bits={bits} arc {u}->{v}: H_T={H_T}, H_flip={H_flip}")

    print(f"  n={n}: {passes}/{total} have H(T)=H(T')")
    if passes < total:
        print(f"  *** H is NOT invariant under arc flip! The formula implies the")
        print(f"      POLYNOMIAL difference factors as (x-1)*D(x), not that H is invariant.")
        print(f"      We need: F(T,1)-F(T',1) = 0 iff (x-1)*D(x) evaluated at x=1 = 0.")
        print(f"      Since (x-1)*D(x) is a polynomial, evaluating at x=1 gives 0 trivially.")
        print(f"      So H(T) = H(T'). But let's double-check the counterexamples...")
