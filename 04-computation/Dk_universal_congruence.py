#!/usr/bin/env python3
"""
Investigate: which D_k mod 2^{n-1-k} are universal (tournament-independent)?

D_k = sum_P C(forward(P), k) where forward(P) = #{i : A[P_i][P_{i+1}]=1}

Equivalently: D_k = sum_{k-subsets S of positions} sum_P prod_{i in S} A[P_i][P_{i+1}]

For D_2: sum over pairs of positions (i,j) of sum_P A[P_i][P_{i+1}]*A[P_j][P_{j+1}]

For non-adjacent positions i < j with j > i+1:
  sum_P A[P_i][P_{i+1}]*A[P_j][P_{j+1}] 
  = sum over 4-tuples (u,v,w,x) of distinct vertices of A[u][v]*A[w][x] * (factor)

For adjacent positions j = i+1:
  sum_P A[P_i][P_{i+1}]*A[P_{i+1}][P_{i+2}]
  = sum over 3-tuples (u,v,w) of A[u][v]*A[v][w] * (n-3)! * (n-3 slots for remaining)
  Wait, need to be more careful.

Let me just compute D_2 for many tournaments and check the mod.
"""
from itertools import permutations, combinations
from collections import Counter
import random
import math

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

def compute_Dk_fast(A, n):
    """Compute D_k via forward edge count distribution."""
    D = [0] * n
    binom_cache = {}
    def binom(n, k):
        if k < 0 or k > n: return 0
        if (n,k) not in binom_cache:
            binom_cache[(n,k)] = math.comb(n, k)
        return binom_cache[(n,k)]
    
    for P in permutations(range(n)):
        fwd = sum(A[P[i]][P[i+1]] for i in range(n-1))
        for k in range(n):
            D[k] += binom(fwd, k)
    return D

# Check D_2 mod 2^{n-3} for n = 3, 5, 7
for n in [3, 5, 7]:
    print(f"\n=== n={n}: D_k congruences ===")
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    m = len(tiles)
    
    num_tests = min(100, 1 << m)
    Dk_mod_data = [Counter() for _ in range(n)]
    
    for trial in range(num_tests):
        if m <= 10:
            bits = trial % (1 << m)
        else:
            bits = random.getrandbits(m)
        A = tiling_to_tournament(bits, n)
        D = compute_Dk_fast(A, n)
        
        for k in range(n):
            req = max(1, 2**(n-1-k))
            Dk_mod_data[k][D[k] % req] += 1
    
    for k in range(n):
        req = max(1, 2**(n-1-k))
        vals = sorted(Dk_mod_data[k].items())
        is_universal = (len(vals) == 1)
        print(f"  D_{k} mod {req}: {dict(vals)} {'UNIVERSAL' if is_universal else 'varies'}")

# Now try to find D_2 as a formula
print("\n=== D_2 formula exploration ===")
for n in [3, 5, 7]:
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    m = len(tiles)
    
    # D_2 = sum over pairs of positions × perm counts
    # = sum_{i<j, j>i+1} sum_P A[P_i][P_{i+1}]*A[P_j][P_{j+1}]  (non-adjacent)
    # + sum_{i} sum_P A[P_i][P_{i+1}]*A[P_{i+1}][P_{i+2}]          (adjacent)
    
    # Non-adjacent pair (i,j) with j > i+1:
    # sum_P A[P_i][P_{i+1}]*A[P_j][P_{j+1}] = 
    #   sum_{(u,v) u!=v} sum_{(w,x) w!=x, w,x != u or v at most} A[u][v]*A[w][x] * count
    
    # For non-adjacent: P_i=u, P_{i+1}=v, P_j=w, P_{j+1}=x with u,v,w,x all distinct
    # Count of such perms: (n-4)! * P(n-4, n-4-remainder) = (n-4)! * 1 for non-overlapping
    # Wait, positions i,i+1,j,j+1 are 4 distinct positions. The 4 values u,v,w,x at these
    # positions must be distinct. The remaining n-4 values go to remaining n-4 positions.
    # Count = (n-4)! for each valid (u,v,w,x) quadruple.
    
    # Number of non-adjacent pairs (i,j): C(n-1,2) - (n-2) = C(n-1,2) - (n-2)
    n_nonadj = math.comb(n-1, 2) - (n-2)
    
    # For non-adjacent: sum = n_nonadj * (n-4)! * sum_{u≠v, w≠x, all distinct} A[u][v]*A[w][x]
    # sum_{u≠v, w≠x, all distinct} A[u][v]*A[w][x] = E*(E-1) - ... hmm
    # Actually it's sum_{u≠v} A[u][v] * (sum_{w≠x, w,x not in {u,v}} A[w][x])
    # = sum_{u≠v} A[u][v] * (E - deg_out(u) - deg_in(v) + A[u][v]... no this gets messy)
    # where E = C(n,2) edges total
    
    # Adjacent pair at position i: P_i=u, P_{i+1}=v, P_{i+2}=w with u,v,w distinct
    # sum = (n-2) * (n-3)! * sum_{u≠v≠w, distinct} A[u][v]*A[v][w]
    n_adj = n - 2
    
    # sum_{u≠v≠w distinct} A[u][v]*A[v][w] = sum_v (out_deg(v)) * (in_deg(v))
    # Wait: sum_{u: u→v} sum_{w: v→w} = in_deg(v) * out_deg(v)
    # Hmm no: A[u][v] means u→v, so sum_{u≠v} A[u][v] = in_deg(v).
    # A[v][w] means v→w, so sum_{w≠v} A[v][w] = out_deg(v).
    # But need u,v,w all distinct: u ≠ w. So subtract cases where u = w.
    # sum_{u≠v, w≠v, u≠w} A[u][v]*A[v][w] = sum_v (in(v) * out(v) - #{u : u→v and v→u})
    # But in a tournament, u→v implies not v→u. So A[u][v]*A[v][u] = 0 always.
    # Wait, A[u][v]=1 and A[v][u]=0 in a tournament. So A[u][v]*A[v][u]=0.
    # So u=w gives A[u][v]*A[v][u] = 0 for all u,v.
    # Therefore sum = sum_v in_deg(v) * out_deg(v) = sum_v d_in(v) * d_out(v).
    # Since d_in(v) + d_out(v) = n-1, this is sum_v d_in(v) * (n-1-d_in(v)).
    
    print(f"\nn={n}: D_2 decomposition")
    print(f"  Non-adjacent pairs: {n_nonadj}, adjacent pairs: {n_adj}")
    
    # For a specific tournament, compute
    bits_test = 0
    A = tiling_to_tournament(bits_test, n)
    
    # Verify
    D = compute_Dk_fast(A, n)
    
    # Direct D_2 computation via formula
    E = sum(A[i][j] for i in range(n) for j in range(n) if i != j)  # = C(n,2)
    
    # Score sequence
    scores = [sum(A[i][j] for j in range(n)) for i in range(n)]
    sum_din_dout = sum(s * (n-1-s) for s in scores)
    
    # Non-adjacent contribution
    # sum_{u≠v,w≠x, all distinct} A[u][v]*A[w][x] = 
    # (sum_{u≠v} A[u][v])^2 - sum_{u≠v} A[u][v]^2 - 2*sum_{u≠v,w: w∈{u,v}} A[u][v]*A[w][x']
    # This is getting complicated. Let me just verify numerically.
    
    # Direct computation
    adj_sum = 0
    for i_pos in range(n-2):
        for P in permutations(range(n)):
            adj_sum += A[P[i_pos]][P[i_pos+1]] * A[P[i_pos+1]][P[i_pos+2]]
    
    nonadj_sum = D[2] - adj_sum
    
    print(f"  D_2 = {D[2]}, adj_contribution = {adj_sum}, nonadj_contribution = {nonadj_sum}")
    print(f"  adj/(n-2)/(n-3)! = {adj_sum / (n-2) / math.factorial(n-3):.1f} = sum(din*dout)?")
    print(f"  sum(din*dout) = {sum_din_dout}")
    
    adj_per_pos = adj_sum // (n-2)
    print(f"  adj per position = {adj_per_pos} = {adj_per_pos / math.factorial(n-3):.1f} * (n-3)!")

