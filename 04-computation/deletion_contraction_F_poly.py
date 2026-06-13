#!/usr/bin/env python3
"""
deletion_contraction_F_poly.py — Deletion-contraction for F(T,x).

Kind-pasteur verified H(T) = H(T\e) + H(T/e) at n=4,5 (100%).
This script tests whether the identity lifts to the F-polynomial level.

Convention (from kind-pasteur): for arc e = u->v:
  T\e = delete arc u->v (NOT flip; result is a digraph, not tournament)
  T/e = contract u,v into w; w gets IN from u (tail), OUT from v (head)

GOAL: Find a formula relating F(T,x), F(T\e,x), F(T/e,x).

Author: opus-2026-03-07-S45
"""
from itertools import permutations
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

def ham_paths_by_fwd(adj, n):
    """Count Ham paths by number of forward edges. Returns dict {fwd_count: path_count}."""
    full = (1 << n) - 1
    # dp[mask][v] = dict: fwd_count -> count
    dp = [[None]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = {0: 1}
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] is None:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    new_mask = mask | (1 << u)
                    if dp[new_mask][u] is None:
                        dp[new_mask][u] = {}
                    for fc, cnt in dp[mask][v].items():
                        dp[new_mask][u][fc+1] = dp[new_mask][u].get(fc+1, 0) + cnt
                elif adj[u][v]:
                    # backward edge u<-v, still a valid step v->u in the path
                    new_mask = mask | (1 << u)
                    if dp[new_mask][u] is None:
                        dp[new_mask][u] = {}
                    for fc, cnt in dp[mask][v].items():
                        dp[new_mask][u][fc] = dp[new_mask][u].get(fc, 0) + cnt
                # else: no edge between v and u (happens in digraphs)
    F = {}
    for v in range(n):
        if dp[full][v] is not None:
            for fc, cnt in dp[full][v].items():
                F[fc] = F.get(fc, 0) + cnt
    return F

def ham_paths_by_fwd_digraph(adj, n):
    """For a general digraph: count Ham paths by fwd edges.
    A step from v to u is valid iff adj[v][u] OR adj[u][v].
    fwd = number of steps where adj[v][u]=1.
    """
    full = (1 << n) - 1
    dp = [[None]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = {0: 1}
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] is None:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u] or adj[u][v]:
                    new_mask = mask | (1 << u)
                    if dp[new_mask][u] is None:
                        dp[new_mask][u] = {}
                    fwd_inc = 1 if adj[v][u] else 0
                    for fc, cnt in dp[mask][v].items():
                        nf = fc + fwd_inc
                        dp[new_mask][u][nf] = dp[new_mask][u].get(nf, 0) + cnt
    F = {}
    for v in range(n):
        if dp[full][v] is not None:
            for fc, cnt in dp[full][v].items():
                F[fc] = F.get(fc, 0) + cnt
    return F

def delete_arc(adj, n, u, v):
    """Delete arc u->v. Result has no edge between u and v."""
    B = [row[:] for row in adj]
    B[u][v] = 0
    return B

def contract_arc(adj, n, u, v):
    """Contract arc u->v. w gets IN from u, OUT from v."""
    others = [x for x in range(n) if x != u and x != v]
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    # vertex 0 = w, vertices 1.. = others
    for i, x in enumerate(others):
        if adj[x][u]:
            B[i+1][0] = 1
        if adj[v][x]:
            B[0][i+1] = 1
    for i, x in enumerate(others):
        for j, y in enumerate(others):
            B[i+1][j+1] = adj[x][y]
    return B, new_n

def dict_to_list(F_dict, degree):
    """Convert {fwd: count} dict to list [F_0, F_1, ..., F_{degree}]."""
    return [F_dict.get(k, 0) for k in range(degree + 1)]

def poly_str(F_list):
    """Pretty-print polynomial."""
    return str(F_list)

# ============================================================
# VERIFY H(T) = H(T\e) + H(T/e) first
# ============================================================
print("=" * 60)
print("VERIFY: H(T) = H(T\\e) + H(T/e) at n=4")
print("=" * 60)

n = 4
m = n*(n-1)//2
total = 0
ok = 0

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F_T = ham_paths_by_fwd(adj, n)
    H_T = sum(F_T.values())

    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue
            del_adj = delete_arc(adj, n, u, v)
            F_del = ham_paths_by_fwd_digraph(del_adj, n)
            H_del = sum(F_del.values())

            con_adj, con_n = contract_arc(adj, n, u, v)
            F_con = ham_paths_by_fwd(con_adj, con_n)
            H_con = sum(F_con.values())

            total += 1
            if H_T == H_del + H_con:
                ok += 1
            elif total - ok <= 3:
                print(f"  FAIL: bits={bits} arc {u}->{v}: H_T={H_T}, H_del={H_del}, H_con={H_con}")

print(f"  Result: {ok}/{total} pass")

# ============================================================
# DECOMPOSE F(T,x) by arc usage
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) DECOMPOSITION BY ARC USAGE")
print("=" * 60)

n = 4
m = n*(n-1)//2

for bits in [0, 3, 7, 15, 31, 63]:
    if bits >= (1 << m):
        continue
    adj = tournament_from_bits(n, bits)
    F_T = dict_to_list(ham_paths_by_fwd(adj, n), n-1)

    # Pick first arc
    u, v = None, None
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                u, v = i, j
                break
        if u is not None:
            break

    print(f"\nbits={bits}, arc {u}->{v}:")
    print(f"  F(T) = {F_T}")

    # Decompose: for each Hamiltonian path, classify by whether u,v consecutive
    F_uv_fwd = {}  # paths with ...u,v... (arc goes forward)
    F_vu_bwd = {}  # paths with ...v,u... (arc goes backward)
    F_none = {}    # paths where u,v not consecutive

    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
        uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
        if uses_uv:
            F_uv_fwd[fwd] = F_uv_fwd.get(fwd, 0) + 1
        elif uses_vu:
            F_vu_bwd[fwd] = F_vu_bwd.get(fwd, 0) + 1
        else:
            F_none[fwd] = F_none.get(fwd, 0) + 1

    F_uv_list = dict_to_list(F_uv_fwd, n-1)
    F_vu_list = dict_to_list(F_vu_bwd, n-1)
    F_none_list = dict_to_list(F_none, n-1)

    print(f"  F_uv (uses u->v, fwd) = {F_uv_list}")
    print(f"  F_vu (uses v->u, bwd) = {F_vu_list}")
    print(f"  F_none (no u,v)       = {F_none_list}")
    print(f"  Check: sum = {[F_uv_list[k]+F_vu_list[k]+F_none_list[k] for k in range(n)]}")

    # T\e: no edge between u,v. Only F_none paths survive.
    del_adj = delete_arc(adj, n, u, v)
    F_del = dict_to_list(ham_paths_by_fwd_digraph(del_adj, n), n-1)
    print(f"  F(T\\e)   = {F_del}")
    print(f"  == F_none? {F_del == F_none_list}")

    # T/e: contract. The contracted paths come from F_uv and F_vu.
    # For F_uv paths: ...u,v... contributes +1 to fwd. After contraction, this step
    # disappears, so the contracted path has fwd_count - 1 forward edges.
    # For F_vu paths: ...v,u... contributes 0 to fwd. After contraction, this step
    # disappears, so the contracted path has the same fwd_count.
    #
    # BUT: the contraction changes the vertex set, so the adjacency may differ.
    # Specifically, w's edges are NOT the same as u's or v's in general.
    # So the fwd count in the contracted tournament may differ from what we'd expect.

    con_adj, con_n = contract_arc(adj, n, u, v)
    F_con = dict_to_list(ham_paths_by_fwd(con_adj, con_n), con_n-1)
    print(f"  F(T/e)   = {F_con}")

    # Naive expectation: F_con[k] = F_uv[k+1] + F_vu[k]
    naive = [(F_uv_fwd.get(k+1, 0) + F_vu_bwd.get(k, 0)) for k in range(con_n)]
    print(f"  F_uv[k+1]+F_vu[k] = {naive}")
    print(f"  Naive match? {naive == F_con}")

    # If naive doesn't match, the issue is that contraction changes which edges
    # are "forward" for the OTHER steps in the path.

# ============================================================
# TEST: F(T,x) = F(T\e,x) + x * G(T/e, x) for some G?
# ============================================================
print("\n" + "=" * 60)
print("SEARCHING FOR F(T,x) DELETION-CONTRACTION")
print("=" * 60)

n = 4
m = n*(n-1)//2

# For each tournament and each arc, compute F(T), F(T\e), F(T/e)
# and look for polynomial relationships.

# Collect all (F_T, F_del, F_con) triples
triples = []
for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F_T = dict_to_list(ham_paths_by_fwd(adj, n), n-1)

    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue
            del_adj = delete_arc(adj, n, u, v)
            F_del = dict_to_list(ham_paths_by_fwd_digraph(del_adj, n), n-1)
            con_adj, con_n = contract_arc(adj, n, u, v)
            F_con = dict_to_list(ham_paths_by_fwd(con_adj, con_n), con_n-1)
            triples.append((F_T, F_del, F_con, bits, u, v))

# Test: F_T[k] = F_del[k] + F_con[k-1]  (shift by x)
test_shift = sum(1 for (ft, fd, fc, *_) in triples
                 if all(ft[k] == fd[k] + (fc[k-1] if k > 0 else 0) for k in range(n)))
print(f"  F_T[k] = F_del[k] + F_con[k-1]: {test_shift}/{len(triples)}")

# Test: F_T[k] = F_del[k] + F_con[k]  (no shift)
test_noshift = sum(1 for (ft, fd, fc, *_) in triples
                   if all(ft[k] == fd[k] + (fc[k] if k < len(fc) else 0) for k in range(n)))
print(f"  F_T[k] = F_del[k] + F_con[k]: {test_noshift}/{len(triples)}")

# Test: F_T[k] = F_del[k] + F_con[k-1] + F_con[k]  (x+1 multiplier)
test_xp1 = sum(1 for (ft, fd, fc, *_) in triples
               if all(ft[k] == fd[k] + (fc[k-1] if k > 0 else 0) + (fc[k] if k < len(fc) else 0)
                      for k in range(n)))
print(f"  F_T[k] = F_del[k] + (x+1)*F_con[k]: {test_xp1}/{len(triples)}")

# Show some examples of the actual relationship
print("\nSample relationships:")
for i in range(min(5, len(triples))):
    ft, fd, fc, bits, u, v = triples[i]
    diff = [ft[k] - fd[k] for k in range(n)]
    print(f"  bits={bits} arc {u}->{v}: F_T-F_del = {diff}, F_con = {fc}")

# Hmm. Let's look at F(T,x) - F(T\e,x) more carefully.
# F(T) - F(T\e) = F_uv + F_vu (the paths using edge (u,v) in either direction)
# In F(T\e), those paths don't exist.
# So F(T) - F(T\e) counts paths using (u,v) consecutively.

# The question is: how does this relate to F(T/e)?
# In T/e, a Ham path on n-1 vertices corresponds to a Ham path in T
# where u,v are consecutive AND the step from u to v (or v to u) is "contracted".

# But the key subtlety: in T/e, the "forward edges" of the contracted path
# are computed using T/e's adjacency, which may not match T's adjacency
# for the non-contracted edges.

# Actually wait: for edges among "others" (neither u nor v), the adjacency
# is UNCHANGED by contraction. The only changes are for edges involving w.
# And w's edges come from u (for incoming) and v (for outgoing).

# So for a path in T/e that visits w at position j:
# - The edge from predecessor to w uses adj[predecessor][u] (since IN from u)
# - The edge from w to successor uses adj[v][successor] (since OUT from v)
# These are EXACTLY the edges that would be used in T for the subsequence
# ...predecessor, u, v, successor... !
# The fwd count in T/e for this path = (fwd count of original T path) - (1 if u->v was forward)
# because the u->v step is removed but all other edges match.

# Wait, that's exactly the naive formula! So why doesn't it match?
# Let me check more carefully...

print("\n" + "=" * 60)
print("DETAILED PATH-BY-PATH COMPARISON")
print("=" * 60)

n = 4
adj = tournament_from_bits(n, 0)  # transitive tournament
u, v = 1, 0  # first arc

print(f"T (bits=0), arc {u}->{v}")
print(f"Adjacency:")
for i in range(n):
    print(f"  {[adj[i][j] for j in range(n)]}")

con_adj, con_n = contract_arc(adj, n, u, v)
print(f"\nT/e adjacency (w=merged {u},{v}):")
for i in range(con_n):
    print(f"  {[con_adj[i][j] for j in range(con_n)]}")

# Others = vertices not u or v
others = [x for x in range(n) if x != u and x != v]
print(f"Others: {others}")
print(f"In T/e: vertex 0 = w (merged {u},{v}), vertices {list(range(1, con_n))} = {others}")

# For each T path using u,v consecutively, show the contracted path and compare fwd counts
print(f"\nPath-by-path:")
for P in permutations(range(n)):
    fwd_T = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
    # Does this path use u,v consecutively?
    uv_pos = None
    for i in range(n-1):
        if (P[i] == u and P[i+1] == v) or (P[i] == v and P[i+1] == u):
            uv_pos = i
            break
    if uv_pos is None:
        continue

    # Contract: remove the u-v step, replace u and v with w
    # The contracted path visits vertices {w} union (others that appear in P)
    # which is all of T/e's vertices.
    contracted = []
    for x in P:
        if x == u or x == v:
            if not contracted or contracted[-1] != 0:
                contracted.append(0)  # w
        else:
            contracted.append(others.index(x) + 1)

    fwd_con = sum(1 for i in range(len(contracted)-1) if con_adj[contracted[i]][contracted[i+1]])
    uv_fwd = 1 if adj[P[uv_pos]][P[uv_pos+1]] else 0

    print(f"  T path {P}: fwd={fwd_T}, uv_step fwd={uv_fwd}, "
          f"contracted={contracted}: fwd_con={fwd_con}, "
          f"naive expect={fwd_T - uv_fwd}")
