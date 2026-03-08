#!/usr/bin/env python3
"""
deletion_contraction_F_poly_v2.py — F(T,x) deletion-contraction, corrected.

Key insight from v1: the naive fwd-count formula works PATH-BY-PATH
  fwd_contracted = fwd_T - (1 if u->v step was forward)
This means F(T/e) encodes the fwd count EXACTLY.

The H identity H(T) = H(T\e) + H(T/e) holds (kind-pasteur verified)
where T\e is a digraph and Ham paths follow DIRECTED arcs only.

For F(T,x): in T\e (arc deleted), Ham paths can only use existing arcs.
Since we deleted u->v, paths through u,v consecutively can only go v->u
(IF v->u exists... but in a tournament with u->v, v->u does NOT exist).
So T\e has NO arc between u and v, meaning paths cannot use u,v consecutively.
F(T\e, x) = F_none(x) where F_none counts paths NOT using (u,v) edge.

BUT: "forward edge" in T\e means different things! In T, every consecutive
pair has an arc. In T\e, there's no arc u-v in either direction.
So F(T\e, x) is well-defined only for paths that DON'T visit u,v consecutively.

For these paths, the fwd count is the SAME as in T (same arcs for all other pairs).
So F(T\e, x) = F_none(x) exactly. VERIFIED in v1.

For F(T/e, x): the naive formula says
  F(T/e)[k] = F_uv[k+1] + F_vu[k]
where F_uv[k] = # paths with ...u,v... having fwd_T = k.
This was VERIFIED path-by-path in v1.

BUT the v1 test showed "Naive match? False" because I was computing
F(T/e) using the wrong Hamiltonian path counter. Let me use the
DIRECTED-only counter for T/e (since T/e may not be a tournament).

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

def ham_paths_dp_directed(adj, n):
    """Count Ham paths following DIRECTED arcs only. Returns total count."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def ham_paths_by_fwd_directed(adj, n):
    """Count Ham paths following DIRECTED arcs, tracking fwd edges.
    For a tournament: all consecutive pairs have an arc.
    For a digraph: only follow existing arcs. fwd = 1 for every step (all are forward by construction).
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
                if adj[v][u]:  # only follow directed arcs
                    new_mask = mask | (1 << u)
                    if dp[new_mask][u] is None:
                        dp[new_mask][u] = {}
                    for fc, cnt in dp[mask][v].items():
                        dp[new_mask][u][fc+1] = dp[new_mask][u].get(fc+1, 0) + cnt
    F = {}
    for v in range(n):
        if dp[full][v] is not None:
            for fc, cnt in dp[full][v].items():
                F[fc] = F.get(fc, 0) + cnt
    return F

def compute_F_tournament(adj, n):
    """F(T,x) for a tournament: fwd(P) = # steps where arc goes forward."""
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
                new_mask = mask | (1 << u)
                fwd_inc = adj[v][u]
                if dp[new_mask][u] is None:
                    dp[new_mask][u] = {}
                for fc, cnt in dp[mask][v].items():
                    nf = fc + fwd_inc
                    dp[new_mask][u][nf] = dp[new_mask][u].get(nf, 0) + cnt
    F = [0] * n
    for v in range(n):
        if dp[full][v] is not None:
            for fc, cnt in dp[full][v].items():
                F[fc] += cnt
    return F

def delete_arc(adj, n, u, v):
    B = [row[:] for row in adj]
    B[u][v] = 0
    return B

def contract_arc(adj, n, u, v):
    """Contract u->v. w gets IN from u, OUT from v."""
    others = sorted([x for x in range(n) if x != u and x != v])
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    # vertex 0 = w
    for i, x in enumerate(others):
        if adj[x][u]:  # x -> u => x -> w
            B[i+1][0] = 1
        if adj[v][x]:  # v -> x => w -> x
            B[0][i+1] = 1
        # What about adj[u][x] or adj[x][v]? These are NOT inherited by w.
        # This is the convention.
    for i, x in enumerate(others):
        for j, y in enumerate(others):
            B[i+1][j+1] = adj[x][y]
    return B, new_n

# ============================================================
# VERIFY H(T) = H(T\e) + H(T/e)  (directed-only counting)
# ============================================================
print("=" * 60)
print("H(T) = H(T\\e) + H(T/e) — directed-only counting")
print("=" * 60)

n = 4
m = n*(n-1)//2
ok = 0
total = 0

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    H_T = sum(compute_F_tournament(adj, n))

    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue
            del_adj = delete_arc(adj, n, u, v)
            H_del = ham_paths_dp_directed(del_adj, n)
            con_adj, con_n = contract_arc(adj, n, u, v)
            H_con = ham_paths_dp_directed(con_adj, con_n)
            total += 1
            if H_T == H_del + H_con:
                ok += 1

print(f"  n=4: {ok}/{total} pass")

# ============================================================
# Now: F(T,x) decomposition
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) = F_none(x) + F_through(x)")
print("F_through(x) = x * F_uv_con(x) + F_vu_con(x)")
print("where F_uv_con[k] = F_uv[k+1], F_vu_con[k] = F_vu[k]")
print("=" * 60)

# KEY INSIGHT:
# F(T,x) = F_none(x) + F_uv(x) + F_vu(x)
# F(T\e, x) should equal the F_none part (only paths not using u-v edge)
# But in a directed digraph, "F_none" means: paths not visiting u,v consecutively,
# counted by the number of FORWARD arcs. Since all arcs among non-{u,v} pairs
# are unchanged, this is well-defined.

# However! In T\e, the forward direction for every step is the ONLY direction
# (since each step follows a directed arc). So F(T\e, x) has all coefficients
# at position n-1 (every step is "forward" by construction).

# Wait, that's not right either. In T\e, for vertices a,b that are not {u,v},
# there's exactly one arc between them (same as T). So the "forward/backward"
# distinction for those pairs is the same as in T.

# For T\e, paths cannot use u,v consecutively (no arc). So F(T\e) tracks
# the same fwd count as T for those paths.

# Actually, in a DIGRAPH, not every pair has an arc. A "forward edge" in a path
# means the arc goes in the direction of travel. But if the only existing arc
# goes backward, we can't travel in that direction in a directed-only model.

# I think the issue is: in a TOURNAMENT, every step in a Ham path has an arc.
# Some arcs are "forward" (direction of travel) and some "backward".
# F(T,x) counts x^{# forward steps}.

# In T\e (digraph), some pairs have no arc at all (u,v pair). A Ham path
# in T\e CANNOT use u,v consecutively. For all other pairs, there IS an arc
# in the direction of T. So a Ham path in T\e has every step following a
# directed arc. The notion of "forward" is then: does the arc go in the
# direction of travel? Since we only follow directed arcs, EVERY step is forward!

# So in T\e, F(T\e, x) = H(T\e) * x^{n-1} ??
# No wait, that's not right. Let me reconsider.

# In a tournament, a Ham path v_0, v_1, ..., v_{n-1} visits all vertices.
# For step i, there's an arc v_i -> v_{i+1} (forward) or v_{i+1} -> v_i (backward).
# fwd(P) = # of forward steps.

# In T\e, the pair (u,v) has NO arc. So we can't have u,v consecutive in any direction.
# For all other pairs, there's exactly one arc direction (same as T).
# A "Hamiltonian path" in the DIGRAPH T\e means: visit all n vertices v_0,...,v_{n-1}
# such that for each step i, there exists arc v_i -> v_{i+1} in T\e.
# Since T\e has the arc a->b iff T has arc a->b AND (a,b) != (u,v),
# this means: every step must use a directed arc from T, EXCEPT we cannot use u->v.
# But we CAN use v->u if that arc exists in T (which it doesn't, since T has u->v).
# Wait, in T we have u->v, so v->u does NOT exist. In T\e, we remove u->v, leaving NO arc.
# So neither direction works for the u,v pair.

# Conclusion: paths in T\e are exactly paths in T that do NOT use u,v consecutively.
# And the fwd count: since every step in such a path follows a directed arc from T,
# each step IS forward. So fwd count = n-1 for EVERY Ham path in T\e.

# Hmm wait, that contradicts the earlier result where F(T\e) = F_none = [0,6,6,0].
# The F_none counts include backward steps (v_{i+1} -> v_i direction).

# I think the confusion is: in a TOURNAMENT Ham path, a step from v_i to v_{i+1}
# always exists as one of the two arcs. The "forward edge count" tracks how many
# steps go WITH the arc direction. But in the DIGRAPH T\e model (directed-only),
# we can ONLY travel along arcs, so we can't even have backward steps.

# So the two models are different:
# MODEL 1 (tournament): Any permutation is a Ham path. fwd = # arcs going forward.
# MODEL 2 (digraph): Only permutations where every step has a directed arc are Ham paths.

# For MODEL 1, F(T\e) doesn't quite make sense because T\e is not a tournament.
# For MODEL 2, every step is forward, so there's no notion of "partial fwd count".

# THE RIGHT THING: Stay in MODEL 1 (tournament model).
# F(T,x) = F_none(x) + F_uv(x) + F_vu(x)
# F_none is the same in T and T\e (paths not using u,v consecutively have the same fwd count).
# F_uv and F_vu are the paths using u,v consecutively.

# For CONTRACTION: F(T/e) in MODEL 1 on n-1 vertices.
# Each path in T/e is a permutation of {w, others}. For each step, there's
# an arc in one direction. w->x if v->x, x->w if x->u.

# The naive formula: a T-path ...u,v... with fwd=k maps to T/e-path (contract u,v to w)
# with fwd = k - (1 if u->v step is forward) = k - 1 (since u->v IS an arc = forward).
# A T-path ...v,u... with fwd=k maps to fwd = k - 0 = k (since v->u = backward in T).

# So F(T/e)[k] should be F_uv[k+1] + F_vu[k].

# But v1 showed this is FALSE. Why?

# The issue: when we contract u,v to w, the edges involving w may be
# DIFFERENT from those in T for the adjacent vertices.
# Specifically, for vertex x:
#   In T: if path has ...x,u,v,y..., the x->u step is forward iff adj[x][u]=1.
#   In T/e: the path has ...x,w,y..., the x->w step is forward iff con_adj[x_new][w]=1.
#   con_adj[x_new][0] = adj[x][u]. So x->w forward iff x->u forward. SAME.
#   And w->y forward iff v->y forward. Also SAME.

# So the non-(u,v) edges have the SAME forward/backward status!
# And the (u,v) step itself is removed.
# So the naive formula SHOULD work.

# Let me recheck why v1 showed False...

print("\nDetailed recheck at n=4, bits=0:")
adj = tournament_from_bits(4, 0)
n = 4
u, v = 1, 0  # arc 1->0

print(f"Arc {u}->{v}, adj[{u}][{v}]={adj[u][v]}")
# This is the transitive tournament with all arcs going from higher to lower index.

# Compute F_uv and F_vu from T paths
F_uv = {}
F_vu = {}
for P in permutations(range(n)):
    fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
    uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
    uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
    if uses_uv:
        F_uv[fwd] = F_uv.get(fwd, 0) + 1
    elif uses_vu:
        F_vu[fwd] = F_vu.get(fwd, 0) + 1

print(f"F_uv = {dict(sorted(F_uv.items()))}")
print(f"F_vu = {dict(sorted(F_vu.items()))}")

# Expected F(T/e)[k] = F_uv[k+1] + F_vu[k]
max_k = n - 2  # degree of F(T/e) is n-2
expected = {k: F_uv.get(k+1, 0) + F_vu.get(k, 0) for k in range(n-1)}
print(f"Expected F(T/e) = {dict(sorted(expected.items()))}")

# Actual F(T/e) using MODEL 1 (every permutation is a path, fwd counts arcs)
con_adj, con_n = contract_arc(adj, n, u, v)
print(f"\nT/e adjacency:")
for i in range(con_n):
    print(f"  {[con_adj[i][j] for j in range(con_n)]}")

F_con_actual = compute_F_tournament(con_adj, con_n)
print(f"Actual F(T/e) = {F_con_actual}")

# Check: is T/e a tournament?
is_tournament = True
for i in range(con_n):
    for j in range(i+1, con_n):
        if con_adj[i][j] + con_adj[j][i] != 1:
            is_tournament = False
            print(f"  NOT tournament: con_adj[{i}][{j}]={con_adj[i][j]}, con_adj[{j}][{i}]={con_adj[j][i]}")

print(f"T/e is tournament: {is_tournament}")

# ============================================================
# THE REAL ISSUE: T/e may NOT be a tournament!
# ============================================================
print("\n" + "=" * 60)
print("CHECK: Is T/e always a tournament?")
print("=" * 60)

n = 4
m = n*(n-1)//2
non_tournament = 0
is_tournament_count = 0

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue
            con_adj, con_n = contract_arc(adj, n, u, v)
            ok = True
            for i in range(con_n):
                for j in range(i+1, con_n):
                    if con_adj[i][j] + con_adj[j][i] != 1:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                is_tournament_count += 1
            else:
                non_tournament += 1

total = is_tournament_count + non_tournament
print(f"  T/e is tournament: {is_tournament_count}/{total}")
print(f"  T/e is NOT tournament: {non_tournament}/{total}")

# ============================================================
# WHEN T/e IS a tournament, does the naive formula work?
# ============================================================
print("\n" + "=" * 60)
print("NAIVE FORMULA: F(T/e)[k] = F_uv[k+1] + F_vu[k]")
print("Only testing when T/e IS a tournament")
print("=" * 60)

n = 4
m = n*(n-1)//2
ok_count = 0
fail_count = 0
skip_count = 0

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue

            con_adj, con_n = contract_arc(adj, n, u, v)
            is_t = all(con_adj[i][j] + con_adj[j][i] == 1
                       for i in range(con_n) for j in range(i+1, con_n))
            if not is_t:
                skip_count += 1
                continue

            # Compute F_uv, F_vu from T
            F_uv_dict = {}
            F_vu_dict = {}
            for P in permutations(range(n)):
                fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
                uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
                uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
                if uses_uv:
                    F_uv_dict[fwd] = F_uv_dict.get(fwd, 0) + 1
                elif uses_vu:
                    F_vu_dict[fwd] = F_vu_dict.get(fwd, 0) + 1

            expected = [F_uv_dict.get(k+1, 0) + F_vu_dict.get(k, 0) for k in range(con_n)]
            actual = compute_F_tournament(con_adj, con_n)

            if expected == actual:
                ok_count += 1
            else:
                if fail_count < 5:
                    print(f"  FAIL: bits={bits} arc {u}->{v}")
                    print(f"    expected = {expected}")
                    print(f"    actual   = {actual}")
                fail_count += 1

print(f"\n  Tournament T/e: ok={ok_count}, fail={fail_count}, skip(non-tournament)={skip_count}")

# ============================================================
# KEY TEST: path-by-path fwd count match (should be 100%)
# ============================================================
print("\n" + "=" * 60)
print("PATH-BY-PATH FWD COUNT MATCH")
print("=" * 60)

n = 4
m = n*(n-1)//2
mismatch = 0
total_paths = 0

for bits in range(min(16, 1 << m)):
    adj = tournament_from_bits(n, bits)
    u, v = None, None
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                u, v = i, j
                break
        if u is not None:
            break

    con_adj, con_n = contract_arc(adj, n, u, v)
    others = sorted([x for x in range(n) if x != u and x != v])

    for P in permutations(range(n)):
        # Check if u,v consecutive
        uv_pos = None
        for i in range(n-1):
            if (P[i] == u and P[i+1] == v) or (P[i] == v and P[i+1] == u):
                uv_pos = i
                break
        if uv_pos is None:
            continue

        fwd_T = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        uv_fwd = 1 if adj[P[uv_pos]][P[uv_pos+1]] else 0

        # Build contracted path
        contracted = []
        for x in P:
            if x == u or x == v:
                if not contracted or contracted[-1] != 0:
                    contracted.append(0)
            else:
                contracted.append(others.index(x) + 1)

        fwd_con = sum(1 for i in range(len(contracted)-1) if con_adj[contracted[i]][contracted[i+1]])

        total_paths += 1
        if fwd_con != fwd_T - uv_fwd:
            mismatch += 1
            if mismatch <= 3:
                print(f"  MISMATCH: bits={bits} P={P} fwd_T={fwd_T} uv_fwd={uv_fwd} fwd_con={fwd_con}")
                # Debug: check each non-uv step
                for i in range(n-1):
                    if i == uv_pos:
                        continue
                    a, b = P[i], P[i+1]
                    # Map to contracted
                    a_c = 0 if (a == u or a == v) else others.index(a) + 1
                    b_c = 0 if (b == u or b == v) else others.index(b) + 1
                    fwd_orig = adj[a][b]
                    fwd_new = con_adj[a_c][b_c]
                    if fwd_orig != fwd_new:
                        print(f"    Step {i}: {a}->{b} (mapped to {a_c}->{b_c}): orig_fwd={fwd_orig}, con_fwd={fwd_new}")
                        print(f"    In T: adj[{a}][{b}]={adj[a][b]}")
                        print(f"    In T/e: con_adj[{a_c}][{b_c}]={con_adj[a_c][b_c]}")
                        # If a or b is one of u,v, the mapping changes things
                        if a == u:
                            print(f"    a={a}=u, mapping to w. But w's OUT comes from v, not u!")
                            print(f"    adj[u][b]={adj[u][b]}, adj[v][b]={adj[v][b]}")
                        elif a == v:
                            print(f"    a={a}=v, mapping to w. w's OUT comes from v.")
                            print(f"    adj[v][b]={adj[v][b]}")
                        if b == u:
                            print(f"    b={b}=u, mapping to w. w's IN comes from u.")
                            print(f"    adj[a][u]={adj[a][u]}")
                        elif b == v:
                            print(f"    b={b}=v, mapping to w. But w's IN comes from u, not v!")
                            print(f"    adj[a][v]={adj[a][v]}, adj[a][u]={adj[a][u]}")

print(f"\n  Total paths checked: {total_paths}")
print(f"  Mismatches: {mismatch}")
if mismatch > 0:
    print(f"  => The naive formula FAILS because contraction changes edge directions!")
    print(f"     When a path has ...x,u,v,y..., the x->u step maps to x->w.")
    print(f"     In T: adj[x][u] determines forward/backward.")
    print(f"     In T/e: w's IN comes from u, so con_adj[x_map][0] = adj[x][u]. SAME.")
    print(f"     But for ...x,v,u,y..., the x->v step maps to x->w.")
    print(f"     In T: adj[x][v] determines forward/backward.")
    print(f"     In T/e: con_adj[x_map][0] = adj[x][u], which may differ from adj[x][v]!")
else:
    print(f"  => The naive formula works perfectly!")
