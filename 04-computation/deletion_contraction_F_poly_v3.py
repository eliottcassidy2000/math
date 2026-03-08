#!/usr/bin/env python3
"""
deletion_contraction_F_poly_v3.py — F(T,x) and deletion-contraction.

H(T) = H(T\e) + H(T/e) VERIFIED (kind-pasteur, 100% at n=4,5).
T\e = delete arc u->v (digraph, no arc between u,v).
T/e = contract u,v into w (IN from u, OUT from v).

Question: How does this lift to F(T,x)?

Key finding from v2: contraction changes edge directions for v->u paths,
so naive formula F(T/e)[k] = F_uv[k+1] + F_vu[k] FAILS.

This script: find the CORRECT F-polynomial relationship.

Author: opus-2026-03-07-S45
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
    """F(T,x) for tournament T (all permutations valid, count forward arcs)."""
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def ham_paths_dp(adj, n):
    """Ham path count in DIRECTED graph (only following arcs)."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def delete_arc(adj, n, u, v):
    B = [row[:] for row in adj]
    B[u][v] = 0
    return B

def contract_arc(adj, n, u, v):
    others = sorted([x for x in range(n) if x != u and x != v])
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    for i, x in enumerate(others):
        if adj[x][u]:
            B[i+1][0] = 1
        if adj[v][x]:
            B[0][i+1] = 1
    for i, x in enumerate(others):
        for j, y in enumerate(others):
            B[i+1][j+1] = adj[x][y]
    return B, new_n

# ============================================================
# VERIFY H IDENTITY FIRST
# ============================================================
print("=" * 60)
print("VERIFY: H(T) = H(T\\e) + H(T/e)")
print("=" * 60)

for n in [4, 5]:
    m = n*(n-1)//2
    ok = 0
    total = 0
    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)
        for u in range(n):
            for v in range(n):
                if u != v and adj[u][v]:
                    del_adj = delete_arc(adj, n, u, v)
                    H_del = ham_paths_dp(del_adj, n)
                    con_adj, con_n = contract_arc(adj, n, u, v)
                    H_con = ham_paths_dp(con_adj, con_n)
                    total += 1
                    if H_T == H_del + H_con:
                        ok += 1
    print(f"  n={n}: {ok}/{total} pass")

# ============================================================
# F(T,x) DECOMPOSITION BY ARC USAGE
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) vs CONTRACTION — searching for correct formula")
print("=" * 60)

# For each tournament T and arc u->v:
# F(T,x) = F_none(x) + F_uv(x) + F_vu(x)
# where:
#   F_none[k] = # perms with fwd=k not using u,v consecutively
#   F_uv[k] = # perms with fwd=k and ...u,v... (forward arc step)
#   F_vu[k] = # perms with fwd=k and ...v,u... (backward arc step)
#
# F(T\e) = F_none (paths not using u,v edge, directed-only)
# But F_none counts all perms, not just directed ones...
#
# Actually, H(T\e) counts only paths following DIRECTED arcs.
# F(T\e) in the directed sense: every step follows an arc.
# Since there's no arc between u,v, these paths never have u,v consecutive.
#
# But F for a digraph is tricky: "forward" = following a directed arc.
# In a tournament, every step has exactly one arc direction.
# In T\e, every non-(u,v) step has exactly one arc. But u,v steps are impossible.
# So every step in a T\e Hamiltonian path IS forward.
# Thus H(T\e) = F(T\e, x) evaluated at... well, F(T\e) is just H(T\e)*x^{n-1}.
# No, that makes no sense either.

# I think the right way to think about this:
# In a TOURNAMENT, for any permutation P, we define fwd(P) = # forward arcs.
# This requires that every pair has an arc.
# In a DIGRAPH, not every pair has an arc. So "F(T\e, x)" with this definition
# doesn't naturally make sense.

# ALTERNATIVE: Define F_digraph(D, x) = sum over directed Ham paths of x^{fwd(P)}
# where fwd(P) = number of steps. But that's just H(D) * x^{n-1} since
# every step follows an arc.

# So deletion-contraction for F(T,x) requires a DIFFERENT approach.

# APPROACH: Work entirely in the tournament world.
# F(T,x) = F_none(x) + F_uv(x) + F_vu(x)
# F_uv counts paths with ...u,v... step (arc goes forward, contributing x)
# F_vu counts paths with ...v,u... step (arc goes backward, contributing 1)
# Now: what is F_uv + F_vu as a polynomial?

# For CONTRACTION: each T-path with u,v consecutive maps to a path in T/e
# (a tournament on n-1 vertices if T/e is a tournament).
# But the fwd count in T/e DIFFERS from the naive formula because
# for v->u paths, the neighboring edges change direction.

# Let me compute F(T/e) directly and see how it relates to F_uv + F_vu.

n = 4
m = n*(n-1)//2

print(f"\nn={n}: examining F(T), F(T\\e), F(T/e) for all arcs")

# For each arc, record the "delta polynomial" F(T) - F_none(T,e) - x*F(T/e)
# and see if there's a pattern.

all_deltas = []

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)

    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue

            # F_none: paths not using u,v consecutively (with T's fwd counting)
            F_none = [0]*n
            F_uv_poly = [0]*n
            F_vu_poly = [0]*n
            for P in permutations(range(n)):
                fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
                uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
                uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
                if uses_uv:
                    F_uv_poly[fwd] += 1
                elif uses_vu:
                    F_vu_poly[fwd] += 1
                else:
                    F_none[fwd] += 1

            # H(T\e) counts DIRECTED Ham paths. F_none counts tournament perms not using u,v.
            # H(T\e) = sum(F_none) + F_vu (paths using v->u which is NOT an arc in T\e)
            # Wait: F_none = perms NOT using u,v consecutively.
            # H(T\e) = # perms where every step follows a directed arc in T\e.
            # T\e has no arc u-v. So any perm using u,v consecutively is NOT a Ham path in T\e.
            # For perms NOT using u,v consecutively: every step has an arc (same as T).
            # But wait: in T, both v->u and u->v have an arc. In T\e, NEITHER has an arc.
            # So: a perm not using u,v consecutively has all steps following T arcs.
            # Is every such step a directed arc in T\e? YES: T\e only removes u->v.
            # For steps not involving u,v pair: same arc as T.
            # So H(T\e) = sum(F_none) -- this should hold!
            del_adj = delete_arc(adj, n, u, v)
            H_del = ham_paths_dp(del_adj, n)
            # Actually F_none counts perms where u,v DON'T appear consecutively,
            # with fwd counting based on T. The TOTAL count sum(F_none) should equal H(T\e).
            # Debug:
            if H_del != sum(F_none):
                # Count explicitly
                H_del_check = 0
                for P in permutations(range(n)):
                    valid = True
                    for i in range(n-1):
                        if not del_adj[P[i]][P[i+1]] and not del_adj[P[i+1]][P[i]]:
                            # No arc in either direction between consecutive vertices
                            valid = False
                            break
                    if valid:
                        H_del_check += 1
                # Hmm: H_del counts paths following ONLY directed arcs.
                # A path ...a,b... requires del_adj[a][b]=1 (not del_adj[b][a]=1).
                # In a tournament, every step has one direction.
                # For steps not involving u,v: exactly one of adj[a][b], adj[b][a] is 1.
                # In T\e, same arcs except u->v removed.
                # So F_none perms: every step not involving u,v has an arc.
                # But these are PERMUTATIONS, not directed paths.
                # A permutation P visits v_0,v_1,...,v_{n-1}. Some steps go "with" the arc
                # (forward) and some "against" (backward).
                # H(T\e) only counts perms where EVERY step goes WITH the arc.
                # F_none perms include ones where some steps go AGAINST the arc.
                # So sum(F_none) != H(T\e) in general!
                pass

            # F(T/e) as tournament polynomial
            con_adj, con_n = contract_arc(adj, n, u, v)
            # Check if T/e is a tournament
            is_t = all(con_adj[i][j] + con_adj[j][i] == 1
                       for i in range(con_n) for j in range(i+1, con_n))
            if not is_t:
                # Skip non-tournament contractions for now
                continue

            F_con = compute_F(con_adj, con_n)

            # H identity check
            H_T = sum(F_T)
            H_con = sum(F_con)
            # H(T) = H(T\e) + H(T/e)
            assert H_T == H_del + H_con, f"H identity failed"

            # Now: F_uv + F_vu counts paths using u,v consecutively.
            # sum(F_uv) + sum(F_vu) should equal H_con.
            H_through = sum(F_uv_poly) + sum(F_vu_poly)
            assert H_through == H_con, f"H_through mismatch: {H_through} vs {H_con}"

            # F(T,x) = F_none(x) + F_uv(x) + F_vu(x)
            # We want: F_uv(x) + F_vu(x) = some function of F(T/e, x)?

            # Direct comparison:
            through_poly = [F_uv_poly[k] + F_vu_poly[k] for k in range(n)]
            # F_con has degree n-2
            # through_poly has degree n-1

            # Maybe: through_poly(x) = (1 + x) * F_con(x) / something?
            # Check: through_poly(x) = sum_k through[k]*x^k
            # (1+x)*F_con(x) = sum_k (F_con[k] + F_con[k-1])*x^k

            prod_1px = [(F_con[k] if k < con_n else 0) + (F_con[k-1] if k > 0 and k-1 < con_n else 0) for k in range(n)]
            if through_poly == prod_1px:
                all_deltas.append(('match_1px', bits, u, v))
            else:
                # Try: through_poly = x * F_con(x) + F_con(x)?
                # That's the same as (1+x)*F_con.
                # Try: through_poly = x * G(x) + G(x) for some other G?
                all_deltas.append(('nomatch', bits, u, v, through_poly, F_con, prod_1px))

# Count
matches = sum(1 for d in all_deltas if d[0] == 'match_1px')
nomatches = sum(1 for d in all_deltas if d[0] == 'nomatch')
print(f"  through(x) = (1+x)*F(T/e,x): {matches}/{matches+nomatches}")

# Show some nomatches
for d in all_deltas:
    if d[0] == 'nomatch':
        _, bits, u, v, through, fcon, prod = d
        print(f"  bits={bits} arc {u}->{v}: through={through}, (1+x)*F_con={prod}")
        break

# Try other relationships
print(f"\n  Trying other polynomial relationships:")

# F_uv(x) = x * something?
# F_uv[k] = # paths with ...u,v... and fwd=k.
# The u->v step contributes 1 to fwd. So F_uv[k] > 0 only for k >= 1.
# Factor out x: F_uv(x) = x * G_uv(x) where G_uv[k] = F_uv[k+1].

# Similarly F_vu[k] = # paths with ...v,u... and fwd=k.
# The v->u step contributes 0. So F_vu can have F_vu[0] > 0.

# So through(x) = x*G_uv(x) + F_vu(x).

# F(T,x) = F_none(x) + x*G_uv(x) + F_vu(x)
# H: sum F_T = sum F_none + sum G_uv + sum F_vu

# Can we express G_uv and F_vu in terms of F(T/e)?
# If the fwd count SHIFTED perfectly: G_uv[k] = F_uv[k+1] counts paths
# where the non-(u,v) steps have fwd_count = k.
# And F_vu[k] counts paths where the non-(v,u) steps have fwd_count = k.
# The contracted path should have the SAME non-(u,v) fwd count...
# EXCEPT it doesn't, because for v->u paths, the edges involving w differ.

# Let me just compute, for each contracted path in T/e, how many T-paths
# it corresponds to in the F_uv and F_vu categories, with what fwd counts.

print(f"\n  Path-by-path mapping at n=4:")
adj = tournament_from_bits(4, 0)
n = 4
u, v = 1, 0
others = sorted([x for x in range(n) if x != u and x != v])
con_adj, con_n = contract_arc(adj, n, u, v)

# For each contracted path, find all T-paths that map to it
for Q in permutations(range(con_n)):
    fwd_con = sum(1 for i in range(con_n-1) if con_adj[Q[i]][Q[i+1]])

    # Map Q back to T-paths: w (=vertex 0) expands to either (u,v) or (v,u)
    w_pos = list(Q).index(0)  # position of w in Q

    # Expand w to u,v or v,u
    for order in [(u, v), (v, u)]:
        P_list = list(Q)
        P_list[w_pos:w_pos+1] = list(order)
        # Map back to original vertex indices
        P = tuple(others[x-1] if x > 0 else P_list[P_list.index(x)] for x in P_list)
        # Actually, let me be more careful
        P = []
        for x in P_list:
            if x == 0:
                # This should have been expanded already... wait
                pass
            else:
                P.append(x)
        # Hmm, the expansion is already done. P_list has u,v or v,u at w_pos.
        # But the other entries are 1, 2 (mapped from others = [2, 3]).
        # Need to map 1 -> others[0] = 2, 2 -> others[1] = 3.
        P_mapped = []
        for x in P_list:
            if x == u or x == v:
                P_mapped.append(x)
            elif isinstance(x, int) and x > 0:
                P_mapped.append(others[x-1])
            else:
                P_mapped.append(x)
        P_tuple = tuple(P_mapped)
        fwd_T = sum(1 for i in range(n-1) if adj[P_tuple[i]][P_tuple[i+1]])
        is_uv = (P_tuple[w_pos] == u and P_tuple[w_pos+1 if w_pos+1 < n else w_pos] == v) if w_pos < n-1 else False
        # Simpler check
        uv_step = any(P_tuple[i]==u and P_tuple[i+1]==v for i in range(n-1))

        step_type = "u->v(fwd)" if uv_step else "v->u(bwd)"
        print(f"    Q={Q} fwd_con={fwd_con} -> P={P_tuple} fwd_T={fwd_T} ({step_type})")

# ============================================================
# CLEANER APPROACH: Compute the "through polynomial" directly
# ============================================================
print("\n" + "=" * 60)
print("RELATIONSHIP: through(T,e,x) vs F(T/e,x)")
print("through[k] = F_uv[k] + F_vu[k]")
print("=" * 60)

n = 4
m = n*(n-1)//2

# For each (T, e) where T/e is a tournament, compute through(x) and F(T/e,x)
# and find the transformation matrix.

examples = []
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
                continue

            F_con = compute_F(con_adj, con_n)

            through = [0]*n
            for P in permutations(range(n)):
                fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
                uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
                uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
                if uses_uv or uses_vu:
                    through[fwd] += 1

            examples.append((through, F_con, bits, u, v))

# Check: is through always (1+x)*F_con?
match_1px = 0
for through, fcon, *_ in examples:
    prod = [(fcon[k] if k < len(fcon) else 0) + (fcon[k-1] if k > 0 and k-1 < len(fcon) else 0) for k in range(n)]
    if through == prod:
        match_1px += 1

print(f"  through(x) = (1+x)*F(T/e,x): {match_1px}/{len(examples)}")

if match_1px < len(examples):
    # Find the actual relationship
    # through has degree n-1 = 3, F_con has degree n-2 = 2.
    # through(x) = A(x)*F_con(x) for some A(x)?
    # through has n coefficients, F_con has n-1.
    # If A(x) = a + bx, then through[k] = a*F_con[k] + b*F_con[k-1]
    # This gives n equations in 2 unknowns. Likely overdetermined.

    # Let's just look at the actual values
    print("\n  Examples:")
    for through, fcon, bits, u, v in examples[:10]:
        prod = [(fcon[k] if k < len(fcon) else 0) + (fcon[k-1] if k > 0 and k-1 < len(fcon) else 0) for k in range(n)]
        print(f"    bits={bits} arc {u}->{v}: through={through}, F_con={fcon}, (1+x)*F_con={prod}")

    # Maybe: through(x) = (ax+b) * F_con(x) + R(x)?
    # Or: F(T,x) = F_none(x) + (1+x)*F_con(x) + correction?
    # Where correction depends on the "edge disagreement" between v and u's neighbors.

    # The disagreement count: for how many x in others do adj[x][u] != adj[x][v]
    # and adj[u][x] != adj[v][x]?
    print("\n  Disagreement analysis:")
    for through, fcon, bits, u, v in examples[:10]:
        adj = tournament_from_bits(n, bits)
        others_list = sorted([x for x in range(n) if x != u and x != v])
        # IN disagreement: adj[x][u] != adj[x][v]
        in_disagree = sum(1 for x in others_list if adj[x][u] != adj[x][v])
        # OUT disagreement: adj[u][x] != adj[v][x]
        out_disagree = sum(1 for x in others_list if adj[u][x] != adj[v][x])
        prod = [(fcon[k] if k < len(fcon) else 0) + (fcon[k-1] if k > 0 and k-1 < len(fcon) else 0) for k in range(n)]
        diff = [through[k] - prod[k] for k in range(n)]
        print(f"    bits={bits} arc {u}->{v}: in_dis={in_disagree}, out_dis={out_disagree}, diff={diff}")

# ============================================================
# APPROACH: F(T,x) = F(T_flip, x) + (x-1) * G(e, x)
# where T_flip flips arc u->v to v->u
# ============================================================
print("\n" + "=" * 60)
print("FLIP APPROACH: F(T,x) - F(T_flip,x)")
print("=" * 60)

n = 4
m = n*(n-1)//2

for bits in range(min(8, 1 << m)):
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)

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

    diff = [F_T[k] - F_flip[k] for k in range(n)]
    # Check if (x-1) divides diff(x)
    # diff evaluated at x=1: sum(diff) = H(T) - H(T_flip)
    H_diff = sum(diff)
    # diff evaluated at x=0: diff[0]
    # diff[k] = F_T[k] - F_flip[k]

    con_adj, con_n = contract_arc(adj, n, u, v)
    is_t = all(con_adj[i][j] + con_adj[j][i] == 1
               for i in range(con_n) for j in range(i+1, con_n))

    if is_t:
        F_con = compute_F(con_adj, con_n)
    else:
        F_con = None

    print(f"  bits={bits} arc {u}->{v}: F_T={F_T}, F_flip={F_flip}")
    print(f"    diff = {diff}, H_diff = {H_diff}")
    if F_con is not None:
        # Try: diff(x) = (x-1) * F_con(x)?
        test = [-(F_con[k] if k < len(F_con) else 0) + (F_con[k-1] if k > 0 and k-1 < len(F_con) else 0) for k in range(n)]
        print(f"    (x-1)*F_con = {test}")
        print(f"    match: {diff == test}")

        # Try: diff(x) = 2*(x-1) * F_con(x)?
        test2 = [2*t for t in test]
        print(f"    2*(x-1)*F_con = {test2}")
        print(f"    match: {diff == test2}")
