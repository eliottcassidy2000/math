#!/usr/bin/env python3
"""
Deep analysis of the two regular SC tournaments at n=8:
  - T_A: |Aut|=9, H=621, Omega PERFECT, mu={9:18, 3:2}
  - T_B: |Aut|=3, H=621, Omega IMPERFECT (C5!), mu={3:2, 7:12, 11:6}

Both have scores (4,4,4,4,3,3,3,3), Fix(beta)=39, 20 3-cycles.
What distinguishes them? Is T_A the Paley tournament?

Instance: opus-2026-03-05-S8
"""

from itertools import combinations, permutations
from collections import Counter


def count_ham_paths_dp(T, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1][v] for v in range(n))


def count_concordant(T, n, alpha):
    half = n // 2
    orbits = []
    seen = set()
    for v in range(n):
        if v not in seen:
            orbits.append((v, alpha[v]))
            seen.add(v)
            seen.add(alpha[v])
    count = 0
    for perm in permutations(range(len(orbits))):
        for bits in range(1 << half):
            path = [0] * n
            for i in range(half):
                path[i] = orbits[perm[i]][(bits >> i) & 1]
                path[n - 1 - i] = alpha[path[i]]
            if all(T[path[i]][path[i+1]] for i in range(n - 1)):
                count += 1
    return count


def count_automorphisms(T, n):
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    sg = {}
    for v, s in enumerate(scores):
        sg.setdefault(s, []).append(v)
    possible = [sg[scores[v]] for v in range(n)]
    count = [0]
    auts = []

    def bt(v, perm, used):
        if v == n:
            count[0] += 1
            auts.append(tuple(perm))
            return
        for t in possible[v]:
            if t in used:
                continue
            perm[v] = t
            ok = all(T[v][u] == T[t][perm[u]] and T[u][v] == T[perm[u]][t] for u in range(v))
            if ok:
                used.add(t)
                bt(v + 1, perm, used)
                used.remove(t)
            perm[v] = -1

    bt(0, [-1]*n, set())
    return count[0], auts


def build_sc(n, alpha):
    free_arcs = []
    determined = set()
    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in determined:
                continue
            ai, aj = alpha[i], alpha[j]
            u, v = min(aj, ai), max(aj, ai)
            if (u, v) == (i, j):
                free_arcs.append(((i, j), None))
                determined.add((i, j))
            else:
                free_arcs.append(((i, j), (aj, ai)))
                determined.add((i, j))
                determined.add((u, v))
    num_free = len(free_arcs)
    tournaments = []
    for bits in range(1 << num_free):
        T = [[0]*n for _ in range(n)]
        for k, (primary, linked) in enumerate(free_arcs):
            val = (bits >> k) & 1
            i, j = primary
            T[i][j] = val
            T[j][i] = 1 - val
            if linked is not None:
                a, b = linked
                T[a][b] = val
                T[b][a] = 1 - val
        tournaments.append(T)
    return tournaments


def find_3cycles(T, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if T[a][b] and T[b][c] and T[c][a]:
            cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            cycles.append((a, c, b))
    return cycles


def independence_poly(adj, m):
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        verts = [i for i in range(m) if mask & (1 << i)]
        ok = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            coeffs[len(verts)] += 1
    return coeffs


def mu_value(T, n, cycle_verts):
    remaining = [v for v in range(n) if v not in cycle_verts]
    if len(remaining) < 3:
        return 1
    sub_cycles = []
    for a, b, c in combinations(remaining, 3):
        if T[a][b] and T[b][c] and T[c][a]:
            sub_cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            sub_cycles.append((a, c, b))
    if not sub_cycles:
        return 1
    sm = len(sub_cycles)
    cs = [set(c) for c in sub_cycles]
    adj = [[0]*sm for _ in range(sm)]
    for i in range(sm):
        for j in range(i+1, sm):
            if cs[i] & cs[j]:
                adj[i][j] = adj[j][i] = 1
    coeffs = independence_poly(adj, sm)
    return sum(c * (2**k) for k, c in enumerate(coeffs))


def build_omega(T, n):
    cycles = find_3cycles(T, n)
    m = len(cycles)
    cs = [set(c) for c in cycles]
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cs[i] & cs[j]:
                adj[i][j] = adj[j][i] = 1
    return cycles, adj


def chromatic_number_bound(adj, m):
    """Greedy coloring upper bound."""
    colors = [-1] * m
    for v in range(m):
        used = set()
        for u in range(m):
            if adj[v][u] and colors[u] >= 0:
                used.add(colors[u])
        c = 0
        while c in used:
            c += 1
        colors[v] = c
    return max(colors) + 1 if colors else 0


def clique_number(adj, m):
    """Find maximum clique size."""
    adj_sets = [set(j for j in range(m) if adj[i][j]) for i in range(m)]
    best = [0]

    def bt(clique, candidates):
        if not candidates:
            best[0] = max(best[0], len(clique))
            return
        if len(clique) + len(candidates) <= best[0]:
            return
        for i, v in enumerate(candidates):
            new_cands = [u for u in candidates[i+1:] if u in adj_sets[v]]
            bt(clique + [v], new_cands)

    bt([], list(range(m)))
    return best[0]


def has_induced_C5(adj_sets, m):
    for v0 in range(m):
        for v1 in adj_sets[v0]:
            if v1 <= v0:
                continue
            for v2 in adj_sets[v1]:
                if v2 <= v0 or v2 in adj_sets[v0]:
                    continue
                for v3 in adj_sets[v2]:
                    if v3 <= v0 or v3 == v1 or v3 in adj_sets[v0] or v3 in adj_sets[v1]:
                        continue
                    for v4 in adj_sets[v3]:
                        if v4 <= v0 or v4 == v1 or v4 == v2:
                            continue
                        if v4 in adj_sets[v0] and v4 not in adj_sets[v1] and v4 not in adj_sets[v2]:
                            return True, (v0, v1, v2, v3, v4)
    return False, None


def main():
    n = 8
    alpha = [7, 6, 5, 4, 3, 2, 1, 0]

    print("=" * 70)
    print("DEEP ANALYSIS: Two regular SC tournaments at n=8")
    print("=" * 70)

    # Find the two specific tournaments
    tournaments = build_sc(n, alpha)

    T_A = None  # |Aut|=9, H=621
    T_B = None  # |Aut|=3, H=621
    T_657 = None  # H=657

    for T in tournaments:
        out_deg = sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)])
        if out_deg != [3, 3, 3, 3, 4, 4, 4, 4]:
            continue

        H = count_ham_paths_dp(T, n)
        if H != 621 and H != 657:
            continue

        scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
        sg = {}
        for v, s in enumerate(scores):
            sg.setdefault(s, []).append(v)
        if all(len(g) == 1 for g in sg.values()):
            aut_count = 1
        else:
            aut_count, _ = count_automorphisms(T, n)

        if aut_count <= 1:
            continue

        if H == 621 and aut_count == 9 and T_A is None:
            T_A = [row[:] for row in T]
            print(f"Found T_A: H={H}, |Aut|={aut_count}")
        elif H == 621 and aut_count == 3 and T_B is None:
            T_B = [row[:] for row in T]
            print(f"Found T_B: H={H}, |Aut|={aut_count}")
        elif H == 657 and T_657 is None:
            T_657 = [row[:] for row in T]
            print(f"Found T_657: H={H}, |Aut|={aut_count}")

        if T_A and T_B and T_657:
            break

    for label, T in [("T_A (|Aut|=9, H=621)", T_A),
                      ("T_B (|Aut|=3, H=621)", T_B),
                      ("T_657 (H=657)", T_657)]:
        if T is None:
            print(f"\n{label}: NOT FOUND")
            continue

        print(f"\n{'='*60}")
        print(f"{label}")
        print(f"{'='*60}")

        H = count_ham_paths_dp(T, n)
        fix = count_concordant(T, n, alpha)
        aut_size, auts = count_automorphisms(T, n)

        print(f"  H = {H}")
        print(f"  |Aut| = {aut_size}")
        print(f"  Fix(beta) = {fix}")
        print(f"  (H-Fix)/2 = {(H-fix)//2}  ({'even' if (H-fix)//2 % 2 == 0 else 'odd'})")
        print(f"  H/|Aut| = {H//aut_size}  ({'even' if (H//aut_size) % 2 == 0 else 'odd'})")
        print(f"  H mod 4 = {H % 4}")
        print(f"  H mod 8 = {H % 8}")

        out_deg = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
        print(f"\n  Adjacency matrix:")
        for row in T:
            print(f"    {''.join(str(x) for x in row)}")
        print(f"  Out-degrees: {out_deg}")

        # Check doubly regular
        print(f"\n  Regularity analysis:")
        is_reg = len(set(out_deg)) == 1
        if not is_reg:
            # Nearly regular check
            vals = set(out_deg)
            print(f"  Out-degree values: {vals}")
        else:
            print(f"  Regular with degree {out_deg[0]}")
            # Doubly regular: for each oriented pair (i,j), count |{k: T(i,k)=T(j,k)=1}|
            # For i->j: common out-neighbors
            common_adj = {}
            common_non = {}
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    ct = sum(1 for k in range(n) if k != i and k != j and T[i][k] and T[j][k])
                    if T[i][j]:
                        common_adj.setdefault(ct, 0)
                        common_adj[ct] += 1
                    else:
                        common_non.setdefault(ct, 0)
                        common_non[ct] += 1
            print(f"  Common out-neighbors (i->j): {common_adj}")
            print(f"  Common out-neighbors (j->i): {common_non}")
            if len(common_adj) == 1 and len(common_non) == 1:
                print(f"  DOUBLY REGULAR!")

        # Automorphism group structure
        print(f"\n  Automorphism group ({aut_size} elements):")
        cycle_types = Counter()
        for perm in auts:
            # Compute cycle type
            seen = set()
            cycles = []
            for v in range(n):
                if v in seen:
                    continue
                cyc = []
                u = v
                while u not in seen:
                    seen.add(u)
                    cyc.append(u)
                    u = perm[u]
                cycles.append(len(cyc))
            ct = tuple(sorted(cycles, reverse=True))
            cycle_types[ct] += 1
        for ct, cnt in sorted(cycle_types.items()):
            print(f"    {ct}: {cnt} elements")

        # Omega analysis
        cycles_3, omega_adj = build_omega(T, n)
        m = len(cycles_3)
        adj_sets = [set(j for j in range(m) if omega_adj[i][j]) for i in range(m)]
        print(f"\n  Omega(T) analysis:")
        print(f"    Vertices (3-cycles): {m}")
        edges = sum(sum(row) for row in omega_adj) // 2
        print(f"    Edges: {edges}")
        degrees = [sum(omega_adj[i]) for i in range(m)]
        print(f"    Degree sequence: {sorted(degrees)}")

        c5_result = has_induced_C5(adj_sets, m)
        print(f"    Has induced C5: {c5_result[0]}")
        if c5_result[0]:
            c5 = c5_result[1]
            print(f"    C5 vertices: {c5}")
            print(f"    C5 cycles: {[cycles_3[v] for v in c5]}")

        omega_clique = clique_number(omega_adj, m)
        omega_chi = chromatic_number_bound(omega_adj, m)
        print(f"    Clique number (omega): {omega_clique}")
        print(f"    Chromatic number (upper bound): {omega_chi}")
        print(f"    Perfect requires omega = chi: {'CONSISTENT' if omega_clique == omega_chi else 'VIOLATED' if omega_clique < omega_chi else '?'}")

        coeffs = independence_poly(omega_adj, m)
        I2 = sum(c * (2**k) for k, c in enumerate(coeffs))
        print(f"    Independence polynomial: {coeffs}")
        print(f"    I(Omega, 2) = {I2}")
        print(f"    OCF (I=H): {I2 == H}")

        # Check roots
        import numpy as np
        poly_coeffs = list(reversed(coeffs[:omega_clique+2]))  # highest degree first
        if len(poly_coeffs) > 1:
            roots = np.roots(poly_coeffs)
            real_roots = [r for r in roots if abs(r.imag) < 1e-10]
            print(f"    All roots real: {len(real_roots) == len(roots)}")
            if real_roots:
                print(f"    Root range: [{min(r.real for r in real_roots):.4f}, {max(r.real for r in real_roots):.4f}]")
            neg_real = all(r.real < 0 for r in real_roots)
            print(f"    All roots negative: {neg_real}")

        # mu value analysis
        print(f"\n  mu values for 3-cycles:")
        mu_data = []
        for c in cycles_3:
            mu = mu_value(T, n, set(c))
            mu_data.append((c, mu))
        mu_dist = Counter(m for _, m in mu_data)
        print(f"    Distribution: {dict(mu_dist)}")
        print(f"    Sum of mu values: {sum(m for _, m in mu_data)}")

        # Per-vertex 3-cycle analysis
        print(f"\n  Per-vertex D_v = sum of mu over 3-cycles through v:")
        for v in range(n):
            v_cycles = [(c, m) for c, m in mu_data if v in c]
            D_v = sum(m for _, m in v_cycles)
            print(f"    v={v}: {len(v_cycles)} cycles, D_v = {D_v}")

        # 5-cycles
        print(f"\n  5-cycles:")
        cycles_5 = []
        for verts in combinations(range(n), 5):
            for perm in permutations(verts):
                if all(T[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                    canon = min(perm[j:] + perm[:j] for j in range(5))
                    if canon not in cycles_5:
                        cycles_5.append(canon)
                    break
        print(f"    Count: {len(cycles_5)}")
        if cycles_5:
            mu5_dist = Counter()
            for c in cycles_5:
                mu = mu_value(T, n, set(c))
                mu5_dist[mu] += 1
            print(f"    mu distribution: {dict(mu5_dist)}")
            print(f"    Sum of 5-cycle mu: {sum(mu5_dist[k]*k for k in mu5_dist)}")

        # Paley comparison
        print(f"\n  Paley P(7) comparison:")
        P7 = [[0]*7 for _ in range(7)]
        qr = {1, 2, 4}
        for i in range(7):
            for j in range(7):
                if i != j and (j - i) % 7 in qr:
                    P7[i][j] = 1
        P7_H = count_ham_paths_dp(P7, 7)
        P7_aut, _ = count_automorphisms(P7, 7)
        print(f"    P(7): H={P7_H}, |Aut|={P7_aut}")

        # Check if T minus any vertex is isomorphic to P7
        for v in range(n):
            verts = [u for u in range(n) if u != v]
            subT = [[T[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
            sub_H = count_ham_paths_dp(subT, 7)
            sub_out = sorted([sum(subT[i][j] for j in range(7) if j != i) for i in range(7)])
            if sub_out == [3]*7:  # Regular on 7 vertices
                print(f"    T - {v}: regular, H={sub_H}")

        # Signed position identity
        print(f"\n  Signed position identity check:")
        # Enumerate paths for this
        paths = []
        def enum_paths(path, visited):
            if len(path) == n:
                paths.append(tuple(path))
                return
            last = path[-1]
            for v in range(n):
                if v not in visited and T[last][v]:
                    visited.add(v)
                    path.append(v)
                    enum_paths(path, visited)
                    path.pop()
                    visited.remove(v)
        for start in range(n):
            enum_paths([start], {start})

        print(f"    Total paths: {len(paths)}")

        # Check all pairs
        fail_count = 0
        for i in range(n):
            for j in range(i+1, n):
                if T[i][j]:
                    lhs = sum((-1)**P.index(i) for P in paths if P.index(i) < P.index(j))
                    rhs = sum((-1)**P.index(j) for P in paths if P.index(j) < P.index(i))
                    if lhs != rhs:
                        fail_count += 1
                        if fail_count <= 3:
                            print(f"    FAIL: ({i},{j}): LHS={lhs}, RHS={rhs}")
        print(f"    Failures: {fail_count} / {sum(sum(row) for row in T)//2} arcs")

        # Deeper: pos-weighted path analysis
        print(f"\n  Pos-weighted analysis:")
        # For each vertex v, compute sum_{P} (-1)^{pos(v,P)}
        for v in range(min(4, n)):
            total = sum((-1)**P.index(v) for P in paths)
            print(f"    v={v}: sum_P (-1)^pos(v) = {total}")


if __name__ == "__main__":
    main()
