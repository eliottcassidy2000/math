"""
Try to prove: NOT all vertices can be bad.

Key idea: If ALL n vertices are bad (b1(T)=0 but b1(T\v)=1 for all v),
then for each v, the free component of T\v spans n-1 vertices.

This means: for each vertex v, every 3-cycle NOT containing v either
(a) was already free in T (Dom_T = empty), or
(b) was uniquely dominated by v.

If ALL vertices are bad, this places extreme constraints.

Approach: count total unique dominations.
For each cycle C = {a,b,c}, dom(C) = #{external dominators}.
If dom(C) = 0: free. Contributes to f(v) for all v not in C.
If dom(C) = 1, say dom = {d}: uniquely dominated. Contributes to f(d).
If dom(C) >= 2: multi-dominated. Removing any single dominator doesn't free it.

For v to be bad: the freed cycles (dom subset {v}) must form a connected
free component spanning n-1 vertices. CRITICALLY: cycles with dom >= 2
remain dominated in T\v and thus CAN'T be in the free component.

So the free component of T\v consists of:
- Free cycles (dom=0) not containing v
- Cycles with dom={v} (uniquely v-dominated)

These must span ALL n-1 remaining vertices and be connected.

Now: if ALL n vertices are bad, then for each v, these freed cycles
span n-1 vertices. Let's count.

Each free cycle (dom=0) has 3 vertices. It appears in the freed set of
ALL n-3 vertices not in the cycle.
Each unique-dom cycle (dom=1, say dom={d}) appears only in d's freed set.

Total freed-cycle-vertex incidences = sum over all free components
= sum_v (sum over cycles in F_v of 3)

But each free cycle appears in n-3 freed sets, covering 3 vertices each time.
Each unique-dom cycle appears in 1 freed set, covering 3 vertices.

For the spanning requirement: each vertex u != v must appear in some
freed cycle of F_v.

If u appears only in cycles with dom >= 2 (not freed for any v):
then u is NOT in any freed cycle for any v != u.
But we need u in the freed set of every v != u (since F_v spans n-1 vertices).
Contradiction! So u must be in some cycle with dom <= 1.

Actually: for v != u, u must be in some freed cycle of F_v.
A freed cycle containing u either:
(a) is free (dom=0): then u is in a free cycle.
(b) has dom={v}: then u is in a cycle uniquely dominated by v.

If (a): u is in some free cycle. This free cycle contributes to F_w for all w
not in the cycle. Fine.

If only (b): for EACH v != u, u must be in a cycle uniquely dominated by v.
That means: for each v != u, there exists a 3-cycle C_v containing u with
Dom_T(C_v) = {v}.

How many such cycles? For each of the n-1 other vertices v, we need a
cycle through u dominated by v. These cycles are distinct (different
dominators). But u is in at most C(n-1, 2) = (n-1)(n-2)/2 three-cycles.
Each has at most n-3 external vertices (potential dominators).
The uniquely-dominated ones have exactly 1 dominator each.

Key constraint: For EACH v != u, there's a cycle through u with dom = {v}.
This means v dominates some cycle containing u. So either:
v -> u and v -> (two others in cycle), OR
u -> v and (two others in cycle) -> v. But wait: v is the dominator,
so v -> all 3 or all 3 -> v.

If v -> u,a,b (cycle {u,a,b}): v must be external to cycle.
If u,a,b -> v: v must be external.

In BOTH cases: v is NOT in the cycle (since v is the external dominator).

So for each v != u: exists cycle {u,a_v,b_v} with Dom = {v}.
There are n-1 such cycles (one per v). Each uses u and 2 other vertices
from V\{u,v}.

These n-1 cycles all contain u. Each also contains 2 other vertices from
V\{u} (total pool: n-1 vertices minus the dominator v).

Can n-1 distinct cycles through u exist? u participates in at most
C(n-1, 2) = (n-1)(n-2)/2 three-cycles. For n >= 4, we need n-1 such
cycles with distinct dominators. Since C(n-1,2) >= n-1 for n >= 4,
this is in principle possible.

But each cycle {u,a_v,b_v} has dom = {v}. This means v -> u,a_v,b_v
(or u,a_v,b_v -> v) AND no other external vertex dominates the cycle.

THE KEY: the dominator v either beats or loses to all of {u,a_v,b_v}.

Case A: v -> u (v beats u). Then v -> u,a_v,b_v.
Case B: u -> v (u beats v). Then u,a_v,b_v -> v.

For a fixed u: some v's beat u (out-neighbors of u in reverse, i.e.,
vertices in in(u)), and some v's lose to u (out-neighbors of u).

Let S = out(u) (vertices u beats), T_in = in(u) (vertices beating u).
|S| = s(u), |T_in| = n-1-s(u).

For v in T_in (v -> u): v dominates {u,a_v,b_v} from above.
  v -> u (yes), v -> a_v, v -> b_v. So a_v, b_v in out(v).
  Also: a_v, b_v != v and != u. So a_v, b_v in out(v) \ {u}.
  Since v -> u, u in out(v). So a_v, b_v in out(v) \ {u}.

For v in S (u -> v): {u,a_v,b_v} -> v. u -> v (yes), a_v -> v, b_v -> v.
  So a_v, b_v in in(v) \ {u}.
  Since u -> v, u in out(u) and v in in(u)... wait, u -> v means v in out(u) = S.

This gives very specific constraints on the structure.

Can we derive a contradiction from ALL vertices being bad?
Let's try at n=6.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def find_3cycles(A, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def get_dom(A, n, cyc):
    doms = []
    a, b, c = cyc
    for d in range(n):
        if d in {a,b,c}: continue
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            doms.append(d)
    return doms

def main():
    print("PROOF ATTEMPT: Not all vertices can be bad")
    print("=" * 60)

    # For each vertex u, count: how many other vertices v have a cycle
    # through u uniquely dominated by v?
    n = 6
    total = 2**(n*(n-1)//2)

    # Focus on tournaments where many vertices are bad
    max_bad_tours = []

    for bits in range(total):
        A = bits_to_adj(bits, n)

        # Quick check: compute b1
        cycles = find_3cycles(A, n)
        if not cycles:
            continue

        # Check domination
        from itertools import combinations as comb
        def shared_de(c1, c2):
            e1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
            e2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
            return len(e1 & e2) > 0

        nc = len(cycles)
        adj = [[] for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if shared_de(cycles[i], cycles[j]):
                    adj[i].append(j)
                    adj[j].append(i)

        visited = [False]*nc
        components = []
        for start in range(nc):
            if visited[start]: continue
            comp = []
            stack = [start]
            while stack:
                vv = stack.pop()
                if visited[vv]: continue
                visited[vv] = True
                comp.append(vv)
                for uu in adj[vv]:
                    if not visited[uu]: stack.append(uu)
            components.append(comp)

        def is_dom(cyc):
            return len(get_dom(A, n, cyc)) > 0

        b1 = sum(1 for comp in components if all(not is_dom(cycles[ci]) for ci in comp))
        if b1 != 0:
            continue

        # Count bad vertices (expensive but OK for n=6)
        bad_count = 0
        for v in range(n):
            keep = [i for i in range(n) if i != v]
            B = np.zeros((n-1, n-1), dtype=int)
            for i, ki in enumerate(keep):
                for j, kj in enumerate(keep):
                    B[i][j] = A[ki][kj]
            cyc_v = find_3cycles(B, n-1)
            if not cyc_v:
                continue
            nc_v = len(cyc_v)
            adj_v = [[] for _ in range(nc_v)]
            for i in range(nc_v):
                for j in range(i+1, nc_v):
                    if shared_de(cyc_v[i], cyc_v[j]):
                        adj_v[i].append(j)
                        adj_v[j].append(i)
            vis = [False]*nc_v
            comps_v = []
            for start in range(nc_v):
                if vis[start]: continue
                comp = []
                stack = [start]
                while stack:
                    vv = stack.pop()
                    if vis[vv]: continue
                    vis[vv] = True
                    comp.append(vv)
                    for uu in adj_v[vv]:
                        if not vis[uu]: stack.append(uu)
                comps_v.append(comp)

            def is_dom_v(cyc):
                return any(
                    all(B[d][cyc[k]] for k in range(3)) or all(B[cyc[k]][d] for k in range(3))
                    for d in range(n-1) if d not in cyc
                )

            b1_v = sum(1 for comp in comps_v if all(not is_dom_v(cyc_v[ci]) for ci in comp))
            if b1_v > 0:
                bad_count += 1

        if bad_count == 3:
            max_bad_tours.append(bits)

    print(f"n=6: {len(max_bad_tours)} tournaments with exactly 3 bad vertices")

    # For these 3-bad tournaments: analyze the domination structure
    # Question: Can we prove a vertex u always has fewer than n-1
    # distinct uniquely-dominating vertices?
    dom_per_vertex = defaultdict(int)  # u -> #{v: v uniquely dominates a cycle thru u}

    for bits in max_bad_tours[:200]:  # sample
        A = bits_to_adj(bits, n)
        cycles = find_3cycles(A, n)

        # For each vertex u: how many distinct v uniquely dominate a cycle through u?
        for u in range(n):
            unique_dominator_of_u = set()
            for cyc in cycles:
                if u not in cyc:
                    continue
                doms = get_dom(A, n, cyc)
                if len(doms) == 1:
                    unique_dominator_of_u.add(doms[0])
            dom_per_vertex[len(unique_dominator_of_u)] += 1

    print("\nCount of unique dominators of cycles through u:")
    for k in sorted(dom_per_vertex.keys()):
        print(f"  {k}: {dom_per_vertex[k]}")

    # What is the MAXIMUM number of distinct unique dominators for any vertex?
    max_unique_doms = 0
    for bits in max_bad_tours:
        A = bits_to_adj(bits, n)
        cycles = find_3cycles(A, n)
        for u in range(n):
            unique_dominator_of_u = set()
            for cyc in cycles:
                if u not in cyc:
                    continue
                doms = get_dom(A, n, cyc)
                if len(doms) == 1:
                    unique_dominator_of_u.add(doms[0])
            max_unique_doms = max(max_unique_doms, len(unique_dominator_of_u))

    print(f"\nMax unique-dominators across all vertices: {max_unique_doms}")
    print(f"n-1 = {n-1}")
    print(f"If max < n-1, then not all vertices can be bad!")

    # Check: for the GOOD vertices, do they have fewer unique dominators?
    print(f"\n{'='*60}")
    print("UNIQUE DOMINATOR COUNT: bad vs good vertices")
    print(f"{'='*60}")

    bad_udom = defaultdict(int)
    good_udom = defaultdict(int)

    for bits in max_bad_tours:
        A = bits_to_adj(bits, n)
        cycles = find_3cycles(A, n)

        # Identify bad vertices
        bad = set()
        for v in range(n):
            keep = [i for i in range(n) if i != v]
            B = np.zeros((n-1, n-1), dtype=int)
            for i, ki in enumerate(keep):
                for j, kj in enumerate(keep):
                    B[i][j] = A[ki][kj]
            cyc_v = find_3cycles(B, n-1)
            if not cyc_v: continue
            nc_v = len(cyc_v)
            adj_v = [[] for _ in range(nc_v)]
            for i in range(nc_v):
                for j in range(i+1, nc_v):
                    if shared_de(cyc_v[i], cyc_v[j]):
                        adj_v[i].append(j)
                        adj_v[j].append(i)
            vis = [False]*nc_v
            comps_v = []
            for start in range(nc_v):
                if vis[start]: continue
                comp = []
                stack = [start]
                while stack:
                    vv = stack.pop()
                    if vis[vv]: continue
                    vis[vv] = True
                    comp.append(vv)
                    for uu in adj_v[vv]:
                        if not vis[uu]: stack.append(uu)
                comps_v.append(comp)

            def is_dom_v(cyc):
                return any(
                    all(B[d][cyc[k]] for k in range(3)) or all(B[cyc[k]][d] for k in range(3))
                    for d in range(n-1) if d not in cyc
                )

            b1_v = sum(1 for comp in comps_v if all(not is_dom_v(cyc_v[ci]) for ci in comp))
            if b1_v > 0:
                bad.add(v)

        for u in range(n):
            unique_doms = set()
            for cyc in cycles:
                if u not in cyc: continue
                doms = get_dom(A, n, cyc)
                if len(doms) == 1:
                    unique_doms.add(doms[0])
            if u in bad:
                bad_udom[len(unique_doms)] += 1
            else:
                good_udom[len(unique_doms)] += 1

    print("BAD vertices - #{unique dominators of cycles through v}:")
    for k in sorted(bad_udom.keys()):
        print(f"  {k}: {bad_udom[k]}")
    print("GOOD vertices - #{unique dominators of cycles through v}:")
    for k in sorted(good_udom.keys()):
        print(f"  {k}: {good_udom[k]}")

if __name__ == '__main__':
    main()
