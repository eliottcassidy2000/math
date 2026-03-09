"""
Attempt to prove Case 2 at n=5 algebraically.

At n=5 with b1=0 and free cycles present: Lemma B gives n-5=0 guaranteed good
vertices. We need at least 1.

The issue: (F, D) adjacent pair with V(F)={a,b,c}, V(D)={a,b,x}, dom(D)={d}.
At-risk set = {a,b,c,x,d} which could be all 5 vertices.

Question: Is there always a vertex NOT in the at-risk set? No, at n=5 the
at-risk set can be all 5 vertices. But maybe a more careful analysis works.

Alternative approaches:
1. Show that at n=5 with free cycles, there are always MULTIPLE (F,D) pairs
   giving different at-risk sets
2. Show that the "at-risk" vertices are actually still good for structural reasons
3. Use the isolation characterization directly at n=5
4. Show that at n=5 the all-dominated case covers all b1=0 tournaments
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

def shared_directed_edge(c1, c2):
    edges1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
    edges2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
    return len(edges1 & edges2) > 0

def compute_b1(A, n):
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    nc = len(cycles)
    adj_list = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cycles[i], cycles[j]):
                adj_list[i].append(j)
                adj_list[j].append(i)
    visited = [False]*nc
    b1 = 0
    for start in range(nc):
        if visited[start]: continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]: continue
            visited[v] = True
            comp.append(v)
            for u in adj_list[v]:
                if not visited[u]: stack.append(u)
        if all(not get_dom(A, n, cycles[ci]) for ci in comp):
            b1 += 1
    return b1

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def main():
    n = 5

    # Exhaustive analysis at n=5
    total = 2**(n*(n-1)//2)
    print(f"Exhaustive analysis at n={n} ({total} tournaments)")
    print("=" * 70)

    # Approach 1: Check if all b1=0 n=5 tournaments with free cycles have |dom(D)|>=2
    print("\nApproach 1: Does adjacent dominated cycle always have |dom|>=2?")
    n_free_and_b1_0 = 0
    n_has_multidom = 0
    n_only_singledom = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
        if not free_cycs:
            continue

        n_free_and_b1_0 += 1
        dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]

        has_multi = False
        for fc in free_cycs:
            for dc, doms in dom_cycs:
                if shared_directed_edge(fc, dc) and len(doms) >= 2:
                    has_multi = True
                    break
            if has_multi:
                break

        if has_multi:
            n_has_multidom += 1
        else:
            n_only_singledom += 1

    print(f"  b1=0 with free cycles: {n_free_and_b1_0}")
    print(f"  has (F adj D with |dom(D)|>=2): {n_has_multidom}")
    print(f"  only single-dom: {n_only_singledom}")
    print(f"  If has_multidom covers all: guaranteed n-4={n-4} good vertices")

    # Approach 2: When only single-dom, is the unique dominator d in V(F)?
    print("\nApproach 2: When single-dom, is d in V(F)?")
    n_d_in_F = 0
    n_d_not_in_F = 0
    n_d_in_F_good = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
        if not free_cycs:
            continue

        dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]

        for fc in free_cycs:
            fc_set = set(fc)
            for dc, doms in dom_cycs:
                if not shared_directed_edge(fc, dc):
                    continue
                if len(doms) == 1:
                    d = doms[0]
                    if d in fc_set:
                        n_d_in_F += 1
                        # at-risk = V(F) u V(D) u {d} but d in V(F)
                        # at-risk = V(F) u V(D) = 4 vertices, good = n-4 = 1
                        n_d_in_F_good += 1
                    else:
                        n_d_not_in_F += 1

    print(f"  d in V(F): {n_d_in_F} (gives n-4={n-4} good)")
    print(f"  d not in V(F): {n_d_not_in_F} (gives n-5={n-5} good)")

    # Approach 3: For every tournament, find the BEST (F,D) pair
    print("\nApproach 3: Best (F,D) pair analysis")
    min_good_per_tour = []

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
        if not free_cycs:
            continue

        dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]

        best_good = 0
        for fc in free_cycs:
            fc_set = set(fc)
            for dc, doms in dom_cycs:
                if not shared_directed_edge(fc, dc):
                    continue
                dc_set = set(dc)
                at_risk = fc_set | dc_set
                if len(doms) == 1:
                    at_risk.add(doms[0])
                good = n - len(at_risk)
                best_good = max(best_good, good)

        min_good_per_tour.append(best_good)

    print(f"  min best_good: {min(min_good_per_tour)}")
    print(f"  distribution: {dict(sorted(defaultdict(int, ((g, min_good_per_tour.count(g)) for g in set(min_good_per_tour))).items()))}")

    # Approach 4: n=5 has at most C(5,3)=10 vertex triples, giving at most 10 directed 3-cycles
    # How many free cycles can there be at n=5 with b1=0?
    print("\nApproach 4: Structure of n=5 b1=0 with free cycles")

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
        if not free_cycs:
            continue

        dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]

        # Find best pair
        best_good = 0
        best_pair = None
        for fc in free_cycs:
            fc_set = set(fc)
            for dc, doms in dom_cycs:
                if not shared_directed_edge(fc, dc):
                    continue
                dc_set = set(dc)
                at_risk = fc_set | dc_set
                if len(doms) == 1:
                    at_risk.add(doms[0])
                good = n - len(at_risk)
                if good > best_good:
                    best_good = good
                    best_pair = (fc, dc, doms, at_risk)

        if best_good == 0:
            # This tournament has NO (F,D) pair giving a guaranteed good vertex
            # Examine it
            scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
            bad_verts = [v for v in range(n) if compute_b1(delete_vertex(A, n, v), n-1) > 0]
            print(f"\n  bits={bits}, scores={scores}, #free={len(free_cycs)}, #dom={len(dom_cycs)}")
            print(f"  bad vertices: {bad_verts}")
            print(f"  free cycles: {free_cycs}")
            for dc, doms in dom_cycs:
                adj_free = [fc for fc in free_cycs if shared_directed_edge(fc, dc)]
                if adj_free:
                    print(f"  dom cycle {dc}, dom={doms}, adj to free: {adj_free}")

    # Approach 5: Maybe at n=5, the d-in-V(F) case ALWAYS applies?
    print("\n\nApproach 5: Does every n=5 b1=0 tour with free cycles have d in V(F)?")
    n_always_d_in_F = 0
    n_not_always = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
        if not free_cycs:
            continue

        dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]

        # Check: for EVERY (F,D) adj pair, is there always one where d in V(F) or |dom|>=2?
        found_good_pair = False
        for fc in free_cycs:
            fc_set = set(fc)
            for dc, doms in dom_cycs:
                if not shared_directed_edge(fc, dc):
                    continue
                if len(doms) >= 2:
                    found_good_pair = True
                    break
                if doms[0] in fc_set:
                    found_good_pair = True
                    break
            if found_good_pair:
                break

        if found_good_pair:
            n_always_d_in_F += 1
        else:
            n_not_always += 1

    print(f"  Always has good pair (|dom|>=2 or d in V(F)): {n_always_d_in_F}")
    print(f"  No good pair: {n_not_always}")

    if n_not_always == 0:
        print("\n  *** ALGEBRAIC PROOF POSSIBLE: at n=5, every (F,D) pair satisfies")
        print("  *** |dom(D)|>=2 OR d in V(F), giving at-risk <= 4, good >= 1")
    else:
        print(f"\n  Need to handle {n_not_always} cases differently")

    # Approach 6: Can we prove d in V(F) or |dom|>=2 must hold at n=5?
    # At n=5: V(F) = {a,b,c}, V(D) = {a,b,x} (share edge a->b).
    # The remaining vertices: {c, x, d} where d = dom(D).
    # If d not in V(F)={a,b,c} and d not in V(D)={a,b,x}: d is the 5th vertex.
    # But {a,b,c,x,d} = all 5 vertices.
    # V(F) = {a,b,c}, V(D) = {a,b,x}. If |dom(D)|=1 with d not in V(F):
    #   d must be x or the 5th vertex.
    #   If d = x: d in V(D) but not V(F). at_risk = {a,b,c,x,d}={a,b,c,x} but d=x, so {a,b,c,x}. Good = 1.
    #   Hmm wait, d = x means dom = {x} and V(D) = {a,b,x}. But x is IN the cycle D,
    #   so x can't be an external dominator. dom only counts EXTERNAL vertices!
    print("\n\nApproach 6: At n=5, V(F)={a,b,c}, V(D)={a,b,x}, external d")
    print("d must be external to D, so d not in {a,b,x}.")
    print("d is either c or the 5th vertex y.")
    print("If d = c: d in V(F), at-risk = {a,b,c,x} = 4, good = 1.")
    print("If d = y (5th vertex): at-risk = {a,b,c,x,y} = 5, good = 0.")
    print("So we need: for SOME (F,D) pair, d in V(F) OR |dom(D)|>=2.")
    print("Equivalently: NOT all (F,D) pairs have d = y (5th vertex) and |dom|=1.")

    # Check: when d = y (5th vertex), what does this structurally mean?
    n_all_d_is_5th = 0
    n_some_d_is_5th = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
        if not free_cycs:
            continue

        dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]

        all_pairs_d_is_5th = True
        some_d_is_5th = False

        for fc in free_cycs:
            fc_set = set(fc)
            for dc, doms in dom_cycs:
                if not shared_directed_edge(fc, dc):
                    continue
                dc_set = set(dc)
                for d in doms:
                    if d not in fc_set and d not in dc_set:
                        some_d_is_5th = True
                    else:
                        all_pairs_d_is_5th = False

        if all_pairs_d_is_5th:
            n_all_d_is_5th += 1
        elif some_d_is_5th:
            n_some_d_is_5th += 1

    print(f"\n  ALL (F,D) pairs have d outside V(F)uV(D): {n_all_d_is_5th}")
    print(f"  SOME pairs have d outside: {n_some_d_is_5th}")
    print(f"  => {n_all_d_is_5th} tournaments where Lemma B gives 0 guaranteed good")

    if n_all_d_is_5th == 0:
        print("\n  *** EVERY tournament has a (F,D) pair with d in V(F) or V(D)!")
        print("  *** This gives guaranteed good vertex algebraically!")

if __name__ == '__main__':
    main()
