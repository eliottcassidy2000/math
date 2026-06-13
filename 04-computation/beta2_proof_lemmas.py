"""
Verify the key lemmas for the beta_2 = 0 proof.

LEMMA A: If b1(T) = 0 and T has free cycles, then some free cycle F is
directly adjacent (shares directed edge) to some dominated cycle D.

Proof: b1=0 means no all-free connected component. So every component
containing a free cycle also has a dominated cycle. By connectivity,
there's a path from some free cycle to some dominated cycle. On this path,
there must be consecutive cycles where one is free and the other is dominated.
Hence DIRECT adjacency.

LEMMA B (Main): If b1(T) = 0 and T has free cycles, then there exist
at least max(n-3, n-4) good vertices, depending on whether the adjacent
dominated cycle has 1 or 2+ dominators.

LEMMA C: If b1(T) = 0 and T has NO free cycles, does a good vertex still exist?

The full proof structure:
- No 3-cycles: b1=0 for all subgraphs, all good.
- Free cycles exist + b1=0: Lemma A + Lemma B give good vertex for n >= 5.
- All dominated: Lemma C needed.
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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
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
    adj = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cycles[i], cycles[j]):
                adj[i].append(j)
                adj[j].append(i)
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
            for u in adj[v]:
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
    rng = np.random.RandomState(42)

    # LEMMA A: Free cycle always adjacent to dominated cycle when b1=0
    print("=" * 70)
    print("LEMMA A: Free cycle always adjacent to dominated (given b1=0)")
    print("=" * 70)

    for n in [4, 5, 6, 7, 8, 9, 10]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000, 10: 1000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        total_b1_0 = 0
        has_free = 0
        free_adj_dom = 0
        free_not_adj_dom = 0

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue
            total_b1_0 += 1

            cycles = find_3cycles(A, n)
            free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
            if not free_cycs:
                continue
            has_free += 1

            dom_cycs = [cyc for cyc in cycles if get_dom(A, n, cyc)]
            found = any(shared_directed_edge(fc, dc)
                       for fc in free_cycs for dc in dom_cycs)
            if found:
                free_adj_dom += 1
            else:
                free_not_adj_dom += 1

        print(f"n={n} ({method}): b1=0: {total_b1_0}, has free: {has_free}")
        print(f"  free adj dom: {free_adj_dom}, NOT adj: {free_not_adj_dom}")

    # LEMMA B: Given free adj dom, count good vertices
    print("\n" + "=" * 70)
    print("LEMMA B: Good vertex count from free-dom adjacency")
    print("=" * 70)

    for n in [5, 6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        # For each (F, D) pair with F free and D dominated:
        # Vertices not in F and not unique-dom of D are guaranteed good.
        # How many such guaranteed-good vertices?
        min_guaranteed = n  # track minimum across tournaments

        total_with_free_dom = 0
        all_has_multidom = 0  # has F adj D with |dom(D)| >= 2

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue

            cycles = find_3cycles(A, n)
            free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
            if not free_cycs:
                continue

            dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]

            # Find best (F, D) pair: maximize guaranteed good vertices
            best = 0
            has_multi = False
            for fc in free_cycs:
                fc_set = set(fc)
                for dc, doms in dom_cycs:
                    if not shared_directed_edge(fc, dc):
                        continue
                    total_with_free_dom += 1
                    # Vertices guaranteed good: not in F and not unique dominator of D
                    at_risk = fc_set.copy()
                    if len(doms) == 1:
                        at_risk.add(doms[0])
                    guaranteed = n - len(at_risk)
                    best = max(best, guaranteed)
                    if len(doms) >= 2:
                        has_multi = True

            if best > 0:
                min_guaranteed = min(min_guaranteed, best)
            if has_multi:
                all_has_multidom += 1

        print(f"n={n} ({method}): min guaranteed good vertices = {min_guaranteed}")
        print(f"  has multi-dom adj to free: {all_has_multidom}")

    # LEMMA C: All-dominated case analysis
    print("\n" + "=" * 70)
    print("LEMMA C: All-dominated (no free cycles) case")
    print("=" * 70)

    for n in [4, 5, 6, 7, 8, 9]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        all_dom = 0
        all_dom_has_bad = 0
        all_dom_bad_count = defaultdict(int)
        all_dom_all_bad = 0

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            cycles = find_3cycles(A, n)
            if not cycles:
                continue

            free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
            if free_cycs:
                continue  # skip: has free cycles

            # All cycles are dominated
            if compute_b1(A, n) != 0:
                continue  # b1 should be 0 anyway
            all_dom += 1

            bad_count = 0
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad_count += 1

            all_dom_bad_count[bad_count] += 1
            if bad_count > 0:
                all_dom_has_bad += 1
            if bad_count == n:
                all_dom_all_bad += 1

        print(f"\nn={n} ({method}): {all_dom} all-dominated tournaments")
        print(f"  has bad vertex: {all_dom_has_bad}")
        print(f"  all bad: {all_dom_all_bad}")
        print(f"  bad count dist: {dict(sorted(all_dom_bad_count.items()))}")

    # LEMMA C continued: In all-dominated case, prove good vertex via isolation argument
    # For v to be bad in all-dom case: freed(v) = {cycles with dom={v}} isolated + spanning
    # For any freed cycle C (dom={v}) adj to remaining-dom cycle D (dom != {v}):
    # vertices NOT in C that benefit: need v != d for all d in dom(D), which is auto since dom(D)!={v}
    # Wait: dom(D) could CONTAIN v along with others. Then D remains dominated after removing v.
    # So iso_edges(v) counts edges between freed(v) and remaining-dom(v).
    print("\n" + "=" * 70)
    print("LEMMA C': All-dominated case - isolation argument")
    print("=" * 70)

    for n in [5, 6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        n_all_dom_with_bad = 0
        n_bad_has_adj_multidom = 0  # bad v has freed cycle adj to multi-dom remaining

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            cycles = find_3cycles(A, n)
            if not cycles:
                continue
            free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
            if free_cycs:
                continue

            if compute_b1(A, n) != 0:
                continue

            # Check for bad vertices
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) == 0:
                    continue  # v is good
                n_all_dom_with_bad += 1

                # v is bad. Freed(v) = {C: dom(C)={v}, v not in C}
                freed = [cyc for cyc in cycles if v not in set(cyc) and get_dom(A, n, cyc) == [v]]
                # Remaining dom (not through v): dom(C) has elements other than v
                remaining_dom = []
                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    doms = get_dom(A, n, cyc)
                    doms_nv = [d for d in doms if d != v]
                    if doms_nv:
                        remaining_dom.append(cyc)

                # Any freed cycle adj to any remaining dom?
                adj_found = any(shared_directed_edge(fc, rc)
                               for fc in freed for rc in remaining_dom)

                if adj_found:
                    n_bad_has_adj_multidom += 1

        print(f"n={n} ({method}): bad vertices in all-dom tours: {n_all_dom_with_bad}")
        print(f"  of those, freed adj to remaining-dom: {n_bad_has_adj_multidom}")
        print(f"  => always isolated (no adj): {n_all_dom_with_bad - n_bad_has_adj_multidom}")

    # PART 5: Comprehensive test of the proof at all n
    print("\n" + "=" * 70)
    print("COMPREHENSIVE PROOF VERIFICATION")
    print("=" * 70)

    for n in [4, 5, 6, 7, 8, 9, 10]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000, 10: 1000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        total_b1_0 = 0
        good_exists = 0
        good_missing = 0
        proof_covers = 0
        proof_fails = 0

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue
            total_b1_0 += 1

            cycles = find_3cycles(A, n)

            # Check if good vertex exists (ground truth)
            has_good = False
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) == 0:
                    has_good = True
                    break

            if has_good:
                good_exists += 1
            else:
                good_missing += 1

            # Can the proof find a good vertex?
            if not cycles:
                # No 3-cycles: all vertices good
                proof_covers += 1
                continue

            free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]

            if free_cycs:
                # Free cycles exist. Find F adj D.
                dom_cycs = [(cyc, get_dom(A, n, cyc)) for cyc in cycles if get_dom(A, n, cyc)]
                found_good_v = False
                for fc in free_cycs:
                    for dc, doms in dom_cycs:
                        if not shared_directed_edge(fc, dc):
                            continue
                        # Found F adj D. Guaranteed good: not in F and (|dom|>=2 or not unique dom)
                        fc_set = set(fc)
                        if len(doms) >= 2:
                            # Any v not in F is good
                            good_set = set(range(n)) - fc_set
                        else:
                            # dom = {d}. Any v not in F and v != d is good
                            good_set = set(range(n)) - fc_set - set(doms)
                        if good_set:
                            found_good_v = True
                            break
                    if found_good_v:
                        break

                if found_good_v:
                    proof_covers += 1
                else:
                    proof_fails += 1
            else:
                # All dominated case - need separate argument
                # For now: just check if any v has iso_edges > 0
                found_good_v = False
                for v in range(n):
                    freed = [cyc for cyc in cycles if v not in set(cyc) and get_dom(A, n, cyc) == [v]]
                    rem_dom = [cyc for cyc in cycles if v not in set(cyc) and
                              any(d != v for d in get_dom(A, n, cyc))]
                    iso = any(shared_directed_edge(fc, rc) for fc in freed for rc in rem_dom)
                    if not iso and freed:
                        continue  # v might be bad (freed set isolated)
                    # v has no freed cycles or its freed set is not isolated => v is good
                    if not freed or iso:
                        found_good_v = True
                        break

                if found_good_v:
                    proof_covers += 1
                else:
                    proof_fails += 1

        print(f"\nn={n} ({method}): {total_b1_0} with b1=0")
        print(f"  Good vertex exists: {good_exists}/{total_b1_0}")
        print(f"  Proof covers: {proof_covers}/{total_b1_0}")
        print(f"  Proof FAILS: {proof_fails}/{total_b1_0}")
        if proof_fails > 0:
            print(f"  *** PROOF GAP: {proof_fails} tournaments not covered ***")

if __name__ == '__main__':
    main()
