"""
Algebraic proof attempt: In the all-dominated case, at most 1 vertex is bad.

KEY CLAIM: If all 3-cycles of tournament T are dominated, then at most 1 vertex v
satisfies b_1(T\v) > 0.

APPROACH: Suppose v and w are both bad. Show contradiction.

v bad => freed(v) = {C : dom(C) = {v}} is isolated, connected, spans V\{v}
w bad => freed(w) = {C : dom(C) = {w}} is isolated, connected, spans V\{w}

Since all cycles are dominated, freed(v) consists of cycles uniquely dominated by v.
freed(w) consists of cycles uniquely dominated by w.

Key observations to verify:
1. freed(v) and freed(w) are DISJOINT (different unique dominators)
2. Cycles in freed(w) that don't go through v are remaining-dom for v (dom contains w != v)
3. Since freed(v) spans V\{v} (contains all vertices except v),
   and freed(w) spans V\{w}, their vertex sets overlap on V\{v,w}.
4. Since freed(v) is CONNECTED and spans V\{v}, it covers vertices a,b,c
   that form a cycle in freed(w). This cycle in freed(w) shares edges with
   freed(v) cycles. This creates isolation edges for v => contradiction.
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

    print("=" * 70)
    print("PART 1: Two-bad-vertex structure in all-dominated case")
    print("=" * 70)

    for n in [5, 6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 3000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        two_bad = 0
        one_bad = 0
        zero_bad = 0

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
                continue  # want all-dominated

            if compute_b1(A, n) != 0:
                continue

            bad_verts = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad_verts.append(v)

            if len(bad_verts) >= 2:
                two_bad += 1
            elif len(bad_verts) == 1:
                one_bad += 1
            else:
                zero_bad += 1

        print(f"\nn={n} ({method}): 0 bad: {zero_bad}, 1 bad: {one_bad}, 2+ bad: {two_bad}")

    # PART 2: Deep analysis of what happens when ONE vertex is bad in all-dominated
    print("\n" + "=" * 70)
    print("PART 2: Structure of bad vertex in all-dominated case")
    print("=" * 70)

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = 10000
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        stats = {
            'total': 0,
            'freed_count': [],  # |freed(v)| for bad v
            'freed_span': [],   # |vertex span of freed(v)|
            'unique_dom_score': [],  # score of v in T
            'dom_type': defaultdict(int),  # above/below dominator
            'freed_connected': 0,
            'freed_not_connected': 0,
        }

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

            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) == 0:
                    continue

                stats['total'] += 1

                # freed(v) = cycles with dom = {v}
                freed = [cyc for cyc in cycles if v not in set(cyc) and get_dom(A, n, cyc) == [v]]
                stats['freed_count'].append(len(freed))

                span = set()
                for c in freed:
                    span.update(c)
                stats['freed_span'].append(len(span))

                # Is freed set connected?
                if freed:
                    nr = len(freed)
                    adj_r = [set() for _ in range(nr)]
                    for i in range(nr):
                        for j in range(i+1, nr):
                            if shared_directed_edge(freed[i], freed[j]):
                                adj_r[i].add(j)
                                adj_r[j].add(i)
                    visited = set()
                    stack = [0]
                    while stack:
                        x = stack.pop()
                        if x in visited: continue
                        visited.add(x)
                        for y in adj_r[x]:
                            if y not in visited: stack.append(y)
                    if len(visited) == nr:
                        stats['freed_connected'] += 1
                    else:
                        stats['freed_not_connected'] += 1

                # Dom type: does v beat all cycle vertices (above) or lose to all (below)?
                above_count = 0
                below_count = 0
                for cyc in freed:
                    a, b, c = cyc
                    if A[v][a] and A[v][b] and A[v][c]:
                        above_count += 1
                    elif A[a][v] and A[b][v] and A[c][v]:
                        below_count += 1
                if above_count > 0 and below_count == 0:
                    stats['dom_type']['above'] += 1
                elif below_count > 0 and above_count == 0:
                    stats['dom_type']['below'] += 1
                elif above_count > 0 and below_count > 0:
                    stats['dom_type']['mixed'] += 1
                else:
                    stats['dom_type']['none'] += 1

                # Score of v
                score_v = sum(A[v][j] for j in range(n))
                stats['unique_dom_score'].append(score_v)

        print(f"\nn={n} ({method}): {stats['total']} bad vertices in all-dominated case")
        if stats['total'] > 0:
            fc = stats['freed_count']
            fs = stats['freed_span']
            sc = stats['unique_dom_score']
            print(f"  freed count: min={min(fc)}, max={max(fc)}, avg={sum(fc)/len(fc):.1f}")
            print(f"  freed span: min={min(fs)}, max={max(fs)}, avg={sum(fs)/len(fs):.1f}")
            print(f"  freed connected: {stats['freed_connected']}, not connected: {stats['freed_not_connected']}")
            print(f"  dom type: {dict(stats['dom_type'])}")
            print(f"  score: min={min(sc)}, max={max(sc)}, avg={sum(sc)/len(sc):.1f}")
            print(f"  score dist: {dict(sorted(defaultdict(int, ((s, sc.count(s)) for s in set(sc))).items()))}")

    # PART 3: The KEY structural fact for the algebraic proof
    # When v is bad in all-dom case, v uniquely dominates many cycles.
    # How many OTHER vertices also uniquely dominate some cycles?
    # If freed(v) spans V\{v}, and all cycles in freed(v) have dom={v},
    # then for any OTHER cycle not in freed(v) and not through v,
    # it must have dom with some element != v.
    # Question: how many cycles have dom = {w} for some w != v?
    print("\n" + "=" * 70)
    print("PART 3: Unique dominator distribution when v is bad (all-dom)")
    print("=" * 70)

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = 10000
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        n_bad = 0
        # For each bad v: distribution of |{C not through v : dom(C)={w}}| for w != v
        other_unique_counts = []
        # How many cycles not through v have dom that INCLUDES v (but isn't just v)?
        multi_dom_v_counts = []
        # How many other vertices w have dom(C)={w} for some cycle not through v?
        n_other_unique_doms = []

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

            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) == 0:
                    continue
                n_bad += 1

                # Categorize cycles not through v
                freed_v = []  # dom = {v}
                other_unique = defaultdict(list)  # dom = {w}, w != v
                multi_dom_with_v = []  # v in dom but |dom| >= 2
                multi_dom_no_v = []  # v not in dom, |dom| >= 2

                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    doms = sorted(get_dom(A, n, cyc))
                    if doms == [v]:
                        freed_v.append(cyc)
                    elif len(doms) == 1:
                        other_unique[doms[0]].append(cyc)
                    elif v in doms:
                        multi_dom_with_v.append(cyc)
                    else:
                        multi_dom_no_v.append(cyc)

                other_unique_counts.append(sum(len(cl) for cl in other_unique.values()))
                multi_dom_v_counts.append(len(multi_dom_with_v))
                n_other_unique_doms.append(len(other_unique))

        print(f"\nn={n} ({method}): {n_bad} bad vertices")
        if n_bad > 0:
            ouc = other_unique_counts
            mvc = multi_dom_v_counts
            noud = n_other_unique_doms
            print(f"  Cycles with dom={{w}} (w!=v, |dom|=1): min={min(ouc)}, max={max(ouc)}, avg={sum(ouc)/len(ouc):.1f}")
            print(f"  Cycles with v in dom, |dom|>=2: min={min(mvc)}, max={max(mvc)}, avg={sum(mvc)/len(mvc):.1f}")
            print(f"  #Other unique dominators: min={min(noud)}, max={max(noud)}, avg={sum(noud)/len(noud):.1f}")

    # PART 4: THE KEY — freed(v) cycles adjacent to freed(w) cycles
    # If v and w are both bad, their freed sets must not share isolation edges
    # But freed(w) cycles not through v ARE remaining-dom for v
    # So if any freed(v) cycle is adjacent to a freed(w) cycle not through v => contradiction
    print("\n" + "=" * 70)
    print("PART 4: Hypothetical two-bad analysis — freed-freed adjacency")
    print("=" * 70)

    for n in [5, 6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 3000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        # For tournaments with exactly 1 bad vertex v:
        # Pick another vertex w and compute freed(w) = {C : dom(C)={w}}
        # Is freed(v) ever adjacent to freed(w)?
        n_one_bad = 0
        always_adj = 0  # freed(v) is adj to freed(w) for ALL w != v
        sometimes_adj = 0  # freed(v) adj to freed(w) for SOME w
        never_adj = 0  # freed(v) never adj to any freed(w)

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

            bad_verts = []
            for v_cand in range(n):
                Bv = delete_vertex(A, n, v_cand)
                if compute_b1(Bv, n-1) > 0:
                    bad_verts.append(v_cand)

            if len(bad_verts) != 1:
                continue
            n_one_bad += 1
            v = bad_verts[0]

            freed_v = [cyc for cyc in cycles if v not in set(cyc) and get_dom(A, n, cyc) == [v]]

            # For each other vertex w, compute freed(w) cycles not through v
            adj_count = 0
            total_w = 0
            for w in range(n):
                if w == v:
                    continue
                freed_w = [cyc for cyc in cycles if w not in set(cyc) and v not in set(cyc)
                          and get_dom(A, n, cyc) == [w]]
                if not freed_w:
                    continue
                total_w += 1
                # Is any freed(v) cycle adjacent to any freed(w) cycle (both not through v)?
                adj = any(shared_directed_edge(fv, fw) for fv in freed_v for fw in freed_w)
                if adj:
                    adj_count += 1

            if total_w == 0:
                pass  # no other vertex with unique-dom cycles
            elif adj_count == total_w:
                always_adj += 1
            elif adj_count > 0:
                sometimes_adj += 1
            else:
                never_adj += 1

        print(f"\nn={n} ({method}): {n_one_bad} tournaments with exactly 1 bad vertex")
        print(f"  freed(v) adj ALL freed(w): {always_adj}")
        print(f"  freed(v) adj SOME freed(w): {sometimes_adj}")
        print(f"  freed(v) adj NO freed(w): {never_adj}")

    # PART 5: Spanning + connectivity forcing adjacency
    # KEY QUESTION: If freed(v) spans V\{v} and freed(w) spans V\{w},
    # and both are connected, must they be adjacent?
    # Test: for any two CONNECTED SPANNING subsets of cycle graph, are they always adjacent?
    print("\n" + "=" * 70)
    print("PART 5: Two spanning connected freed sets MUST be adjacent")
    print("=" * 70)

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = 10000
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        n_tested = 0
        n_adj = 0
        n_not_adj = 0

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

            # For each pair (v, w): compute freed(v), freed(w)
            # Check if freed(v) restricted to cycles not through v AND not through w
            # overlaps with freed(w) restricted similarly
            for v in range(n):
                for w in range(v+1, n):
                    # freed(v) cycles not through w
                    fv = [cyc for cyc in cycles if v not in set(cyc) and w not in set(cyc)
                          and get_dom(A, n, cyc) == [v]]
                    # freed(w) cycles not through v
                    fw = [cyc for cyc in cycles if v not in set(cyc) and w not in set(cyc)
                          and get_dom(A, n, cyc) == [w]]

                    if not fv or not fw:
                        continue

                    # Check if fv spans V\{v,w} and fw spans V\{v,w}
                    span_fv = set()
                    for c in fv: span_fv.update(c)
                    span_fw = set()
                    for c in fw: span_fw.update(c)

                    target = set(range(n)) - {v, w}
                    if span_fv != target or span_fw != target:
                        continue

                    # Both span V\{v,w}. Are they adjacent?
                    n_tested += 1
                    adj = any(shared_directed_edge(f1, f2) for f1 in fv for f2 in fw)
                    if adj:
                        n_adj += 1
                    else:
                        n_not_adj += 1
                        if n_not_adj <= 3:
                            print(f"  NOT ADJ at n={n}: v={v}, w={w}")
                            print(f"    freed(v) not through w: {fv}")
                            print(f"    freed(w) not through v: {fw}")

        print(f"n={n} ({method}): {n_tested} pairs tested, {n_adj} adjacent, {n_not_adj} NOT adjacent")

    # PART 6: Simplified — Just prove freed(v) adj to freed(w) in restricted graph for v
    # v bad => freed(v) isolated in restricted(v) graph
    # freed(w) cycles not through v ARE remaining-dom for v (dom contains w != v)
    # So freed(v) must NOT be adjacent to any freed(w) cycle not through v
    # Equivalently: freed(v) and freed(w)\{cycles through v} share no edge
    print("\n" + "=" * 70)
    print("PART 6: freed(v) vs freed(w) adjacency in restricted(v) graph")
    print("=" * 70)
    print("If v is bad, freed(v) ISOLATED from remaining-dom including freed(w)\\v")
    print("So freed(v) NOT adj to freed(w)\\v in restricted(v).")
    print("Question: Can this hold for BOTH v and w simultaneously?")

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = 10000
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        n_all_dom = 0
        n_could_both_bad = 0  # freed(v) not adj freed(w)\v AND freed(w) not adj freed(v)\w

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
            n_all_dom += 1

            for v in range(n):
                for w in range(v+1, n):
                    # freed(v) = cycles not through v with dom={v}
                    freed_v = [cyc for cyc in cycles if v not in set(cyc) and get_dom(A, n, cyc) == [v]]
                    # freed(w) = cycles not through w with dom={w}
                    freed_w = [cyc for cyc in cycles if w not in set(cyc) and get_dom(A, n, cyc) == [w]]

                    if not freed_v or not freed_w:
                        continue

                    # freed(w) restricted to not-through-v = freed(w) cycles not containing v
                    freed_w_no_v = [cyc for cyc in freed_w if v not in set(cyc)]
                    freed_v_no_w = [cyc for cyc in freed_v if w not in set(cyc)]

                    # v bad requires: freed(v) not adj freed(w_no_v)
                    v_condition = not any(shared_directed_edge(fv, fw)
                                        for fv in freed_v for fw in freed_w_no_v)
                    # w bad requires: freed(w) not adj freed(v_no_w)
                    w_condition = not any(shared_directed_edge(fw, fv)
                                        for fw in freed_w for fv in freed_v_no_w)

                    if v_condition and w_condition:
                        # Both conditions met. But we also need spanning+connected.
                        # Check spanning
                        span_v = set()
                        for c in freed_v: span_v.update(c)
                        span_w = set()
                        for c in freed_w: span_w.update(c)

                        target_v = set(range(n)) - {v}
                        target_w = set(range(n)) - {w}

                        if span_v == target_v and span_w == target_w:
                            n_could_both_bad += 1
                            if n_could_both_bad <= 5:
                                print(f"  BOTH COULD BE BAD: n={n}, v={v}, w={w}")
                                print(f"    freed(v)={freed_v}, spans={span_v}")
                                print(f"    freed(w)={freed_w}, spans={span_w}")
                                # Check actual b1 values
                                Bv = delete_vertex(A, n, v)
                                Bw = delete_vertex(A, n, w)
                                print(f"    b1(T\\v)={compute_b1(Bv,n-1)}, b1(T\\w)={compute_b1(Bw,n-1)}")

        print(f"\nn={n} ({method}): {n_all_dom} all-dom tours, "
              f"{n_could_both_bad} with BOTH conditions met and spanning")

    # PART 7: For bad v in all-dom, what fraction of cycle graph edges
    # connect freed(v) to remaining-dom?
    # If freed(v) is LARGE relative to total cycles, adjacency is forced.
    print("\n" + "=" * 70)
    print("PART 7: Freed(v) size vs total cycles when v is bad (all-dom)")
    print("=" * 70)

    for n in [5, 6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 3000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        freed_fracs = []
        total_cyc_counts = []

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

            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) == 0:
                    continue
                # v is bad
                restricted = [cyc for cyc in cycles if v not in set(cyc)]
                freed = [cyc for cyc in restricted if get_dom(A, n, cyc) == [v]]
                if restricted:
                    freed_fracs.append(len(freed) / len(restricted))
                    total_cyc_counts.append(len(restricted))

        if freed_fracs:
            print(f"\nn={n} ({method}): {len(freed_fracs)} bad vertices")
            print(f"  freed fraction: min={min(freed_fracs):.3f}, max={max(freed_fracs):.3f}, "
                  f"avg={sum(freed_fracs)/len(freed_fracs):.3f}")
            print(f"  total restricted cycles: min={min(total_cyc_counts)}, max={max(total_cyc_counts)}, "
                  f"avg={sum(total_cyc_counts)/len(total_cyc_counts):.1f}")

if __name__ == '__main__':
    main()
