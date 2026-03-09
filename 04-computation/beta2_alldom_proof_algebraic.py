"""
ALGEBRAIC PROOF of at-most-1-bad in all-dominated case.

THEOREM: If all 3-cycles of T are dominated, at most 1 vertex v has b1(T\v) > 0.

PROOF OUTLINE:
Step 1 (Extreme Score Lemma): If v is bad, score(v) in {0, n-1}.
  - freed(v) must be connected and span V\{v} (by isolation characterization)
  - Let W+ = {u : v->u}, W- = {u : u->v}. Partition V\{v} = W+ u W-.
  - v dominates C from above iff V(C) c W+, from below iff V(C) c W-.
  - Cycles in W+ not adjacent to cycles in W- (disjoint vertex sets => no shared edges).
  - freed(v) connected => all in W+ or all in W-.
  - freed(v) spans V\{v} = W+ u W- => W- = empty or W+ = empty.
  - So score(v) = n-1 or score(v) = 0.

Step 2: score(v) in {0, n-1} => v is in no 3-cycle.
  - score 0: v loses to all, so v->x never, no cycle uses outgoing edge from v.
  - score n-1: v beats all, so x->v never, no cycle uses incoming edge to v.

Step 3: v in no 3-cycle => no other w can be bad.
  - freed(w) cycles don't go through w. They use vertices from V\{w}.
  - v is in no cycle, so no cycle covers vertex v.
  - freed(w) spans at most V\{v,w}, which has n-2 < n-1 vertices.
  - So freed(w) cannot span V\{w} (which needs to cover v).
  - By isolation characterization, w is not bad.

This script verifies each step computationally.
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

    # STEP 1 VERIFICATION: Bad vertex always has extreme score
    print("=" * 70)
    print("STEP 1: Bad vertex in all-dom case has score 0 or n-1")
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

        extreme = 0
        non_extreme = 0

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
                score = sum(A[v][j] for j in range(n))
                if score == 0 or score == n-1:
                    extreme += 1
                else:
                    non_extreme += 1
                    if non_extreme <= 3:
                        print(f"  NON-EXTREME: n={n}, v={v}, score={score}")

        print(f"n={n} ({method}): extreme={extreme}, non-extreme={non_extreme}")

    # STEP 2 VERIFICATION: Extreme-score vertex is in no 3-cycle
    print("\n" + "=" * 70)
    print("STEP 2: Score 0 or n-1 => not in any 3-cycle (trivial)")
    print("=" * 70)
    print("This is trivially true:")
    print("  score 0: v loses to all, so v->x never happens => no cycle through v")
    print("  score n-1: v beats all, so x->v never happens => no cycle through v")
    print("Verification:")

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
        else:
            sample = range(5000)

        violations = 0
        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            cycles = find_3cycles(A, n)
            for v in range(n):
                score = sum(A[v][j] for j in range(n))
                if score not in (0, n-1):
                    continue
                in_cycle = any(v in set(c) for c in cycles)
                if in_cycle:
                    violations += 1

        print(f"  n={n}: violations = {violations}")

    # STEP 3 VERIFICATION: v not in any cycle => no other vertex w can have
    # freed(w) spanning V\{w} (since v is uncoverable by cycles)
    print("\n" + "=" * 70)
    print("STEP 3: v in no cycle => freed(w) can't span V\\{w} for w != v")
    print("=" * 70)
    print("If v is in no 3-cycle, vertex v cannot appear in any cycle.")
    print("freed(w) uses only cycle vertices, so can't include v.")
    print("Thus freed(w) spans at most V\\{v,w}, which has n-2 < n-1 vertices.")

    for n in [5, 6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 3000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        n_tested = 0
        max_span_when_v_absent = []

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

            # Find vertex v not in any cycle (score 0 or n-1)
            for v in range(n):
                score = sum(A[v][j] for j in range(n))
                if score not in (0, n-1):
                    continue
                # v is not in any cycle. Check max span of freed(w) for w != v
                for w in range(n):
                    if w == v:
                        continue
                    # freed(w) = cycles not through w with dom = {w}
                    freed_w = [cyc for cyc in cycles if w not in set(cyc)
                              and get_dom(A, n, cyc) == [w]]
                    if freed_w:
                        span = set()
                        for c in freed_w:
                            span.update(c)
                        n_tested += 1
                        max_span_when_v_absent.append(len(span))
                        if v in span:
                            print(f"  ERROR: v={v} (score {score}) IN span of freed(w={w})!")

        if max_span_when_v_absent:
            print(f"n={n} ({method}): {n_tested} pairs tested")
            print(f"  max span of freed(w): max={max(max_span_when_v_absent)}, "
                  f"need n-1={n-1} but get at most n-2={n-2}")
            above_threshold = sum(1 for s in max_span_when_v_absent if s >= n-1)
            print(f"  span >= n-1: {above_threshold} (should be 0)")
        else:
            print(f"n={n} ({method}): no extreme-score vertices with freed sets")

    # FULL PROOF VERIFICATION: Combine all steps
    print("\n" + "=" * 70)
    print("FULL ALGEBRAIC PROOF VERIFICATION")
    print("=" * 70)
    print("Claim: all-dominated + b1=0 => at most 1 bad vertex")
    print("Proof: bad => extreme score => not in cycle => other can't be bad")

    for n in [4, 5, 6, 7, 8, 9, 10]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000, 10: 1000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        total_alldom = 0
        proof_holds = 0  # at most 1 bad
        proof_fails = 0  # 2+ bad
        step1_verified = 0
        step1_fails = 0

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
            total_alldom += 1

            bad = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.append(v)

            if len(bad) <= 1:
                proof_holds += 1
            else:
                proof_fails += 1

            # Verify step 1 for bad vertices
            for v in bad:
                score = sum(A[v][j] for j in range(n))
                if score in (0, n-1):
                    step1_verified += 1
                else:
                    step1_fails += 1

        print(f"\nn={n} ({method}): {total_alldom} all-dom tournaments")
        print(f"  at most 1 bad: {proof_holds}, 2+ bad: {proof_fails}")
        print(f"  Step 1 (extreme score): verified={step1_verified}, fails={step1_fails}")

    # BONUS: Show that the Extreme Score Lemma follows from the proof
    # The algebraic argument is:
    # freed(v) connected + spanning => all freed in W+ or all in W- => W- or W+ empty
    print("\n" + "=" * 70)
    print("BONUS: Why freed(v) connected + spanning => extreme score")
    print("=" * 70)
    print("freed(v) has cycles in W+ = {u:v->u} and W- = {u:u->v}")
    print("Cycles in W+ share no vertices with cycles in W- (partition)")
    print("No shared vertices => no shared edges => not adjacent")
    print("freed(v) connected => all in one side")
    print("freed(v) spanning V\\{v} = W+ u W- => other side empty")
    print("=> score = 0 or n-1")
    print()

    # Verify the partition argument directly
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
        all_in_wplus = 0
        all_in_wminus = 0
        mixed = 0

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

                W_plus = set(j for j in range(n) if j != v and A[v][j])
                W_minus = set(j for j in range(n) if j != v and A[j][v])

                freed = [cyc for cyc in cycles if v not in set(cyc) and get_dom(A, n, cyc) == [v]]

                in_plus = 0
                in_minus = 0
                for c in freed:
                    if set(c) <= W_plus:
                        in_plus += 1
                    elif set(c) <= W_minus:
                        in_minus += 1
                    else:
                        mixed += 1

                if in_plus > 0 and in_minus == 0:
                    all_in_wplus += 1
                elif in_minus > 0 and in_plus == 0:
                    all_in_wminus += 1

        print(f"n={n} ({method}): {n_bad} bad vertices")
        print(f"  all freed in W+: {all_in_wplus}, all in W-: {all_in_wminus}, mixed: {mixed}")

if __name__ == '__main__':
    main()
