"""
ecyc_formula.py — Find formula for e_cyc (edges in 3-cycles) in tournaments.

We proved: dim(Omega_2) = C(n,3) + 2*c3 - e_cyc

e_cyc = #{directed edges participating in at least one 3-cycle}

e_cyc is NOT a function of c3 alone (verified at n=5,6).
Can we express e_cyc in terms of the score sequence?

Key idea: An edge u->v is in a 3-cycle iff there exists w with v->w->u.
That is, N^+(v) ∩ N^-(u) ≠ ∅, where the intersection excludes u,v.
N^+(v) has s_v elements, N^-(u) has n-1-s_u elements.
Their intersection (excluding the edge u->v itself) has size A^2[u,v].

So edge u->v is NOT in any 3-cycle iff A^2[u,v] = 0,
i.e., N^+(u) → v → N^-(v) with no path back through a 3-cycle.

The number of edges NOT in 3-cycles = n(n-1)/2 - e_cyc.
These are called "bridge edges" or "acyclic edges" in the literature.

For a transitive tournament: e_cyc = 0 (no 3-cycles).
For a regular tournament (odd n): e_cyc = n(n-1)/2 (all edges in 3-cycles).

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

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

def count_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            c3 += 1
    return c3

def edges_in_3cycles(A, n):
    in_cycle = set()
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            in_cycle.add((i,j)); in_cycle.add((j,k)); in_cycle.add((k,i))
        elif A[i][k] and A[k][j] and A[j][i]:
            in_cycle.add((i,k)); in_cycle.add((k,j)); in_cycle.add((j,i))
    return len(in_cycle)

def edges_not_in_3cycles(A, n):
    """Count edges (directed arcs) NOT in any 3-cycle.
    An arc u->v is acyclic iff N^+(v) ∩ N^-(u) \ {u,v} = ∅,
    i.e., A^2[u,v] = 0."""
    A2 = A @ A
    count = 0
    for u in range(n):
        for v in range(n):
            if A[u][v] == 1 and A2[u][v] == 0:
                count += 1
    return count


def main():
    print("=" * 70)
    print("E_CYC FORMULA INVESTIGATION")
    print("=" * 70)

    # Part 1: e_cyc vs score sequence
    print("\n--- Part 1: e_cyc as function of score sequence ---")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        score_ecyc = defaultdict(list)
        score_c3 = defaultdict(list)

        for bits in range(N):
            A = bits_to_adj(bits, n)
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            c3 = count_3cycles(A, n)
            e_cyc = edges_in_3cycles(A, n)
            score_ecyc[scores].append(e_cyc)
            score_c3[scores].append(c3)

        print(f"\n  n={n}: e_cyc by score sequence:")
        for sc in sorted(score_ecyc.keys()):
            vals = score_ecyc[sc]
            c3_vals = score_c3[sc]
            unique_ecyc = sorted(set(vals))
            unique_c3 = sorted(set(c3_vals))
            is_const = len(unique_ecyc) == 1
            print(f"    scores={sc}: c3={unique_c3}, e_cyc={unique_ecyc}, "
                  f"{'CONSTANT' if is_const else 'VARIES'} ({len(vals)} tournaments)")

    # Part 2: Formula e_cyc = n(n-1)/2 - acyclic_edges
    # Acyclic edge u->v means: every w beaten by v also beats u (no "back path").
    # For u->v acyclic: s_v vertices beaten by v, and ALL of them beat u.
    # So u has in-degree >= s_v from N^+(v). That means s_u <= n-1-s_v.
    # More precisely: u->v is acyclic iff A^2[u,v] = 0.
    print("\n--- Part 2: Acyclic edge characterization ---")
    n = 5
    N = 2**(n*(n-1)//2)

    for bits in range(min(N, 20)):
        A = bits_to_adj(bits, n)
        scores = [int(sum(A[i])) for i in range(n)]
        A2 = A @ A
        acyclic = edges_not_in_3cycles(A, n)
        e_cyc = edges_in_3cycles(A, n)

        if bits < 5:
            print(f"\n  bits={bits}: scores={scores}")
            for u in range(n):
                for v in range(n):
                    if A[u][v]:
                        status = "acyclic" if A2[u][v] == 0 else f"cyclic(A^2={A2[u][v]})"
                        print(f"    {u}->{v}: {status}")
            print(f"    e_cyc={e_cyc}, acyclic={acyclic}, total arcs={n*(n-1)//2}")
            print(f"    CHECK: e_cyc + acyclic = {e_cyc + acyclic} = {n*(n-1)//2}? "
                  f"{'YES' if e_cyc + acyclic == n*(n-1)//2 else 'NO'}")

    # Part 3: Can we express acyclic edges in terms of score sequence?
    # An edge u->v is acyclic iff N^+(v) ⊆ N^+(u) ∪ {u}
    # (everyone v beats is also beaten by u, except possibly u itself)
    # This is a DOMINANCE condition: u "covers" all of v's out-neighborhood
    print("\n--- Part 3: Acyclic = dominance edges ---")
    n = 5
    N = 2**(n*(n-1)//2)

    # For transitive tournament: all edges are acyclic (dominance = total order)
    # For regular: no edges are acyclic
    # General: acyclic iff local dominance

    # Let's count acyclic edges as function of (s_u, s_v) pair
    uv_pairs = defaultdict(list)  # (s_u, s_v) -> is_acyclic counts
    for bits in range(N):
        A = bits_to_adj(bits, n)
        scores = [int(sum(A[i])) for i in range(n)]
        A2 = A @ A
        for u in range(n):
            for v in range(n):
                if A[u][v]:
                    is_acyc = 1 if A2[u][v] == 0 else 0
                    uv_pairs[(scores[u], scores[v])].append(is_acyc)

    print(f"\n  n={n}: P(acyclic | s_u, s_v) for edge u->v:")
    for key in sorted(uv_pairs.keys()):
        vals = uv_pairs[key]
        prob = np.mean(vals)
        print(f"    (s_u={key[0]}, s_v={key[1]}): P(acyclic) = {prob:.4f} ({sum(vals)}/{len(vals)})")

    # Part 4: sum A^2[u,v] over edges u->v equals 3*c3 (triple counting)
    print("\n--- Part 4: sum A^2[u,v] over edges = 3*c3 ---")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        violations = 0
        for bits in range(N):
            A = bits_to_adj(bits, n)
            c3 = count_3cycles(A, n)
            A2 = A @ A
            total_A2 = sum(A2[u][v] for u in range(n) for v in range(n) if A[u][v])
            if total_A2 != 3 * c3:
                violations += 1
                if violations <= 3:
                    print(f"  FAIL bits={bits}: sum(A^2)={total_A2}, 3*c3={3*c3}")
        if violations == 0:
            print(f"  n={n}: CONFIRMED sum(A^2[u,v] for u->v) = 3*c3 for ALL tournaments")

    # Part 5: Edge multiplicity in 3-cycles
    # Each 3-cycle contributes 3 edges. If an edge is in k 3-cycles,
    # it's counted k times in the sum of A^2 values.
    # So: sum(A^2[u,v] for u->v) = sum_e (number of 3-cycles through e) = 3*c3.
    # And: e_cyc = #{e : A^2[u,v] > 0} = #{edges in at least one 3-cycle}
    # We want: e_cyc = f(score sequence, ?)

    # Let's try: e_cyc = n(n-1)/2 - sum_v max(0, s_v*(s_v-1)/2 - something) ...
    # Actually let's try a formula based on sum s_v^2 (a classic tournament invariant)
    print("\n--- Part 5: e_cyc vs sum(s_v^2) ---")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        data = defaultdict(list)

        for bits in range(N):
            A = bits_to_adj(bits, n)
            scores = [int(sum(A[i])) for i in range(n)]
            ss2 = sum(s*s for s in scores)
            e_cyc = edges_in_3cycles(A, n)
            # c3 = (sum s_i(s_i-1)/2 - sum_{i<j} A[i][j]*A[j][i]*(n-2))/...
            # Actually c3 = (n(n-1)(2n-1)/12 - sum s_i^2 / 2) ... no
            # c3 = C(n,3) - sum C(s_i, 2) + ...
            # Known: c3 = (C(n,2)^2 - sum s_i^2) / ... no
            # c3 = (sum s_i*(s_i - 1)) / 2 counts the number of 2-paths
            # and c3 = (n choose 3) - #{transitive triples}
            # #{transitive triples} = sum C(s_i, 2)
            # So c3 = C(n,3) - sum C(s_i, 2)
            c3 = comb(n,3) - sum(comb(s,2) for s in scores)
            data[(ss2, c3)].append(e_cyc)

        print(f"\n  n={n}: e_cyc by (sum_s^2, c3):")
        for key in sorted(data.keys()):
            vals = data[key]
            unique = sorted(set(vals))
            ss2, c3 = key
            is_const = len(unique) == 1
            print(f"    ss2={ss2}, c3={c3}: e_cyc={unique}, {'CONST' if is_const else 'VARIES'}")

    # Part 6: Does e_cyc depend on MORE than just c3?
    # At n=6, c3=2 has e_cyc in {5,6}. What distinguishes them?
    print("\n--- Part 6: What distinguishes e_cyc values within same c3? ---")
    n = 6
    N = 2**(n*(n-1)//2)

    c3_2_groups = defaultdict(list)
    for bits in range(N):
        A = bits_to_adj(bits, n)
        c3 = count_3cycles(A, n)
        if c3 != 2: continue
        e_cyc = edges_in_3cycles(A, n)
        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        c3_2_groups[e_cyc].append((bits, scores))

    for ec in sorted(c3_2_groups.keys()):
        group = c3_2_groups[ec]
        score_dist = defaultdict(int)
        for bits, scores in group:
            score_dist[scores] += 1
        print(f"\n  c3=2, e_cyc={ec}: {len(group)} tournaments")
        for sc, cnt in sorted(score_dist.items()):
            print(f"    scores={sc}: {cnt}")

        # Check: do the two 3-cycles SHARE an edge?
        share_count = 0
        for bits, scores in group[:5]:
            A = bits_to_adj(bits, n)
            cycles = []
            for i,j,k in combinations(range(n), 3):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i,j,k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append((i,k,j))
            if len(cycles) == 2:
                c1_edges = {(cycles[0][i], cycles[0][(i+1)%3]) for i in range(3)}
                c2_edges = {(cycles[1][i], cycles[1][(i+1)%3]) for i in range(3)}
                shared = c1_edges & c2_edges
                if shared:
                    share_count += 1
                    if len(c3_2_groups[ec]) <= 100:
                        print(f"      bits={bits}: cycles={cycles}, shared edges={shared}")
        if share_count > 0:
            print(f"    {share_count}/5 samples share an edge between the 2 cycles")

    print("\n--- Part 7: e_cyc = 3*c3 - 3*(shared edges between cycle pairs) ---")
    # If cycles are edge-disjoint: e_cyc = 3*c3 (each cycle adds 3 new edges)
    # But can exceed n(n-1)/2, so there must be sharing for large c3.
    # Let overlap = #{edge pairs shared between distinct 3-cycles}
    # Then e_cyc = min(n(n-1)/2, something involving overlap)

    # Actually: e_cyc = #{edges in the UNION of all 3-cycles}
    # By inclusion-exclusion on edge-level, this is complex.
    # But: e_cyc <= 3*c3 (since each cycle contributes at most 3 edges)
    # And: e_cyc <= n(n-1)/2 (total arcs)
    # Let's check: is e_cyc = min(3*c3, n(n-1)/2)?
    print("\n  Checking e_cyc vs min(3*c3, total_arcs):")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        total_arcs = n*(n-1)//2
        match = 0
        total = 0
        counter_examples = []

        for bits in range(N):
            A = bits_to_adj(bits, n)
            c3 = count_3cycles(A, n)
            e_cyc = edges_in_3cycles(A, n)
            candidate = min(3*c3, total_arcs)
            if e_cyc == candidate:
                match += 1
            else:
                counter_examples.append((bits, c3, e_cyc, candidate))
            total += 1

        print(f"\n  n={n}: e_cyc = min(3*c3, {total_arcs}) in {match}/{total} "
              f"({100*match/total:.1f}%)")
        if counter_examples:
            for bits, c3, ec, cand in counter_examples[:5]:
                print(f"    FAIL bits={bits}: c3={c3}, e_cyc={ec}, min(3*c3,arcs)={cand}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
