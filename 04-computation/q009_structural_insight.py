"""
Q-009: Verify key structural insight about destroyed/created cycles.

Claim: A cycle C through v is destroyed by flipping i->j iff i->j is an arc
of C, which requires {i,j} ⊆ V(C). Similarly for created cycles.

If true, the complement V\V(C) of destroyed/created cycles never contains
i or j, so H(T[complement]) = H(T'[complement]) — the mu-weight of
destroyed/created cycles is the SAME in T and T'.

This simplifies the sum-equality: destroyed cycles contribute +mu(C) and
created cycles contribute -mu(C'), where mu values depend only on the
complement sub-tournament (unchanged by the flip).

Combined with the complement adjacency formula for persist-change cycles,
this gives a complete structural picture of delta(sum mu).

Author: opus-2026-03-05-S2
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, mu
)


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


if __name__ == "__main__":
    for n in [5, 6]:
        print(f"\n{'='*60}")
        print(f"n={n}: Do all destroyed/created cycles contain {{i,j}}?")
        print(f"{'='*60}")

        destroyed_with_ij = 0
        destroyed_without_ij = 0
        created_with_ij = 0
        created_without_ij = 0
        total_flips = 0

        limit = 1024 if n == 5 else 300
        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= limit:
                break
            cycles_T = find_odd_cycles(T)
            for v in range(n):
                cycles_T_v = [c for c in cycles_T if v in set(c)]
                set_T_v = set(tuple(c) for c in cycles_T_v)
                others = [u for u in range(n) if u != v]
                for idx_i in range(len(others)):
                    i = others[idx_i]
                    for idx_j in range(idx_i + 1, len(others)):
                        j = others[idx_j]
                        if T[i][j] == 0:
                            continue
                        total_flips += 1

                        T2 = flip_arc(T, i, j)
                        cycles_T2 = find_odd_cycles(T2)
                        cycles_T2_v = [c for c in cycles_T2 if v in set(c)]
                        set_T2_v = set(tuple(c) for c in cycles_T2_v)

                        destroyed = set_T_v - set_T2_v
                        created = set_T2_v - set_T_v

                        for c in destroyed:
                            vset = set(c)
                            if i in vset and j in vset:
                                destroyed_with_ij += 1
                            else:
                                destroyed_without_ij += 1
                                if destroyed_without_ij <= 3:
                                    print(f"  UNEXPECTED: destroyed C={c} without "
                                          f"{{i={i},j={j}}} ⊆ V(C)={vset}")
                                    # Check if it uses arc i->j
                                    cl = list(c)
                                    uses = any(cl[k]==i and cl[(k+1)%len(cl)]==j
                                               for k in range(len(cl)))
                                    print(f"    Uses arc i->j? {uses}")

                        for c in created:
                            vset = set(c)
                            if i in vset and j in vset:
                                created_with_ij += 1
                            else:
                                created_without_ij += 1
                                if created_without_ij <= 3:
                                    print(f"  UNEXPECTED: created C={c} without "
                                          f"{{i={i},j={j}}} ⊆ V(C)={vset}")

        print(f"\n  Total flips: {total_flips}")
        print(f"  Destroyed with {{i,j}} ⊆ V(C): {destroyed_with_ij}")
        print(f"  Destroyed WITHOUT {{i,j}} ⊆ V(C): {destroyed_without_ij}")
        print(f"  Created with {{i,j}} ⊆ V(C): {created_with_ij}")
        print(f"  Created WITHOUT {{i,j}} ⊆ V(C): {created_without_ij}")

        if destroyed_without_ij == 0 and created_without_ij == 0:
            print(f"\n  *** CONFIRMED: All destroyed/created cycles contain {{i,j}} ***")
            print(f"  This means complement V\\V(C) never contains i or j,")
            print(f"  so mu(C) = H(T[V\\V(C)]) = H(T'[V\\V(C)]).")
        else:
            print(f"\n  INSIGHT FAILS: some destroyed/created cycles don't contain {{i,j}}")

    # At n=6, further decompose: how many destroyed/created are 3-cycles vs 5-cycles?
    print(f"\n{'='*60}")
    print(f"n=6: Destroyed/created cycle length distribution")
    print(f"{'='*60}")

    from collections import Counter
    destroyed_len = Counter()
    created_len = Counter()

    for t_idx, T in enumerate(all_tournaments(6)):
        if t_idx >= 300:
            break
        cycles_T = find_odd_cycles(T)
        for v in range(6):
            cycles_T_v = [c for c in cycles_T if v in set(c)]
            set_T_v = set(tuple(c) for c in cycles_T_v)
            others = [u for u in range(6) if u != v]
            for idx_i in range(len(others)):
                i = others[idx_i]
                for idx_j in range(idx_i + 1, len(others)):
                    j = others[idx_j]
                    if T[i][j] == 0:
                        continue
                    T2 = flip_arc(T, i, j)
                    cycles_T2 = find_odd_cycles(T2)
                    cycles_T2_v = [c for c in cycles_T2 if v in set(c)]
                    set_T2_v = set(tuple(c) for c in cycles_T2_v)

                    for c in (set_T_v - set_T2_v):
                        destroyed_len[len(c)] += 1
                    for c in (set_T2_v - set_T_v):
                        created_len[len(c)] += 1

    print(f"  Destroyed by length: {dict(sorted(destroyed_len.items()))}")
    print(f"  Created by length: {dict(sorted(created_len.items()))}")

    # For n=6: 3-cycle complement has 3 verts (mu ∈ {1,3})
    # 5-cycle complement has 1 vert (mu = 1 always)
    print(f"\n  At n=6:")
    print(f"  - 3-cycle: complement has 3 vertices, mu ∈ {{1, 3}}")
    print(f"  - 5-cycle: complement has 1 vertex, mu = 1 always")

    # Verify: for destroyed/created 5-cycles, mu always = 1
    print(f"\n  Verifying mu=1 for all destroyed/created 5-cycles...")
    mu_not_1 = 0
    for t_idx, T in enumerate(all_tournaments(6)):
        if t_idx >= 100:
            break
        cycles_T = find_odd_cycles(T)
        for v in range(6):
            cycles_T_v = [c for c in cycles_T if v in set(c)]
            set_T_v = set(tuple(c) for c in cycles_T_v)
            Tv, ol = delete_vertex(T, v)
            tc = find_odd_cycles(Tv)
            cache = (Tv, ol, tc)

            others = [u for u in range(6) if u != v]
            for idx_i in range(len(others)):
                i = others[idx_i]
                for idx_j in range(idx_i + 1, len(others)):
                    j = others[idx_j]
                    if T[i][j] == 0:
                        continue
                    T2 = flip_arc(T, i, j)
                    cycles_T2 = find_odd_cycles(T2)
                    cycles_T2_v = [c for c in cycles_T2 if v in set(c)]
                    set_T2_v = set(tuple(c) for c in cycles_T2_v)

                    for c in (set_T_v - set_T2_v):
                        if len(c) == 5:
                            m = mu(T, v, c, cache)
                            if m != 1:
                                mu_not_1 += 1
    print(f"  5-cycles with mu != 1: {mu_not_1}")
