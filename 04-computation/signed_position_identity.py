"""
Signed Position Identity interpretation of Even-Odd Split Lemma.

Rewriting the alternating sum:
  sum_S (-1)^|S| a(S)*b(M\S) = sum_P (-1)^{pos(i)-1}

where the sum is over all Ham paths P containing i->j as a consecutive pair,
and pos(i) is the 0-indexed position of i in P.

So the Even-Odd Split Lemma says:
  sum_{P: i->j in P} (-1)^{pos_P(i)} = sum_{P': j->i in P'} (-1)^{pos_{P'}(j)}

where pos gives the position (0-indexed) of the first vertex in the ij edge.

This is a signed counting identity relating adjacencies in T and T'.

Instance: opus-2026-03-05-S4
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import all_tournaments
from itertools import permutations


def verify_signed_position(n, max_tournaments=None):
    print(f"=== Signed Position Identity at n={n} ===")
    total = 0
    passes = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if max_tournaments and t_idx >= max_tournaments:
            break
        for i in range(n):
            for j in range(n):
                if i == j or T[i][j] != 1:
                    continue
                total += 1

                # LHS: sum over Ham paths P with i->j consecutive
                # Weight: (-1)^{position of i in P}
                lhs = 0
                for perm in permutations(range(n)):
                    valid = True
                    for k in range(n-1):
                        if T[perm[k]][perm[k+1]] != 1:
                            valid = False
                            break
                    if not valid:
                        continue
                    for k in range(n-1):
                        if perm[k] == i and perm[k+1] == j:
                            lhs += (-1)**k

                # RHS: sum over Ham paths P' with j->i consecutive (in T')
                # But T' only differs from T in the i<->j arc, and
                # the paths with j->i consecutive don't use i->j, so
                # they're valid in T iff they're valid in T' (arcs among others same).
                # Wait -- paths with j->i use j->i which is NOT in T (T has i->j).
                # So we need T' for these paths.
                T2 = [row[:] for row in T]
                T2[i][j] = 0
                T2[j][i] = 1

                rhs = 0
                for perm in permutations(range(n)):
                    valid = True
                    for k in range(n-1):
                        if T2[perm[k]][perm[k+1]] != 1:
                            valid = False
                            break
                    if not valid:
                        continue
                    for k in range(n-1):
                        if perm[k] == j and perm[k+1] == i:
                            rhs += (-1)**k

                if lhs == rhs:
                    passes += 1
                else:
                    print(f"  FAIL: T#{t_idx} {i}->{j}: LHS={lhs}, RHS={rhs}")
                    if total - passes > 5:
                        return

    print(f"Result: {passes}/{total}")

    # Also display some example values
    print("\nExample signed adjacency values:")
    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 5:
            break
        for i in range(min(n, 3)):
            for j in range(min(n, 3)):
                if i == j or T[i][j] != 1:
                    continue
                val = 0
                for perm in permutations(range(n)):
                    valid = True
                    for k in range(n-1):
                        if T[perm[k]][perm[k+1]] != 1:
                            valid = False
                            break
                    if not valid:
                        continue
                    for k in range(n-1):
                        if perm[k] == i and perm[k+1] == j:
                            val += (-1)**k
                print(f"  T#{t_idx} adj_signed({i},{j}) = {val}", end="")

                # Also show unsigned adj for comparison
                adj = sum(1 for perm in permutations(range(n))
                         if all(T[perm[k]][perm[k+1]] == 1 for k in range(n-1))
                         and any(perm[k]==i and perm[k+1]==j for k in range(n-1)))
                print(f"  (unsigned adj = {adj})")
                break
            break


def check_stronger_identity(n, max_tournaments=10):
    """
    Check if an even stronger identity holds:
    Does the signed adj decompose by position?

    sum_{P: i->j at pos k} 1 on T-side vs T'-side for each k?
    """
    print(f"\n=== Per-position decomposition at n={n} ===")
    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= max_tournaments:
            break

        i, j = 0, 1
        if T[i][j] != 1:
            i, j = 1, 0
            if T[i][j] != 1:
                continue

        T2 = [row[:] for row in T]
        T2[i][j] = 0
        T2[j][i] = 1

        by_pos_T = [0] * (n-1)
        by_pos_T2 = [0] * (n-1)

        for perm in permutations(range(n)):
            # Check T
            valid_T = True
            for k in range(n-1):
                if T[perm[k]][perm[k+1]] != 1:
                    valid_T = False
                    break
            if valid_T:
                for k in range(n-1):
                    if perm[k] == i and perm[k+1] == j:
                        by_pos_T[k] += 1

            # Check T'
            valid_T2 = True
            for k in range(n-1):
                if T2[perm[k]][perm[k+1]] != 1:
                    valid_T2 = False
                    break
            if valid_T2:
                for k in range(n-1):
                    if perm[k] == j and perm[k+1] == i:
                        by_pos_T2[k] += 1

        print(f"  T#{t_idx} ({i}->{j}):")
        print(f"    adj(i,j) by pos: {by_pos_T}  total={sum(by_pos_T)}")
        print(f"    adj'(j,i) by pos: {by_pos_T2}  total={sum(by_pos_T2)}")

        signed_T = sum((-1)**k * by_pos_T[k] for k in range(n-1))
        signed_T2 = sum((-1)**k * by_pos_T2[k] for k in range(n-1))
        print(f"    signed: T={signed_T}, T'={signed_T2}, equal={signed_T==signed_T2}")

        # Check: does per-position match? (probably not)
        per_pos_match = all(by_pos_T[k] == by_pos_T2[k] for k in range(n-1))
        print(f"    per-position match: {per_pos_match}")


if __name__ == "__main__":
    for n in [4, 5]:
        verify_signed_position(n, max_tournaments=None if n <= 4 else 1024)

    for n in [4, 5, 6]:
        check_stronger_identity(n, max_tournaments=5)
