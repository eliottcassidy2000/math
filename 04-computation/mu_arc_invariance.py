"""
Key structural observation: mu(C) is invariant under arc flips not involving
any vertex of C.

If C is a cycle through v and we flip arc i->j to j->i where i,j not in V(C),
then mu(C) = mu'(C) because:
1. T-v changes to T'-v (i<->j flipped in T-v)
2. C\{v} is unchanged
3. The "available" vertices for Omega(T-v)|_{avoid C\{v}} include i and j
4. The odd cycles in the available set might change...

Wait, but if i and j are NOT in C, they ARE in the available set.
The flip changes the tournament on the available set, which changes
the odd cycles there, which changes mu.

But the computation shows mu doesn't change at n=5!

Hypothesis: at n=5, the available set is too small for the flip to matter.
For 3-cycles at n=5: available has 2 vertices (< 3 for cycles), so mu=1 trivially.
For 5-cycles: available has 0 vertices, mu=1 trivially.

So at n=5, mu invariance is trivial (all mu=1 by THM-008).
The real test is n=6, where mu can be 1 or 3.

If i,j are both in the available set AND the available set has 3 vertices,
flipping i<->j could change whether those 3 vertices form a 3-cycle,
changing mu from 1 to 3 or vice versa.

Let me test at n=6.

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, find_odd_cycles, delete_vertex, mu
)


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


if __name__ == "__main__":
    print("="*70)
    print("Does mu(C) change under arc flips not in V(C)?")
    print("="*70)

    for n in [5, 6]:
        print(f"\n--- n={n} ---")
        mu_changed = 0
        mu_same = 0
        mu_changed_details = []

        limit = 200 if n == 6 else 1024
        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= limit:
                break
            all_cycles = find_odd_cycles(T)
            for v in range(n):
                cycles_v = [c for c in all_cycles if v in set(c)]
                others = [u for u in range(n) if u != v]
                for i in others:
                    for j in others:
                        if i == j or T[i][j] == 0:
                            continue
                        T2 = flip_arc(T, i, j)
                        for c in cycles_v:
                            vset = set(c)
                            # Only look at cycles NOT using the arc
                            uses_arc = any(
                                c[k] == i and c[(k+1) % len(c)] == j
                                for k in range(len(c))
                            )
                            if uses_arc:
                                continue
                            # Also skip if i or j is in V(C) but arc not used
                            # We want: i,j NOT in V(C) at all
                            if i in vset or j in vset:
                                continue
                            # Verify cycle still valid in T'
                            still_valid = all(
                                T2[c[k]][c[(k+1) % len(c)]] == 1
                                for k in range(len(c))
                            )
                            if not still_valid:
                                continue

                            Tv, ol = delete_vertex(T, v)
                            tc = find_odd_cycles(Tv)
                            mu_T = mu(T, v, c, (Tv, ol, tc))

                            Tv2, ol2 = delete_vertex(T2, v)
                            tc2 = find_odd_cycles(Tv2)
                            mu_T2 = mu(T2, v, c, (Tv2, ol2, tc2))

                            if mu_T != mu_T2:
                                mu_changed += 1
                                if len(mu_changed_details) < 5:
                                    mu_changed_details.append(
                                        (t_idx, v, c, i, j, mu_T, mu_T2))
                            else:
                                mu_same += 1

        print(f"  i,j NOT in V(C): mu changed={mu_changed}, same={mu_same}")
        if mu_changed > 0:
            print(f"  mu DOES change when i,j not in V(C)!")
            for d in mu_changed_details:
                print(f"    T#{d[0]}, v={d[1]}, C={d[2]}, flip {d[3]}->{d[4]}: "
                      f"mu {d[5]} -> {d[6]}")
        else:
            print(f"  mu does NOT change when i,j not in V(C)")

        # Also check: i in V(C) but j not in V(C) (or vice versa)
        print(f"\n  Now checking: i in V(C), j NOT in V(C) (one endpoint in cycle):")
        mu_changed_partial = 0
        mu_same_partial = 0
        details_partial = []

        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= limit:
                break
            all_cycles = find_odd_cycles(T)
            for v in range(n):
                cycles_v = [c for c in all_cycles if v in set(c)]
                others = [u for u in range(n) if u != v]
                for i in others:
                    for j in others:
                        if i == j or T[i][j] == 0:
                            continue
                        T2 = flip_arc(T, i, j)
                        for c in cycles_v:
                            vset = set(c)
                            # Check: exactly one of i,j in V(C) \ {v}
                            # and the arc i->j not used in cycle
                            uses_arc = any(
                                c[k] == i and c[(k+1) % len(c)] == j
                                for k in range(len(c))
                            )
                            if uses_arc:
                                continue
                            i_in = i in vset
                            j_in = j in vset
                            if not (i_in ^ j_in):  # exactly one
                                continue
                            still_valid = all(
                                T2[c[k]][c[(k+1) % len(c)]] == 1
                                for k in range(len(c))
                            )
                            if not still_valid:
                                continue

                            Tv, ol = delete_vertex(T, v)
                            tc = find_odd_cycles(Tv)
                            mu_T = mu(T, v, c, (Tv, ol, tc))

                            Tv2, ol2 = delete_vertex(T2, v)
                            tc2 = find_odd_cycles(Tv2)
                            mu_T2 = mu(T2, v, c, (Tv2, ol2, tc2))

                            if mu_T != mu_T2:
                                mu_changed_partial += 1
                                if len(details_partial) < 5:
                                    details_partial.append(
                                        (t_idx, v, c, i, j, mu_T, mu_T2,
                                         i_in, j_in))
                            else:
                                mu_same_partial += 1

        print(f"  mu changed={mu_changed_partial}, same={mu_same_partial}")
        if mu_changed_partial > 0:
            print(f"  mu DOES change with one endpoint in cycle!")
            for d in details_partial[:5]:
                who = f"i={d[3]} in C" if d[7] else f"j={d[4]} in C"
                print(f"    T#{d[0]}, v={d[1]}, C={d[2]}, flip {d[3]}->{d[4]} "
                      f"({who}): mu {d[5]} -> {d[6]}")

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("""
Key question: When does mu(C) change under an arc flip i<->j?

Case 1: i,j both OUTSIDE V(C) — mu changes because the flip affects
  the available vertex set's tournament structure.

Case 2: One of i,j in V(C)\{v}, arc not used in cycle — the flip affects
  T-v but not the cycle itself. The available set is V(T-v) \ C\{v}.
  If the one in C is in C\{v}, it's excluded from available set.
  If only j is NOT in V(C), then j IS in the available set and the
  flip changes the tournament there.

The arc-reversal approach needs to track ALL these changes simultaneously.
""")
