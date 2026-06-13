#!/usr/bin/env python3
"""
Claim A via inclusion-exclusion on cycles through v.

We proved: H(T) - H(T-v) = sum_{S: S ni some C ni v} 2^|S|

The RHS sums over independent sets containing >= 1 cycle through v.
By inclusion-exclusion on which cycles through v are used:

  RHS = sum_{C ni v} [sum_{S: C in S, S indep} 2^|S|]
      - sum_{C1,C2 ni v, C1 != C2} [sum_{S: C1,C2 in S, S indep} 2^|S|]
      + ...

Let f(C) = sum_{S: C in S, S indep} 2^|S| = 2 * sum_{S': S' indep in Omega(T)|_{avoid N(C)}} 2^|S'|
         = 2 * I(Omega(T)|_{non-conflicting with C}, 2)

Wait, more carefully:
  f(C) = sum_{S indep in Omega: C in S} 2^|S|
        = 2 * [sum_{S' indep in Omega: S' subset of {cycles not conflicting with C}} 2^|S'|]
        = 2 * I(Omega(T) - N[C], 2)

where N[C] = closed neighborhood of C in Omega(T) = {C} union {cycles sharing vertex with C}.

And mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2).

So Claim A says:
  sum_{inc-exc over cycles through v} = 2 * sum_{C ni v} I(Omega(T-v)|_{avoid C\{v}}, 2)

The LHS can also be written as:
  sum_{C ni v} 2*I(Omega-N[C], 2) - sum_{C1<C2} 4*I(Omega-N[C1]-N[C2], 2) + ...

QUESTION: Does the inclusion-exclusion collapse to give exactly
  2 * sum_{C ni v} mu(C)?

Let's test this at n=5 in detail.

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import defaultdict

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def find_ALL_directed_odd_cycles(A, n, vertices=None):
    if vertices is None:
        vertices = list(range(n))
    all_cycles = []
    for length in range(3, len(vertices)+1, 2):
        for subset in combinations(vertices, length):
            for perm in permutations(subset):
                is_cycle = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1)%length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = list(perm).index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    all_cycles.append(canon)
    return list(set(all_cycles))

def I_subgraph(cycles, full_adj, allowed_indices, x=2):
    """Compute I(subgraph on allowed_indices, x)."""
    allowed = list(allowed_indices)
    na = len(allowed)
    total = 0
    for mask in range(2**na):
        nodes = [allowed[i] for i in range(na) if (mask >> i) & 1]
        ok = True
        for a in range(len(nodes)):
            for b in range(a+1, len(nodes)):
                if full_adj[nodes[a]][nodes[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += x**len(nodes)
    return total

def main():
    n = 5
    m = n*(n-1)//2
    print(f"--- Claim A inclusion-exclusion structure at n={n} ---\n")

    # Pick a tournament with interesting structure
    for bits_int in [100, 200, 920, 40]:
        if bits_int >= 2**m:
            continue
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        cycles = find_ALL_directed_odd_cycles(A, n)
        nc = len(cycles)
        if nc == 0:
            continue

        # Build conflict graph
        adj = [[False]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if set(cycles[i]) & set(cycles[j]):
                    adj[i][j] = adj[j][i] = True

        H = I_subgraph(cycles, adj, range(nc))

        print(f"\n{'='*50}")
        print(f"bits={bits_int}, H={H}, #cycles={nc}")
        for i, c in enumerate(cycles):
            print(f"  C{i}: {c} (len={len(c)})")

        for v in range(n):
            cycles_v = [i for i in range(nc) if v in cycles[i]]
            cycles_not_v = [i for i in range(nc) if v not in cycles[i]]

            remaining = [u for u in range(n) if u != v]
            cycles_Tv = find_ALL_directed_odd_cycles(A, n, remaining)
            nc_v = len(cycles_Tv)

            # H(T-v) via subgraph of Omega on cycles not through v
            H_Tv = I_subgraph(cycles, adj, cycles_not_v)

            diff = H - H_Tv

            # RHS of Claim A: 2 * sum_C mu(C)
            rhs = 0
            mu_vals = {}
            for ci in cycles_v:
                # mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
                forbidden = set(cycles[ci]) - {v}
                # Cycles in T-v not conflicting with C\{v}
                allowed = []
                for j in range(nc_v):
                    if not (set(cycles_Tv[j]) & forbidden):
                        allowed.append(j)
                # Build sub-adjacency
                adj_tv = [[False]*nc_v for _ in range(nc_v)]
                for ii in range(nc_v):
                    for jj in range(ii+1, nc_v):
                        if set(cycles_Tv[ii]) & set(cycles_Tv[jj]):
                            adj_tv[ii][jj] = adj_tv[jj][ii] = True
                mu = I_subgraph(cycles_Tv, adj_tv, allowed)
                mu_vals[ci] = mu
                rhs += mu

            rhs *= 2

            # Also compute: f(C) = sum_{S: C in S, S indep} 2^|S|
            # = 2 * I(Omega - N[C], 2) where N[C] = {C} + neighbors of C
            f_vals = {}
            for ci in cycles_v:
                # Non-conflicting with C (excluding C itself)
                non_conflict = [j for j in range(nc) if j != ci and not adj[ci][j]]
                f_vals[ci] = 2 * I_subgraph(cycles, adj, non_conflict)

            # Inclusion-exclusion: LHS via ie
            # For subsets of cycles_v
            ie_sum = 0
            for r in range(1, len(cycles_v)+1):
                for subset in combinations(cycles_v, r):
                    # Independent sets containing ALL cycles in subset
                    # These cycles must be mutually non-adjacent
                    mutual_ok = True
                    for a in range(len(subset)):
                        for b_idx in range(a+1, len(subset)):
                            if adj[subset[a]][subset[b_idx]]:
                                mutual_ok = False
                                break
                        if not mutual_ok:
                            break
                    if not mutual_ok:
                        continue  # can't have both in same indep set

                    # Non-conflicting with all cycles in subset
                    non_conflict = set(range(nc))
                    for ci in subset:
                        non_conflict -= {ci}
                        non_conflict -= {j for j in range(nc) if adj[ci][j]}

                    contrib = ((-1)**(r+1)) * (2**r) * I_subgraph(cycles, adj, list(non_conflict))
                    ie_sum += contrib

            if len(cycles_v) <= 5:
                print(f"\n  v={v}: diff={diff}, RHS(Claim A)={rhs}, IE_sum={ie_sum}")
                print(f"    Cycles through v: {[cycles[i] for i in cycles_v]}")
                for ci in cycles_v:
                    print(f"    C{ci}: f(C)={f_vals[ci]}, mu(C)={mu_vals[ci]}")
                if diff == rhs:
                    print(f"    Claim A: VERIFIED")
                else:
                    print(f"    Claim A: FAIL (diff={diff}, rhs={rhs})")

                # KEY: is f(C) = 2*mu(C)?
                for ci in cycles_v:
                    if f_vals[ci] == 2 * mu_vals[ci]:
                        print(f"    f(C{ci}) = 2*mu(C{ci}): YES")
                    else:
                        print(f"    f(C{ci}) = 2*mu(C{ci}): NO ({f_vals[ci]} vs {2*mu_vals[ci]})")

    print("\n\n" + "="*50)
    print("ANALYSIS")
    print("="*50)
    print()
    print("f(C) = 2 * I(Omega(T) - N[C], 2)")
    print("mu(C) = I(Omega(T-v) - N'[C], 2)")
    print()
    print("where N[C] = closed nbhd of C in Omega(T)")
    print("  and N'[C] = cycles in T-v conflicting with C\\{v}")
    print()
    print("If f(C) = 2*mu(C) for all C, then inclusion-exclusion is trivial")
    print("and Claim A reduces to: each cycle contributes independently.")
    print("But this is unlikely for cycles that share vertices...")

if __name__ == "__main__":
    main()
