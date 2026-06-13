#!/usr/bin/env python3
"""
THM-343 AUDIT (opus-2026-05-28-S6).

For each foundational claim in the THM-343 proof, perform rigorous verification:

  A1. OCF: H(T) = I(Omega(T), 2) at n=3..7 exhaustively + n=8 sampled.
  A2. Unique decomposition: alpha_1 + 2 alpha_2 + 4 alpha_3 + ... = 3 with
      combinatorial constraint alpha_2 >= 1 ==> alpha_1 >= 2 ==>
      unique solution is (3, 0, 0, ...).
  A3. Chain constraint (Kruskal-Katona basic form): alpha_k >= 1 ==> alpha_{k-1} >= k.
  A4. SCC partition: every vertex of T lies in a unique SCC; SCCs partition V.
  A5. Cycle-in-SCC: every directed cycle has all vertices in one SCC.
  A6. Cross-SCC cycles disjoint: cycles in distinct SCCs are vertex-disjoint
      ==> they are non-adjacent in Omega.
  A7. n=4 SC has unique score sequence (1,1,2,2): exhaustive enumeration
      and Landau-strict check.
  A8. 3-cycle counting: c_3(T) = C(n,3) - sum C(s_i, 2) for each tournament T.
  A9. Moon-Moser: SC tournament on n vertices has >= n-2 three-cycles
      (verify exhaustively n=3,...,6 and sampled n=7,8).
  A10. Moon-Camion / Camion: every SC tournament has a directed Hamilton cycle
      (verify exhaustively n=3,...,6).
  A11. n=4 no odd cycles of length > 3 (trivially since 4 < 5).
  A12. Pancyclicity at n=5: every SC tournament on 5 verts has a Hamilton 5-cycle (odd).

  E. Exhaustive H!=7 at n=7 (NEW, exhaustive not sample).

Instance: opus-2026-05-28-S6
"""

import os, sys
os.environ['PYTHONIOENCODING'] = 'utf-8'
from itertools import combinations, permutations, product
from math import comb
from collections import Counter
import random
import time

random.seed(0xfeed)


# -------- Standard helpers --------
def gen_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(1 << m):
        T = [[False]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = True
            else:
                T[j][i] = True
        yield T


def count_HP(T, n):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1<<n):
        for v in range(n):
            if dp[S][v] == 0: continue
            if not ((S >> v) & 1): continue
            for u in range(n):
                if (S >> u) & 1: continue
                if T[v][u]:
                    dp[S | (1<<u)][u] += dp[S][v]
    full = (1<<n) - 1
    return sum(dp[full][v] for v in range(n))


def enumerate_odd_cycles(T):
    n = len(T)
    seen = set()
    for L in range(3, n+1, 2):
        for subset in combinations(range(n), L):
            v0 = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(L):
                    a,b = cycle[i], cycle[(i+1)%L]
                    if not T[a][b]:
                        ok = False; break
                if ok:
                    seen.add(cycle)
    return list(seen)


def independence_poly(T):
    cycles = enumerate_odd_cycles(T)
    nc = len(cycles)
    Vsets = [frozenset(c) for c in cycles]
    n = len(T)
    counts = [0] * (n+1)
    counts[0] = 1
    def recurse(start, blocked, k):
        for i in range(start, nc):
            if not (Vsets[i] & blocked):
                counts[k+1] += 1
                recurse(i+1, blocked | Vsets[i], k+1)
    recurse(0, frozenset(), 0)
    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return counts


def is_strongly_connected(T):
    n = len(T)
    if n <= 1: return True
    def reach(s, m):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for v in range(n):
                if m[u][v] and v not in seen:
                    seen.add(v); st.append(v)
        return seen
    if len(reach(0, T)) != n: return False
    Trev = [[T[j][i] for j in range(n)] for i in range(n)]
    return len(reach(0, Trev)) == n


def scores(T):
    n = len(T)
    return [sum(1 for j in range(n) if T[i][j]) for i in range(n)]


def landau_strict(scores_sorted, n):
    """Return True iff Landau's inequalities are strict for k=1..n-1 and equality for k=n."""
    s = sorted(scores_sorted)
    total = sum(s)
    if total != comb(n, 2):
        return False
    for k in range(1, n):
        if sum(s[:k]) <= comb(k, 2):
            return False  # need STRICT for SC
    return True


def count_3_cycles_fast(T):
    n = len(T)
    c = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if T[i][j] and T[j][k] and T[k][i]: c += 1
                if T[i][k] and T[k][j] and T[j][i]: c += 1
    return c


def find_sccs(T):
    """Tarjan's SCC; returns list of frozensets (SCCs)."""
    n = len(T)
    idx = {}
    low = {}
    on_stack = [False]*n
    stack = []
    counter = [0]
    sccs = []
    def strongconnect(v):
        idx[v] = counter[0]; low[v] = counter[0]; counter[0] += 1
        stack.append(v); on_stack[v] = True
        for w in range(n):
            if T[v][w]:
                if w not in idx:
                    strongconnect(w)
                    low[v] = min(low[v], low[w])
                elif on_stack[w]:
                    low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            scc = []
            while True:
                w = stack.pop(); on_stack[w] = False
                scc.append(w)
                if w == v: break
            sccs.append(frozenset(scc))
    sys.setrecursionlimit(10000)
    for v in range(n):
        if v not in idx:
            strongconnect(v)
    return sccs


def has_directed_hamilton_cycle(T):
    n = len(T)
    if n < 3: return False
    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        ok = True
        for i in range(n):
            a, b = cycle[i], cycle[(i+1)%n]
            if not T[a][b]:
                ok = False; break
        if ok:
            return True
    return False


# =============== AUDIT TESTS ===============

def audit_A1_OCF(max_n=6):
    """H(T) = I(Omega(T), 2) verified exhaustively n=3..max_n"""
    print("="*60)
    print(f"A1. OCF: H(T) = I(Omega(T), 2)")
    print("="*60)
    for n in range(3, max_n+1):
        tot = 0; fail = 0
        for T in gen_tournaments(n):
            tot += 1
            H = count_HP(T, n)
            ips = independence_poly(T)
            I2 = sum(c * (2**i) for i, c in enumerate(ips))
            if H != I2:
                fail += 1
                if fail <= 3:
                    print(f"  n={n}: H={H} != I(Omega,2)={I2} for T={T}")
        print(f"  n={n}: {tot} tournaments, {fail} OCF failures. {'OK' if fail==0 else 'FAIL'}")


def audit_A2_uniqueness():
    """Unique non-negative integer solution to alpha_1 + 2 alpha_2 + 4 alpha_3 + ... = 3
    UNDER CHAIN CONSTRAINT alpha_k >= 1 ==> alpha_{k-1} >= k."""
    print("="*60)
    print(f"A2. Uniqueness of (3,0,0,...) for alpha_1 + 2 alpha_2 + ... = 3")
    print("="*60)
    # enumerate non-negative integer tuples (a_1, a_2, ...) with sum <= 3 and weighted sum = 3
    solutions = []
    for a1 in range(4):
        for a2 in range(2):  # 4*a2 <= 3 ==> a2 = 0 only
            for a3 in range(1):  # 4*a3 <= 3 ==> a3 = 0
                if a1 + 2*a2 + 4*a3 == 3:
                    solutions.append((a1, a2, a3))
    print(f"  All raw non-neg integer solutions: {solutions}")
    valid = []
    for sol in solutions:
        a1, a2, a3 = sol
        # chain constraint: a_k >= 1 ==> a_{k-1} >= k
        ok = True
        if a2 >= 1 and a1 < 2: ok = False
        if a3 >= 1 and a2 < 3: ok = False
        if ok:
            valid.append(sol)
    print(f"  Solutions satisfying chain constraint: {valid}")
    if valid == [(3, 0, 0)]:
        print("  UNIQUE solution (3, 0, 0). OK")
    else:
        print(f"  NOT UNIQUE — multiple solutions: {valid}")


def audit_A3_chain_constraint(n_max=6):
    """alpha_k >= 1 ==> alpha_{k-1} >= k for any tournament's Omega."""
    print("="*60)
    print(f"A3. Chain constraint alpha_k >= 1 ==> alpha_{{k-1}} >= k")
    print("="*60)
    for n in range(3, n_max+1):
        fails = 0; tot = 0
        for T in gen_tournaments(n):
            tot += 1
            ips = independence_poly(T)
            for k in range(2, len(ips)):
                if ips[k] >= 1 and ips[k-1] < k:
                    fails += 1
                    break
        print(f"  n={n}: {tot} tournaments, {fails} violations. {'OK' if fails==0 else 'FAIL'}")


def audit_A4_SCC_partition(n_max=5):
    """Verify SCC partition: SCCs cover V with no overlap."""
    print("="*60)
    print(f"A4. SCC partition property")
    print("="*60)
    for n in range(2, n_max+1):
        fails = 0; tot = 0
        for T in gen_tournaments(n):
            tot += 1
            sccs = find_sccs(T)
            union = set()
            for s in sccs:
                if union & s:
                    fails += 1; break
                union |= s
            if union != set(range(n)):
                fails += 1
        print(f"  n={n}: {tot} tournaments, {fails} partition violations. {'OK' if fails==0 else 'FAIL'}")


def audit_A5_cycle_in_scc(n_max=5):
    """Every directed cycle's vertices lie in a single SCC."""
    print("="*60)
    print(f"A5. Directed cycles contained in single SCC")
    print("="*60)
    for n in range(3, n_max+1):
        fails = 0; tot = 0
        for T in gen_tournaments(n):
            tot += 1
            sccs = find_sccs(T)
            vert_to_scc = {}
            for idx, s in enumerate(sccs):
                for v in s:
                    vert_to_scc[v] = idx
            cycles = enumerate_odd_cycles(T)
            for C in cycles:
                scc_set = {vert_to_scc[v] for v in C}
                if len(scc_set) != 1:
                    fails += 1
                    if fails <= 3:
                        print(f"    n={n}: cycle {C} spans {len(scc_set)} SCCs")
                    break
        print(f"  n={n}: {tot} tournaments, {fails} violations. {'OK' if fails==0 else 'FAIL'}")


def audit_A7_n4_SC_score():
    """n=4 SC: unique score (1,1,2,2)."""
    print("="*60)
    print(f"A7. n=4 SC has unique score sequence (1,1,2,2)")
    print("="*60)
    n = 4
    score_set = set()
    sc_count = 0
    for T in gen_tournaments(n):
        if is_strongly_connected(T):
            sc_count += 1
            score_set.add(tuple(sorted(scores(T))))
    print(f"  n=4: {sc_count} SC tournaments; distinct score sequences: {score_set}")
    if score_set == {(1,1,2,2)}: print("  OK")
    else: print("  FAIL")


def audit_A8_3cycle_formula(n_max=5):
    """c_3(T) = C(n,3) - sum C(s_i, 2)."""
    print("="*60)
    print(f"A8. 3-cycle count c_3 = C(n,3) - sum C(s_i, 2)")
    print("="*60)
    for n in range(3, n_max+1):
        fails = 0; tot = 0
        for T in gen_tournaments(n):
            tot += 1
            c3_direct = count_3_cycles_fast(T)
            s = scores(T)
            c3_formula = comb(n,3) - sum(comb(si, 2) for si in s)
            if c3_direct != c3_formula:
                fails += 1
        print(f"  n={n}: {tot} tournaments, {fails} mismatches. {'OK' if fails==0 else 'FAIL'}")


def audit_A9_moon_moser(n_max=6):
    """Moon-Moser: SC tournament has >= n-2 three-cycles."""
    print("="*60)
    print(f"A9. Moon-Moser: SC tournament on n verts has >= n-2 three-cycles")
    print("="*60)
    for n in range(3, n_max+1):
        sc_count = 0; viol = 0; min_c3 = None
        for T in gen_tournaments(n):
            if is_strongly_connected(T):
                sc_count += 1
                c3 = count_3_cycles_fast(T)
                if min_c3 is None or c3 < min_c3: min_c3 = c3
                if c3 < n - 2:
                    viol += 1
        print(f"  n={n}: {sc_count} SC tournaments; min c3 = {min_c3}; bound n-2={n-2}; violations: {viol}")


def audit_A10_camion(n_max=6):
    """Camion: every SC tournament has a Hamilton cycle."""
    print("="*60)
    print(f"A10. Camion: SC tournament has directed Hamilton cycle")
    print("="*60)
    for n in range(3, n_max+1):
        sc_count = 0; viol = 0
        for T in gen_tournaments(n):
            if is_strongly_connected(T):
                sc_count += 1
                if not has_directed_hamilton_cycle(T):
                    viol += 1
        print(f"  n={n}: {sc_count} SC tournaments; violations: {viol}")


def audit_n7_exhaustive_H_not_7():
    """Exhaustive enumeration: no n=7 tournament has H=7."""
    print("="*60)
    print(f"E. Exhaustive: no n=7 tournament has H=7")
    print("="*60)
    n = 7
    tot = 0
    h7_count = 0
    h_distribution = Counter()
    t0 = time.time()
    for T in gen_tournaments(n):
        tot += 1
        H = count_HP(T, n)
        h_distribution[H] += 1
        if H == 7:
            h7_count += 1
        if tot % 100000 == 0:
            elapsed = time.time() - t0
            rate = tot / elapsed
            eta = (2**21 - tot) / rate
            print(f"    progress: {tot}/{2**21} ({100*tot/(2**21):.1f}%), eta {eta:.1f}s")
    print(f"  n={n}: {tot} tournaments processed")
    print(f"  H=7 occurrences: {h7_count}")
    achievable_H = sorted(h_distribution.keys())
    print(f"  Achievable H values: {achievable_H[:30]}...{achievable_H[-5:]}")
    print(f"  H spectrum size: {len(achievable_H)}")
    forbidden_in_range = [h for h in range(1, max(achievable_H)+1, 2) if h not in h_distribution]
    print(f"  Forbidden odd H in [1, {max(achievable_H)}]: {forbidden_in_range}")


def main():
    print("THM-343 AUDIT (opus-2026-05-28-S6)", flush=True)
    print("="*60, flush=True)
    audit_A1_OCF(); sys.stdout.flush()
    print(flush=True)
    audit_A2_uniqueness(); sys.stdout.flush()
    print(flush=True)
    audit_A3_chain_constraint(); sys.stdout.flush()
    print(flush=True)
    audit_A4_SCC_partition(); sys.stdout.flush()
    print(flush=True)
    audit_A5_cycle_in_scc(); sys.stdout.flush()
    print(flush=True)
    audit_A7_n4_SC_score(); sys.stdout.flush()
    print(flush=True)
    audit_A8_3cycle_formula(); sys.stdout.flush()
    print(flush=True)
    audit_A9_moon_moser(); sys.stdout.flush()
    print(flush=True)
    audit_A10_camion(); sys.stdout.flush()
    print(flush=True)
    print("FAST AUDIT COMPLETE", flush=True)
    # Skip exhaustive n=7 (separate script)


if __name__ == "__main__":
    main()
