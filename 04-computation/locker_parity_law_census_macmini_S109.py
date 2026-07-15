#!/usr/bin/env python3
"""locker_parity_law_census_macmini_S109.py -- mac-mini-2026-07-15-S109.
THM-865 exploration: the locker tournament D_n (u<v: v->u iff u|v or v=u+1, else u->v).
(1) census of directed odd cycles alpha_1(D_n), n=3..13; check EVEN (the law, via THM-466).
(2) stratify by top vertex: t_m = #odd cycles through m in D_m (alpha_1(D_n) = sum_{m<=n} t_m).
(3) within t_m: parity tables by exit vertex (m -> b), by entry vertex (u -> m), by (exit,entry),
    and by length -- to locate the involution.
"""
import sys
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

def build_Dn(n):
    """adj[u] = set of out-neighbours, vertices 1..n."""
    adj = {u: set() for u in range(1, n + 1)}
    for u in range(1, n + 1):
        for v in range(u + 1, n + 1):
            if v % u == 0 or v == u + 1:
                adj[v].add(u)      # bigger beats smaller (down-arc)
            else:
                adj[u].add(v)      # smaller beats bigger (up-arc)
    return adj

def cycles_through_top(n):
    """All simple directed cycles through vertex n in D_n. Returns list of tuples
    (length, exit_vertex b, entry_vertex u, cycle_as_tuple starting at n)."""
    adj = build_Dn(n)
    out = []
    # DFS from n; cycle recorded when we return to n
    stack = [(n, (n,), frozenset((n,)))]
    while stack:
        v, path, seen = stack.pop()
        for w in adj[v]:
            if w == n and len(path) >= 3:
                out.append((len(path), path[1], path[-1], path))
            elif w != n and w not in seen:
                stack.append((w, path + (w,), seen | {w}))
    return out

def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    total = 0
    print("m : t_m (odd cycles through top vertex m) ; cumulative alpha_1(D_m) ; parity")
    per_m = {}
    for m in range(3, NMAX + 1):
        cyc = cycles_through_top(m)
        odd = [c for c in cyc if c[0] % 2 == 1]
        per_m[m] = odd
        t = len(odd)
        total += t
        # exit/entry parity tables
        from collections import Counter
        by_exit = Counter(c[1] for c in odd)
        by_entry = Counter(c[2] for c in odd)
        by_len = Counter(c[0] for c in odd)
        divisors = sorted(d for d in range(2, m) if m % d == 0)
        print(f"m={m:2d}: t_m={t:8d} (par {t%2})  alpha1(D_{m})={total:9d} (par {total%2})"
              f"  ALL={len(cyc):8d}")
        print(f"      by_len(odd): {dict(sorted(by_len.items()))}")
        print(f"      exits: divisors(m)={divisors} + {{{m-1}}}")
        print(f"      by_exit  (odd): {dict(sorted(by_exit.items()))}")
        print(f"      by_exit  PARITY: { {k: v % 2 for k, v in sorted(by_exit.items())} }")
        print(f"      by_entry PARITY: { {k: v % 2 for k, v in sorted(by_entry.items())} }")
        # (exit,entry) parity matrix, only odd-parity cells
        by_ee = Counter((c[1], c[2]) for c in odd)
        oddcells = sorted(k for k, v in by_ee.items() if v % 2 == 1)
        print(f"      odd-parity (exit,entry) cells: {oddcells}")
    print()
    print("SUMMARY alpha_1(D_n) parity:", "ALL EVEN" if all(
        sum(len(per_m[m]) for m in range(3, n + 1)) % 2 == 0 for n in range(3, NMAX + 1))
        else "PARITY FAILURE -- LAW REFUTED")

if __name__ == "__main__":
    main()
