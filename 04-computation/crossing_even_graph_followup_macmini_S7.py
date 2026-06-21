#!/usr/bin/env python3
"""
FOLLOW-UP: nail the even-graph <-> 2-page-crossing correspondence and push n=7.
mac-mini-2026-06-20-S7.

Findings from part 1: brute 2-page = Guy Z(n) for n=3..6 (reproduces Abrego et al).
Optimizer is NOT length-parity. Z/C(n,4) = 0,0,1/5,1/5,9/35,9/35,2/7,2/7 (paired,
no clean closed form vs E[c3]).

This part:
  (A) Correct even-page test. For a 2-coloring of E(K_n), "page d even" means each
      vertex has even degree on page d. Total degree of each vertex is n-1.
      - n ODD  => n-1 even => page0 even  <=>  page1 even. Both-even possible.
      - n EVEN => n-1 odd  => exactly one page is even at each vertex; cannot have
        BOTH pages globally even. So "both even" is impossible for n even.
      Correct definition for the EVEN-GRAPH correspondence: a 2-page assignment
      where PAGE 0 (the "above" page) is an EVEN GRAPH. (Its complement-on-spine is
      page 1; for n odd page 1 is then also even.)
  (B) Among assignments with page-0 an even graph, what is the MIN crossing number?
      Does it still equal Z(n)? (H2)
  (C) Smarter (cyclic-construction) check at n=7,8,9: evaluate the KNOWN optimal
      2-page construction (Blazek-Koman / cylindrical) and confirm it equals Z(n),
      and check whether that optimal drawing's page-0 is an even graph.
  (D) The cylindrical construction page-0: which edges go above? Standard optimal:
      put vertices on a circle 0..n-1; page assignment by edge "type". The
      Guy/Blazek-Koman optimum: split the n vertices into two near-equal arcs is
      the CYLINDRICAL drawing, distinct from book. We test the BOOK optimum's
      structure directly via local search at n=7.
"""

from itertools import combinations
from math import comb
import random

def guy(n):
    return (n//2)*((n-1)//2)*((n-2)//2)*((n-3)//2)//4

def edges_of(n):
    return [(i, j) for i in range(n) for j in range(i+1, n)]

def interleave(e, f):
    (a, b), (c, d) = e, f
    return (a < c < b < d) or (c < a < d < b)

def crossing_pairs(n):
    edges = edges_of(n)
    cps = []
    for x in range(len(edges)):
        for y in range(x+1, len(edges)):
            e, f = edges[x], edges[y]
            if len(set(e) | set(f)) == 4 and interleave(e, f):
                cps.append((x, y))
    return edges, cps

def cross_count(n, assign, cps):
    return sum(1 for (x, y) in cps if assign[x] == assign[y])

def page_degrees(n, edges, assign, page):
    deg = [0]*n
    for idx, (i, j) in enumerate(edges):
        if assign[idx] == page:
            deg[i] += 1; deg[j] += 1
    return deg

def is_even_page(n, edges, assign, page):
    return all(d % 2 == 0 for d in page_degrees(n, edges, assign, page))

# ----- exact brute for n<=6 with corrected even definition -----
def brute(n, constraint=None):
    edges, cps = crossing_pairs(n)
    m = len(edges)
    if m > 18:
        return None, None
    best, ba = None, None
    for mask in range(1 << m):
        assign = [(mask >> k) & 1 for k in range(m)]
        if constraint and not constraint(n, edges, assign):
            continue
        c = cross_count(n, assign, cps)
        if best is None or c < best:
            best, ba = c, assign
    return best, ba

def page0_even(n, edges, assign):
    return is_even_page(n, edges, assign, 0)

# ----- local search (simulated annealing-lite) for n=7,8,9 -----
def local_search(n, restarts=40, iters=8000, seed=0, constraint=None):
    rng = random.Random(seed)
    edges, cps = crossing_pairs(n)
    m = len(edges)
    # incidence of each edge in cps for fast delta
    involved = [[] for _ in range(m)]
    for pi, (x, y) in enumerate(cps):
        involved[x].append(pi); involved[y].append(pi)
    best_overall = None; best_assign = None
    for r in range(restarts):
        assign = [rng.randint(0, 1) for _ in range(m)]
        if constraint:
            # repair toward constraint not implemented for even; only run constraint-free LS
            pass
        cur = cross_count(n, assign, cps)
        T = 2.0
        for it in range(iters):
            k = rng.randrange(m)
            # delta if we flip edge k
            delta = 0
            for pi in involved[k]:
                x, y = cps[pi]
                other = y if x == k else x
                same_now = (assign[x] == assign[y])
                # after flip, assign[k] changes
                newk = 1 - assign[k]
                if x == k:
                    same_new = (newk == assign[y])
                else:
                    same_new = (assign[x] == newk)
                delta += (1 if same_new else 0) - (1 if same_now else 0)
            if delta <= 0 or rng.random() < pow(2.718281828, -delta / max(T, 1e-9)):
                assign[k] = 1 - assign[k]
                cur += delta
            T *= 0.9995
        if best_overall is None or cur < best_overall:
            best_overall = cur; best_assign = assign[:]
    return best_overall, best_assign

def describe_optimizer(n, edges, assign):
    p0 = sum(1 for a in assign if a == 0)
    p1 = len(assign) - p0
    even0 = is_even_page(n, edges, assign, 0)
    even1 = is_even_page(n, edges, assign, 1)
    deg0 = page_degrees(n, edges, assign, 0)
    return p0, p1, even0, even1, deg0

def main():
    print("="*78)
    print("(A/B) Corrected even-page constraint: min crossings with page-0 even")
    print("="*78)
    print(f"{'n':>2} {'Z(n)':>5} {'nu2 (free)':>10} {'nu2 (page0 even)':>16} {'both=Z?':>8}")
    for n in range(3, 7):
        Z = guy(n)
        free, fa = brute(n)
        ev, ea = brute(n, constraint=page0_even)
        both = (free == Z and (ev is None or ev == Z))
        print(f"{n:>2} {Z:>5} {str(free):>10} {str(ev):>16} {str(both):>8}")
        if ea is not None:
            edges = edges_of(n)
            p0,p1,e0,e1,deg0 = describe_optimizer(n, edges, ea)
            print(f"      page0-even optimizer: |page0|={p0} |page1|={p1} "
                  f"even0={e0} even1={e1} deg0={deg0}")

    print()
    print("="*78)
    print("(C) Local search for n=7,8,9 (free), compare to Z(n)")
    print("="*78)
    for n in [7, 8, 9]:
        Z = guy(n)
        best, ba = local_search(n, restarts=60, iters=12000, seed=12345+n)
        edges = edges_of(n)
        p0, p1, e0, e1, deg0 = describe_optimizer(n, edges, ba)
        print(f" n={n}: Z={Z}  LS_best={best}  match={best==Z}  "
              f"|page0|={p0}|page1|={p1} optimizer_page0_even={e0} page1_even={e1}")
        print(f"        page0 degree seq: {sorted(deg0)}")

    print()
    print("="*78)
    print("(D) Does an OPTIMAL 2-page drawing ALWAYS admit an even page-0?")
    print("    Search over all optimal assignments (n<=6) for one with page0 even.")
    print("="*78)
    for n in range(3, 7):
        edges, cps = crossing_pairs(n)
        m = len(edges)
        Z = guy(n)
        found_even = False
        n_opt = 0
        for mask in range(1 << m):
            assign = [(mask >> k) & 1 for k in range(m)]
            if cross_count(n, assign, cps) == Z:
                n_opt += 1
                if is_even_page(n, edges, assign, 0):
                    found_even = True
        print(f" n={n}: #optimal assignments={n_opt}  "
              f"some-optimal-has-even-page0={found_even}")

if __name__ == "__main__":
    main()
