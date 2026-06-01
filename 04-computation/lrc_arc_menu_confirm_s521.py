#!/usr/bin/env python3
"""
lrc_arc_menu_confirm_s521.py   claudebox-2026-06-01-S521

Confirm two claims from lrc_arc_menu_geometry_s521.py with an exact
refinement-based canonical form (fast enough for m=9, no networkx needed):

 (A) menu(n) = 2*F(n-3) for n>=5  (Fibonacci), menu(4)=1   -> confirm menu(10)=26
 (B) the menu and #feasible flip-sets are CONSTANT for all L in (1/2,1)
     (contradicting HYP-1987's "menu grows with L"): test at L near 1/2 and near 1.

Canonical form: color-refine vertices by (color, sorted out-nbr colors, sorted
in-nbr colors) to a stable partition, then take the lex-min adjacency string
over labelings that respect the color order (permute only within equal cells).
Verified against the brute permutation canon for m<=7 inside this script.
"""
from __future__ import annotations
from fractions import Fraction
from itertools import permutations
from functools import lru_cache

HALF = Fraction(1, 2)

# ---- reuse feasibility + flip-set machinery (copied from geometry script) ----
def w_add(a, b): return (a[0] + b[0], a[1] + b[1])
def w_lt(a, b):  return (a[0], a[1]) < (b[0], b[1])

def feasible(m, constraints):
    N = m + 1
    edges = [(b, a, w) for (a, b, w) in constraints] + [(m, v, (Fraction(0), 0)) for v in range(m)]
    dist = [None] * N
    dist[m] = (Fraction(0), 0)
    for _ in range(N - 1):
        changed = False
        for (u, v, w) in edges:
            if dist[u] is None: continue
            nd = w_add(dist[u], w)
            if dist[v] is None or w_lt(nd, dist[v]):
                dist[v] = nd; changed = True
        if not changed: break
    for (u, v, w) in edges:
        if dist[u] is None: continue
        nd = w_add(dist[u], w)
        if dist[v] is None or w_lt(nd, dist[v]): return False
    return True

def realizable(m, S, L):
    cons = []
    for i in range(1, m):
        cons.append((i - 1, i, (Fraction(0), -1)))
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            if (i, j) in S: cons.append((i - 1, j - 1, (-HALF, -1)))
            else:           cons.append((j - 1, i - 1, (HALF, -1)))
    cons.append((m - 1, 0, (L, 0)))
    return feasible(m, cons)

def upsets(m):
    res = []
    def rec(i, prev, cur):
        if i == m: res.append(cur); return
        for tau in range(max(i + 1, prev), m + 2):
            S = set(cur)
            for j in range(tau, m + 1): S.add((i, j))
            rec(i + 1, tau, S)
    rec(1, 2, set())
    return res

def build_adj(m, S):
    adj = [[0] * m for _ in range(m)]
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            if (i, j) in S: adj[i - 1][j - 1] = 1
            else:           adj[j - 1][i - 1] = 1
    return tuple(tuple(r) for r in adj)

def scores(adj): return tuple(sorted(sum(r) for r in adj))
def is_transitive(adj): return scores(adj) == tuple(range(len(adj)))

def Hcount(adj):
    m = len(adj); full = (1 << m) - 1
    @lru_cache(None)
    def dp(mask, last):
        if mask == full: return 1
        return sum(dp(mask | (1 << x), x) for x in range(m) if not (mask >> x) & 1 and adj[last][x])
    return sum(dp(1 << s, s) for s in range(m))

# ---- brute canonical (reference, m<=7) ----
def canon_brute(adj):
    m = len(adj); best = None
    for p in permutations(range(m)):
        flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or flat < best: best = flat
    return best

# ---- refinement-based canonical ----
def refine(adj):
    m = len(adj)
    color = [sum(adj[v]) for v in range(m)]   # out-degree seed
    while True:
        sig = []
        for v in range(m):
            outc = tuple(sorted(color[w] for w in range(m) if adj[v][w]))
            inc  = tuple(sorted(color[w] for w in range(m) if adj[w][v]))
            sig.append((color[v], outc, inc))
        order = sorted(set(sig))
        rank = {s: i for i, s in enumerate(order)}
        newcolor = [rank[sig[v]] for v in range(m)]
        if newcolor == color: break
        color = newcolor
    return color

def canon_fast(adj):
    m = len(adj)
    color = refine(adj)
    cells = {}
    for v in range(m): cells.setdefault(color[v], []).append(v)
    cell_order = [cells[c] for c in sorted(cells)]
    best = None
    # permute within each cell; cell order fixed by color
    def gen(idx, perm):
        nonlocal best
        if idx == len(cell_order):
            p = perm  # p[position] = vertex
            flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
            if best is None or flat < best: best = flat
            return
        for q in permutations(cell_order[idx]):
            gen(idx + 1, perm + list(q))
    gen(0, [])
    return best

def menu(m, L, canon):
    feas = [S for S in upsets(m) if realizable(m, S, L)]
    raw = set(build_adj(m, S) for S in feas)
    classes = {}
    for adj in raw:
        c = canon(adj)
        if c not in classes:
            classes[c] = (Hcount(adj), scores(adj), is_transitive(adj))
    return len(feas), classes

def fib(k):  # F(1)=F(2)=1
    a, b = 1, 1
    for _ in range(k - 1): a, b = b, a + b
    return a

def main():
    print("claudebox-S521 confirmation\n")
    # verify canon_fast == canon_brute on menu counts for m<=7
    print("VALIDATION (fast canon vs brute canon), L=(n-2)/n:")
    for n in range(4, 8):
        m = n - 1; L = Fraction(n - 2, n)
        fb = menu(m, L, canon_brute)[1]
        ff = menu(m, L, canon_fast)[1]
        ok = "OK" if len(fb) == len(ff) else "MISMATCH"
        print(f"  n={n} m={m}: brute={len(fb)} fast={len(ff)}  {ok}")
    print()

    print("(A) menu(n) vs 2*F(n-3) [n>=5], at LRC arc L=(n-2)/n:")
    print(f"  {'n':>2} {'m':>2} {'#feasS':>7} {'menu':>5} {'2F(n-3)':>8} {'2^(m-1)':>8}")
    seq = []
    for n in range(4, 11):
        m = n - 1; L = Fraction(n - 2, n)
        nf, classes = menu(m, L, canon_fast)
        pred = "" if n < 5 else str(2 * fib(n - 3))
        seq.append(len(classes))
        print(f"  {n:>2} {m:>2} {nf:>7} {len(classes):>5} {pred:>8} {2**(m-1):>8}")
    print(f"  menu sequence n=4..10: {seq}")
    print(f"  2*F(n-3) n=5..10:      {[2*fib(n-3) for n in range(5,11)]}")
    print()

    print("(B) L-independence on (1/2,1): menu & #feasS at L near 1/2 and near 1")
    print(f"  {'n':>2} {'m':>2} {'L':>8} {'#feasS':>7} {'menu':>5}")
    for n in (6, 7, 8, 9):
        m = n - 1
        for L in (Fraction(101, 200), Fraction(199, 200)):  # 0.505, 0.995
            nf, classes = menu(m, L, canon_fast)
            print(f"  {n:>2} {m:>2} {float(L):>8.3f} {nf:>7} {len(classes):>5}")

if __name__ == "__main__":
    main()
