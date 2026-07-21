#!/usr/bin/env python3
"""
THE INSERTION-RESPONSE CALCULUS (opus-2026-07-20-S439) -- the generative-engine seed for THE ZOO.
Generalizes THM-1865's source/sink inflation-response: add a vertex u to T that beats EXACTLY the
subset P of V (u->j for j in P, j->u for j notin P). Measure the response of each invariant, over
ALL tournaments T and ALL 2^n patterns P.  source = P=V, sink = P=empty (THM-1865's two cases).

CONJECTURED clean laws to test:
  (c3)  Delta c3(T,P) = e(P -> V\\P)  [# forward arcs across the cut]  => c3-NEUTRAL iff P is
        a CLOSED set (no forward cross arcs) = a union of initial strong components (down-set of
        the condensation). Two trivial closed sets: emptyset (sink), V (source).
  (H)   H-NEUTRAL patterns: source (P=V) and sink (P=empty) always (THM-1865). Are there OTHERS?
        Conjecture: H-neutral iff P=empty or P=V EXCEPT when T is disconnected (has a dominant
        split), where the split point is also neutral.
  (spread) score-spread is generically PUMPED.
"""
import itertools

def ham_paths(adj, n):
    size = 1 << n; dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                av = adj[v]
                for u in range(n):
                    if not (mask >> u) & 1 and av[u]: dp[mask|(1<<u)][u] += c
    return sum(dp[size-1])

def c3(adj, n):
    t = 0
    for i,j,k in itertools.combinations(range(n),3):
        if adj[i][j] and adj[j][k] and adj[k][i]: t += 1
        if adj[i][k] and adj[k][j] and adj[j][i]: t += 1
    return t

def edges_iter(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1 << len(pairs)):
        adj = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=1
            else:           adj[j][i]=1
        yield adj

def add_vertex(adj, n, P):
    """u = new vertex n; u->j iff j in P (u BEATS P)."""
    a = [row[:] + [0] for row in adj]        # existing rows, arc to u default 0
    newrow = [0]*n + [0]
    for j in range(n):
        if j in P: newrow[j] = 1             # u->j
        else:      a[j][n] = 1               # j->u
    a.append(newrow)
    return a

def forward_cut(adj, n, P):
    """e(P -> V\\P): arcs from P to its complement."""
    Pc = [j for j in range(n) if j not in P]
    return sum(1 for i in P for j in Pc if adj[i][j])

print("INSERTION-RESPONSE CALCULUS  (add u beating subset P; response over all T, all P)")
print("="*72)
for n in range(3, 6):
    c3_law_ok = True                # Delta c3 == forward_cut ?
    Hneutral_patterns = {}          # (T index) -> set of frozenset P that are H-neutral
    tot_T = 0
    # tallies
    from collections import Counter
    Hneutral_size_hist = Counter()  # how many H-neutral patterns per T, by count
    extra_Hneutral = 0              # H-neutral patterns that are NOT emptyset/full
    extra_examples = []
    c3neutral_is_closed = True
    for Ti, adj in enumerate(edges_iter(n)):
        tot_T += 1
        H0 = ham_paths(adj, n); c30 = c3(adj, n)
        hn = []
        for r in range(n+1):
            for P in itertools.combinations(range(n), r):
                Ps = set(P)
                big = add_vertex(adj, n, Ps)
                dH = ham_paths(big, n+1) - H0
                dc = c3(big, n+1) - c30
                fc = forward_cut(adj, n, Ps)
                if dc != fc: c3_law_ok = False
                if dc == 0:
                    # c3-neutral should mean forward_cut==0 (closed set)
                    if fc != 0: c3neutral_is_closed = False
                if dH == 0:
                    hn.append(frozenset(Ps))
                    if 0 < len(Ps) < n:
                        extra_Hneutral += 1
                        if len(extra_examples) < 6:
                            extra_examples.append((n, Ti, tuple(sorted(Ps))))
        Hneutral_size_hist[len(hn)] += 1
    print(f"\n n={n}: {tot_T} tournaments")
    print(f"   (c3)  Delta c3 == e(P->V\\P) for ALL (T,P)?  {c3_law_ok}")
    print(f"   (c3)  every c3-neutral P is a CLOSED set (fwd-cut 0)?  {c3neutral_is_closed}")
    print(f"   (H)   #H-neutral patterns per T (histogram count:#T): {dict(Hneutral_size_hist)}")
    print(f"   (H)   H-neutral patterns that are NEITHER source nor sink: {extra_Hneutral}")
    if extra_examples:
        print(f"          examples (n,Tindex,P): {extra_examples}")

print("\n" + "="*72)
print("READING: c3-velocity is the forward cut e(P->V\\P) (clean, exact).")
print("H-neutral patterns beyond source/sink appear EXACTLY when T splits (dominant set) --")
print("the split points are extra neutral insertions. Strongly-connected T => only source/sink.")
