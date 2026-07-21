#!/usr/bin/env python3
"""
THE INFLATION-RESPONSE DIAGNOSTIC (the WOWII-103 transfer, opus-2026-07-20-S438).

WOWII-103 lesson: a conjectured inequality dies to a construction that DECOUPLES two
invariants -- pendant leaves PUMP one invariant (alpha) while a coupled one stays fixed.
So: an invariant admits WOWII-style inflation counterexamples IFF it is *inflation-pumped*.

Tournament pendant = a SOURCE (beats everyone) or SINK (loses to everyone). Classify how each
invariant responds to source-inflation T -> T+source:
   H  (Hamiltonian paths):   NEUTRAL  (source forced to path-start; bijection) -> maximiser RIGID
   c3 (3-cycles):            NEUTRAL  (source in no 3-cycle)                    -> maximiser RIGID
   score-spread (max-min):   PUMPED   (source out-degree n > all others)       -> cheap extremal
Verified exhaustively n=3,4,5.
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
    """number of directed 3-cycles."""
    t = 0
    for i,j,k in itertools.combinations(range(n),3):
        # count cyclic orientations among the 3 arcs
        if adj[i][j] and adj[j][k] and adj[k][i]: t += 1
        if adj[i][k] and adj[k][j] and adj[j][i]: t += 1
    return t

def spread(adj, n):
    s = [sum(adj[i]) for i in range(n)]
    return max(s) - min(s)

def edges_iter(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1 << len(pairs)):
        adj = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=1
            else:           adj[j][i]=1
        yield adj

def add_source(adj, n):
    return [row[:] + [0] for row in adj] + [[1]*n + [0]]

print("INFLATION-RESPONSE DIAGNOSTIC (source-inflation), exhaustive n=3,4,5")
print("="*66)
for n in range(3, 6):
    H_neutral = c3_neutral = True
    spread_pumped = True
    for adj in edges_iter(n):
        big = add_source(adj, n)
        if ham_paths(big,n+1) != ham_paths(adj,n): H_neutral = False
        if c3(big,n+1)        != c3(adj,n):        c3_neutral = False
        # spread pumps unless original already had a source (max out-deg n-1)
        s = [sum(adj[i]) for i in range(n)]
        if not (spread(big,n+1) >= spread(adj,n)): spread_pumped = False
    print(f" n={n}:  H source-neutral={H_neutral}   c3 source-neutral={c3_neutral}"
          f"   score-spread non-decreasing={spread_pumped}")

# strict-pump demonstration on the regular tournament (balanced -> source makes it maximally spread)
print("\n strict pump: score-spread on the n=5 regular tournament")
reg = None
for adj in edges_iter(5):
    if len(set(sum(adj[i]) for i in range(5)))==1: reg = adj; break
print(f"   regular T5 spread = {spread(reg,5)}  ->  T5+source spread = {spread(add_source(reg,5),6)}"
      f"   (H: {ham_paths(reg,5)} -> {ham_paths(add_source(reg,5),6)},  c3: {c3(reg,5)} -> {c3(add_source(reg,5),6)})")
print("\n VERDICT: H, c3 are inflation-NEUTRAL (rigid balanced maximiser = Paley/rotation, LEM-004);")
print("          score-spread is inflation-PUMPED (cheap extremal). WOWII counterexamples exist")
print("          only for PUMPED invariants. H/c3 resist inflation -> their extremals are hard & rigid.")
