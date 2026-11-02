#!/usr/bin/env python3
"""prime_squaresum.py (klein-2026-11-01) -- the square-sum graph induced on primes: p ~ q iff p + q is a perfect square.

(1) the 25 primes below 100, layer by layer (Ogg's supersingular layers Fi22, B\\Fi22, M\\B and the two size-halves of the
    ten non-supersingular primes below 100, exactly as the owner's 2026-11-01 directive lists them: 67 present, 71 once);
(2) the two congruence mechanisms, checked by brute force on all primes below 10^5:
    (a) bipartite: odd p + odd q = even square = 0 mod 4 forces {p, q} = {1, 3} mod 4, and 2 + q = odd square = 1 mod 8
        forces q = 7 mod 8, so A = {2} u {p = 1 mod 4}, B = {p = 3 mod 4} is a bipartition;
    (b) p = q = 1 mod 3 is never an edge (p + q = 2 mod 3 is a non-residue);
(3) null model for the owner's split of the ten non-supersingular primes into two edge-free pentads;
(4) isolated primes below 100 and their first prime partners; growth of the prime graph to 10^5;
(5) the prime square-sum path problem: a Hamiltonian path on the primes below N needs no isolated prime, at most two
    leaves and (bipartite) | |A| - |B| | <= 1, i.e. Chebyshev's bias pi(N;3,4) - pi(N;1,4) in {0, 1, 2}; a Hamiltonian
    cycle needs bias = 1 and min degree >= 2.  Candidates are searched with a node budget.  Standard library only."""
import math, sys, itertools
sys.setrecursionlimit(10000)
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
def sq_graph(P):
    Ps = set(P); adj = {p: set() for p in P}; top = max(P)
    for p in P:
        s = math.isqrt(p + 2)
        while s * s <= p + top:
            q = s * s - p
            if q > p and q in Ps: adj[p].add(q); adj[q].add(p)
            s += 1
    return adj
def components(adj):
    seen = set(); comps = []
    for v in adj:
        if v in seen: continue
        stack = [v]; comp = set()
        while stack:
            x = stack.pop()
            if x in comp: continue
            comp.add(x); stack.extend(adj[x] - comp)
        seen |= comp; comps.append(sorted(comp))
    return comps
def side(p): return 'A' if p == 2 or p % 4 == 1 else 'B'

print('(1) the 25 primes below 100')
P = primes_below(100); adj = sq_graph(P)
E = sorted((p, q) for p in P for q in adj[p] if p < q)
print('    %d edges; by square:' % len(E), {s: [e for e in E if sum(e) == s] for s in sorted({sum(e) for e in E})})
print('    isolated:', [p for p in P if not adj[p]], '; leaves:', [p for p in P if len(adj[p]) == 1],
      '; components:', [len(c) for c in components(adj)])
layers = {'Fi22': [2,3,5,7,11,13], 'B\\Fi22': [17,19,23,31,47], 'M\\B': [29,41,59,71], 'nonSS-1': [37,43,53,61,67], 'nonSS-2': [73,79,83,89,97]}
assert sorted(sum(layers.values(), [])) == P
for name, L in layers.items():
    Ls = set(L); inner = [(e, sum(e)) for e in E if e[0] in Ls and e[1] in Ls]
    print('    %-8s %-24s sides %s internal edges %s' % (name, L, ''.join(side(p) for p in L), inner))
A = [p for p in P if side(p) == 'A']; B = [p for p in P if side(p) == 'B']
print('    bipartition A = {2} u {1 mod 4} =', A, '(%d) ; B = {3 mod 4} =' % len(A), B, '(%d)' % len(B))

print('(2) congruence mechanisms, brute force on primes below 10^5')
P5 = primes_below(100000); adj5 = sq_graph(P5)
E5 = [(p, q) for p in P5 for q in adj5[p] if p < q]
bip = all(side(p) != side(q) for p, q in E5)
two = all(q % 8 == 7 for q in adj5[2])
mod3 = not any(p % 3 == 1 and q % 3 == 1 for p, q in E5)
print('    edges %d; bipartite A|B: %s; neighbours of 2 all = 7 mod 8: %s; no edge with p = q = 1 mod 3: %s' % (len(E5), bip, two, mod3))
assert bip and two and mod3

print('(3) null model: splits of the ten non-supersingular primes below 100 into two pentads with no internal edge on either side')
NS = [37,43,53,61,67,73,79,83,89,97]; NSs = set(NS)
innerNS = [e for e in E if e[0] in NSs and e[1] in NSs]
good = sum(1 for half in itertools.combinations(NS, 5) if 37 in half and not any(set(e) <= set(half) or not (set(e) & set(half)) for e in innerNS))
print('    internal edges among the ten: %s; edge-free splits: %d of %d' % (innerNS, good, math.comb(9, 4)))

print('(4) isolated primes below 100: first prime partner; growth')
for p in (37, 43, 67):
    s = math.isqrt(p) + 1
    while True:
        q = s * s - p
        if q > 1 and all(q % d for d in range(2, math.isqrt(q) + 1)): print('    %d: first partner %d (square %d^2 = %d)' % (p, q, s, s * s)); break
        s += 1
for N in (100, 1000, 10000, 100000):
    Pn = primes_below(N); a = sq_graph(Pn); comps = components(a); degs = [len(a[p]) for p in Pn]
    iso = [p for p in Pn if not a[p]]
    bias = sum(1 for p in Pn if p % 4 == 3) - sum(1 for p in Pn if p % 4 == 1)
    print('    N=%6d: primes %5d, edges %6d, isolated %d %s, components %d (largest %d), max degree %d, |A|-|B| = %d' % (
        N, len(Pn), sum(degs) // 2, len(iso), iso[:7], len(comps), max(len(c) for c in comps), max(degs), 1 - bias))

print('(5) the prime square-sum path problem: primes below N, N over the primes up to 2000 (N = next prime after the last vertex)')
def ham_path(adj, budget):
    """Warnsdorff-ordered DFS for a Hamiltonian path; returns (path or None, nodes, exhausted?)."""
    n = len(adj); best = [None]; nodes = [0]
    verts = sorted(adj, key=lambda v: len(adj[v]))
    starts = [v for v in verts if len(adj[v]) == 1] or verts
    def dfs(path, used):
        nodes[0] += 1
        if nodes[0] > budget: return 'budget'
        if len(path) == n: best[0] = list(path); return 'found'
        v = path[-1]
        cand = [w for w in adj[v] if w not in used]
        # a vertex of remaining degree 1 that is not reachable as the final vertex kills the branch
        ends = [w for w in cand if len(adj[w] - used) <= 1]
        if len(ends) > 1: return None
        if ends: cand = ends
        cand.sort(key=lambda w: len(adj[w] - used))
        for w in cand:
            path.append(w); used.add(w)
            r = dfs(path, used)
            if r: return r
            path.pop(); used.discard(w)
        return None
    for s in starts:
        r = dfs([s], {s})
        if r == 'found': return best[0], nodes[0], False
        if r == 'budget': return None, nodes[0], False
    return None, nodes[0], True
allp = primes_below(2100)
found_first = None; rows = []
for i in range(1, len(allp)):
    N = allp[i]; Pn = allp[:i]; a = sq_graph(Pn)
    if len(Pn) < 3 or any(not a[p] for p in Pn): continue
    leaves = [p for p in Pn if len(a[p]) == 1]
    nA = sum(1 for p in Pn if side(p) == 'A'); nB = len(Pn) - nA
    if len(leaves) > 2 or abs(nA - nB) > 1: continue
    path, nodes, exhausted = ham_path(a, 3_000_000)
    rows.append((N, len(Pn), len(leaves), nA - nB, 'PATH' if path else ('NONE' if exhausted else 'UNKNOWN'), nodes))
    if path and found_first is None:
        found_first = (N, path)
print('    candidates passing (no isolated, <= 2 leaves, | |A|-|B| | <= 1):')
for r in rows: print('      N=%5d (%3d primes, %d leaves, |A|-|B| = %+d): %s after %d nodes' % r)
if found_first:
    N, path = found_first
    assert all(math.isqrt(path[j] + path[j+1]) ** 2 == path[j] + path[j+1] for j in range(len(path) - 1)) and len(set(path)) == len(path)
    print('    first Hamiltonian path on the primes below N = %d:' % N, path)
else:
    print('    no Hamiltonian path found below 2100 within the budget')
