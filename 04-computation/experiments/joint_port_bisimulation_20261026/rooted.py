# Rooted two-port tournament states and the re-mark/complement dynamics of literal lines (rooted-lines conventions):
# a state h = [C; s, t] is a tournament with an ordered, jointly marked entrance s and exit t (s != t);
# its successors are [C^op xor E(P); s, t] over all Hamilton paths P of C from s to t ("legal re-marking with fixed
# ports, then complement").  mac-mini-2026-10-25-pyth (2026-10-26).
import itertools, subprocess
from functools import lru_cache

def tournaments(n):
    """all tournament classes on n vertices from nauty's gentourng (upper-triangle bit strings, i<j, bit 1 = i->j)"""
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    res = []
    for s in out:
        adj = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == '1': adj[i] |= 1 << j
                else: adj[j] |= 1 << i
                k += 1
        res.append(tuple(adj))
    return res

def arc(adj, u, v): return (adj[u] >> v) & 1

def relabel(adj, perm):
    """perm[old] = new"""
    n = len(adj); new = [0] * n
    for u in range(n):
        row = adj[u]
        for v in range(n):
            if (row >> v) & 1: new[perm[u]] |= 1 << perm[v]
    return tuple(new)

_perm_cache = {}
def middle_perms(n):
    if n not in _perm_cache:
        _perm_cache[n] = list(itertools.permutations(range(1, n - 1)))
    return _perm_cache[n]

def canon_rooted(adj, s, t):
    """canonical form of (adj; s, t): move s to 0 and t to n-1, then minimise over permutations of the rest"""
    n = len(adj)
    rest = [v for v in range(n) if v not in (s, t)]
    best = None
    for mp in middle_perms(n):
        perm = [0] * n; perm[s] = 0; perm[t] = n - 1
        for v, w in zip(rest, mp): perm[v] = w
        r = relabel(adj, perm)
        if best is None or r < best: best = r
    return best

def canon_free(adj):
    n = len(adj)
    best = None
    for perm in itertools.permutations(range(n)):
        r = relabel(adj, perm)
        if best is None or r < best: best = r
    return best

def ham_paths(adj, s=None, t=None):
    n = len(adj); full = (1 << n) - 1
    out = []
    starts = [s] if s is not None else range(n)
    def dfs(path, used):
        if used == full:
            if t is None or path[-1] == t: out.append(tuple(path))
            return
        u = path[-1]
        cand = adj[u] & ~used
        while cand:
            v = (cand & -cand).bit_length() - 1; cand &= cand - 1
            if t is not None and v == t and (used | (1 << v)) != full: continue
            path.append(v); dfs(path, used | (1 << v)); path.pop()
    for a in starts:
        dfs([a], 1 << a)
    return out

def partner(adj, P):
    """C^op xor E(P): keep path arcs, reverse every other arc"""
    n = len(adj)
    onpath = set()
    for a, b in zip(P, P[1:]): onpath.add((a, b))
    new = [0] * n
    for u in range(n):
        for v in range(n):
            if u == v: continue
            if (u, v) in onpath: new[u] |= 1 << v
            elif (v, u) in onpath: pass
            elif arc(adj, v, u): new[u] |= 1 << v      # reversed
    return tuple(new)

def is_strong(adj):
    n = len(adj)
    # strong iff from vertex 0 every vertex reachable and 0 reachable from every vertex
    def reach(rows):
        seen = 1; frontier = 1
        while frontier:
            nxt = 0
            f = frontier
            while f:
                v = (f & -f).bit_length() - 1; f &= f - 1
                nxt |= rows[v]
            frontier = nxt & ~seen; seen |= nxt
        return seen
    full = (1 << n) - 1
    rev = [0] * n
    for u in range(n):
        for v in range(n):
            if arc(adj, u, v): rev[v] |= 1 << u
    return reach(adj) == full and reach(rev) == full

def H(adj): return len(ham_paths(adj))
