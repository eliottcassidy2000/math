"""natural_matchings_switching_classes_20261002_lib.py -- helpers for the OPEN-Q-060 natural-matching runner.

Tournaments are adjacency lists of lists (A[i][j] = 1 iff i -> j); graphs are edge lists of pairs (i, j)
with i < j.  nauty is called through subprocess (geng, gentourng, labelg, dreadnaut, with or without the
'nauty-' prefix).
"""
import itertools
import re
import shutil
import subprocess

from sympy import factorint, isprime, primitive_root

A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536]
A002854 = [1, 1, 1, 2, 3, 7, 16, 54, 243, 2038, 33120, 1182004, 87723296, 12886193064]
A049313 = [1, 1, 1, 1, 2, 2, 6, 12, 79, 792, 19576, 886288, 75369960, 11856006240]
# index n; A002854 / A049313 values for n = 11..13 as recorded in
# 05-knowledge/results/switching_classes_level_burnside_cbx2.out


def nauty(name):
    for cand in (name, 'nauty-' + name):
        if shutil.which(cand):
            return cand
    raise RuntimeError('nauty program %s not found (tried %s and nauty-%s)' % (name, name, name))


def run(cmd, inp=None):
    return subprocess.run(cmd, input=inp, capture_output=True, text=True, check=True).stdout


# ----------------------------------------------------------------------------------------------- formats
def g6_of_edges(E, n):
    S = {(min(e), max(e)) for e in E}
    bits = [1 if (i, j) in S else 0 for j in range(1, n) for i in range(j)]
    while len(bits) % 6:
        bits.append(0)
    s = chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        s += chr(v + 63)
    return s


def edges_of_g6(s):
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        bits += [(v >> (5 - i)) & 1 for i in range(6)]
    E = []
    k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                E.append((i, j))
            k += 1
    return n, E


def d6(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6:
        bits.append(0)
    s = '&' + (chr(n + 63) if n <= 62 else '~' + ''.join(chr(((n >> t) & 63) + 63) for t in (12, 6, 0)))
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        s += chr(v + 63)
    return s


def canon_many(strings):
    """nauty canonical labels (labelg) of graph6 / digraph6 strings, in order."""
    if not strings:
        return []
    out = run([nauty('labelg'), '-q'], '\n'.join(strings) + '\n').split()
    assert len(out) == len(strings), (len(out), len(strings))
    return out


# ------------------------------------------------------------------------------------------- tournaments
def gentourng(n):
    """all tournaments on n vertices up to isomorphism (gentourng upper-triangle format)."""
    res = []
    for s in run([nauty('gentourng'), '-q', str(n)]).split():
        A = [[0] * n for _ in range(n)]
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == '1':
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                k += 1
        res.append(A)
    return res


def scores(A):
    return [sum(r) for r in A]


def switch(A, U):
    n = len(A)
    U = set(U)
    return [[(A[j][i] if ((i in U) != (j in U)) else A[i][j]) if i != j else 0 for j in range(n)]
            for i in range(n)]


def is_tournament(A):
    n = len(A)
    return all(A[i][i] == 0 for i in range(n)) and all(A[i][j] + A[j][i] == 1
                                                       for i in range(n) for j in range(i + 1, n))


def is_aut(A, g):
    n = len(A)
    return sorted(g) == list(range(n)) and all(A[g[i]][g[j]] == A[i][j] for i in range(n) for j in range(n))


def aut_tournament(A):
    """all automorphisms (backtracking; small n)."""
    n = len(A)
    sc = scores(A)
    res = []
    g = [-1] * n
    used = [False] * n

    def bt(k):
        if k == n:
            res.append(tuple(g))
            return
        for v in range(n):
            if used[v] or sc[v] != sc[k]:
                continue
            if all(A[a][k] == A[g[a]][v] for a in range(k)):
                g[k] = v
                used[v] = True
                bt(k + 1)
                used[v] = False
                g[k] = -1
    bt(0)
    return res


def class_stabilizer(A):
    """all g in S_n mapping the switching class of A to itself (backtracking on the switching
    invariant pi(x,y,z) = parity of arcs of A agreeing with the cyclic order x -> y -> z -> x)."""
    n = len(A)

    def val(x, y, z):
        return (A[x][y] + A[y][z] + A[z][x]) & 1
    res = []
    g = [-1] * n
    used = [False] * n

    def bt(k):
        if k == n:
            res.append(tuple(g))
            return
        for v in range(n):
            if used[v]:
                continue
            ok = True
            for a in range(k):
                for b in range(a + 1, k):
                    if val(a, b, k) != val(g[a], g[b], v):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                g[k] = v
                used[v] = True
                bt(k + 1)
                used[v] = False
                g[k] = -1
    bt(0)
    return res


def descendant_d6s(A):
    """for each v: switch the in-neighbours of v (v becomes a source), delete v; digraph6 strings."""
    n = len(A)
    out = []
    for v in range(n):
        U = {x for x in range(n) if x != v and A[x][v]}
        B = switch(A, U)
        assert all(B[v][x] for x in range(n) if x != v)
        idx = [x for x in range(n) if x != v]
        out.append(d6([[B[i][j] for j in idx] for i in idx]))
    return out


def class_key(A):
    """sorted canonical forms of the n descendants: an isomorphism invariant of the switching class
    (two classes are isomorphic iff their descendant multisets share an element, iff the keys agree)."""
    return tuple(sorted(canon_many(descendant_d6s(A))))


# ------------------------------------------------------------------------------------------------ graphs
def pair_orbits(n, gens):
    """orbits of the group <gens> on unordered pairs (union-find)."""
    par = {}

    def f(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for p in pairs:
        f(p)
    for g in gens:
        for (i, j) in pairs:
            a, b = g[i], g[j]
            q = (a, b) if a < b else (b, a)
            ra, rb = f((i, j)), f(q)
            if ra != rb:
                par[ra] = rb
    orb = {}
    for p in pairs:
        orb.setdefault(f(p), []).append(p)
    return list(orb.values())


def dreadnaut_gens(graphs):
    """graphs: list of (n, edge list).  Automorphism-group generators of each graph (dreadnaut)."""
    if not graphs:
        return []
    lines = []
    for n, E in graphs:
        adj = [[] for _ in range(n)]
        for (i, j) in E:
            adj[i].append(j)
            adj[j].append(i)
        lines.append('n=%d g' % n)
        for i in range(n):
            lines.append(' '.join(map(str, adj[i])) + (';' if i < n - 1 else '.'))
        lines.append('+a x')
    lines.append('q')
    out = run([nauty('dreadnaut')], '\n'.join(lines) + '\n')
    res, cur, gen = [], [], None
    for line in out.splitlines():
        if line.startswith('('):
            if gen is not None:
                cur.append(gen)
            gen = line.strip()
        elif gen is not None and line[:1] == ' ' and line.strip()[:1] in '0123456789(':
            gen += ' ' + line.strip()
        else:
            if gen is not None:
                cur.append(gen)
                gen = None
            if 'grpsize' in line:
                res.append(cur)
                cur = []
    assert len(res) == len(graphs), (len(res), len(graphs))
    perms = []
    for (n, E), gl in zip(graphs, res):
        P = []
        for gtxt in gl:
            g = list(range(n))
            for cyc in re.findall(r'\(([^)]*)\)', gtxt):
                xs = list(map(int, cyc.split()))
                for k in range(len(xs)):
                    g[xs[k]] = xs[(k + 1) % len(xs)]
            P.append(g)
        ES = {(min(e), max(e)) for e in E}
        for g in P:
            assert all((min(g[i], g[j]), max(g[i], g[j])) in ES for (i, j) in E)
        perms.append(P)
    return perms


def eps(E, g):
    """(-1)^(number of edges {i<j} with g(i) > g(j)): the orientation sign of g on E."""
    c = sum(1 for (i, j) in E if g[i] > g[j])
    return -1 if c & 1 else 1


def twisted_flags(graphs):
    """for each (n, E): True iff some automorphism reverses an odd number of edges (eps is a
    homomorphism on Aut, so it suffices to test generators)."""
    return [any(eps(E, g) == -1 for g in gl) for (n, E), gl in zip(graphs, dreadnaut_gens(graphs))]


def is_euler(E, n):
    deg = [0] * n
    for (i, j) in E:
        deg[i] += 1
        deg[j] += 1
    return all(d % 2 == 0 for d in deg)


def all_graphs(n):
    return [edges_of_g6(s)[1] for s in run([nauty('geng'), '-q', str(n)]).split()]


# ----------------------------------------------------------------------------------------- rigid blocks
def paley(p):
    sq = {(x * x) % p for x in range(1, p)}
    A = [[1 if (j - i) % p in sq else 0 for j in range(p)] for i in range(p)]
    g = primitive_root(p)
    return A, [[(x + 1) % p for x in range(p)], [(g * g * x) % p for x in range(p)]]


def rquart(p):
    """p = 5 mod 8: arcs x -> y iff y - x in M u gM (M = fourth powers); H = <x+1, g^4 x>."""
    g = primitive_root(p)
    M = {pow(g, 4 * k, p) for k in range((p - 1) // 4)}
    D = M | {(g * x) % p for x in M}
    A = [[1 if (j - i) % p in D else 0 for j in range(p)] for i in range(p)]
    return A, [[(x + 1) % p for x in range(p)], [(pow(g, 4, p) * x) % p for x in range(p)]]


def single():
    return [[0]], []


def lex(TA, TB):
    """A[B]: vertex (a, b) -> a * |B| + b; generators of H_B wr H_A."""
    (A, ga), (B, gb) = TA, TB
    na, nb = len(A), len(B)
    n = na * nb
    C = [[0] * n for _ in range(n)]
    for a in range(na):
        for b in range(nb):
            for a2 in range(na):
                for b2 in range(nb):
                    if a != a2:
                        C[a * nb + b][a2 * nb + b2] = A[a][a2]
                    elif b != b2:
                        C[a * nb + b][a2 * nb + b2] = B[b][b2]
    gens = [[s[a] * nb + b for a in range(na) for b in range(nb)] for s in ga]
    gens += [[(t[b] if a == 0 else b) + a * nb for a in range(na) for b in range(nb)] for t in gb]
    return C, gens


def prime_block(p):
    if p % 4 == 3:
        return paley(p)
    if p % 8 == 5:
        return rquart(p)
    raise ValueError(p)


def in_sigma(w):
    return w % 2 == 1 and all(p % 8 != 1 for p in factorint(w))


def rigid(w, order=None):
    """rigid block of order w in Sigma: lexicographic product of prime blocks (prime factors in the
    given order, default increasing)."""
    if w == 1:
        return single()
    ps = order if order is not None else [p for p, e in sorted(factorint(w).items()) for _ in range(e)]
    T = prime_block(ps[0])
    for p in ps[1:]:
        T = lex(T, prime_block(p))
    return T


def compose(blocks, Q):
    """Q[B_1, ..., B_m] with generators of H_1 x ... x H_m."""
    sizes = [len(B) for B, _ in blocks]
    off = [sum(sizes[:i]) for i in range(len(sizes))]
    n = sum(sizes)
    C = [[0] * n for _ in range(n)]
    gens = []
    for i, (B, gb) in enumerate(blocks):
        for x in range(sizes[i]):
            for j in range(len(blocks)):
                for y in range(sizes[j]):
                    if i == j:
                        if x != y:
                            C[off[i] + x][off[j] + y] = B[x][y]
                    else:
                        C[off[i] + x][off[j] + y] = Q[i][j]
        for t in gb:
            g = list(range(n))
            for x in range(sizes[i]):
                g[off[i] + x] = off[i] + t[x]
            gens.append(g)
    return C, gens


def transitive(m):
    return [[1 if i < j else 0 for j in range(m)] for i in range(m)]


def is_transitive_group(n, gens):
    seen = {0}
    fr = [0]
    while fr:
        x = fr.pop()
        for g in gens:
            if g[x] not in seen:
                seen.add(g[x])
                fr.append(g[x])
    return len(seen) == n


def invariant_graphs(n, gens, euler_only=False, limit=18):
    """all nonempty <gens>-invariant graphs (optionally only Euler ones), as edge lists."""
    orbs = pair_orbits(n, gens)
    k = len(orbs)
    assert k <= limit, k
    res = []
    for mask in range(1, 1 << k):
        E = [p for i in range(k) if mask >> i & 1 for p in orbs[i]]
        if euler_only and not is_euler(E, n):
            continue
        res.append(E)
    return res, k


def representations(n):
    """the forced-class constructions of Theorem 6 (all two-block splits for even n)."""
    reps = []
    if n % 2 == 0:
        for w1 in range(1, n // 2 + 1, 2):
            if in_sigma(w1) and in_sigma(n - w1):
                reps.append((w1, n - w1))
    else:
        for w in range(1, (n - 1) // 2 + 1, 2):
            if in_sigma(w) and in_sigma(n - 2 * w):
                reps.append((w, w, n - 2 * w))
        if isprime(n) and in_sigma(n):
            reps.append((n,))
    return reps


def build(rep):
    if len(rep) == 1:
        return rigid(rep[0])
    return compose([rigid(w) for w in rep], transitive(len(rep)))


def prime_tournament(A):
    """True iff A has no module other than singletons and V (closure of every pair is V)."""
    n = len(A)
    for x in range(n):
        for y in range(x + 1, n):
            S = {x, y}
            changed = True
            while changed:
                changed = False
                for z in range(n):
                    if z not in S and len({A[z][a] for a in S}) > 1:
                        S.add(z)
                        changed = True
            if len(S) < n:
                return False
    return True


# --------------------------------------------------------------------------------------- matching
def max_matching(adj, nright):
    """adj[u] = list of right vertices; returns (size, match_r) -- augmenting paths with a greedy
    free-vertex scan first."""
    match_r = {}

    def aug(u, seen):
        for v in adj[u]:
            if v not in match_r:
                match_r[v] = u
                return True
        for v in adj[u]:
            if v in seen:
                continue
            seen.add(v)
            if aug(match_r[v], seen):
                match_r[v] = u
                return True
        return False
    for u in sorted(range(len(adj)), key=lambda u: len(adj[u])):
        aug(u, set())
    return len(match_r), match_r


def hall_violator(adj, match_r):
    matched = set(match_r.values())
    Ls = {u for u in range(len(adj)) if u not in matched}
    Rs = set()
    fr = list(Ls)
    while fr:
        u = fr.pop()
        for v in adj[u]:
            if v not in Rs:
                Rs.add(v)
                w = match_r.get(v)
                if w is not None and w not in Ls:
                    Ls.add(w)
                    fr.append(w)
    return Ls, Rs


def perms_all(n):
    return itertools.permutations(range(n))
