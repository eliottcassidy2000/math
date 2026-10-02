#!/usr/bin/env python3
"""procgen_tbij_20261001_lib.py -- library for the tbij lane (OPEN-Q-060, bijective form).

Conventions
  * vertices 0..n-1; a tournament is a tuple `out` of n bitmasks, out[i] = set of j with i -> j.
  * a graph is a tuple `adj` of n bitmasks (symmetric, no loops).
  * edges of K_n are indexed by EIDX[n][(i,j)] (i<j) in lexicographic order; an edge set is an int mask.
  * permutations are tuples p with p[i] = image of i.
  * the reference tournament T0 has i -> j for i < j.
All functions are self-contained (only nauty binaries gentourng, geng, labelg are called as subprocesses).
"""
import itertools
import math
import subprocess
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache

# ----------------------------------------------------------------------------------------------
# basic combinatorics
# ----------------------------------------------------------------------------------------------


def popcount(x):
    return bin(x).count('1')


@lru_cache(maxsize=None)
def edge_list(n):
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


@lru_cache(maxsize=None)
def edge_index(n):
    return {e: k for k, e in enumerate(edge_list(n))}


def eid(n, i, j):
    if i > j:
        i, j = j, i
    return edge_index(n)[(i, j)]


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def zee(mu):
    c = Counter(mu)
    z = 1
    for k, m in c.items():
        z *= k ** m * math.factorial(m)
    return z


def v2(x):
    return (x & -x).bit_length() - 1


def is_level_type(mu):
    return len({v2(l) for l in mu}) == 1


# ----------------------------------------------------------------------------------------------
# closed forms (independent implementations)
# ----------------------------------------------------------------------------------------------


def pair_orbits(mu):
    """number of orbits on unordered pairs of a permutation of cycle type mu"""
    o = sum(l // 2 for l in mu)
    for a in range(len(mu)):
        for b in range(a + 1, len(mu)):
            o += math.gcd(mu[a], mu[b])
    return o


def a049313_closed(n):
    """Babai-Cameron Thm 7.2 / THM-479: sum over level permutations."""
    tot = Fraction(0)
    for mu in partitions(n):
        if not is_level_type(mu):
            continue
        k = len(mu)
        o = pair_orbits(mu)
        if mu[0] % 2 == 1:  # all odd
            tot += Fraction(2 ** (o - k + 1), zee(mu))
        else:
            tot += Fraction(2 ** (o - k), zee(mu))
    assert tot.denominator == 1
    return int(tot)


def euler_fixed_count(mu):
    """#Euler graphs fixed by a permutation of type mu = |cycle space^g| = 2^(o - k + [some cycle odd])...
    computed directly as 2^(dim ker(g-1) on Z): by Brauer this equals #switching classes of graphs fixed.
    Mallows-Sloane/Robinson: 2^(o2 - k + 1) if some cycle is odd, 2^(o2 - k) if all even (n>=1)."""
    k = len(mu)
    o = pair_orbits(mu)
    if any(l % 2 == 1 for l in mu):
        return 2 ** (o - k + 1)
    return 2 ** (o - k)


def a002854_closed(n):
    tot = Fraction(0)
    for mu in partitions(n):
        tot += Fraction(euler_fixed_count(mu), zee(mu))
    assert tot.denominator == 1
    return int(tot)


A049313 = [None, 1, 1, 1, 2, 2, 6, 12, 79, 792, 19576, 886288, 75369960, 11856006240,
           3467430423264, 1893448825054528, 1938818712501985736]
A002854 = [None, 1, 1, 2, 3, 7, 16, 54, 243, 2038, 33120, 1182004, 87723296, 12886193064,
           3633057074584, 1944000150734320, 1967881448329407496]

# ----------------------------------------------------------------------------------------------
# nauty interface
# ----------------------------------------------------------------------------------------------


def gentourng(n):
    """all tournaments on n vertices up to isomorphism, as out-mask tuples"""
    if n == 1:
        return [(0,)]
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True, check=True).stdout
    res = []
    for line in out.split():
        res.append(parse_upper(line, n))
    return res


def parse_upper(s, n):
    """gentourng ascii: upper triangle row by row, '1' at (i,j), i<j, means i -> j"""
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1':
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            k += 1
    return tuple(out)


def geng(n, extra=()):
    if n == 1:
        return [(0,)]
    out = subprocess.run(['geng', '-q'] + list(extra) + [str(n)], capture_output=True, text=True,
                         check=True).stdout
    return [g6_decode(l) for l in out.split()]


def g6_decode(s):
    data = [ord(c) - 63 for c in s]
    n = data[0]
    bits = []
    for d in data[1:]:
        for b in range(5, -1, -1):
            bits.append((d >> b) & 1)
    adj = [0] * n
    k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            k += 1
    return tuple(adj)


def _pack6(bits):
    while len(bits) % 6:
        bits.append(0)
    s = []
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        s.append(chr(v + 63))
    return ''.join(s)


def _size6(n):
    """graph6 / digraph6 size field N(n)"""
    if n <= 62:
        return chr(n + 63)
    assert n <= 258047
    return '~' + ''.join(chr(((n >> s) & 63) + 63) for s in (12, 6, 0))


def g6_encode(adj):
    n = len(adj)
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append((adj[i] >> j) & 1)
    return _size6(n) + _pack6(bits)


def d6_encode(out):
    n = len(out)
    bits = []
    for i in range(n):
        for j in range(n):
            bits.append((out[i] >> j) & 1)
    return '&' + _size6(n) + _pack6(bits)


def labelg_batch(strings, digraph=False):
    """canonical forms (strings) for a list of graph6 / digraph6 strings via nauty labelg"""
    if not strings:
        return []
    inp = '\n'.join(strings) + '\n'
    args = ['labelg', '-q']
    out = subprocess.run(args, input=inp, capture_output=True, text=True, check=True).stdout.split()
    assert len(out) == len(strings)
    return out


# ----------------------------------------------------------------------------------------------
# tournaments
# ----------------------------------------------------------------------------------------------


def scores(out):
    return tuple(popcount(x) for x in out)


def switch(out, U):
    """reverse every arc with exactly one end in U (U a bitmask)"""
    n = len(out)
    full = (1 << n) - 1
    new = list(out)
    for i in range(n):
        if (U >> i) & 1:
            other = full & ~U
        else:
            other = U
        # arcs between i and `other` reverse: out-neighbours in `other` become in-neighbours and vice versa
        inn = full & ~out[i] & ~(1 << i)
        new[i] = (out[i] & ~other) | (inn & other)
    return tuple(new)


def apply_perm_t(out, p):
    """relabel: vertex i becomes p[i]; arc i->j becomes p[i]->p[j]"""
    n = len(out)
    new = [0] * n
    for i in range(n):
        m = out[i]
        img = 0
        while m:
            b = m & -m
            j = b.bit_length() - 1
            img |= 1 << p[j]
            m ^= b
        new[p[i]] = img
    return tuple(new)


def flip_vector(out, n=None):
    """edge mask of pairs {i<j} where the tournament has j -> i (i.e. differs from T0)"""
    n = len(out)
    x = 0
    for (i, j), k in edge_index(n).items():
        if (out[j] >> i) & 1:
            x |= 1 << k
    return x


def from_flip(x, n):
    out = [0] * n
    for (i, j), k in edge_index(n).items():
        if (x >> k) & 1:
            out[j] |= 1 << i
        else:
            out[i] |= 1 << j
    return tuple(out)


def cut_mask(W, n):
    m = 0
    for (i, j), k in edge_index(n).items():
        if ((W >> i) & 1) != ((W >> j) & 1):
            m |= 1 << k
    return m


def is_cut(x, n):
    """is the edge set x a cut delta(W)?  (W normalised to contain vertex 0 or not)"""
    s = [0] * n
    for j in range(1, n):
        s[j] = (x >> eid(n, 0, j)) & 1
    for (i, j), k in edge_index(n).items():
        if ((x >> k) & 1) != (s[i] ^ s[j]):
            return False
    return True


def mod4_euler(out):
    """out-degree == in-degree (mod 4) at every vertex  (only possible for n odd)"""
    n = len(out)
    return all((2 * popcount(o) - (n - 1)) % 4 == 0 for o in out)


def mod4_member(out):
    """for n odd: the unique member of the switching class of `out` with out == in (mod 4) everywhere"""
    n = len(out)
    assert n % 2 == 1
    target = ((n - 1) // 2) & 1
    p = [popcount(o) & 1 for o in out]
    W = 0
    for i in range(n):
        if p[i] != target:
            W |= 1 << i
    # switching at U flips parity exactly on the even one of U, U^c; W has even size (n odd, both weights
    # congruent to C(n,2) mod 2), so switch at U = W
    assert popcount(W) % 2 == 0
    T = switch(out, W)
    assert mod4_euler(T)
    return T


def skew_seidel(out):
    n = len(out)
    return [[0 if i == j else (1 if (out[i] >> j) & 1 else -1) for j in range(n)] for i in range(n)]


def cyclic_triangles(out):
    n = len(out)
    return math.comb(n, 3) - sum(math.comb(popcount(o), 2) for o in out)


# ----------------------------------------------------------------------------------------------
# graphs
# ----------------------------------------------------------------------------------------------


def degrees(adj):
    return tuple(popcount(a) for a in adj)


def is_euler(adj):
    return all(popcount(a) % 2 == 0 for a in adj)


def nedges(adj):
    return sum(popcount(a) for a in adj) // 2


def adj_from_mask(x, n):
    adj = [0] * n
    for (i, j), k in edge_index(n).items():
        if (x >> k) & 1:
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return tuple(adj)


def mask_from_adj(adj):
    n = len(adj)
    x = 0
    for (i, j), k in edge_index(n).items():
        if (adj[i] >> j) & 1:
            x |= 1 << k
    return x


def apply_perm_g(adj, p):
    n = len(adj)
    new = [0] * n
    for i in range(n):
        m = adj[i]
        img = 0
        while m:
            b = m & -m
            j = b.bit_length() - 1
            img |= 1 << p[j]
            m ^= b
        new[p[i]] = img
    return tuple(new)


def components(adj):
    n = len(adj)
    seen = 0
    comps = []
    for s in range(n):
        if (seen >> s) & 1:
            continue
        comp = 1 << s
        frontier = 1 << s
        while frontier:
            b = frontier & -frontier
            v = b.bit_length() - 1
            frontier ^= b
            new = adj[v] & ~comp
            comp |= new
            frontier |= new
        seen |= comp
        comps.append(comp)
    return comps


# ----------------------------------------------------------------------------------------------
# automorphism groups by backtracking (full element lists; fine for n <= 9)
# ----------------------------------------------------------------------------------------------


def _refine_colors(n, nbr_masks_list, init):
    """equitable-ish refinement: colour = (old colour, multiset of neighbour colours for each relation)"""
    col = list(init)
    while True:
        sig = []
        for v in range(n):
            parts = [col[v]]
            for nb in nbr_masks_list:
                parts.append(tuple(sorted(col[u] for u in range(n) if (nb[v] >> u) & 1)))
            sig.append(tuple(parts))
        keys = sorted(set(sig))
        newcol = [keys.index(s) for s in sig]
        if len(set(newcol)) == len(set(col)):
            return newcol
        col = newcol


def automorphisms_graph(adj, limit=None):
    """all automorphisms of an undirected graph (list of tuples)"""
    n = len(adj)
    col = _refine_colors(n, [adj], [0] * n)
    return _backtrack(n, col, lambda i, j, a, b: (((adj[i] >> j) & 1) == ((adj[a] >> b) & 1)), limit)


def automorphisms_tournament(out, limit=None):
    n = len(out)
    inn = [((1 << n) - 1) & ~out[i] & ~(1 << i) for i in range(n)]
    col = _refine_colors(n, [out, inn], [0] * n)
    return _backtrack(n, col, lambda i, j, a, b: (((out[i] >> j) & 1) == ((out[a] >> b) & 1)), limit)


def _backtrack(n, col, compat, limit):
    order = sorted(range(n), key=lambda v: (sum(1 for u in range(n) if col[u] == col[v]), v))
    res = []
    img = [-1] * n
    used = [False] * n

    def rec(t):
        if limit is not None and len(res) >= limit:
            return
        if t == n:
            res.append(tuple(img))
            return
        v = order[t]
        for a in range(n):
            if used[a] or col[a] != col[v]:
                continue
            ok = True
            for s in range(t):
                u = order[s]
                if not compat(v, u, a, img[u]) or not compat(u, v, img[u], a):
                    ok = False
                    break
            if ok:
                img[v] = a
                used[a] = True
                rec(t + 1)
                used[a] = False
                img[v] = -1

    rec(0)
    return res


def class_stabilizer(out):
    """all g in S_n with g(T) in the switching class of T  (g(T): arc i->j becomes g(i)->g(j)).
    Condition: the flip pattern d(i,j) = [T(g i, g j) != T(i, j)] is a cut, i.e. d(i,j) = s_i + s_j."""
    n = len(out)
    A = [[(out[i] >> j) & 1 for j in range(n)] for i in range(n)]
    res = []
    img = [-1] * n
    used = [False] * n
    s = [0] * n

    def rec(t):
        if t == n:
            res.append(tuple(img))
            return
        for a in range(n):
            if used[a]:
                continue
            # s[t] determined by pair (0, t) if t > 0
            if t == 0:
                st = 0
            else:
                st = s[0] ^ (A[img[0]][a] != A[0][t])
            ok = True
            for u in range(1, t):
                if (A[img[u]][a] != A[u][t]) != (s[u] ^ st):
                    ok = False
                    break
            if ok:
                img[t] = a
                used[a] = True
                s[t] = st
                rec(t + 1)
                used[a] = False
                img[t] = -1

    rec(0)
    return res


# ----------------------------------------------------------------------------------------------
# permutations
# ----------------------------------------------------------------------------------------------


def cycle_type(p):
    n = len(p)
    seen = [False] * n
    ct = []
    for i in range(n):
        if not seen[i]:
            l = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = p[j]
                l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))


def cycles(p):
    n = len(p)
    seen = [False] * n
    cs = []
    for i in range(n):
        if not seen[i]:
            c = []
            j = i
            while not seen[j]:
                seen[j] = True
                c.append(j)
                j = p[j]
            cs.append(c)
    return cs


def perm_sign(p):
    return (-1) ** sum(len(c) - 1 for c in cycles(p))


def compose(p, q):
    """(p o q)(i) = p[q[i]]"""
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    r = [0] * len(p)
    for i, a in enumerate(p):
        r[a] = i
    return tuple(r)


def eps_orient(adj, g):
    """(-1)^(number of edges of the graph whose T0-orientation g reverses), g in Aut(adj)"""
    n = len(adj)
    c = 0
    for i in range(n):
        m = adj[i] >> (i + 1)
        j = i + 1
        while m:
            if m & 1:
                if g[i] > g[j]:
                    c += 1
            m >>= 1
            j += 1
    return -1 if c % 2 else 1


def eps_cycle(adj, g):
    """cycle form: (-1)^(# even cycles of g, length 2k, whose antipodal pairs {x, g^k x} are edges)"""
    c = 0
    for cyc in cycles(g):
        L = len(cyc)
        if L % 2 == 0:
            x, y = cyc[0], cyc[L // 2]
            if (adj[x] >> y) & 1:
                c += 1
    return -1 if c % 2 else 1


def eps_arcs(adj, g):
    """sign of the permutation g induces on the 2|E| arcs (directed edges) of the graph"""
    n = len(adj)
    arcs = [(i, j) for i in range(n) for j in range(n) if (adj[i] >> j) & 1]
    pos = {a: k for k, a in enumerate(arcs)}
    perm = [pos[(g[i], g[j])] for (i, j) in arcs]
    return perm_sign(perm) if arcs else 1
