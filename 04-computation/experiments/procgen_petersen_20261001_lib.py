#!/usr/bin/env python3
"""Helpers for procgen_petersen_20261001_run.py (Petersen/Heawood families in Paley coordinates,
arc-HP parities, the Gale tournament of linear K6, and the Collatz-side parity-bridge checks).

Only stdout is used by the runner; nothing here writes files except the compiled HP engine,
which is built into a scratch directory given by the caller."""
import itertools
import math
import os
import subprocess
from fractions import Fraction

import networkx as nx

# ----------------------------------------------------------------------------- graph encodings

def g6(G):
    """graph6 string of a simple graph (own encoder; vertices taken in sorted order)."""
    nodes = sorted(G.nodes())
    n = len(nodes)
    assert n < 63
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if G.has_edge(nodes[i], nodes[j]) else 0)
    while len(bits) % 6:
        bits.append(0)
    s = chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = 2 * v + b
        s += chr(v + 63)
    return s


def digraph6(n, arcs):
    """digraph6 string (nauty format, '&' prefix)."""
    A = set(arcs)
    bits = [1 if (i, j) in A else 0 for i in range(n) for j in range(n)]
    while len(bits) % 6:
        bits.append(0)
    s = "&" + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = 2 * v + b
        s += chr(v + 63)
    return s


def labelg(strings):
    """canonical forms by nauty's labelg (graph6 or digraph6 input)."""
    if not strings:
        return []
    data = "\n".join(strings) + "\n"
    out = subprocess.run(["labelg", "-q"], input=data.encode(), capture_output=True, check=True).stdout.decode().split()
    assert len(out) == len(strings)
    return out


def canon(G):
    return labelg([g6(nx.convert_node_labels_to_integers(G))])[0]


def canon_many(graphs):
    return labelg([g6(nx.convert_node_labels_to_integers(G)) for G in graphs])

# ----------------------------------------------------------------------------- Delta-Y moves

def triangles(G):
    out = []
    for a, b, c in itertools.combinations(sorted(G.nodes()), 3):
        if G.has_edge(a, b) and G.has_edge(b, c) and G.has_edge(a, c):
            out.append((a, b, c))
    return out


def delta_y(G, tri, newv):
    a, b, c = tri
    H = G.copy()
    H.remove_edges_from([(a, b), (b, c), (a, c)])
    H.add_node(newv)
    H.add_edges_from([(newv, a), (newv, b), (newv, c)])
    return H


def y_delta(G, v):
    """returns None if the move would create a multiple edge"""
    a, b, c = sorted(G.neighbors(v))
    if G.has_edge(a, b) or G.has_edge(b, c) or G.has_edge(a, c):
        return None
    H = G.copy()
    H.remove_node(v)
    H.add_edges_from([(a, b), (b, c), (a, c)])
    return H


def closure(G0):
    """Delta-Y / Y-Delta closure. Returns (family dict canon->graph, arrows {(from,to,kind)}, #multi-edge skips)."""
    G0 = nx.convert_node_labels_to_integers(G0)
    c0 = canon(G0)
    fam = {c0: G0}
    arrows = set()
    queue = [c0]
    skips = 0
    while queue:
        c = queue.pop()
        G = fam[c]
        kids, kinds = [], []
        nv = max(G.nodes()) + 1
        for tri in triangles(G):
            kids.append(delta_y(G, tri, nv))
            kinds.append('DY')
        for v in list(G.nodes()):
            if G.degree(v) == 3:
                H = y_delta(G, v)
                if H is None:
                    skips += 1
                    continue
                kids.append(H)
                kinds.append('YD')
        if kids:
            for H, k, ch in zip(kids, kinds, canon_many(kids)):
                arrows.add((c, ch, k))
                if ch not in fam:
                    fam[ch] = nx.convert_node_labels_to_integers(H)
                    queue.append(ch)
    return fam, arrows, skips


def dy_descendants(c0, arrows):
    seen = {c0}
    st = [c0]
    while st:
        c = st.pop()
        for (a, b, k) in arrows:
            if a == c and k == 'DY' and b not in seen:
                seen.add(b)
                st.append(b)
    return seen


def nx_class_count(graphs):
    """independent isomorphism-class count (WL hash buckets + VF2)."""
    reps = []
    for G in graphs:
        h = nx.weisfeiler_lehman_graph_hash(G, iterations=4)
        if not any(h2 == h and nx.is_isomorphic(G, G2) for h2, G2 in reps):
            reps.append((h, G))
    return len(reps)


def aut_order(G):
    return sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())


def describe(G):
    degs = tuple(sorted((d for _, d in G.degree()), reverse=True))
    try:
        girth = nx.girth(G)
    except Exception:
        girth = None
    return dict(n=G.number_of_nodes(), m=G.number_of_edges(), degs=degs, tri=len(triangles(G)),
                bip=nx.is_bipartite(G), girth=girth, aut=aut_order(G))


def smooth(G):
    """suppress degree-2 vertices (unless a multi-edge would arise) and delete degree <= 1 vertices, repeatedly."""
    H = G.copy()
    changed = True
    while changed:
        changed = False
        for v in sorted(H.nodes(), key=str):
            d = H.degree(v)
            if d <= 1:
                H.remove_node(v)
                changed = True
                break
            if d == 2:
                a, b = list(H.neighbors(v))
                if not H.has_edge(a, b):
                    H.remove_node(v)
                    H.add_edge(a, b)
                    changed = True
                    break
    return H


def count_hc_undirected(G):
    nodes = sorted(G.nodes())
    n = len(nodes)
    first = nodes[0]
    cnt = 0
    for p in itertools.permutations(nodes[1:]):
        if p[0] > p[-1]:
            continue
        cyc = (first,) + p
        if all(G.has_edge(cyc[i], cyc[(i + 1) % n]) for i in range(n)):
            cnt += 1
    return cnt

# ----------------------------------------------------------------------------- Fano / Paley coordinates

D = (1, 2, 4)
QR7 = frozenset(D)
NQR7 = frozenset({3, 5, 6})
LINES = {t: frozenset((d + t) % 7 for d in D) for t in range(7)}


def paley(q):
    QR = {(x * x) % q for x in range(1, q)}
    return [[1 if i != j and (j - i) % q in QR else 0 for j in range(q)] for i in range(q)]


def collineations():
    lines = set(LINES.values())
    return [p for p in itertools.permutations(range(7))
            if all(frozenset(p[x] for x in ln) in lines for ln in lines)]


def dy_on_lines(G, lines):
    H = G.copy()
    nv = max(H.nodes()) + 1
    for ln in lines:
        a, b, c = sorted(ln)
        assert H.has_edge(a, b) and H.has_edge(b, c) and H.has_edge(a, c)
        H = delta_y(H, (a, b, c), nv)
        nv += 1
    return H


def split_digraph(A):
    """bipartite double ('split') of a digraph: (s,'out') ~ (t,'in') for every arc s -> t."""
    S = nx.Graph()
    n = len(A)
    for s in range(n):
        for t in range(n):
            if A[s][t]:
                S.add_edge((s, 'out'), (t, 'in'))
    return S

# ----------------------------------------------------------------------------- HP engines

def build_engine(src, workdir):
    os.makedirs(workdir, exist_ok=True)
    exe = os.path.join(workdir, "procgen_petersen_hp")
    subprocess.run(["gcc", "-O2", "-o", exe, src], check=True)
    return exe


def engine_all(exe, n, edges):
    inp = "%d %d\n" % (n, len(edges)) + "".join("%d %d\n" % e for e in edges) + "A\n"
    out = subprocess.run([exe], input=inp.encode(), capture_output=True, check=True).stdout.decode()
    stats, allodd = {}, []
    for line in out.strip().split("\n"):
        w = line.split()
        if w[0] == "ALLODD":
            allodd.append((int(w[1]), int(w[2])))
        elif w[0] == "TOTAL":
            for k in range(0, len(w), 2):
                stats[w[k]] = int(w[k + 1])
        elif w[0] == "NODD_HIST":
            stats["NODD_HIST"] = list(map(int, w[1:]))
    return stats, allodd


def engine_list(exe, n, edges, masks):
    inp = "%d %d\n" % (n, len(edges)) + "".join("%d %d\n" % e for e in edges) + "L %d\n" % len(masks)
    inp += "".join("%d\n" % x for x in masks)
    out = subprocess.run([exe], input=inp.encode(), capture_output=True, check=True).stdout.decode()
    res = []
    for line in out.strip().split("\n"):
        w = list(map(int, line.split()))
        res.append((w[1], w[2:]))
    return res


def py_hp_counts(n, arcs):
    """independent pure-Python count: returns (H, {arc: #directed HPs through it}) by DP over (set, end)
    combined with an explicit path-extension recursion for the arc counts (different method from the C engine)."""
    out = [[] for _ in range(n)]
    for (a, b) in arcs:
        out[a].append(b)
    cnt = {a: 0 for a in arcs}
    H = 0
    # depth-first enumeration of all directed Hamiltonian paths (n <= 10 here)
    path = []
    used = [False] * n

    def rec(v, depth):
        nonlocal H
        if depth == n:
            H += 1
            for i in range(n - 1):
                cnt[(path[i], path[i + 1])] += 1
            return
        for w in out[v]:
            if not used[w]:
                used[w] = True
                path.append(w)
                rec(w, depth + 1)
                path.pop()
                used[w] = False
    for s in range(n):
        used[s] = True
        path.append(s)
        rec(s, 1)
        path.pop()
        used[s] = False
    return H, cnt


def orient(edges, mask):
    return [(v, u) if (mask >> k) & 1 else (u, v) for k, (u, v) in enumerate(edges)]

# ----------------------------------------------------------------------------- tournaments (n = 6)

def parse_tourn(s, n):
    """gentourng's default output: upper-triangle row-wise 0/1 string, '1' means i -> j for i < j."""
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1':
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def tour_arcs(A):
    n = len(A)
    return [(i, j) for i in range(n) for j in range(n) if A[i][j]]


def is_locally_transitive(A):
    n = len(A)
    for v in range(n):
        for nb in ([w for w in range(n) if A[v][w]], [w for w in range(n) if A[w][v]]):
            for a, b, c in itertools.combinations(nb, 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[b][a] and A[c][b] and A[a][c]):
                    return False
    return True


def is_transitive(A):
    return sorted(sum(r) for r in A) == list(range(len(A)))


def switch(A, W):
    n = len(A)
    return [[(A[i][j] if ((i in W) == (j in W)) else 1 - A[i][j]) if i != j else 0 for j in range(n)] for i in range(n)]

# ----------------------------------------------------------------------------- linear K6 (exact)

def det3(u, v, w):
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]))


def orient4(a, b, c, d):
    return det3([b[i] - a[i] for i in range(3)], [c[i] - a[i] for i in range(3)], [d[i] - a[i] for i in range(3)])


def sgn(x):
    return (x > 0) - (x < 0)


def seg_pierces(P, d, e, a, b, c):
    """signed piercing of the open segment de through the triangle abc (exact; general position asserted)."""
    s1 = sgn(orient4(P[a], P[b], P[c], P[d]))
    s2 = sgn(orient4(P[a], P[b], P[c], P[e]))
    assert s1 and s2
    if s1 == s2:
        return 0
    t1 = sgn(orient4(P[d], P[e], P[a], P[b]))
    t2 = sgn(orient4(P[d], P[e], P[b], P[c]))
    t3 = sgn(orient4(P[d], P[e], P[c], P[a]))
    assert t1 and t2 and t3
    return s2 if t1 == t2 == t3 else 0


def lk_triangles(P, A, B):
    a, b, c = A
    return sum(seg_pierces(P, B[i], B[(i + 1) % 3], a, b, c) for i in range(3))


def positive_weights(g):
    """positive rational w with sum w_i g_i = 0 (exists iff the configuration is totally cyclic)."""
    for free in itertools.combinations(range(6), 4):
        rest = [i for i in range(6) if i not in free]
        for base in itertools.product((1, 2, 3), repeat=4):
            sx = sum(base[t] * g[free[t]][0] for t in range(4))
            sy = sum(base[t] * g[free[t]][1] for t in range(4))
            a, b = g[rest[0]]
            c, d = g[rest[1]]
            det = a * d - b * c
            if det == 0:
                continue
            w0 = Fraction(-sx * d + sy * c, det)
            w1 = Fraction(-a * sy + b * sx, det)
            if w0 > 0 and w1 > 0:
                w = [Fraction(0)] * 6
                for t in range(4):
                    w[free[t]] = Fraction(base[t])
                w[rest[0]] = w0
                w[rest[1]] = w1
                return w
    raise RuntimeError("not totally cyclic")


def affine_dual(gw):
    """six points in Q^3 whose space of affine dependences is spanned by the coordinate rows of gw (sum gw = 0)."""
    M = [[Fraction(gw[i][r]) for i in range(6)] for r in range(2)]
    a, b = M[0][4], M[0][5]
    c, d = M[1][4], M[1][5]
    det = a * d - b * c
    assert det != 0
    basis = []
    for k in range(4):
        x = [Fraction(0)] * 6
        x[k] = Fraction(1)
        r0, r1 = -M[0][k], -M[1][k]
        x[4] = (d * r0 - b * r1) / det
        x[5] = (-c * r0 + a * r1) / det
        basis.append(x)
    ones = [Fraction(1)] * 6
    for drop in range(4):
        rest = [basis[i] for i in range(4) if i != drop]
        if rank([ones] + rest) == 4:
            return [(rest[0][i], rest[1][i], rest[2][i]) for i in range(6)]
    raise RuntimeError


def rank(rows):
    Mm = [r[:] for r in rows]
    rk = 0
    ncol = len(Mm[0])
    for col in range(ncol):
        piv = next((r for r in range(rk, len(Mm)) if Mm[r][col] != 0), None)
        if piv is None:
            continue
        Mm[rk], Mm[piv] = Mm[piv], Mm[rk]
        for r in range(len(Mm)):
            if r != rk and Mm[r][col] != 0:
                f = Mm[r][col] / Mm[rk][col]
                Mm[r] = [Mm[r][i] - f * Mm[rk][i] for i in range(ncol)]
        rk += 1
    return rk


def chirotope(P):
    """signs of orient4 on increasing 4-tuples (exact)."""
    return {q: sgn(orient4(*(P[i] for i in q))) for q in itertools.combinations(range(6), 4)}


def o4(chi, a, b, c, d):
    t = [a, b, c, d]
    inv = sum(1 for i in range(4) for j in range(i + 1, 4) if t[i] > t[j])
    s = chi[tuple(sorted(t))]
    return -s if inv % 2 else s


def seg_pierces_chi(chi, d, e, a, b, c):
    """same as seg_pierces, from the chirotope only."""
    s1 = o4(chi, a, b, c, d)
    s2 = o4(chi, a, b, c, e)
    assert s1 and s2
    if s1 == s2:
        return 0
    t1 = o4(chi, d, e, a, b)
    t2 = o4(chi, d, e, b, c)
    t3 = o4(chi, d, e, c, a)
    return s2 if t1 == t2 == t3 else 0


TRIANGLE_PAIRS = []
for _A in itertools.combinations(range(6), 3):
    _B = tuple(sorted(set(range(6)) - set(_A)))
    if _A < _B:
        TRIANGLE_PAIRS.append((_A, _B))


def pierce_rule(T, A, B):
    """number of b in B that are 'uniform-opposite' to A in the tournament T."""
    cnt = 0
    for b in B:
        rest = [x for x in B if x != b]
        if all(T[b][a] for a in A) and all(T[x][b] for x in rest):
            cnt += 1
        elif all(T[a][b] for a in A) and all(T[b][x] for x in rest):
            cnt += 1
    return cnt

# ----------------------------------------------------------------------------- PL linking (floating point)

def lk_polygons(P1, P2):
    """linking number of two closed polygons in R^3 (numpy arrays), via a projection to the xy-plane."""
    s = 0.0
    n1, n2 = len(P1), len(P2)
    for i in range(n1):
        a, b = P1[i], P1[(i + 1) % n1]
        for j in range(n2):
            c, d = P2[j], P2[(j + 1) % n2]
            r = b[:2] - a[:2]
            q = d[:2] - c[:2]
            den = r[0] * q[1] - r[1] * q[0]
            if abs(den) < 1e-14:
                continue
            w = c[:2] - a[:2]
            t = (w[0] * q[1] - w[1] * q[0]) / den
            u = (w[0] * r[1] - w[1] * r[0]) / den
            if 0 < t < 1 and 0 < u < 1:
                z1 = a[2] + t * (b[2] - a[2])
                z2 = c[2] + u * (d[2] - c[2])
                sg = (1 if den > 0 else -1) * (1 if z1 > z2 else -1)
                s += sg
    return s / 2

# ----------------------------------------------------------------------------- Collatz side

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def necklaces(L):
    out = set()
    for bits in itertools.product((0, 1), repeat=L):
        rots = [bits[i:] + bits[:i] for i in range(L)]
        if len(set(rots)) < L:
            continue
        out.add(min(rots))
    return sorted(out)


def cycle_point(w, q, shortcut):
    """fixed point of the affine composite along w for x -> x/2 (bit 0) and (q x + 1)/2 [shortcut] or q x + 1 (bit 1)."""
    a, b = Fraction(1), Fraction(0)
    for bit in w:
        if bit == 0:
            a, b = a / 2, b / 2
        elif shortcut:
            a, b = q * a / 2, (q * b + 1) / 2
        else:
            a, b = q * a, q * b + 1
    return b / (1 - a)


def follows_word(x, w, q, shortcut):
    y = x
    for bit in w:
        if y.denominator % 2 == 0 or (y.numerator % 2) != bit:
            return False
        if bit == 0:
            y = y / 2
        else:
            y = (q * y + 1) / 2 if shortcut else q * y + 1
    return y == x


def code_real(w):
    L = len(w)
    return sum(w[i] << (L - 1 - i) for i in range(L))


def code_2adic(w):
    return sum(w[i] << i for i in range(len(w)))


def cw_poly(w):
    """coefficients of c_w(q) (shortcut map): c_w(q) = sum_{t: w_t=1} q^{#ones after t} 2^t."""
    k = sum(w)
    coeff = [0] * max(k, 1)
    for t, bit in enumerate(w):
        if bit:
            after = sum(w[t + 1:])
            coeff[after] += 2 ** t
    return coeff


def positive_cycles(q, X=20000, maxsteps=3000, cap=10 ** 15):
    found = {}
    for x0 in range(1, X + 1):
        seen = set()
        x = x0
        k = 0
        while k < maxsteps and x < cap:
            if x in seen:
                cyc = []
                y = x
                while True:
                    cyc.append(y)
                    y = (q * y + 1) // 2 if y % 2 else y // 2
                    if y == x:
                        break
                found[min(cyc)] = len(cyc)
                break
            seen.add(x)
            x = (q * x + 1) // 2 if x % 2 else x // 2
            k += 1
    return found

# ----------------------------------------------------------------------------- anti-circulant tournaments

def anti_full_s(m, free):
    """extend free signs s(1..m) to s on (Z/2m) - {0} by s(-d) = -(-1)^d s(d)."""
    n = 2 * m
    s = {d: free[d - 1] for d in range(1, m + 1)}
    for d in range(1, m):
        s[n - d] = -((-1) ** d) * s[d]
    return s


def anti_tournament(m, s):
    """x -> y iff s(y - x) = (-1)^x on Z/2m."""
    n = 2 * m
    return [[1 if i != j and s[(j - i) % n] == (1 if i % 2 == 0 else -1) else 0 for j in range(n)] for i in range(n)]


def anti_pattern_key(m, s):
    """canonical representative of s under multipliers u in (Z/2m)^* and global negation (both give isomorphic tournaments)."""
    n = 2 * m
    best = None
    for u in range(1, n, 2):
        if math.gcd(u, n) != 1:
            continue
        for sign in (1, -1):
            t = tuple(sign * s[(u * d) % n] for d in range(1, n))
            if best is None or t < best:
                best = t
    return best


def prime_factors(n):
    out = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def mu_tournament(p, m):
    """QR_p restricted to the subgroup mu_2m of F_p^* (p = 3 mod 4, 2m | p - 1), labelled by exponents of a generator."""
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // r, p) != 1 for r in set(prime_factors(p - 1))))
    zeta = pow(g, (p - 1) // (2 * m), p)
    pts = [pow(zeta, i, p) for i in range(2 * m)]
    QR = {(x * x) % p for x in range(1, p)}
    return [[1 if i != j and (pts[j] - pts[i]) % p in QR else 0 for j in range(2 * m)] for i in range(2 * m)]


def par_engine(exe, A):
    """H mod 2 and per-arc parities from the bitset engine procgen_petersen_20261001_par.c"""
    n = len(A)
    inp = "%d\n" % n + "".join("".join(str(x) for x in r) + "\n" for r in A)
    out = subprocess.run([exe], input=inp.encode(), capture_output=True, check=True).stdout.decode().split("\n")
    H2 = int(out[0].split()[1])
    par = {}
    for line in out[1:]:
        if line.strip():
            a, b, p = map(int, line.split())
            par[(a, b)] = p
    return H2, par


def exact_hp_python(A):
    """independent exact Hamiltonian-path and arc counts (pure Python subset DP with big integers)."""
    n = len(A)
    N = 1 << n
    f = [[0] * n for _ in range(N)]
    g = [[0] * n for _ in range(N)]
    for v in range(n):
        f[1 << v][v] = 1
        g[1 << v][v] = 1
    for S in range(1, N):
        fs, gs = f[S], g[S]
        for v in range(n):
            if not (S >> v) & 1:
                continue
            fv, gv = fs[v], gs[v]
            if fv:
                for w in range(n):
                    if not (S >> w) & 1 and A[v][w]:
                        f[S | (1 << w)][w] += fv
            if gv:
                for w in range(n):
                    if not (S >> w) & 1 and A[w][v]:
                        g[S | (1 << w)][w] += gv
    H = sum(f[N - 1])
    c = {}
    for a in range(n):
        for b in range(n):
            if A[a][b]:
                c[(a, b)] = sum(f[S][a] * g[(N - 1) ^ S][b] for S in range(N) if (S >> a) & 1 and not (S >> b) & 1)
    return H, c


def anti_engine(exe, A):
    """H mod 2 and the parity of c on each difference class d = 1..m (arc between 0 and d), from
    procgen_petersen_20261001_anti.c (anti-circulant tournaments only; N <= 26)."""
    n = len(A)
    inp = "%d\n" % n + "".join("".join(str(x) for x in r) + "\n" for r in A)
    r = subprocess.run([exe], input=inp.encode(), capture_output=True)
    if r.returncode != 0:
        raise RuntimeError("anti engine exit code %d" % r.returncode)
    out = r.stdout.decode().split("\n")
    H2 = int(out[0].split()[1])
    par = {}
    for line in out[1:]:
        if line.strip():
            d, a, b, p = map(int, line.split())
            par[d] = (a, b, p)
    return H2, par


def f27_tournament():
    """QR_27 - 0 as an anti-circulant on Z/26: F_27 = F_3[x]/(x^3 - x - 1), vertices = powers of a generator."""
    def mul(a, b):
        c = [0] * 5
        for i in range(3):
            for j in range(3):
                c[i + j] = (c[i + j] + a[i] * b[j]) % 3
        c[2] = (c[2] + c[4]) % 3
        c[1] = (c[1] + c[4]) % 3
        c[1] = (c[1] + c[3]) % 3
        c[0] = (c[0] + c[3]) % 3
        return (c[0], c[1], c[2])
    elems = [e for e in itertools.product(range(3), repeat=3) if e != (0, 0, 0)]
    one = (1, 0, 0)

    def order(g):
        x, k = g, 1
        while x != one:
            x = mul(x, g)
            k += 1
        return k
    gen = next(g for g in elems if order(g) == 26)
    pts = [one]
    for _ in range(25):
        pts.append(mul(pts[-1], gen))
    squares = {mul(e, e) for e in elems}
    assert len(squares) == 13
    sub = lambda a, b: tuple((a[i] - b[i]) % 3 for i in range(3))
    return [[1 if i != j and sub(pts[j], pts[i]) in squares else 0 for j in range(26)] for i in range(26)]


# ----------------------------------------------------------------------------- linear K_{n+3} in R^n (floating point)

def gale_float(P):
    """Gale vectors (N x 2) of N points in R^n with N = n + 3 (basis of the affine dependences)."""
    import numpy as np
    N = len(P)
    M = np.vstack([np.ones(N), P.T])
    _, _, vt = np.linalg.svd(M)
    return vt[M.shape[0]:].T


def lk_simplices(P, A, B):
    """signed intersection number of the boundary of conv(B) with the simplex conv(A) in R^n (|A| = |B| = (n+3)/2)."""
    import numpy as np
    n = P.shape[1]
    tot = 0
    for i in range(len(B)):
        F = [B[j] for j in range(len(B)) if j != i]
        M = np.zeros((n + 2, len(A) + len(F)))
        for c, a in enumerate(A):
            M[:n, c] = P[a]
        for c, f in enumerate(F):
            M[:n, len(A) + c] = -P[f]
        M[n, :len(A)] = 1
        M[n + 1, len(A):] = 1
        rhs = np.zeros(n + 2)
        rhs[n] = rhs[n + 1] = 1
        x = np.linalg.solve(M, rhs)
        if np.all(x > 0):
            V = [P[F[j]] - P[F[0]] for j in range(1, len(F))] + [P[A[j]] - P[A[0]] for j in range(1, len(A))]
            tot += ((-1) ** i) * (1 if np.linalg.det(np.array(V)) > 0 else -1)
    return tot


def walk_linked(G):
    """splits predicted linked by the Gale walk: windows of N consecutive rays among the 2N rays +-g_i (sorted by
    angle) containing exactly N/2 negative rays, at which the walk crosses (eps_{k-1} = eps_k)."""
    N = len(G)
    rays = sorted((math.atan2(s * G[i][1], s * G[i][0]), i, s) for i in range(N) for s in (1, -1))
    eps = [s for (_, _, s) in rays]
    lab = [i for (_, i, _) in rays]
    R = 2 * N
    out = set()
    for k in range(R):
        W = frozenset(lab[(k + j) % R] for j in range(N) if eps[(k + j) % R] == -1)
        if len(W) == N // 2 and eps[(k - 1) % R] == eps[k]:
            out.add(min(W, frozenset(range(N)) - W, key=lambda s: tuple(sorted(s))))
    return out
