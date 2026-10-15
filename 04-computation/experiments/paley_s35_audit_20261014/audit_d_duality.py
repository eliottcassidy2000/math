#!/usr/bin/env python3
"""Audit (d) of S35: planar Fourier duality U(G) = chi(-1)^{|E|/2} q^{(|V|-|V*|)/2} U(G*), and U(Q_3) = q U(K_{2,2,2}).

Independent of networkx's planarity embedding: plane embeddings are found by brute force over rotation systems
(genus 0 <=> #faces = |E| - |V| + 2).  Faces are traced with the face on the LEFT of each half-edge (next half-edge
after u->v is v->w, w the clockwise neighbour of u at v); the dual arc of u->v runs from the face left of u->v to the
face left of v->u (i.e. left -> right), as in the S35 statement.  U by brute force (translation x_0 = 0, factor q).
Extra tests: graphs with |E|/2 odd at q = 3 mod 4 (sign), a 2-connected graph with several inequivalent embeddings,
and the cube identity at all primes 5 <= q <= 160 with an O(q^3) cube evaluator written here.
"""
import itertools, math
import numpy as np
import numba as nb
import networkx as nx
import sympy as sp


def chi_table(q):
    c = np.full(q, -1, dtype=np.int64)
    c[0] = 0
    for x in range(1, q):
        c[(x * x) % q] = 1
    return c


@nb.njit(cache=True)
def U_brute(q, chi, nv, eu, ev):
    tot = 0
    v = np.zeros(nv, dtype=np.int64)
    N = q ** (nv - 1)
    for code in range(N):
        c = code
        for i in range(1, nv):
            v[i] = c % q
            c //= q
        p = 1
        for k in range(eu.shape[0]):
            p *= chi[(v[ev[k]] - v[eu[k]]) % q]
            if p == 0:
                break
        tot += p
    return tot * q


def U(q, nv, arcs):
    eu = np.array([a for a, b in arcs], dtype=np.int64)
    ev = np.array([b for a, b in arcs], dtype=np.int64)
    return int(U_brute(q, chi_table(q), nv, eu, ev))


def faces_of(rot):
    """rot[v] = cyclic (ccw) list of neighbours.  Returns list of faces (lists of half-edges), face on the left."""
    pos = {v: {w: i for i, w in enumerate(r)} for v, r in rot.items()}
    seen = set()
    faces = []
    for u in rot:
        for v in rot[u]:
            if (u, v) in seen:
                continue
            f = []
            a, b = u, v
            while (a, b) not in seen:
                seen.add((a, b))
                f.append((a, b))
                r = rot[b]
                w = r[(pos[b][a] - 1) % len(r)]       # clockwise neighbour of a at b
                a, b = b, w
            faces.append(f)
    return faces


def plane_embeddings(nv, edges):
    nbrs = {v: sorted({b for a, b in edges if a == v} | {a for a, b in edges if b == v}) for v in range(nv)}
    target = len(edges) - nv + 2
    choices = []
    for v in range(nv):
        first, rest = nbrs[v][0], nbrs[v][1:]
        choices.append([[first] + list(p) for p in itertools.permutations(rest)])
    for combo in itertools.product(*choices):
        rot = {v: combo[v] for v in range(nv)}
        F = faces_of(rot)
        if len(F) == target:
            yield rot, F


def dual(arcs, F):
    fid = {}
    for i, f in enumerate(F):
        for he in f:
            fid[he] = i
    return len(F), [(fid[(u, v)], fid[(v, u)]) for (u, v) in arcs]


def check(name, nv, arcs, primes, max_emb=1):
    edges = [(min(a, b), max(a, b)) for a, b in arcs]
    out = []
    n_emb = 0
    duals_iso = set()
    for rot, F in plane_embeddings(nv, edges):
        n_emb += 1
        if n_emb > max_emb:
            continue
        nf, darcs = dual(arcs, F)
        Dg = nx.MultiGraph([(a, b) for a, b in darcs])
        duals_iso.add(nx.weisfeiler_lehman_graph_hash(nx.Graph(Dg)) if not any(a == b for a, b in darcs) else 'loop')
        E = len(arcs)
        for q in primes:
            a = U(q, nv, arcs)
            b = U(q, nf, darcs)
            cm1 = 1 if q % 4 == 1 else -1
            pred = cm1 ** (E // 2) * q ** ((nv - nf) // 2) * b if (nv - nf) % 2 == 0 else None
            out.append((n_emb, q, a, b, pred == a))
    print(f'{name}: |V| = {nv}, |E| = {len(arcs)}, plane embeddings (rotation systems) found: {n_emb}; '
          f'(embedding, q, U(G), U(G*), duality holds): {out}', flush=True)
    return all(x[-1] for x in out)


@nb.njit(cache=True)
def cube_U(q, chi):
    # Q_3, vertices 0..7, arcs u -> u | 2^k (bit orientation).  Even vertices 0,3,5,6; odd 1,2,4,7 have all
    # neighbours even and all signs +1 (each odd vertex has an even number of larger neighbours).
    T = np.zeros((q, q), dtype=np.int64)
    for b in range(q):
        for c in range(q):
            s = 0
            for y in range(q):
                s += chi[y] * chi[(y - b) % q] * chi[(y - c) % q]
            T[b, c] = s
    tot = 0
    for x3 in range(q):
        for x5 in range(q):
            t1 = T[x3, x5]                     # vertex 1: neighbours 0, 3, 5
            if t1 == 0:
                continue
            for x6 in range(q):
                t2 = T[x3, x6]                 # vertex 2: neighbours 0, 3, 6
                t4 = T[x5, x6]                 # vertex 4: neighbours 0, 5, 6
                t7 = T[(x5 - x3) % q, (x6 - x3) % q]   # vertex 7: neighbours 3, 5, 6
                tot += t1 * t2 * t4 * t7
    return tot * q


def beukers(N):
    P = [0] * N
    P[0] = 1
    for m in range(1, N):
        for step in (2 * m, 4 * m):
            if step >= N:
                continue
            for _ in range(4):
                for i in range(N - 1, step - 1, -1):
                    P[i] -= P[i - step]
    return [0] + P[:N - 1]


if __name__ == '__main__':
    allok = True
    cube = [(u, u | (1 << k)) for u in range(8) for k in range(3) if not u >> k & 1]
    allok &= check('cube Q_3 (bit orientation; dual = octahedron)', 8, cube, [5, 7, 11])
    k6p6 = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4)]
    allok &= check('K6 - P6 (S35 orientation; |E|/2 = 5 odd)', 6, k6p6, [5, 7, 11, 13])
    # is K6 - P6 self-dual as a plane graph?
    for rot, F in plane_embeddings(6, [(min(a, b), max(a, b)) for a, b in k6p6]):
        nf, darcs = dual(k6p6, F)
        print('K6 - P6 dual isomorphic to K6 - P6 (undirected):',
              nx.is_isomorphic(nx.Graph([(a, b) for a, b in darcs]), nx.Graph([(a, b) for a, b in k6p6])))
        break
    w4 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (4, 1)]
    allok &= check('W_4 (|E|/2 = 4)', 5, w4, [7, 11])
    w5 = [(0, i) for i in range(1, 6)] + [(i, i % 5 + 1) for i in range(1, 6)]
    allok &= check('W_5 (|E|/2 = 5 odd)', 6, w5, [7, 11])
    w6 = [(0, i) for i in range(1, 7)] + [(i, i % 6 + 1) for i in range(1, 7)]
    allok &= check('W_6 (|E|/2 = 6)', 7, w6, [5, 7])
    # a 2-connected, not 3-connected graph with several embeddings, random orientations, |E| even
    rng = np.random.default_rng(5)
    g2 = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 4), (4, 3), (4, 5), (5, 1), (5, 3)]
    g2 = [(a, b) if rng.random() < 0.5 else (b, a) for a, b in g2]
    allok &= check('2-connected graph, several embeddings, random orientation', 6, g2, [7, 11], max_emb=4)
    g3 = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2), (0, 4), (1, 3)]
    g3 = [(a, b) if rng.random() < 0.5 else (b, a) for a, b in g3]
    allok &= check('random-orientation 5-vertex graph', 5, g3, [7, 11, 13], max_emb=3)
    print('duality checks all hold:', allok)
    # cube identity at all primes 5..160 vs q * octahedron closed form q(q-1)((5-q) a_q + 3q)
    a = beukers(200)
    bad = []
    n = 0
    for q in sp.primerange(5, 161):
        Uc = int(cube_U(q, chi_table(q)))
        pred = q * q * (q - 1) * ((5 - q) * a[q] + 3 * q)
        n += 1
        if Uc != pred:
            bad.append((q, Uc, pred))
    print(f'U(Q_3) (O(q^3) evaluator) = q * q(q-1)((5-q) a_q + 3q) at {n} primes 5..157: failures {bad}')
    # brute-force sanity of cube_U
    for q in (5, 7):
        print(f'cube_U vs brute force at q = {q}:', int(cube_U(q, chi_table(q))) == U(q, 8, cube))
