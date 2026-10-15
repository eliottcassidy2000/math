#!/usr/bin/env python3
"""Planar Fourier (Kramers-Wannier) duality for Paley pattern sums, checked numerically.

PROVED in the S35 note (section 4d): for a connected plane graph G with an orientation omega and |E| even,
    U(G, omega) = chi(-1)^{|E|/2} q^{(|V| - |V*|)/2} U(G*, omega*),
where G* is the plane dual and omega* the induced orientation (the dual edge of u -> v crosses it from the face on
its left to the face on its right).  Proof: Fourier expansion turns the vertex sum into a sum over flows (q^{|V|-|E|}),
flows on a plane graph are face potentials modulo constants (1/q), and chi-hat = chi(-1) G chi with G^2 = chi(-1) q.
Checks: the cube Q_3 against the octahedron (U(Q_3) = q U(K_{2,2,2})), self-duality of K_6 - P_6 and of the wheels
W_4, W_5, and random plane graphs with |E| even, at several primes in both classes mod 4.
Session opus-2026-10-14-S35.
"""
import itertools, random
import numpy as np
import networkx as nx


def chi_table(q):
    chi = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        chi[x] = 1 if pow(x, (q - 1) // 2, q) == 1 else -1
    return chi


def U(q, n, arcs):
    """unrestricted sum over x: V -> F_q of prod chi(x_v - x_u), by vertex elimination order brute force (small n)"""
    chi = chi_table(q)
    # fix nothing; brute force over q^n (n <= 7 at small q)
    tot = 0
    for x in itertools.product(range(q), repeat=n):
        p = 1
        for (u, v) in arcs:
            p *= chi[(x[v] - x[u]) % q]
            if p == 0:
                break
        tot += p
    return tot


def dual_with_orientation(G, arcs):
    ok, emb = nx.check_planarity(G)
    assert ok
    seen = set(); faces = []
    for u, v in emb.edges():
        if (u, v) in seen:
            continue
        faces.append(emb.traverse_face(u, v, mark_half_edges=seen))
    # half-edge (u, v) belongs to the face traversed with u -> v; that face lies to the LEFT of u -> v in networkx's
    # counterclockwise traversal convention (traverse_face follows the face to the left... we fix the convention
    # empirically by the cube test and keep it for all checks)
    fid = {}
    for i, f in enumerate(faces):
        for j in range(len(f)):
            fid[(f[j], f[(j + 1) % len(f)])] = i
    darcs = [(fid[(u, v)], fid[(v, u)]) for (u, v) in arcs]
    return len(faces), darcs


def check(name, n, arcs, primes):
    G = nx.Graph([(u, v) for u, v in arcs])
    nf, darcs = dual_with_orientation(G, arcs)
    E = len(arcs)
    res = []
    for q in primes:
        a = U(q, n, arcs); b = U(q, nf, darcs)
        sgn = (1 if q % 4 == 1 else -1) ** (E // 2)
        pred = sgn * q ** ((n - nf) // 2) * b if (n - nf) % 2 == 0 else None
        res.append((q, a, b, pred == a if pred is not None else 'odd exponent'))
    print(f'{name}: |V| = {n}, |E| = {E}, |V*| = {nf};  (q, U(G), U(G*), duality holds): {res}', flush=True)


if __name__ == '__main__':
    cube = [(u, u | (1 << k)) for u in range(8) for k in range(3) if not u >> k & 1]
    check('cube Q_3 (dual: octahedron)', 8, cube, [5, 7])
    k6p6 = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4)]
    check('K_6 - P_6 (self-dual)', 6, k6p6, [5, 7, 11, 13])
    w4 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (4, 1)]
    check('wheel W_4 (self-dual)', 5, w4, [5, 7, 11, 13])
    w5 = [(0, i) for i in range(1, 6)] + [(i, i % 5 + 1) for i in range(1, 6)]
    check('wheel W_5 (self-dual)', 6, w5, [5, 7, 13])
    random.seed(7)
    tries = 0; done = 0
    while done < 3 and tries < 500:
        tries += 1
        n = random.choice([5, 6])
        H = nx.gnm_random_graph(n, random.choice([8, 10]), seed=random.randrange(10 ** 6))
        if not nx.is_connected(H) or not nx.check_planarity(H)[0] or H.number_of_edges() % 2:
            continue
        if min(dict(H.degree()).values()) < 2:
            continue
        arcs = [(u, v) if random.random() < 0.5 else (v, u) for u, v in H.edges()]
        check(f'random plane graph #{done + 1} {sorted(H.edges())}', n, arcs, [5, 7, 11])
        done += 1
