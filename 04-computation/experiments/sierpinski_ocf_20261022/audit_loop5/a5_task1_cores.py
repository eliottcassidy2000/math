# Audit (HYP-9371 section 8), task 1: independent check of the core lists cores_e1..e4.json and cores_e5_v5.json.
# Own code only (no import from cores_enum.py / loop_engine.py):
#   * validity: connected, minimum degree >= 3 (a loop counts 2), E = V + excess;
#   * canonical form = lexicographically smallest upper-triangular multiplicity vector over ALL V! vertex
#     permutations (brute force, numpy-vectorised); pairwise non-isomorphism = distinct canonical forms;
#   * |Aut| = #(vertex permutations fixing the multiplicity matrix) * prod_(u<v) m_uv! * prod_v m_vv! 2^(m_vv),
#     compared with the stored 'aut';
#   * networkx MultiGraph isomorphism (GraphMatcher counts parallel edges and loops) on the 71 cubic 8-vertex cores:
#     all pairs non-isomorphic, and the vertex-automorphism counts agree with the brute force;
#   * mass formula (configuration model): with C(x, y) = sum over the listed connected cores of x^V y^excess/|Aut|,
#     [x^V y^e] exp(C) must equal (1/V!) sum over labelled degree sequences d_v >= 3, sum d = 2E, of
#     (2E-1)!!/prod d_v!  (E = V + e), at every point up to excess 4 (all V) and at excess 5 for V <= 5.
#     Every connected core has positive weight 1/|Aut|, so valid + pairwise non-isomorphic + exact mass match
#     at every point => the list is complete.
import os, sys, json, itertools, time
from fractions import Fraction
from math import factorial
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PERMS = {}
def perms(V):
    if V not in _PERMS: _PERMS[V] = np.array(list(itertools.permutations(range(V))), dtype=np.int8)
    return _PERMS[V]

def mult_matrix(V, edges):
    M = np.zeros((V, V), dtype=np.int64)
    for a, b in edges:
        if a == b: M[a, a] += 1
        else: M[a, b] += 1; M[b, a] += 1
    return M

def canon_and_vaut(V, M):
    P = perms(V).astype(np.int64)
    MP = M[P[:, :, None], P[:, None, :]]                 # (V!, V, V): MP[k, i, j] = M[p_k(i), p_k(j)]
    iu = np.triu_indices(V)
    vecs = MP[:, iu[0], iu[1]]                           # (V!, V(V+1)/2)
    order = np.lexsort(vecs.T[::-1])
    canon = tuple(int(x) for x in vecs[order[0]])
    vaut = int(np.sum(np.all(vecs == M[iu], axis=1)))
    return canon, vaut

def edge_factor(M):
    V = M.shape[0]; f = 1
    for u in range(V):
        f *= factorial(int(M[u, u])) * 2 ** int(M[u, u])
        for v in range(u + 1, V): f *= factorial(int(M[u, v]))
    return f

def valid(V, edges, excess):
    if len(edges) != V + excess: return False
    deg = [0] * V
    G = nx.MultiGraph(); G.add_nodes_from(range(V))
    for a, b in edges:
        deg[a] += 1; deg[b] += 1; G.add_edge(a, b)
    return min(deg) >= 3 and nx.is_connected(G)

def check_list(fn, excess, nx_check_V=None):
    cs = json.load(open(os.path.join(D, fn)))
    byV = {}; bad_valid = 0; bad_aut = 0; canons = {}; dup = 0; mass = {}
    t0 = time.time()
    for k, c in enumerate(cs):
        V, edges = c['V'], [tuple(e) for e in c['edges']]
        byV[V] = byV.get(V, 0) + 1
        if not valid(V, edges, excess): bad_valid += 1
        M = mult_matrix(V, edges)
        can, vaut = canon_and_vaut(V, M)
        aut = vaut * edge_factor(M)
        if aut != c['aut']: bad_aut += 1
        key = (V, can)
        if key in canons: dup += 1
        canons[key] = k
        mass[V] = mass.get(V, Fraction(0)) + Fraction(1, aut)
    print(f"{fn}: {len(cs)} cores, by V {dict(sorted(byV.items()))}; invalid {bad_valid}; |Aut| mismatches {bad_aut}; "
          f"isomorphic duplicates {dup}; sum 1/|Aut| = {sum(mass.values())} ({time.time() - t0:.0f}s)", flush=True)
    if nx_check_V is not None:
        sub = [c for c in cs if c['V'] == nx_check_V]
        Gs = []
        for c in sub:
            G = nx.MultiGraph(); G.add_nodes_from(range(c['V']))
            for a, b in c['edges']: G.add_edge(a, b)
            Gs.append(G)
        iso_pairs = sum(1 for i in range(len(Gs)) for j in range(i + 1, len(Gs)) if nx.is_isomorphic(Gs[i], Gs[j]))
        aut_ok = 0; cubic = 0
        for c, G in zip(sub, Gs):
            nv = sum(1 for _ in GraphMatcher(G, G).isomorphisms_iter())
            M = mult_matrix(c['V'], c['edges'])
            aut_ok += (nv * edge_factor(M) == c['aut'])
            cubic += all(d == 3 for _, d in G.degree())
        print(f"   networkx on the {len(sub)} cores with V = {nx_check_V}: {iso_pairs} isomorphic pairs among "
              f"{len(Gs) * (len(Gs) - 1) // 2}; all cubic: {cubic == len(sub)}; |Aut| (networkx vertex automorphisms x "
              f"edge factor) equal to stored for {aut_ok}/{len(sub)}", flush=True)
    return mass

def dfact(m):
    r = 1
    for k in range(m, 0, -2): r *= k
    return r

def config_total(V, E):
    """(1/V!) sum over labelled degree sequences (d_v >= 3, sum = 2E) of (2E-1)!!/prod d_v!"""
    # sum over compositions = coefficient of x^(2E) in (sum_(d>=3) x^d/d!)^V, exact
    poly = [Fraction(0)] * (2 * E + 1); poly[0] = Fraction(1)
    single = [Fraction(0)] * (2 * E + 1)
    for d in range(3, 2 * E + 1): single[d] = Fraction(1, factorial(d))
    for _ in range(V):
        new = [Fraction(0)] * (2 * E + 1)
        for i, a in enumerate(poly):
            if a:
                for j in range(3, 2 * E + 1 - i): new[i + j] += a * single[j]
        poly = new
    return poly[2 * E] * dfact(2 * E - 1) / factorial(V)

if __name__ == '__main__':
    # sanity test of networkx multigraph handling: 2 vertices, double edge + a loop at each, vs quadruple edge
    G1 = nx.MultiGraph([(0, 1), (0, 1), (0, 0), (1, 1)]); G2 = nx.MultiGraph([(0, 1)] * 4)
    print(f"networkx sanity: double edge + 2 loops vs quadruple edge isomorphic? {nx.is_isomorphic(G1, G2)} (must be False)")
    masses = {}
    for e, fn, vnx in ((1, 'cores_e1.json', None), (2, 'cores_e2.json', None), (3, 'cores_e3.json', None),
                       (4, 'cores_e4.json', 8), (5, 'cores_e5_v5.json', None)):
        masses[e] = check_list(fn, e, vnx)
    # exp of the connected generating function, truncated at excess <= 5 and V <= 8
    C = {(V, e): m for e, d in masses.items() for V, m in d.items()}
    A = {(0, 0): Fraction(1)}; term = {(0, 0): Fraction(1)}
    for k in range(1, 6):
        new = {}
        for (v1, e1), a in term.items():
            for (v2, e2), b in C.items():
                if e1 + e2 <= 5 and v1 + v2 <= 10: new[(v1 + v2, e1 + e2)] = new.get((v1 + v2, e1 + e2), 0) + a * b
        term = {key: v / k for key, v in new.items()}
        for key, v in term.items(): A[key] = A.get(key, 0) + v
    npts = 0; nbad = 0
    for e in range(1, 6):
        for V in range(1, 2 * e + 1):
            if e == 5 and V > 5: continue
            lhs = A.get((V, e), Fraction(0)); rhs = config_total(V, V + e)
            npts += 1; nbad += lhs != rhs
            if e >= 4: print(f"   mass at (V, E) = ({V}, {V + e}): cores+products {lhs} vs configuration model {rhs}: {'ok' if lhs == rhs else 'MISMATCH'}")
    print(f"mass formula: {npts} (V, E) points (excess 1-4 all V, excess 5 with V <= 5), {nbad} mismatches")
