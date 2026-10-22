# General loop engine for the linear-forest cluster expansion (THM-4757), any excess.
#
# A core is a connected multigraph (loops allowed) with minimum degree >= 3; its edges are core paths with a reference
# orientation tail -> head.  A run decomposition of a subdivided core is summed in three layers:
#   * core vertices: a set of disjoint half-edge pairs (pass-throughs); block size s = #half-edges - #pairs,
#     Moebius weight (-1)^(s-1) (s-1)!;
#   * each core path: type m1 (one segment, i.e. one run covers the whole path) or m2 (>= 2 segments; the internal
#     segments are free runs and the junctions have weight -1); an m2 path's end segments belong to runs that end
#     at the first/last internal junction;
#   * each run touching the core (a component of the half-edge link graph) gets an orientation sign; closed runs
#     (components without an end) are excluded.
# Path propagators (all functions of Shat, so they commute):
#   m1, traversal sign s:            R_s = sum_l s^l Shat^l = s Shat (I - s Shat)^-1
#   m2, end-segment signs (s_t, s_h): M_(s_t,s_h) = -A(s_t Shat) A(s_h Shat) (I - Shat^2) (I + Shat^2)^-1,  A(w) = w (I - w)^-1
# The cluster sum of a core is (1/|G|) sum_(labelled decompositions) mu * Hom(propagators), |G| = automorphisms of the
# core counting parallel-edge permutations and loop flips.  mac-mini-2026-10-22-sierp.
import itertools
from math import factorial
from collections import defaultdict
import numpy as np

def mob(s): return (-1) ** (s - 1) * factorial(s - 1)

def pairings(items):
    """all sets of disjoint pairs of items (partial matchings)"""
    items = list(items)
    if not items: yield []; return
    h, rest = items[0], items[1:]
    for p in pairings(rest): yield p
    for i, g in enumerate(rest):
        for p in pairings(rest[:i] + rest[i + 1:]): yield [(h, g)] + p

def automorphisms(V, edges):
    """|G| = sum over vertex permutations preserving all multiplicities of prod m_uv! * prod_v (m_vv! 2^m_vv)"""
    mult = defaultdict(int)
    for a, b in edges: mult[(min(a, b), max(a, b))] += 1
    tot = 0
    for p in itertools.permutations(range(V)):
        ok = all(mult.get((min(p[a], p[b]), max(p[a], p[b])), 0) == m for (a, b), m in mult.items())
        if not ok: continue
        w = 1
        for (a, b), m in mult.items():
            w *= factorial(m) * (2 ** m if a == b else 1)
        tot += w
    return tot

def decompositions(V, edges):
    """Yield (weight, props) with props a tuple per edge of ('m1', s) or ('m2', s_t, s_h); weight includes Moebius
    weights and the junction sign is inside the m2 propagator.  Summed over core configurations, path types and run
    orientations; closed runs excluded."""
    E = len(edges)
    halves_at = defaultdict(list)
    for e, (a, b) in enumerate(edges):
        halves_at[a].append(('t', e)); halves_at[b].append(('h', e))
    per_vertex = [list(pairings(halves_at[v])) for v in range(V)]
    out = defaultdict(int)
    for cfg in itertools.product(*per_vertex):
        mu = 1; link = {}
        for v in range(V):
            s = len(halves_at[v]) - len(cfg[v]); mu *= mob(s)
            for x, y in cfg[v]: link[x] = y; link[y] = x
        for types in itertools.product(('m1', 'm2'), repeat=E):
            # half-edge graph: path link for m1 edges, core link for pairs
            def nbrs(hf):
                r = []
                k, e = hf
                if types[e] == 'm1': r.append(('h' if k == 't' else 't', e))
                if hf in link: r.append(link[hf])
                return r
            seen = set(); comps = []; bad = False
            for e in range(E):
                for k in ('t', 'h'):
                    h0 = (k, e)
                    if h0 in seen: continue
                    # collect component
                    comp = []; stack = [h0]; seen.add(h0)
                    while stack:
                        x = stack.pop(); comp.append(x)
                        for y in nbrs(x):
                            if y not in seen: seen.add(y); stack.append(y)
                    # ends: half-edges with fewer than 2 links (each half-edge has <= 2 links)
                    ends = [x for x in comp if len(nbrs(x)) < 2]
                    if not ends: bad = True; break
                    comps.append((comp, ends))
                if bad: break
            if bad: continue
            # orient each component: walk from ends[0]; record traversal sign of each path segment touched
            # a walk moves: at a half-edge, either cross its path (m1) or cross the core vertex (pair link).
            # sign convention: traversing path e from tail to head = +1.
            orient = []   # per component: dict half-edge -> sign of the segment at that half-edge, for sigma = +1
            for comp, ends in comps:
                start = ends[0]
                sg = {}
                # Determine the first move from the start: the start half-edge is an end, so it has exactly one link
                # or none (isolated: an unpaired half-edge of an m2 path = a run consisting of that end segment only).
                prev = None; cur = start
                # the segment at cur is traversed "from the end-side": if cur is an m2 half-edge with no core link it is
                # an end at the internal junction and the run proceeds out through the core vertex; if it is unpaired at
                # the core vertex the run starts at the core vertex and enters the path.
                k, e = cur
                if cur in link and types[e] == 'm2':
                    # run starts at the internal junction, travels through this end segment toward the core vertex
                    entering = False
                else:
                    entering = True   # run starts at the core vertex and enters path e through cur
                while True:
                    k, e = cur
                    # sign of traversal of the segment at cur
                    if entering: s = +1 if k == 't' else -1      # entering at tail -> along reference
                    else: s = -1 if k == 't' else +1             # leaving at tail -> against reference
                    sg[cur] = s
                    if entering:
                        # move through the path
                        if types[e] == 'm1':
                            other = ('h' if k == 't' else 't', e)
                            sg[other] = s                            # same segment, same sign
                            if other in link:
                                cur = link[other]; entering = True; continue
                            break                                    # run ends at the core vertex of 'other'
                        break                                        # m2: run ends at the internal junction
                    else:
                        # leaving the path at cur through the core vertex
                        if cur in link:
                            cur = link[cur]; entering = True; continue
                        break
                orient.append(sg)
            # sum over component orientations
            for sigmas in itertools.product((1, -1), repeat=len(comps)):
                sgn = {}
                for (comp, ends), sg, sig in zip(comps, orient, sigmas):
                    for x, s in sg.items(): sgn[x] = s * sig
                props = []
                for e in range(E):
                    if types[e] == 'm1': props.append(('m1', sgn[('t', e)]))
                    else: props.append(('m2', sgn[('t', e)], sgn[('h', e)]))
                out[tuple(props)] += mu
    return out

def propagator(desc, Sh):
    n = Sh.shape[0]; I = np.eye(n)
    if desc[0] == 'm1':
        s = desc[1]; return s * Sh @ np.linalg.inv(I - s * Sh)
    s1, s2 = desc[1], desc[2]
    A1 = s1 * Sh @ np.linalg.inv(I - s1 * Sh); A2 = s2 * Sh @ np.linalg.inv(I - s2 * Sh)
    S2 = Sh @ Sh
    return -A1 @ A2 @ (I - S2) @ np.linalg.inv(I + S2)

def hom(V, edges, mats):
    letters = 'abcdefgh'
    terms = []
    for (a, b), M in zip(edges, mats):
        if a == b: terms.append((letters[a], np.diag(M)))
        else: terms.append((letters[a] + letters[b], M))
    spec = ','.join(t[0] for t in terms) + '->'
    return float(np.einsum(spec, *[t[1] for t in terms], optimize='greedy').real)

def core_sum(V, edges, Sh, cache=None):
    dec = decompositions(V, edges)
    G = automorphisms(V, edges)
    pc = {} if cache is None else cache
    tot = 0.0
    for props, w in dec.items():
        if w == 0: continue
        mats = []
        for d in props:
            if d not in pc: pc[d] = propagator(d, Sh)
            mats.append(pc[d])
        tot += w * hom(V, edges, mats)
    return tot / G

TWO_LOOP_CORES = {
    'figure-eight': (1, [(0, 0), (0, 0)]),
    'theta': (2, [(0, 1), (0, 1), (0, 1)]),
    'dumbbell': (2, [(0, 0), (1, 1), (0, 1)]),
}

if __name__ == '__main__':
    import sys
    sys.path.insert(0, '.')
    from two_loop_regular import two_loop
    rng = np.random.default_rng(5)
    # a random regular tournament on 11 vertices via 3-cycle reversals from R_11 (not vertex-transitive)
    n = 11; d = 5
    A = np.array([[1 if 0 < (j - i) % n <= d else 0 for j in range(n)] for i in range(n)])
    for _ in range(4000):
        u, v, w = rng.choice(n, 3, replace=False)
        if A[u, v] and A[v, w] and A[w, u]: A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
    S = (A - A.T).astype(float); Sh = S / n
    I = np.eye(n); S2 = Sh @ Sh
    G = -S2 @ np.linalg.inv(I + S2); H = S2 @ np.linalg.inv(I - S2); O = Sh @ np.linalg.inv(I - S2)
    gd, hd = np.diag(G), np.diag(H)
    ref = {'figure-eight': np.sum(2 * gd * hd - 1.5 * gd ** 2),
           'theta': np.sum(G * (H * H - O * O - (2 / 3) * G * G)),
           'dumbbell': -(hd - gd) @ G @ (hd - gd)}
    cache = {}
    for name, (V, edges) in TWO_LOOP_CORES.items():
        val = core_sum(V, edges, Sh, cache)
        print(f"{name:13s}: engine = {val:+.12e}   closed form = {ref[name]:+.12e}   |G| = {automorphisms(V, edges)}")

# ---------------------------------------------------------------------------------------------------------------
# Faster decomposition: runs that never pass through a core vertex (an m1 path with both ends free, or an unpaired
# end segment of an m2 path) have a free orientation that is summed inside the propagator ('*' signs), so only runs
# through core pass-throughs are enumerated explicitly.
def decompositions_fast(V, edges):
    E = len(edges)
    halves_at = defaultdict(list)
    for e, (a, b) in enumerate(edges):
        halves_at[a].append(('t', e)); halves_at[b].append(('h', e))
    per_vertex = [list(pairings(halves_at[v])) for v in range(V)]
    out = defaultdict(int)
    for cfg in itertools.product(*per_vertex):
        mu = 1; link = {}
        for v in range(V):
            s = len(halves_at[v]) - len(cfg[v]); mu *= mob(s)
            for x, y in cfg[v]: link[x] = y; link[y] = x
        for types in itertools.product(('m1', 'm2'), repeat=E):
            def nbrs(hf):
                r = []
                k, e = hf
                if types[e] == 'm1': r.append(('h' if k == 't' else 't', e))
                if hf in link: r.append(link[hf])
                return r
            seen = set(); comps = []; bad = False
            for e in range(E):
                for k in ('t', 'h'):
                    h0 = (k, e)
                    if h0 in seen: continue
                    comp = []; stack = [h0]; seen.add(h0)
                    while stack:
                        x = stack.pop(); comp.append(x)
                        for y in nbrs(x):
                            if y not in seen: seen.add(y); stack.append(y)
                    ends = [x for x in comp if len(nbrs(x)) < 2]
                    if not ends: bad = True; break
                    comps.append((comp, ends))
                if bad: break
            if bad: continue
            free = []; fixed = []
            for comp, ends in comps:
                if not any(x in link for x in comp): free.append(comp)
                else: fixed.append((comp, ends))
            orient = []
            for comp, ends in fixed:
                start = ends[0]; sg = {}; cur = start
                k, e = cur
                entering = not (cur in link and types[e] == 'm2')
                while True:
                    k, e = cur
                    s = (+1 if k == 't' else -1) if entering else (-1 if k == 't' else +1)
                    sg[cur] = s
                    if entering:
                        if types[e] == 'm1':
                            other = ('h' if k == 't' else 't', e); sg[other] = s
                            if other in link: cur = link[other]; entering = True; continue
                        break
                    else:
                        if cur in link: cur = link[cur]; entering = True; continue
                        break
                orient.append(sg)
            freeset = set(x for comp in free for x in comp)
            for sigmas in itertools.product((1, -1), repeat=len(fixed)):
                sgn = {}
                for sg, sig in zip(orient, sigmas):
                    for x, s in sg.items(): sgn[x] = s * sig
                props = []
                for e in range(E):
                    t, h = ('t', e), ('h', e)
                    if types[e] == 'm1':
                        props.append(('m1', '*') if t in freeset else ('m1', sgn[t]))
                    else:
                        props.append(('m2', '*' if t in freeset else sgn[t], '*' if h in freeset else sgn[h]))
                out[tuple(props)] += mu
    return out

def propagator2(desc, Sh, z=1.0):
    n = Sh.shape[0]; I = np.eye(n); X = z * Sh
    def R(s): return s * X @ np.linalg.inv(I - s * X)
    def Aw(s): return R(s)
    if desc[0] == 'm1':
        return R(1) + R(-1) if desc[1] == '*' else R(desc[1])
    def Asum(s): return Aw(1) + Aw(-1) if s == '*' else Aw(s)
    X2 = X @ X
    return -Asum(desc[1]) @ Asum(desc[2]) @ (I - X2) @ np.linalg.inv(I + X2)

def core_sum_fast(V, edges, Sh, z=1.0, dec=None):
    if dec is None: dec = decompositions_fast(V, edges)
    G = automorphisms(V, edges)
    pc = {}; tot = 0.0
    for props, w in dec.items():
        if w == 0: continue
        mats = []
        for d in props:
            if d not in pc: pc[d] = propagator2(d, Sh, z)
            mats.append(pc[d])
        tot += w * hom(V, edges, mats)
    return tot / G

def cores(excess):
    """connected multigraphs with loops, min degree >= 3, E - V = excess, up to isomorphism"""
    found = []
    for V in range(1, 2 * excess + 1):
        E = V + excess
        pairs = [(a, b) for a in range(V) for b in range(a, V)]
        seen = set()
        for multiset in itertools.combinations_with_replacement(range(len(pairs)), E):
            edges = [pairs[i] for i in multiset]
            deg = [0] * V
            for a, b in edges: deg[a] += 1; deg[b] += 1
            if min(deg) < 3: continue
            par = list(range(V))
            def f(x):
                while par[x] != x: par[x] = par[par[x]]; x = par[x]
                return x
            for a, b in edges: par[f(a)] = f(b)
            if len({f(x) for x in range(V)}) != 1: continue
            canon = min(tuple(sorted((min(p[a], p[b]), max(p[a], p[b])) for a, b in edges)) for p in itertools.permutations(range(V)))
            if canon in seen: continue
            seen.add(canon); found.append((V, list(canon)))
    return found

def hom_c(V, edges, mats):
    letters = 'abcdefgh'
    terms = []
    for (a, b), M in zip(edges, mats):
        if a == b: terms.append((letters[a], np.diag(M)))
        else: terms.append((letters[a] + letters[b], M))
    spec = ','.join(t[0] for t in terms) + '->'
    return complex(np.einsum(spec, *[t[1] for t in terms], optimize='greedy'))

def excess_series(Sh, excess, Lmax, K=64, r=0.6, decs=None):
    """coefficients c_L, L <= Lmax, of X_excess(z) (Shat -> z Shat), by Cauchy sums on |z| = r"""
    cs = cores(excess)
    if decs is None: decs = [decompositions_fast(V, edges) for V, edges in cs]
    vals = np.zeros(K, dtype=complex)
    for kk in range(K):
        z = r * np.exp(2j * np.pi * kk / K)
        tot = 0
        for (V, edges), dec in zip(cs, decs):
            G = automorphisms(V, edges); pc = {}; sub = 0
            for props, w in dec.items():
                if w == 0: continue
                mats = []
                for d in props:
                    if d not in pc: pc[d] = propagator2(d, Sh.astype(complex), z)
                    mats.append(pc[d])
                sub += w * hom_c(V, edges, mats)
            tot += sub / G
        vals[kk] = tot
    return [(np.sum(vals * np.exp(-2j * np.pi * np.arange(K) * L / K)) / K / r ** L).real for L in range(Lmax + 1)]
