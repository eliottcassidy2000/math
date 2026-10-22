# Sign-network engine for the linear-forest loop expansion (THM-4757), any excess; built for five loops (excess 4).
#
# The run-decomposition sum of a core (loop_engine.decompositions_fast) is a tensor network over half-edge signs:
# o(x) = +1 if the run through half-edge x moves away from x's core vertex, -1 if it moves towards it.
#   * vertex v (degree d, k pass-through pairs already fixed): V_v[o] = sum over partial matchings pi of its free
#     half-edges of mob(d - k - |pi|) prod_(x,y in pi) [o(x) = -o(y)]  (a pass-through enters through one half-edge
#     and leaves through the other; an unmatched half-edge is a run end, both directions allowed);
#   * edge e (tail half-edge t, head h): m1 (one segment, traversal sign s = o(t) = -o(h)) carries R_s; m2 (>= 2
#     segments) carries Q_(s_t, s_h) with s_t = o(t), s_h = -o(h).  Free runs are summed automatically ('*').
# Contracting the network sums every configuration, including those with closed runs (cycles of m1 edges joined by
# pass-throughs), which the cluster expansion excludes.  Inclusion-exclusion removes them:
#     Valid = sum over sets A of pairwise half-edge-disjoint closed runs of (-1)^|A| N(A),
# N(A) = the network with A's edges forced to m1 with one orientation sign per closed run, A's pairs fixed (k_v),
# and everything else free (closed runs allowed).
# Output modes:
#   'basis5' : weight tensor over the 5-matrix basis of loop_engine_fast (integers; equals decomp.c's tensors);
#   'general': the core's value for a regular tournament (matrices R_s, Q, contracted over [n]^V), batched over z;
#   'drt'    : doubly regular T: tensor over {P, Shat} labels of the edges, batched over z; every propagator is
#              alpha P + gamma Shat, with alpha, gamma from a = u^2/(1-u^2), b = u^2/(1+u^2), g = z/(1-u^2),
#              h = -2u^2 z/((1+u^2)(1-u^2)), u^2 = -z^2/n:
#              m1: P a, S g s;   m2: P ((b - a) + (b + a) o_t o_h)/2, S (h/2)(o_t - o_h).
# mac-mini-2026-10-22-loop4.
import os, sys, json
from math import factorial
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from loop_engine import pairings

LET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

def mob(s): return (-1) ** (s - 1) * factorial(s - 1)

_VT = {}
def vertex_tensor(d, dk):
    """V[o_1..o_d] (index 0 = +1, 1 = -1) for d free half-edges at a vertex with (total half-edges - fixed pairs) =
    dk; d = 0 gives the scalar mob(dk)"""
    key = (d, dk)
    if key not in _VT:
        import itertools
        idx = np.array(list(itertools.product((0, 1), repeat=d)), dtype=int).T.reshape(d, 2 ** d)
        flat = np.zeros(2 ** d, dtype=np.int64)
        for pi in pairings(list(range(d))):
            mask = np.ones(2 ** d, bool)
            for x, y in pi: mask &= idx[x] != idx[y]
            flat[mask] += mob(dk - len(pi))
        _VT[key] = flat.reshape((2,) * d)
    return _VT[key]

class SignCore:
    def __init__(self, V, edges, aut):
        self.V, self.edges, self.aut = V, [tuple(e) for e in edges], aut
        self.E = E = len(edges)
        self.hv = [self.edges[x // 2][x % 2] for x in range(2 * E)]       # vertex of half-edge x (2e tail, 2e+1 head)
        self.H = [[x for x in range(2 * E) if self.hv[x] == v] for v in range(V)]
        self.runs = self._closed_runs()
        self.sets = self._compatible_sets()

    def _closed_runs(self):
        """closed runs = cycles of the link graph through m1 edges and pass-through pairs, i.e. closed trails with a
        pairing at every visit; returns list of (frozenset of pairs, {edge: eps}) with eps = +1 if traversed tail ->
        head in the run's reference orientation"""
        hv, E = self.hv, self.E
        found = {}
        for e0 in range(E):
            for d0 in (0, 1):
                entry0 = 2 * e0 + d0
                def dfs(exit_he, used, pairs, seq):
                    v = hv[exit_he]
                    if v == hv[entry0]:
                        key = frozenset(pairs + [frozenset((exit_he, entry0))])
                        if key not in found: found[key] = {e: (1 if d == 0 else -1) for e, d in seq}
                    for e in range(e0 + 1, E):
                        if e in used: continue
                        for d in (0, 1):
                            ent = 2 * e + d
                            if hv[ent] != v: continue
                            dfs(2 * e + 1 - d, used | {e}, pairs + [frozenset((exit_he, ent))], seq + [(e, d)])
                dfs(2 * e0 + 1 - d0, {e0}, [], [(e0, d0)])
        return [(k, v) for k, v in found.items()]

    def _compatible_sets(self):
        """all sets of pairwise half-edge-disjoint closed runs (including the empty set); recursion depth = |A|"""
        masks = []
        for pairs, _ in self.runs:
            m = 0
            for p in pairs:
                for x in p: m |= 1 << x
            masks.append(m)
        R = len(masks); out = []
        def rec(start, chosen, used):
            out.append(tuple(chosen))
            for i in range(start, R):
                if not (masks[i] & used): rec(i + 1, chosen + [i], used | masks[i])
        rec(0, [], 0)
        return out

    def network(self, A):
        """structure of N(A): vertex tensors (free half-edges per vertex, total-pairs-adjusted Moebius argument),
        free edges, and the closed runs of A with their edges and signs"""
        used = set(); kv = [0] * self.V; cyc = []
        for r in A:
            pairs, eps = self.runs[r]
            for p in pairs:
                x, y = tuple(p); used |= {x, y}; kv[self.hv[x]] += 1
            cyc.append(eps)
        aedges = set(e for eps in cyc for e in eps)
        verts = []
        for v in range(self.V):
            free = [x for x in self.H[v] if x not in used]
            verts.append((free, len(self.H[v]) - kv[v]))
        fedges = [e for e in range(self.E) if e not in aedges]
        return verts, fedges, cyc

    def contract(self, mode, ops, z_batch=None):
        """sum over compatible sets A of (-1)^|A| N(A).  ops: mode-specific tensors (see build_ops)."""
        total = None
        for A in self.sets:
            verts, fedges, cyc = self.network(A)
            spec_in = []; operands = []
            ol = {}                                              # half-edge -> letter
            pos = [0]
            def new():
                c = LET[pos[0]]; pos[0] += 1; return c
            K = new() if mode in ('general', 'drt') else ''
            for v, (free, dk) in enumerate(verts):
                for x in free: ol[x] = new()
            tl = {e: new() for e in range(self.E)} if mode in ('basis5', 'drt') else {}
            fl = {v: new() for v in range(self.V)} if mode == 'general' else {}
            for v, (free, dk) in enumerate(verts):
                spec_in.append(''.join(ol[x] for x in free)); operands.append(vertex_tensor(len(free), dk))
            for e in fedges:
                a, b = self.edges[e]
                if mode == 'basis5':
                    spec_in.append(ol[2 * e] + ol[2 * e + 1] + tl[e]); operands.append(ops['E'])
                elif mode == 'drt':
                    spec_in.append(K + ol[2 * e] + ol[2 * e + 1] + tl[e]); operands.append(ops['E'])
                else:
                    if a == b:
                        spec_in.append(K + ol[2 * e] + ol[2 * e + 1] + fl[a]); operands.append(ops['Ediag'])
                    else:
                        spec_in.append(K + ol[2 * e] + ol[2 * e + 1] + fl[a] + fl[b]); operands.append(ops['E'])
            for eps in cyc:
                s = new()
                for e, ep in eps.items():
                    a, b = self.edges[e]
                    key = 'F+' if ep == 1 else 'F-'
                    if mode == 'basis5':
                        spec_in.append(s + tl[e]); operands.append(ops[key])
                    elif mode == 'drt':
                        spec_in.append(K + s + tl[e]); operands.append(ops[key])
                    else:
                        if a == b:
                            spec_in.append(K + s + fl[a]); operands.append(ops[key + 'diag'])
                        else:
                            spec_in.append(K + s + fl[a] + fl[b]); operands.append(ops[key])
            if mode == 'basis5': out = ''.join(tl[e] for e in range(self.E))
            elif mode == 'drt': out = K + ''.join(tl[e] for e in range(self.E))
            else: out = K
            if mode == 'general':
                groups = [set([fl[v]] + [ol[x] for x in verts[v][0]]) for v in range(self.V)]
                val = eliminate(spec_in, operands, out, groups)
            else:
                spec = ','.join(spec_in) + '->' + out
                path = np.einsum_path(spec, *operands, optimize='greedy')[0]
                val = np.einsum(spec, *operands, optimize=path)
            val = val * (-1) ** len(A)
            total = val if total is None else total + val
        return total

def eliminate(specs, operands, out, groups):
    """contract a tensor network by eliminating groups of letters (one core vertex: its f index and its half-edge
    signs) one at a time, always the group whose merged tensor is smallest; each step is one small einsum"""
    specs = list(specs); ops = list(operands); groups = [set(g) for g in groups]
    dims = {}
    for s, o in zip(specs, ops):
        for c, d in zip(s, np.shape(o)): dims[c] = d
    while groups:
        live = set(c for s in specs for c in s)
        groups = [g for g in groups if g & live]                   # groups already summed in an earlier merge
        if not groups: break
        best = None
        for gi, g in enumerate(groups):
            idx = [i for i, s in enumerate(specs) if set(s) & g]
            others = set(c for i, s in enumerate(specs) if i not in idx for c in s) | set(out)
            keep = sorted(set(c for i in idx for c in specs[i]) & others)
            size = 1
            for c in keep: size *= dims[c]
            if best is None or size < best[0]: best = (size, gi, idx, keep)
        _, gi, idx, keep = best
        sub = ','.join(specs[i] for i in idx) + '->' + ''.join(keep)
        t = np.einsum(sub, *[ops[i] for i in idx], optimize='greedy')
        specs = [s for i, s in enumerate(specs) if i not in idx] + [''.join(keep)]
        ops = [o for i, o in enumerate(ops) if i not in idx] + [t]
        groups.pop(gi)
    return np.einsum(','.join(specs) + '->' + out, *ops, optimize='greedy')

# ---------------------------------------------------------------------------------------------------------------
# mode-specific edge tensors
def ops_basis5():
    """E[o_t, o_h, b], F+/-[sigma, b]: indices 0 = +1, 1 = -1; basis R_+, R_-, Q_++, Q_+-, Q_--"""
    E = np.zeros((2, 2, 5), dtype=np.int64)
    sg = (1, -1)
    for i in range(2):
        for j in range(2):
            ot, oh = sg[i], sg[j]
            if ot == -oh: E[i, j, 0 if ot == 1 else 1] += 1                # m1, R_(o_t)
            st, sh = ot, -oh                                               # m2, Q_(s_t, s_h)
            E[i, j, 2 if (st, sh) == (1, 1) else 4 if (st, sh) == (-1, -1) else 3] += 1
    Fp = np.zeros((2, 5), dtype=np.int64); Fm = np.zeros((2, 5), dtype=np.int64)
    for k, s in enumerate(sg):
        Fp[k, 0 if s == 1 else 1] = 1                                      # eps = +1: traversal sign sigma
        Fm[k, 0 if -s == 1 else 1] = 1                                     # eps = -1: traversal sign -sigma
    return {'E': E, 'F+': Fp, 'F-': Fm}

def drt_funcs(n, z):
    z = np.asarray(z, dtype=complex)
    u2 = -z * z / n
    a = u2 / (1 - u2); b = u2 / (1 + u2); g = z / (1 - u2); h = -2 * u2 * z / ((1 + u2) * (1 - u2))
    return a, b, g, h

def ops_drt(n, zs):
    """E[K, o_t, o_h, label] and F+/-[K, sigma, label]; label 0 = P, 1 = Shat"""
    a, b, g, h = drt_funcs(n, zs)
    K = len(a); sg = (1, -1)
    E = np.zeros((K, 2, 2, 2), dtype=complex)
    for i in range(2):
        for j in range(2):
            ot, oh = sg[i], sg[j]
            if ot == -oh:
                E[:, i, j, 0] += a; E[:, i, j, 1] += g * ot                # m1
            E[:, i, j, 0] += ((b - a) + (b + a) * ot * oh) / 2             # m2
            E[:, i, j, 1] += h * (ot - oh) / 2
    Fp = np.zeros((K, 2, 2), dtype=complex); Fm = np.zeros((K, 2, 2), dtype=complex)
    for k, s in enumerate(sg):
        Fp[:, k, 0] = a; Fp[:, k, 1] = g * s
        Fm[:, k, 0] = a; Fm[:, k, 1] = -g * s
    return {'E': E, 'F+': Fp, 'F-': Fm}

def ops_general(Sh, zs):
    """E[K, o_t, o_h, i, j], diagonal versions for loops, F+/-[K, sigma, i, j]"""
    n = Sh.shape[0]; I = np.eye(n); sg = (1, -1)
    K = len(zs)
    E = np.zeros((K, 2, 2, n, n), dtype=complex); Fp = np.zeros((K, 2, n, n), dtype=complex); Fm = np.zeros_like(Fp)
    for k, z in enumerate(zs):
        X = z * Sh
        R = {1: X @ np.linalg.inv(I - X), -1: -X @ np.linalg.inv(I + X)}
        Cm = (I - X @ X) @ np.linalg.inv(I + X @ X)
        Q = {(s1, s2): -R[s1] @ R[s2] @ Cm for s1 in sg for s2 in sg}
        for i in range(2):
            for j in range(2):
                ot, oh = sg[i], sg[j]
                if ot == -oh: E[k, i, j] += R[ot]
                E[k, i, j] += Q[(ot, -oh)]
        for m, s in enumerate(sg):
            Fp[k, m] = R[s]; Fm[k, m] = R[-s]
    return {'E': E, 'Ediag': np.einsum('kabii->kabi', E), 'F+': Fp, 'F-': Fm,
            'F+diag': np.einsum('ksii->ksi', Fp), 'F-diag': np.einsum('ksii->ksi', Fm)}

# ---------------------------------------------------------------------------------------------------------------
_CORES = {}
def load_cores(excess, vmax=None):
    key = (excess, vmax)
    if key not in _CORES:
        fn = os.path.join(HERE, f'cores_e{excess}.json' if vmax is None else f'cores_e{excess}_v{vmax}.json')
        _CORES[key] = [(c['V'], [tuple(e) for e in c['edges']], c['aut']) for c in json.load(open(fn))]
    return _CORES[key]

def X_general(Sh, excess, zs, vmax=None, cores=None):
    """X_excess(z) for a regular tournament at the points zs (array), summed over cores (or the listed indices)"""
    ops = ops_general(Sh.astype(complex), np.asarray(zs))
    tot = np.zeros(len(zs), dtype=complex)
    for i, (V, edges, aut) in enumerate(load_cores(excess, vmax)):
        if cores is not None and i not in cores: continue
        tot += SignCore(V, edges, aut).contract('general', ops) / aut
    return tot

def series_general(Sh, excess, Lmax, K=32, r=0.5, vmax=None):
    """[z^L] X_excess by a Cauchy sum on |z| = r using X(-z) = X(z), X(conj z) = conj X(z) (K/4 + 1 evaluations)"""
    assert K % 4 == 0
    zs = r * np.exp(2j * np.pi * np.arange(K // 4 + 1) / K)
    v = X_general(Sh, excess, zs, vmax)
    vals = np.array([v[k % (K // 2)] if k % (K // 2) <= K // 4 else np.conj(v[K // 2 - k % (K // 2)]) for k in range(K)])
    return np.array([(np.mean(vals * np.exp(-2j * np.pi * np.arange(K) * L / K)) / r ** L).real for L in range(Lmax + 1)])
