#!/usr/bin/env python3
"""Partition-lattice compiler for the Paley signed self-avoiding walk sums A_L(q).

A_L(q) = sum over injective v_0..v_L in F_q of prod_{i=1..L} chi(v_i - v_(i-1))   (q = 3 mod 4, chi quadratic)
       = sum over set partitions s of {0..L} with no two consecutive points in a block of mu(s) U(path/s),
mu(s) = prod_B (-1)^(|B|-1) (|B|-1)!  (Moebius function of the partition lattice), and U(G) = unrestricted sum over all
maps V(G) -> F_q of prod_edges chi(v_head - v_tail).
U is evaluated by series-parallel reduction in the edge-function algebra F = {a chi + g 1 + b delta}:
  convolution (series):  chi*chi = 1 - q delta, chi*1 = 0, 1*1 = q 1, delta*f = f
  product (parallel):    chi.chi = 1 - delta,   chi.delta = 0, 1.f = f, delta.delta = delta
  pendant vertex: sum_d f(d) = g q + b;  isolated vertex: q;  loop: f(0) = g + b;  reversal: chi(-d) = -chi(d).
When no reduction applies, an edge with non-chi parts is split (chi / delete / contract).  Graphs whose every edge is
pure chi and min degree >= 3 are CORES: their U is computed by brute force per q and recorded separately, so that
A_L = (exact reducible part) + sum_core coefficient x U_core.
All arithmetic is exact (Python integers) at a list of primes q = 3 mod 4.

Session opus-2026-10-14-S35.  Note: 05-knowledge/results/paley_least_resistance_runs_one_modular_form_20261014.md
Reproduce: python3 04-computation/experiments/paley_run_partition_compiler_20261014.py 2,4,6,8,10   (L = 12: about 10 minutes)
Validation: reproduces A_2..A_8 exactly (closed forms of HYP-9274 (C)); A_10 and A_12 vanish at the too-few-points primes.
"""
import sys, math, itertools, time
from collections import defaultdict
import numpy as np
import numba as nb

QS = [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107]
K = len(QS)


# ---- coefficient vectors: tuples of Python ints, one per q ----
def vconst(c):
    return [c] * K


QV = [q for q in QS]
ZERO = vconst(0)
ONE = vconst(1)


def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vneg(a):
    return [-x for x in a]


def vmul(a, b):
    return [x * y for x, y in zip(a, b)]


def vmulq(a):
    return [x * q for x, q in zip(a, QV)]


def is_zero(a):
    return all(x == 0 for x in a)


# edge function f = (alpha, gamma, beta): alpha chi + gamma 1 + beta delta
def f_rev(f):
    a, g, b = f
    return (vneg(a), g, b)


def f_conv(f1, f2):
    a1, g1, b1 = f1
    a2, g2, b2 = f2
    alpha = vadd(vmul(a1, b2), vmul(b1, a2))
    gamma = vadd(vadd(vmul(a1, a2), vmulq(vmul(g1, g2))), vadd(vmul(g1, b2), vmul(b1, g2)))
    beta = vadd(vneg(vmulq(vmul(a1, a2))), vmul(b1, b2))
    return (alpha, gamma, beta)


def f_prod(f1, f2):
    a1, g1, b1 = f1
    a2, g2, b2 = f2
    alpha = vadd(vmul(a1, g2), vmul(g1, a2))
    gamma = vadd(vmul(a1, a2), vmul(g1, g2))
    beta = vadd(vadd(vneg(vmul(a1, a2)), vmul(g1, b2)), vadd(vmul(b1, g2), vmul(b1, b2)))
    return (alpha, gamma, beta)


CHI = (ONE, ZERO, ZERO)


# ---- graph reduction ----
def reduce_graph(vset, edges, core_acc, scale):
    """nv vertices 0..nv-1 (some may be absent), edges: list of [u, v, f]; returns reducible contribution (vector)
    and adds core terms into core_acc[canonical_core] += scale * coefficient"""
    verts = set(vset)
    edges = [list(e) for e in edges]
    factor = ONE
    changed = True
    while True:
        # loops
        newe = []
        for u, v, f in edges:
            if u == v:
                factor = vmul(factor, vadd(f[1], f[2]))
            else:
                newe.append([u, v, f])
        edges = newe
        if is_zero(factor):
            return ZERO
        # merge parallel edges
        bucket = {}
        for u, v, f in edges:
            key = (min(u, v), max(u, v))
            ff = f if u == key[0] else f_rev(f)
            if key in bucket:
                bucket[key] = f_prod(bucket[key], ff)
            else:
                bucket[key] = ff
        edges = []
        for (u, v), f in bucket.items():
            if is_zero(f[0]) and is_zero(f[1]) and is_zero(f[2]):
                return ZERO
            edges.append([u, v, f])
        deg = defaultdict(int)
        for u, v, f in edges:
            deg[u] += 1
            deg[v] += 1
        acted = False
        for x in list(verts):
            if deg[x] == 0:
                factor = vmulq(factor)
                verts.discard(x)
                acted = True
        if acted:
            continue
        # pendant
        for x in list(verts):
            if deg[x] == 1:
                for i, (u, v, f) in enumerate(edges):
                    if u == x or v == x:
                        a, g, b = f
                        factor = vmul(factor, vadd(vmulq(g), b))
                        del edges[i]
                        break
                verts.discard(x)
                acted = True
                break
        if acted:
            if is_zero(factor):
                return ZERO
            continue
        # series
        for x in list(verts):
            if deg[x] == 2:
                inc = [(i, e) for i, e in enumerate(edges) if e[0] == x or e[1] == x]
                (i1, e1), (i2, e2) = inc
                # orient e1 as u -> x, e2 as x -> w
                if e1[1] == x:
                    u, f1 = e1[0], e1[2]
                else:
                    u, f1 = e1[1], f_rev(e1[2])
                if e2[0] == x:
                    w, f2 = e2[1], e2[2]
                else:
                    w, f2 = e2[0], f_rev(e2[2])
                h = f_conv(f1, f2)
                for i in sorted([i1, i2], reverse=True):
                    del edges[i]
                edges.append([u, w, h])
                verts.discard(x)
                acted = True
                break
        if acted:
            continue
        break
    if not edges:
        # remaining vertices are isolated (handled) -> pure scalar
        return vmul(factor, [q ** len(verts) for q in QV]) if verts else factor
    # split a non-pure edge
    for i, (u, v, f) in enumerate(edges):
        a, g, b = f
        if not (is_zero(g) and is_zero(b)):
            rest = edges[:i] + edges[i + 1:]
            tot = ZERO
            if not is_zero(a):
                tot = vadd(tot, reduce_graph(verts, rest + [[u, v, (a, ZERO, ZERO)]], core_acc, vmul(scale, factor)))
            if not is_zero(g):
                tot = vadd(tot, vmul(g, reduce_graph(verts, rest, core_acc, vmul(vmul(scale, factor), g))))
            if not is_zero(b):
                # contract u, v
                rest2 = [[(u if x == v else x), (u if y == v else y), ff] for x, y, ff in rest]
                tot = vadd(tot, vmul(b, reduce_graph(verts - {v}, rest2, core_acc, vmul(vmul(scale, factor), b))))
            return vmul(factor, tot)
    # pure chi core (min degree >= 3): coefficient = prod alpha
    coef = ONE
    simple = []
    for u, v, f in edges:
        coef = vmul(coef, f[0])
        simple.append((u, v))
    key = canonical_core(sorted(verts), simple)
    tot_scale = vmul(vmul(scale, factor), coef)
    core_acc[key[0]] = vadd(core_acc.get(key[0], ZERO), vmul(tot_scale, vconst(key[1])))
    return ZERO


def canonical_core(verts, edges):
    """canonical form of an oriented simple graph up to relabelling, with the orientation sign:
    returns (key, sign) such that U(G) = sign * U(canonical)"""
    n = len(verts)
    idx = {v: i for i, v in enumerate(verts)}
    E = [(idx[u], idx[v]) for u, v in edges]
    best = None
    for perm in itertools.permutations(range(n)):
        # undirected edge set under perm, orientation sign: count edges oriented high->low after relabelling
        lab = []
        rev = 0
        for u, v in E:
            a, b = perm[u], perm[v]
            if a < b:
                lab.append((a, b))
            else:
                lab.append((b, a))
                rev += 1
        cand = tuple(sorted(lab))
        if best is None or cand < best[0]:
            best = (cand, (-1) ** rev)
    return (n, best[0]), best[1]


def set_partitions_no_consecutive(n):
    """restricted growth strings a[0..n-1], a[0] = 0, a[i] <= max+1, a[i] != a[i-1]"""
    a = [0] * n

    def rec(i, m):
        if i == n:
            yield tuple(a)
            return
        for c in range(m + 2):
            if c == a[i - 1]:
                continue
            a[i] = c
            yield from rec(i + 1, max(m, c))
    yield from rec(1, 0)


def mobius(rgs):
    sizes = defaultdict(int)
    for c in rgs:
        sizes[c] += 1
    mu = 1
    for s in sizes.values():
        mu *= (-1) ** (s - 1) * math.factorial(s - 1)
    return mu


def compile_run(L):
    n = L + 1
    total = ZERO
    cores = {}
    count = 0
    for rgs in set_partitions_no_consecutive(n):
        count += 1
        mu = mobius(rgs)
        nv = max(rgs) + 1
        edges = [[rgs[i - 1], rgs[i], CHI] for i in range(1, n)]
        acc = {}
        val = reduce_graph(set(range(nv)), edges, acc, vconst(mu))
        total = vadd(total, vmul(val, vconst(mu)))
        for k, v in acc.items():
            cores[k] = vadd(cores.get(k, ZERO), v)
    return total, cores, count


# ---- brute-force U for pure-chi cores ----
@nb.njit(cache=True)
def core_U(q, chi, nvert, eu, ev):
    """unrestricted sum over all maps of nvert vertices (vertex 0 fixed at 0, times q) of prod chi(v[ev]-v[eu])"""
    tot = 0
    vals = np.zeros(nvert, dtype=np.int64)
    m = nvert - 1
    N = q ** m
    for code in range(N):
        c = code
        for i in range(1, nvert):
            vals[i] = c % q
            c //= q
        p = 1
        for k in range(eu.shape[0]):
            p *= chi[(vals[ev[k]] - vals[eu[k]]) % q]
            if p == 0:
                break
        tot += p
    return tot * q


def chi_table(q):
    chi = np.zeros(q, dtype=np.int64)
    sq = set((x * x) % q for x in range(1, q))
    for x in range(1, q):
        chi[x] = 1 if x in sq else -1
    return chi


if __name__ == '__main__':
    Ls = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [2, 4, 6, 8]
    for L in Ls:
        t0 = time.time()
        total, cores, count = compile_run(L)
        print(f'L = {L}: {count} partitions, {time.time() - t0:.1f} s; reducible part at q = {QS[:6]}...: {total[:6]}')
        print(f'   cores with nonzero coefficient: {sum(1 for v in cores.values() if not is_zero(v))}')
        for k, v in cores.items():
            if is_zero(v):
                continue
            (n, E) = k
            print(f'   core n = {n}, edges = {E}: coefficient at q = {QS[:4]}: {v[:4]}')
