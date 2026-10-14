#!/usr/bin/env python3
"""AUDIT of HYP-9274 (C)/(D): does a SINGLE run already carry a non-series-parallel (W_4) term?

Extends audit_run_polys_symbolic.py: irreducible pure-chi cores are no longer an error; each is canonicalised
(undirected graph up to relabelling, orientation sign relative to the canonical low->high orientation; cores with an
odd-reversal automorphism or complete cores are 0 by THM-4652 part 4) and its coefficient (a polynomial in q) is
collected.  Result: Sigma*(paths) = P(q) + sum_cores C_core(q) * Sigma(core).
Then: A_10 evaluated with Sigma(W_4) = -q(q-1)(a_q + q) (unrestricted wheel sum, referee) at several primes,
and checked against an independent brute force of A_10 at q = 11 (all 11! orderings... via DFS) .
"""
import sys, itertools, time
from math import factorial
import numpy as np
import numba as nb

sys.setrecursionlimit(100000)
sys.path.insert(0, '.')
from audit_run_polys_symbolic import (padd, pmul, pscale, ONE, Q, ZERO, CHI, rev, conv, pprod, is_zero_label,
                                      partitions_no_adjacent, peval, pstr)


def radd(A, B):
    R = dict(A)
    for k, v in B.items():
        R[k] = padd(R.get(k, ZERO), v)
        if not R[k]:
            del R[k]
    return R


def rmul(A, p):
    R = {}
    for k, v in A.items():
        w = pmul(v, p)
        if w:
            R[k] = w
    return R


def canon_core(verts, edges):
    """edges (u, v, alpha) pure chi.  Returns (key, sign, coeff) or None if the core vanishes by symmetry."""
    vl = sorted(verts)
    n = len(vl)
    coeff = ONE
    oriented = []
    for (u, v, T) in edges:
        coeff = pmul(coeff, T[0])
        oriented.append((u, v))
    und = {(min(u, v), max(u, v)) for (u, v) in oriented}
    if len(und) == n * (n - 1) // 2:
        return None
    best = None
    for perm in itertools.permutations(range(n)):
        mp = {vl[i]: perm[i] for i in range(n)}
        form = tuple(sorted((min(mp[u], mp[v]), max(mp[u], mp[v])) for (u, v) in oriented))
        if best is None or form < best[0]:
            best = (form, mp)
    form, mp = best
    # odd-reversal automorphism check on the canonical form
    fs = set(form)
    for perm in itertools.permutations(range(n)):
        if all((min(perm[a], perm[b]), max(perm[a], perm[b])) in fs for (a, b) in form):
            rv = sum(1 for (a, b) in form if perm[a] > perm[b])
            if rv % 2 == 1:
                return None
    sign = 1
    for (u, v) in oriented:
        if mp[u] > mp[v]:
            sign = -sign
    return form, sign, coeff


def unrestricted(verts, edges):
    verts = set(verts)
    edges = list(edges)
    factor = ONE
    while True:
        changed = False
        ne = []
        for (u, v, T) in edges:
            if u == v:
                factor = pmul(factor, padd(T[1], T[2]))
                changed = True
            else:
                ne.append((u, v, T))
        edges = ne
        if not factor:
            return {}
        groups = {}
        for (u, v, T) in edges:
            key = (min(u, v), max(u, v))
            TT = T if u == key[0] else rev(T)
            if key in groups:
                groups[key] = pprod(groups[key], TT)
                changed = True
            else:
                groups[key] = TT
        if any(is_zero_label(T) for T in groups.values()):
            return {}
        edges = [(k[0], k[1], T) for k, T in groups.items()]
        deg = {v: 0 for v in verts}
        for (u, v, T) in edges:
            deg[u] += 1
            deg[v] += 1
        for v in list(verts):
            if deg[v] == 0:
                factor = pmul(factor, Q)
                verts.discard(v)
                changed = True
        if changed:
            continue
        done = False
        for v in list(verts):
            if deg[v] == 1:
                for i, (a, b, T) in enumerate(edges):
                    if a == v or b == v:
                        factor = pmul(factor, padd(pmul(T[1], Q), T[2]))
                        edges.pop(i)
                        break
                verts.discard(v)
                done = True
                break
        if done:
            if not factor:
                return {}
            continue
        for w in list(verts):
            if deg[w] == 2:
                inc = [(i, e) for i, e in enumerate(edges) if e[0] == w or e[1] == w]
                (i1, e1), (i2, e2) = inc
                u = e1[0] if e1[1] == w else e1[1]
                T1 = e1[2] if e1[1] == w else rev(e1[2])
                v = e2[1] if e2[0] == w else e2[0]
                T2 = e2[2] if e2[0] == w else rev(e2[2])
                edges = [e for j, e in enumerate(edges) if j not in (i1, i2)]
                edges.append((u, v, conv(T1, T2)))
                verts.discard(w)
                done = True
                break
        if done:
            continue
        break
    if not edges:
        return {'P': pmul(factor, {len(verts): 1})}
    for i, (u, v, T) in enumerate(edges):
        a, g, b = T
        if g or b:
            rest = edges[:i] + edges[i + 1:]
            tot = {}
            if a:
                tot = radd(tot, unrestricted(verts, rest + [(u, v, (a, ZERO, ZERO))]))
            if g:
                tot = radd(tot, rmul(unrestricted(verts, rest), g))
            if b:
                rest2 = [((u if x == v else x), (u if y == v else y), S) for (x, y, S) in rest]
                tot = radd(tot, rmul(unrestricted(verts - {v}, rest2), b))
            return rmul(tot, factor)
    c = canon_core(verts, edges)
    if c is None:
        return {}
    form, sign, coeff = c
    return {form: pscale(pmul(coeff, factor), sign)}


def injective_sum_paths(lengths):
    edges, forbid, off = [], {}, 0
    for L in lengths:
        for i in range(L + 1):
            forbid[off + i] = [off + i - 1] if i > 0 else []
        for i in range(L):
            edges.append((off + i, off + i + 1))
        off += L + 1
    total, cnt = {}, 0
    for part in partitions_no_adjacent(off, forbid):
        cnt += 1
        blk, mu = {}, 1
        for bi, b in enumerate(part):
            for x in b:
                blk[x] = bi
            mu *= (-1) ** (len(b) - 1) * factorial(len(b) - 1)
        s = unrestricted(set(range(len(part))), [(blk[a], blk[b], CHI) for a, b in edges])
        if s:
            total = radd(total, {k: pscale(v, mu) for k, v in s.items()})
    return total, cnt


# ---------------------------------------------------------------- numerics
def chi_table(q):
    chi = np.full(q, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, q):
        chi[(x * x) % q] = 1
    return chi


@nb.njit(cache=True)
def core_sum(q, chi, m, tails, heads):
    tot = 0
    v = np.zeros(m, dtype=np.int64)
    total = 1
    for _ in range(m):
        total *= q
    for code in range(total):
        c = code
        for i in range(m):
            v[i] = c % q
            c //= q
        p = 1
        for e in range(tails.shape[0]):
            p *= chi[(v[heads[e]] - v[tails[e]]) % q]
            if p == 0:
                break
        tot += p
    return tot


@nb.njit(cache=True)
def run_sum_dfs(q, chi, L):
    """A_L by DFS over injective sequences v_0 = 0, v_1 = 1 (scaling + translation; L even), times q(q-1)"""
    used = np.zeros(q, dtype=np.int64)
    seq = np.zeros(L + 1, dtype=np.int64)
    nxt = np.zeros(L + 1, dtype=np.int64)
    sgn = np.zeros(L + 1, dtype=np.int64)
    seq[0] = 0
    seq[1] = 1
    used[0] = 1
    used[1] = 1
    sgn[1] = chi[1]
    depth = 2
    nxt[2] = 0
    tot = 0
    while depth >= 2:
        if depth == L + 1:
            tot += sgn[L]
            depth -= 1
            used[seq[depth]] = 0
            nxt[depth] = seq[depth] + 1
            continue
        found = False
        x = nxt[depth]
        while x < q:
            if used[x] == 0:
                found = True
                break
            x += 1
        if not found:
            depth -= 1
            if depth >= 2:
                used[seq[depth]] = 0
                nxt[depth] = seq[depth] + 1
            continue
        seq[depth] = x
        used[x] = 1
        sgn[depth] = sgn[depth - 1] * chi[(x - seq[depth - 1]) % q]
        depth += 1
        if depth <= L:
            nxt[depth] = 0
    return tot * q * (q - 1)


def eta_ap(N):
    P = [0] * (N + 1)
    P[0] = 1
    for n in range(1, N + 1):
        for step in (2 * n, 4 * n):
            if step > N:
                continue
            for _ in range(4):
                for i in range(N, step - 1, -1):
                    P[i] -= P[i - step]
    return [0] + P[:N]


def main():
    t0 = time.time()
    res, cnt = injective_sum_paths([10])
    print(f'A_10: {cnt} partitions, {time.time() - t0:.1f} s')
    print('   polynomial part P(q) =', pstr(res.get('P', {})))
    cores = [k for k in res if k != 'P']
    print(f'   non-vanishing non-series-parallel cores: {len(cores)}')
    for k in cores:
        print(f'      core edges {k} (canonical, low->high): coefficient C(q) = {pstr(res[k])}')
    a = eta_ap(200)
    vals = {}
    mmax = max([1 + max(max(e) for e in k) for k in cores] + [0])
    qlist = (7, 11, 19, 23, 31, 43, 47, 59) if mmax <= 5 else (7, 11, 19, 23)
    for q in qlist:
        chi = chi_table(q)
        tot = peval(res.get('P', {}), q)
        for k in cores:
            m = 1 + max(max(e) for e in k)
            s = int(core_sum(q, chi, m, np.array([e[0] for e in k], dtype=np.int64), np.array([e[1] for e in k], dtype=np.int64)))
            tot += peval(res[k], q) * s
        vals[q] = tot
    print('   A_10(q) = P + sum C*Sigma(core):', vals)
    for q in (11, 19):
        chi = chi_table(q)
        t1 = time.time()
        b = int(run_sum_dfs(q, chi, 10))
        print(f'   independent DFS A_10 at q = {q}: {b}  (agree: {b == vals[q]}; {time.time() - t1:.1f} s)')
    # polynomiality test: is A_10 - P(q) proportional to a_q ?
    for q in [x for x in qlist if x >= 19]:
        print(f'   q = {q}: A_10 - P(q) = {vals[q] - peval(res.get("P", {}), q)};  a_q = {a[q]}')


if __name__ == '__main__':
    main()
