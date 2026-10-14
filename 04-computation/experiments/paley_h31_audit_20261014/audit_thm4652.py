#!/usr/bin/env python3
"""AUDIT of THM-4652 (Paley pattern sums) by brute force.

For oriented multigraphs G on m <= 5 vertices:
  Sigma*(G) = sum over injective v: [m] -> F_q of prod_e chi(v_h(e) - v_t(e))       (q = 7, 11, 19)
  Sigma(G)  = same over all maps.
Tests:
  part 1 (reversal sign), part 2 (odd edge count -> 0), part 3 (single-edge component, odd multiplicity -> 0),
  part 4 (an automorphism of the underlying multigraph reversing an odd number of edges -> 0),
  part 5 (series-parallel => Sigma(G) is ONE polynomial in q for q = 3 mod 4; non-SP control),
  the class-F calculus (alpha chi + gamma + beta delta) symbolically on random SP compositions,
  the 'Reading' (which 5-vertex simple patterns survive parts 2-4 and have a K_4 minor),
  and an exploratory census of zero sums not explained by parts 2-4.
"""
import itertools, random, sys
from fractions import Fraction
import numpy as np
import numba as nb

random.seed(20261014)
QS = [7, 11, 19]


def chi_table(q):
    chi = np.full(q, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, q):
        chi[(x * x) % q] = 1
    return chi


@nb.njit(cache=True)
def pattern_sum(q, chi, m, tails, heads, injective):
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
        if injective:
            ok = True
            for i in range(m):
                for j in range(i + 1, m):
                    if v[i] == v[j]:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue
        p = 1
        for e in range(tails.shape[0]):
            p *= chi[(v[heads[e]] - v[tails[e]]) % q]
            if p == 0:
                break
        tot += p
    return tot


def S(G, m, q, injective=True):
    t = np.array([e[0] for e in G], dtype=np.int64)
    h = np.array([e[1] for e in G], dtype=np.int64)
    return int(pattern_sum(q, chi_table(q), m, t, h, injective))


# ---------------------------------------------------------------- structural predicates
def underlying(G):
    """dict {frozenset(i,j): [count i->j with i<j orientation forward, count backward]}"""
    d = {}
    for (t, h) in G:
        key = (min(t, h), max(t, h))
        f, b = d.get(key, (0, 0))
        if t < h:
            f += 1
        else:
            b += 1
        d[key] = (f, b)
    return d


def odd_reversal_symmetry(G, m):
    U = underlying(G)
    for perm in itertools.permutations(range(m)):
        ok = True
        sign = 0
        for (i, j), (f, b) in U.items():
            pi, pj = perm[i], perm[j]
            key = (min(pi, pj), max(pi, pj))
            if key not in U or sum(U[key]) != f + b:
                ok = False
                break
            # image edges: i->j (f of them) becomes pi->pj ; j->i (b) becomes pj->pi
            # count image edges oriented 'backward' (high -> low) and compare with G's backward count on key
            if pi < pj:
                img_back = b
            else:
                img_back = f
            sign += img_back - U[key][1]
        if ok and sign % 2 == 1:
            return perm
    return None


def components(G, m):
    adj = {i: set() for i in range(m)}
    for t, h in G:
        adj[t].add(h)
        adj[h].add(t)
    seen, comps = set(), []
    for s in range(m):
        if s in seen:
            continue
        st, comp = [s], set()
        while st:
            x = st.pop()
            if x in comp:
                continue
            comp.add(x)
            st.extend(adj[x] - comp)
        seen |= comp
        comps.append(comp)
    return comps


def isolated_edge_component(G, m):
    """component with exactly two vertices joined by an ODD number of edges"""
    for comp in components(G, m):
        if len(comp) == 2:
            i, j = sorted(comp)
            k = sum(1 for (t, h) in G if {t, h} == {i, j})
            if k % 2 == 1:
                return True
    return False


def reduce_mod2(G):
    """parallel classes of odd multiplicity -> single edge with orientation of parity; even -> removed (for Sigma*)"""
    U = underlying(G)
    out = []
    for (i, j), (f, b) in U.items():
        if (f + b) % 2 == 1:
            out.append((i, j) if b % 2 == 0 else (j, i))
    return out


def has_K4_minor(edges, m):
    """simple undirected graph: treewidth <= 2 iff reducible by deleting deg<=1 and suppressing deg-2 vertices"""
    adj = {i: set() for i in range(m)}
    for t, h in edges:
        if t != h:
            adj[t].add(h)
            adj[h].add(t)
    changed = True
    while changed:
        changed = False
        for v in list(adj):
            d = len(adj[v])
            if d <= 1:
                for u in adj[v]:
                    adj[u].discard(v)
                del adj[v]
                changed = True
                break
            if d == 2:
                a, b = list(adj[v])
                adj[a].discard(v)
                adj[b].discard(v)
                adj[a].add(b)
                adj[b].add(a)
                del adj[v]
                changed = True
                break
    return len(adj) > 0


# ---------------------------------------------------------------- class F calculus (part 5)
# element = (alpha, gamma, beta) of polynomials in q (as dict deg->coef);  T(d) = alpha chi(d) + gamma + beta [d=0]
def padd(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}


def pmul(a, b):
    r = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            r[k1 + k2] = r.get(k1 + k2, 0) + v1 * v2
    return {k: v for k, v in r.items() if v}


def pscale(a, c):
    return {k: c * v for k, v in a.items() if c * v}


Q = {1: 1}
ONE = {0: 1}


def series(F, G):
    a1, g1, b1 = F
    a2, g2, b2 = G
    # chi*chi = 1 - q delta ; chi*1 = 0 ; 1*1 = q ; delta*f = f
    alpha = padd(pmul(a1, b2), pmul(b1, a2))
    gamma = padd(padd(pmul(a1, a2), pmul(pmul(g1, g2), Q)), padd(pmul(g1, b2), pmul(b1, g2)))
    beta = padd(pscale(pmul(pmul(a1, a2), Q), -1), pmul(b1, b2))
    return (alpha, gamma, beta)


def parallel(F, G):
    a1, g1, b1 = F
    a2, g2, b2 = G
    # chi.chi = 1 - delta ; chi.delta = 0 ; 1.f = f ; delta.delta = delta
    alpha = padd(pmul(a1, g2), pmul(g1, a2))
    gamma = padd(pmul(a1, a2), pmul(g1, g2))
    beta = padd(padd(pscale(pmul(a1, a2), -1), pmul(g1, b2)), padd(pmul(b1, g2), pmul(b1, b2)))
    return (alpha, gamma, beta)


def peval(a, q):
    return sum(v * q ** k for k, v in a.items())


def random_sp(depth, nxt):
    """returns (edges list on labels, s, t, Fclass, next_label) building series/parallel compositions"""
    if depth == 0 or random.random() < 0.25:
        s, t = nxt, nxt + 1
        if random.random() < 0.5:
            return [(s, t)], s, t, ({0: 1}, {}, {}), nxt + 2
        return [(t, s)], s, t, ({0: -1}, {}, {}), nxt + 2   # reversed edge: chi(s - t) = -chi(t - s)
    E1, s1, t1, F1, nxt = random_sp(depth - 1, nxt)
    E2, s2, t2, F2, nxt = random_sp(depth - 1, nxt)
    if random.random() < 0.5:   # series: identify t1 = s2
        ren = {s2: t1}
        E2 = [(ren.get(a, a), ren.get(b, b)) for a, b in E2]
        t2 = ren.get(t2, t2)
        return E1 + E2, s1, t2, series(F1, F2), nxt
    ren = {s2: s1, t2: t1}      # parallel: identify s2 = s1, t2 = t1
    E2 = [(ren.get(a, a), ren.get(b, b)) for a, b in E2]
    return E1 + E2, s1, t1, parallel(F1, F2), nxt


def relabel(E):
    labs = sorted({x for e in E for x in e})
    mp = {x: i for i, x in enumerate(labs)}
    return [(mp[a], mp[b]) for a, b in E], len(labs)


def main():
    fails = []
    # ---------------- parts 1-4 on random oriented multigraphs, m <= 5
    tested = 0
    zero_unexplained = []
    pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    for trial in range(700):
        m = random.choice([2, 3, 4, 4, 5, 5, 5])
        P = [(i, j) for i in range(m) for j in range(i + 1, m)]
        ne = random.randint(1, min(9, 2 * len(P)))
        G = []
        for _ in range(ne):
            i, j = random.choice(P)
            G.append((i, j) if random.random() < 0.5 else (j, i))
        vals = {q: S(G, m, q) for q in QS}
        tested += 1
        # part 1: reverse one random edge
        k = random.randrange(len(G))
        G2 = list(G)
        G2[k] = (G[k][1], G[k][0])
        for q in QS:
            if S(G2, m, q) != -vals[q]:
                fails.append(('part1', G, q))
        pred0 = False
        if len(G) % 2 == 1:
            pred0 = True
            if any(vals[q] for q in QS):
                fails.append(('part2', G, vals))
        if isolated_edge_component(G, m):
            pred0 = True
            if any(vals[q] for q in QS):
                fails.append(('part3', G, vals))
        perm = odd_reversal_symmetry(G, m)
        if perm is not None:
            pred0 = True
            if any(vals[q] for q in QS):
                fails.append(('part4', G, perm, vals))
        if not pred0 and all(vals[q] == 0 for q in QS):
            zero_unexplained.append((m, G))
    print(f'parts 1-4: {tested} random oriented multigraphs (m <= 5, <= 9 edges) at q = {QS}: failures = {len(fails)}')
    for f in fails[:10]:
        print('   FAIL', f)
    print(f'   zero at all q but not predicted by parts 2-4: {len(zero_unexplained)}; examples: {zero_unexplained[:4]}')
    # explicit part 4 examples
    K5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    K5e = [e for e in K5 if e != (0, 1)]
    K5adj = [e for e in K5 if e not in ((0, 1), (0, 2))]
    K5dis = [e for e in K5 if e not in ((0, 1), (2, 3))]
    for name, G in (('K_5', K5), ('K_5 - e', K5e), ('K_5 - 2 adjacent', K5adj), ('K_5 - 2 disjoint (W_4)', K5dis)):
        print(f'   {name}: odd-reversal symmetry {odd_reversal_symmetry(G, 5)}; Sigma* at {QS} = {[S(G, 5, q) for q in QS]}')
    # ---------------- part 5: class-F calculus vs brute force on random SP graphs (unrestricted sums)
    bad5 = 0
    n5 = 0
    for trial in range(400):
        E, s, t, F, _ = random_sp(4, 0)
        E, m = relabel(E)
        if m > 6:
            continue
        a, g, b = F
        n5 += 1
        for q in (7, 11, 19):
            pred = q * (peval(g, q) * q + peval(b, q))
            got = S(E, m, q, injective=False)
            if pred != got:
                bad5 += 1
                print('   PART5 FAIL', E, q, pred, got)
                break
    print(f'part 5 class-F calculus: {n5} random SP multigraphs (<= 6 vertices), predicted Sigma = q(gamma q + beta) vs brute force at q = 7, 11, 19: failures = {bad5}')
    # polynomiality of unrestricted W_4 sum (non-SP control) and K_4 (vanishes)
    W4 = [(4, 0), (4, 1), (4, 2), (4, 3), (0, 2), (0, 3), (1, 2), (1, 3)]
    print('   unrestricted Sigma(W_4) at q = 7, 11, 19, 23:', [S(W4, 5, q, injective=False) for q in (7, 11, 19, 23)])
    # ---------------- the Reading: all simple graphs on 5 vertices
    seen = {}
    survivors = []
    for mask in range(1 << 10):
        E = [pairs[i] for i in range(10) if mask >> i & 1]
        # canonical form under S_5
        canon = min(tuple(sorted(tuple(sorted((p[a], p[b]))) for a, b in E)) for p in itertools.permutations(range(5)))
        if canon in seen:
            continue
        seen[canon] = E
    print(f'\nReading: {len(seen)} isomorphism classes of simple graphs on 5 vertices')
    for canon, E in sorted(seen.items(), key=lambda kv: len(kv[0])):
        if not has_K4_minor(E, 5):
            continue
        G = [tuple(e) for e in E]           # low -> high orientation
        why = []
        if len(G) % 2:
            why.append('odd size')
        if isolated_edge_component(G, 5):
            why.append('isolated edge')
        if odd_reversal_symmetry(G, 5) is not None:
            why.append('odd-reversal symmetry')
        vals = [S(G, 5, q) for q in QS]
        print(f'   non-SP class with {len(G)} edges {E}: killed by {why if why else "NOTHING"}; Sigma* at {QS} = {vals}')
        if not why:
            survivors.append(E)
    print(f'   survivors (non-SP, not killed): {survivors}')
    print('\nALL THM-4652 CHECKS PASSED' if not fails and not bad5 else '\nSOME CHECKS FAILED')


if __name__ == '__main__':
    main()
