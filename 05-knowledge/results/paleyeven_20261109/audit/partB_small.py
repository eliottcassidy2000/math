#!/usr/bin/env python3
"""Audit P1 (independent, brute force): Part B lemmas on small graphs.

Parity here is computed by BRUTE FORCE over all n! permutations (no search tree, no WL): g is an automorphism iff it
maps the edge set onto itself; it reverses edge {u, v} (u < v) iff g(u) > g(v). A graph is odd iff some automorphism
reverses an odd number of edges.

usage: python3 partB_small.py census_K511.txt trees.txt
Checks:
  (P) brute-force parity = census parity for every connected class with n <= 7;
  (B1) Lemma B1: threshold-graph realisers of every odd k <= 1023 on ceil(log2(k+1)) vertices (direct count);
  (B2.1) every connected class (n <= 11) with adjacent twins is odd (census parity);
  (B2.2) unions of two connected classes, n1 + n2 <= 7: parity(G u H) = parity(G) and parity(H);
  (B2.3) cones over every graph on 2..6 vertices with no universal vertex: parity(K1 v G) = parity(G);
         and K1 v G for G WITH a universal vertex (outside the lemma) for information;
  (B2.4) joins A v G, |A| + |G| <= 7, complements connected, A not isomorphic to G: parity = parity(A) and parity(G);
         A v A: odd iff A odd or |A| odd;
  (B2.5) sgn_G(g) sgn_coG(g) = sgn(g) for random permutations;
  (B2.6) trees (from 'dom_census trees 11'): odd iff bicentral with isomorphic halves (AHU encoding).
"""
import sys
import itertools
import random
from collections import defaultdict


def parse_edges(s):
    return [(int(s[i], 16), int(s[i + 1], 16)) for i in range(0, len(s), 2)]


def adjsets(n, edges):
    A = [set() for _ in range(n)]
    for u, v in edges:
        A[u].add(v)
        A[v].add(u)
    return A


def brute_parity(n, edges):
    """returns (is_odd, |Aut|) by brute force over all permutations."""
    E = set((min(u, v), max(u, v)) for u, v in edges)
    deg = [0] * n
    for u, v in E:
        deg[u] += 1
        deg[v] += 1
    naut = 0
    odd = False
    for g in itertools.permutations(range(n)):
        if any(deg[g[i]] != deg[i] for i in range(n)):
            continue
        ok = True
        rev = 0
        for u, v in E:
            a, b = g[u], g[v]
            if (min(a, b), max(a, b)) not in E:
                ok = False
                break
            if a > b:
                rev += 1
        if ok:
            naut += 1
            if rev % 2:
                odd = True
    return odd, naut


def canon_brute(n, edges):
    E = [(min(u, v), max(u, v)) for u, v in edges]
    best = None
    for g in itertools.permutations(range(n)):
        f = tuple(sorted((min(g[u], g[v]), max(g[u], g[v])) for u, v in E))
        if best is None or f < best:
            best = f
    return best


def connected(n, edges):
    if n == 0:
        return True
    A = adjsets(n, edges)
    seen = {0}
    st = [0]
    while st:
        x = st.pop()
        for y in A[x]:
            if y not in seen:
                seen.add(y)
                st.append(y)
    return len(seen) == n


def complement(n, edges):
    E = set((min(u, v), max(u, v)) for u, v in edges)
    return [(u, v) for u in range(n) for v in range(u + 1, n) if (u, v) not in E]


def count_dom(n, edges):
    A = adjsets(n, edges)
    cl = [(1 << v) | sum(1 << w for w in A[v]) for v in range(n)]
    full = (1 << n) - 1
    c = 0
    for S in range(1 << n):
        dm = 0
        s = S
        while s:
            v = (s & -s).bit_length() - 1
            s &= s - 1
            dm |= cl[v]
        if dm == full:
            c += 1
    return c


def main():
    census, treefile = sys.argv[1], sys.argv[2]
    classes = defaultdict(list)   # n -> list of (edges, parity)
    allc = []
    for line in open(census):
        p = line.split()
        n = int(p[0])
        edges = parse_edges(p[6]) if len(p) > 6 else []
        allc.append((n, edges, p[3], int(p[2])))
        if n <= 7:
            classes[n].append((edges, p[3]))
    # (P)
    bad = 0
    tot = 0
    for n in range(1, 8):
        for edges, par in classes[n]:
            odd, _ = brute_parity(n, edges)
            tot += 1
            if ('O' if odd else 'E') != par:
                bad += 1
    print(f"(P) brute-force parity vs census parity, connected classes n <= 7: {tot} checked, {bad} disagree")
    # (B1)
    def threshold(k):
        """graph with d = k (k odd) on ceil(log2(k+1)) vertices via cones and isolated vertices."""
        N = (k + 1 - 1).bit_length() if k > 0 else 1
        N = max(1, (k).bit_length())   # smallest N with k < 2^N
        def build(k, N):  # returns edge list on vertices 0..N-1 with d = k
            if N == 1:
                assert k == 1
                return []
            if k < 2 ** (N - 1):
                return build(k, N - 1)          # pad with isolated vertex N-1
            x = k - 2 ** (N - 1)
            E = build(x, N - 1)
            return E + [(i, N - 1) for i in range(N - 1)]   # cone vertex N-1
        return N, build(k, N)
    badb1 = 0
    for k in range(1, 1024, 2):
        N, E = threshold(k)
        import math
        if N != math.ceil(math.log2(k + 1)) or count_dom(N, E) != k:
            badb1 += 1
    print(f"(B1) threshold realisers on ceil(log2(k+1)) vertices for odd k <= 1023: {badb1} failures")
    # (B2.1) adjacent twins => odd, over the whole census
    tw = 0
    viol = 0
    for n, edges, par, d in allc:
        A = adjsets(n, edges)
        has = any((A[u] | {u}) == (A[v] | {v}) for u in range(n) for v in A[u] if v > u)
        if has:
            tw += 1
            if par != 'O':
                viol += 1
    print(f"(B2.1) census classes with adjacent twins: {tw}; of these not odd: {viol}")
    # build all graphs (up to iso) on n <= 6 from connected classes and complements of connected classes
    allg = {}
    for n in range(1, 7):
        seen = {}
        for edges, par in classes[n]:
            for E in (edges, complement(n, edges)):
                key = canon_brute(n, E)
                if key not in seen:
                    seen[key] = E
        allg[n] = list(seen.values())
    print("    graphs up to isomorphism on 1..6 vertices (A000088 check: 1, 2, 4, 11, 34, 156):",
          [len(allg[n]) for n in range(1, 7)])
    par_cache = {}
    def par_of(n, E):
        key = (n, canon_brute(n, E)) if n <= 6 else None
        if key and key in par_cache:
            return par_cache[key]
        r = brute_parity(n, E)[0]
        if key:
            par_cache[key] = r
        return r
    # (B2.2) unions
    conn = [(n, edges) for n in range(1, 7) for edges, par in classes[n]]
    nb = 0
    nt = 0
    for (n1, E1), (n2, E2) in itertools.combinations_with_replacement(conn, 2):
        if n1 + n2 > 7:
            continue
        E = E1 + [(u + n1, v + n1) for u, v in E2]
        nt += 1
        if brute_parity(n1 + n2, E)[0] != (par_of(n1, E1) or par_of(n2, E2)):
            nb += 1
    print(f"(B2.2) unions of two connected classes with n1 + n2 <= 7: {nt} tested, {nb} failures")
    # (B2.3) cones
    nb = 0
    nt = 0
    info = defaultdict(int)
    for n in range(2, 7):
        for E in allg[n]:
            A = adjsets(n, E)
            univ = any(len(A[v]) == n - 1 for v in range(n))
            C = E + [(v, n) for v in range(n)]   # apex = vertex n
            pc = brute_parity(n + 1, C)[0]
            pg = par_of(n, E)
            if not univ:
                nt += 1
                if pc != pg:
                    nb += 1
            else:
                info[(pg, pc)] += 1
    print(f"(B2.3) cones K1 v G over graphs on 2..6 vertices without universal vertex: {nt} tested, {nb} failures")
    print(f"       (outside the hypothesis, G with a universal vertex: (G odd?, cone odd?) counts {dict(info)})")
    # (B2.4) joins
    nb = 0
    nt = 0
    nself = 0
    nselfbad = 0
    cocon = [(n, E) for n in range(1, 7) for E in allg[n] if connected(n, complement(n, E))]
    for (n1, E1), (n2, E2) in itertools.combinations_with_replacement(cocon, 2):
        if n1 + n2 > 7:
            continue
        J = E1 + [(u + n1, v + n1) for u, v in E2] + [(u, n1 + v) for u in range(n1) for v in range(n2)]
        pj = brute_parity(n1 + n2, J)[0]
        iso = (n1 == n2 and canon_brute(n1, E1) == canon_brute(n2, E2))
        if not iso:
            nt += 1
            if pj != (par_of(n1, E1) or par_of(n2, E2)):
                nb += 1
        else:
            nself += 1
            if pj != (par_of(n1, E1) or (n1 % 2 == 1)):
                nselfbad += 1
    print(f"(B2.4) joins of co-connected A, G (|A|+|G| <= 7), A not iso G: {nt} tested, {nb} failures; "
          f"A v A: {nself} tested, {nselfbad} failures of 'odd iff A odd or |A| odd'")
    # (B2.5)
    rnd = random.Random(7)
    nb = 0
    for t in range(2000):
        n = rnd.randint(2, 9)
        E = [(u, v) for u in range(n) for v in range(u + 1, n) if rnd.random() < 0.5]
        g = list(range(n))
        rnd.shuffle(g)
        rev = lambda EE: sum(1 for u, v in EE if g[u] > g[v])
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if g[i] > g[j])
        if (rev(E) + rev(complement(n, E)) - inv) % 2:
            nb += 1
    print(f"(B2.5) sgn_G(g) sgn_coG(g) = sgn(g), 2000 random (graph, permutation): {nb} failures")
    # (B2.6) trees
    def centres(n, A):
        deg = [len(A[v]) for v in range(n)]
        alive = set(range(n))
        leaves = [v for v in range(n) if deg[v] <= 1]
        rem = n
        while rem > 2:
            nl = []
            for v in leaves:
                alive.discard(v)
                rem -= 1
                for w in A[v]:
                    if w in alive:
                        deg[w] -= 1
                        if deg[w] == 1:
                            nl.append(w)
            leaves = nl
        return sorted(alive)
    def ahu(v, parent, A):
        return '(' + ''.join(sorted(ahu(w, v, A) for w in A[v] if w != parent)) + ')'
    nb = 0
    nt = 0
    byn = defaultdict(lambda: [0, 0])
    for line in open(treefile):
        p = line.split()
        n = int(p[0])
        E = parse_edges(p[6]) if len(p) > 6 else []
        A = adjsets(n, E)
        c = centres(n, A)
        crit = False
        if len(c) == 2:
            a, b = c
            crit = ahu(a, b, A) == ahu(b, a, A)
        nt += 1
        byn[n][0] += 1
        if p[3] == 'O':
            byn[n][1] += 1
        if crit != (p[3] == 'O'):
            nb += 1
    print(f"(B2.6) trees up to {max(byn)} vertices: {nt} tested, {nb} failures of 'odd iff bicentral with isomorphic halves'")
    print("       trees per order (A000055) / odd:", {n: tuple(v) for n, v in sorted(byn.items())})


if __name__ == '__main__':
    main()
