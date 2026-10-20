#!/usr/bin/env python3
"""Even graphs (Royle-Praeger-Glasby-Freedman-Devillers / von Broemssen) versus tournaments: census for n <= NMAX.
Uses nauty's geng / gentourng for the unlabelled lists and networkx for automorphism groups.
sgn_X(g) = (-1)^{# edges {u<v} of X with g(u) > g(v)}; X is even iff sgn_X = 1 on Aut X.
Reports: counts (even, odd, tournaments), labelled even count, connected even count, rigid counts, |Aut| spectra,
cycle-type symmetry counts t_lambda (tournaments) and e_lambda (even graphs), the stabiliser-containment Hall
deficit (Lemma N1 of THM-4531) with a violator, and the signed even-subgraph count matrix M (n <= 5)."""
import sys, subprocess, itertools
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher, DiGraphMatcher
from collections import Counter, defaultdict
from math import factorial
from sympy import Matrix

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7

def geng(n):
    out = subprocess.run(['geng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    return [nx.from_graph6_bytes(s.encode()) for s in out]

def gentourng(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    Ts = []
    for line in out:
        assert len(line) == n * (n - 1) // 2, (line, n)
        D = nx.DiGraph(); D.add_nodes_from(range(n)); k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[k] == '1': D.add_edge(i, j)
                else: D.add_edge(j, i)
                k += 1
        Ts.append(D)
    return Ts

def cycle_type(g, n):
    seen = [False] * n; ct = []
    for i in range(n):
        if not seen[i]:
            L = 0; j = i
            while not seen[j]:
                seen[j] = True; j = g[j]; L += 1
            ct.append(L)
    return tuple(sorted(ct, reverse=True))

def sign(edges, g):
    s = 0
    for u, v in edges:
        if g[u] > g[v]: s ^= 1
    return -1 if s else 1

def perm_tuple(d, n): return tuple(d[i] for i in range(n))

def conj_contained(AT, AX, n):
    """is some S_n-conjugate of the group AT (set of tuples) contained in AX (set of tuples)?"""
    gens = [g for g in AT if any(x != i for i, x in enumerate(g))]
    for h in itertools.permutations(range(n)):
        hinv = [0] * n
        for i, x in enumerate(h): hinv[x] = i
        ok = True
        for g in gens:
            c = tuple(h[g[hinv[i]]] for i in range(n))
            if c not in AX: ok = False; break
        if ok: return True
    return False

def is_cyclic(A):
    m = len(A)
    for g in A:
        # order of g
        x = g; k = 1; idt = tuple(range(len(g)))
        while x != idt:
            x = tuple(g[i] for i in x); k += 1
        if k == m: return True
    return False

NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 2
for n in range(NMIN, NMAX + 1):
    Gs = geng(n); Ts = gentourng(n)
    ginfo = []
    for G in Gs:
        edges = [(min(u, v), max(u, v)) for u, v in G.edges()]
        auts = [perm_tuple(m, n) for m in GraphMatcher(G, G).isomorphisms_iter()]
        even = all(sign(edges, g) == 1 for g in auts)
        cts = frozenset(cycle_type(g, n) for g in auts)
        ginfo.append(dict(G=G, edges=edges, auts=set(auts), even=even, cts=cts, conn=nx.is_connected(G)))
    tinfo = []
    for T in Ts:
        auts = [perm_tuple(m, n) for m in DiGraphMatcher(T, T).isomorphisms_iter()]
        cts = frozenset(cycle_type(g, n) for g in auts)
        tinfo.append(dict(T=T, auts=set(auts), cts=cts, strong=nx.is_strongly_connected(T)))
    evens = [x for x in ginfo if x['even']]
    odds = [x for x in ginfo if not x['even']]
    print(f"\n===== n = {n}: graphs {len(Gs)}, even {len(evens)}, odd {len(odds)}, tournaments {len(Ts)} =====")
    assert len(evens) == len(Ts)
    lab_even = sum(factorial(n) // len(x['auts']) for x in evens)
    print(f"labelled even graphs {lab_even} of {2**(n*(n-1)//2)}; connected even {sum(1 for x in evens if x['conn'])}; "
          f"connected odd {sum(1 for x in odds if x['conn'])}; strong tournaments {sum(1 for t in tinfo if t['strong'])}")
    print("rigid: even graphs", sum(1 for x in evens if len(x['auts']) == 1), " tournaments", sum(1 for t in tinfo if len(t['auts']) == 1))
    print("|Aut| spectrum even graphs:", sorted(Counter(len(x['auts']) for x in evens).items()))
    print("|Aut| spectrum tournaments :", sorted(Counter(len(t['auts']) for t in tinfo).items()))
    print("even graphs by #edges:", sorted(Counter(len(x['edges']) for x in evens).items()))
    print("sum over even graphs of (-1)^edges =", sum((-1)**len(x['edges']) for x in evens),
          " (THM-4531 G2: should equal #self-converse tournaments)")
    idt = tuple([n])  # placeholder
    lam_t = Counter(); lam_e = Counter(); lam_all = Counter()
    for t in tinfo:
        for ct in t['cts']:
            if ct != tuple([1] * n): lam_t[ct] += 1
    for x in evens:
        for ct in x['cts']:
            if ct != tuple([1] * n): lam_e[ct] += 1
    for x in ginfo:
        for ct in x['cts']:
            if ct != tuple([1] * n): lam_all[ct] += 1
    print("cycle type lambda: t_lambda (tournaments with such an automorphism) / e_lambda (even graphs) / all graphs")
    for ct in sorted(set(lam_t) | set(lam_e), key=lambda c: (-max(c), c)):
        flag = " <-- t > e" if lam_t[ct] > lam_e[ct] else ""
        print(f"  {ct}: {lam_t[ct]} / {lam_e[ct]} / {lam_all[ct]}{flag}")
    # Hall matching (stabiliser containment), nontrivial-Aut tournaments only
    nt = [t for t in tinfo if len(t['auts']) > 1]
    B = nx.Graph()
    for i, t in enumerate(nt):
        B.add_node(('T', i))
        cyc = is_cyclic(t['auts'])
        for j, x in enumerate(evens):
            if len(x['auts']) % len(t['auts']): continue
            if not t['cts'] <= x['cts']: continue
            if cyc or conj_contained(t['auts'], x['auts'], n):
                B.add_edge(('T', i), ('X', j))
    for j in range(len(evens)): B.add_node(('X', j))
    top = [('T', i) for i in range(len(nt))]
    M = nx.bipartite.maximum_matching(B, top_nodes=top) if nt else {}
    matched = sum(1 for u in top if u in M)
    print(f"Hall (Lemma N1): nontrivial-Aut tournaments {len(nt)}, max matched {matched}, deficit {len(nt) - matched}")
    if len(nt) - matched > 0:
        # violator: tournaments reachable from unmatched ones by alternating paths
        unmatched = [u for u in top if u not in M]
        reach = set(unmatched); frontier = list(unmatched); seenX = set()
        while frontier:
            u = frontier.pop()
            for x in B[u]:
                if x in seenX: continue
                seenX.add(x)
                if x in M:
                    v = M[x]
                    if v not in reach: reach.add(v); frontier.append(v)
        print(f"  violator: {len(reach)} tournaments with neighbourhood of size {len(seenX)}")
        print("  violator tournaments (|Aut|, cycle types):",
              sorted(Counter((len(nt[i]['auts']), tuple(sorted(nt[i]['cts'] - {tuple([1]*n)}))) for (_, i) in reach).items()))
        print("  neighbourhood even graphs (|Aut|, edges):",
              sorted(Counter((len(evens[j]['auts']), len(evens[j]['edges'])) for (_, j) in seenX).items()))
    if n <= 5:
        # signed even-subgraph count matrix M[T][X] = Phi_X(T)
        perms = list(itertools.permutations(range(n)))
        rows = []
        for t in tinfo:
            T = t['T']
            tv = {(u, v): (1 if T.has_edge(u, v) else -1) for u in range(n) for v in range(u + 1, n)}
            row = []
            for x in evens:
                E = x['edges']
                copies = {}
                for g in perms:
                    img = frozenset((min(g[u], g[v]), max(g[u], g[v])) for u, v in E)
                    eps = sign(E, g)
                    if img in copies: assert copies[img] == eps, "odd graph?"
                    else: copies[img] = eps
                row.append(sum(eps * (1 if not img else eval('*'.join(str(tv[e]) for e in img))) for img, eps in copies.items()))
            rows.append(row)
        Mx = Matrix(rows)
        print("signed even-subgraph matrix M (rows tournaments, cols even graphs by edges", [len(x['edges']) for x in evens], ")")
        for r in rows: print("  ", r)
        print("  det M =", Mx.det(), " rank", Mx.rank())
