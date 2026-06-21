#!/usr/bin/env python3
"""
half_arc_tournament_search_kpswf5.py -- kind-pasteur-2026-06-21 (Thread C)

THE PRECISE TOURNAMENT STRUCTURE BEHIND HALF-ARC-TRANSITIVITY.

A graph X is HALF-ARC-TRANSITIVE iff Aut(X) is vertex- and edge-transitive but
has exactly TWO arc-orbits. By Tutte's theorem, such X carries an Aut(X)-INVARIANT
ORIENTATION D (an oriented graph / partial tournament): every edge gets a
direction, Aut(X) preserves D, but NO automorphism reverses D. The two arc-orbits
= {arcs of D} and {arcs of D^op}. The converse involution D <-> D^op is the Z_2
that Aut(X) provably CANNOT realize.

This is EXACTLY the tournament converse Z_2 (THM-549/550): the orientation D plays
the role of a tournament, D^op its complement, and "Aut cannot reverse" = "no
relabelling realizes the converse." The half-tiling is the fundamental domain of
that same Z_2 on the FULL tournament; here it appears on a vertex-transitive
oriented graph.

So the right tournament object is NOT a metagraph (those aren't vertex-transitive)
but an Aut-INVARIANT TOURNAMENT/ORIENTATION on a vertex-transitive graph. The
natural family: CIRCULANTS on Z_n (the home of circulant & Paley tournaments).

This script:
  (1) Searches circulant graphs Cay(Z_n, S), S = -S (undirected), for
      half-arc-transitivity, and reads off the invariant orientation.
  (2) For each half-arc circulant, checks: is the invariant orientation a
      vertex-transitive TOURNAMENT (complete) or a proper partial orientation?
  (3) Tests the Paley tournament's place: Paley is ARC-transitive as a digraph
      (too symmetric); its underlying K_p is arc-transitive. Contrast with the
      half-arc case where the orientation is NOT reversible by Aut.

Exact orbit computation; |Aut| capped by using the circulant's own symmetry
(translations + multipliers) plus a VF2 check that this is the FULL Aut for the
small n we test. We verify against the Holt baseline separately (holt_*_kpswf5).

kind-pasteur-2026-06-21-kpswf5
"""
import sys
from itertools import combinations
from math import gcd
import networkx as nx
sys.stdout.reconfigure(line_buffering=True)

def banner(t):
    print("\n" + "="*72 + "\n  " + t + "\n" + "="*72)

# ---- exact orbit machinery (VF2 full Aut; only call on n<=20-ish graphs) ----
def full_aut(G):
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    return [dict(mm) for mm in GM.isomorphisms_iter()]

def transitivity(G, auts=None, cap=200000):
    if auts is None:
        auts = full_aut(G)
    nodes = list(G.nodes())
    # vertex orbits
    rep = {v: v for v in nodes}
    def f(x):
        while rep[x] != x: rep[x]=rep[rep[x]]; x=rep[x]
        return x
    def un(a,b):
        ra,rb=f(a),f(b)
        if ra!=rb: rep[ra]=rb
    for g in auts:
        for v in nodes: un(v,g[v])
    nvorb = len({f(v) for v in nodes})
    # edge orbits
    edges=[frozenset(e) for e in G.edges()]
    repe={e:e for e in edges}
    def fe(x):
        while repe[x]!=x: repe[x]=repe[repe[x]]; x=repe[x]
        return x
    def ue(a,b):
        ra,rb=fe(a),fe(b)
        if ra!=rb: repe[ra]=rb
    for g in auts:
        for e in edges:
            a,b=tuple(e); ue(e,frozenset((g[a],g[b])))
    neorb=len({fe(e) for e in edges})
    # arc orbits
    arcs=[]
    for a,b in G.edges(): arcs+=[(a,b),(b,a)]
    repa={a:a for a in arcs}
    def fa(x):
        while repa[x]!=x: repa[x]=repa[repa[x]]; x=repa[x]
        return x
    def ua(a,b):
        ra,rb=fa(a),fa(b)
        if ra!=rb: repa[ra]=rb
    for g in auts:
        for (a,b) in arcs: ua((a,b),(g[a],g[b]))
    naorb=len({fa(a) for a in arcs})
    return len(auts), nvorb, neorb, naorb

def circulant(n, S):
    G = nx.Graph(); G.add_nodes_from(range(n))
    for x in range(n):
        for s in S:
            G.add_edge(x, (x+s)%n)
    return G

# ---------------------------------------------------------------------------
banner("(1) Circulant search for HALF-ARC-TRANSITIVE graphs on Z_n")
print("  Cay(Z_n,S), S=-S, no 0. Half-arc = VT & ET & arc-orbits=2.")
print("  (Circulants are always vertex-transitive via translation.)\n")
print("  %-4s %-22s %-6s %-4s %-4s %-4s  %s" %
      ("n","connection S (half)","|Aut|","vO","eO","aO","verdict"))
half_arc_circulants=[]
for n in range(5, 28):
    # symmetric connection sets: pick a set of "positive" diffs in 1..n//2,
    # S = those plus their negatives. valency = 2*|pos| (or +1 if n/2 included)
    pos_all = list(range(1, n//2+1))
    # only test connected, 4-regular and 6-regular (Holt-like even valency)
    for k in [2, 3]:               # |pos| -> valency 2k (4-regular, 6-regular)
        for pos in combinations(pos_all, k):
            if n % 2 == 0 and (n//2) in pos:
                continue  # avoid the self-paired diff (gives odd contribution)
            # build S
            S = set()
            for p in pos:
                S.add(p % n); S.add((-p) % n)
            G = circulant(n, S)
            if not nx.is_connected(G): continue
            degs=set(dict(G.degree()).values())
            if len(degs)!=1: continue
            val=next(iter(degs))
            if val not in (4,6): continue
            naut,vo,eo,ao = transitivity(G)
            if vo==1 and eo==1 and ao==2:
                half_arc_circulants.append((n, sorted(pos), val, naut))
                print("  %-4d %-22s %-6d %-4d %-4d %-4d  HALF-ARC-TRANSITIVE <<<" %
                      (n, str(sorted(pos)), naut, vo, eo, ao))
if not half_arc_circulants:
    print("  (no half-arc-transitive circulant found for 4/6-regular up to n=27)")

# ---------------------------------------------------------------------------
banner("(2) Paley tournament / graph: the OPPOSITE extreme (arc-transitive)")
def qr(p): return sorted({(x*x)%p for x in range(1,p)})
for p in [5, 7, 11, 13]:
    Q=qr(p); sym = all(((p-r)%p) in Q for r in Q)
    if sym:  # Paley GRAPH (p=1 mod 4)
        G=circulant(p, set(Q))
        naut,vo,eo,ao = transitivity(G)
        print("  p=%d Paley GRAPH: |Aut|=%d vO=%d eO=%d aO=%d  -> %s"
              % (p, naut, vo, eo, ao,
                 "ARC-TRANSITIVE" if ao==1 else ("half-arc" if ao==2 else "?")))
    else:    # Paley TOURNAMENT (p=3 mod 4): a vertex-transitive DIGRAPH
        # underlying graph is K_p (arc-transitive). The DIGRAPH itself:
        D=nx.DiGraph(); D.add_nodes_from(range(p))
        for x in range(p):
            for s in Q: D.add_edge(x,(x+s)%p)
        # digraph arc-transitivity: Aut of the digraph = AGL-type p(p-1)/2 group
        DM=nx.algorithms.isomorphism.DiGraphMatcher(D,D)
        dauts=[dict(mm) for mm in DM.isomorphisms_iter()]
        arcs=list(D.edges())
        repa={a:a for a in arcs}
        def fa(x):
            while repa[x]!=x: repa[x]=repa[repa[x]]; x=repa[x]
            return x
        def ua(a,b):
            ra,rb=fa(a),fa(b)
            if ra!=rb: repa[ra]=rb
        for g in dauts:
            for (a,b) in arcs: ua((a,b),(g[a],g[b]))
        naorb=len({fa(a) for a in arcs})
        # underlying K_p
        UG=nx.Graph(); UG.add_edges_from((u,v) for u,v in D.edges())
        print("  p=%d Paley TOURNAMENT: |Aut(digraph)|=%d, digraph-arc-orbits=%d "
              "(arc-transitive DIGRAPH=%s); underlying = K_%d (arc-transitive graph)"
              % (p, len(dauts), naorb, naorb==1, p))

print("""
  CONTRAST: Paley is a vertex-transitive tournament whose orientation IS
  reversible by an automorphism family in the AGL sense at the DIGRAPH level it
  is arc-transitive (a single arc-orbit). Its underlying complete graph K_p is
  arc-transitive. So Paley is MAXIMALLY symmetric -- the converse Z_2 IS realized
  (every Paley tournament is self-converse, THM/MEMORY). That is the OPPOSITE of
  the Holt phenomenon, where the converse Z_2 is provably NOT realized by Aut.""")

# ---------------------------------------------------------------------------
banner("(3) The INVARIANT ORIENTATION of a half-arc circulant = a tournament-like object")
print("""
  For each half-arc-transitive circulant found above, the two arc-orbits define
  an Aut-invariant orientation D (pick one orbit as 'forward'). D is an oriented
  Cayley graph Cay(Z_n, S^+) with S^+ a 'positive' half of S. The converse D^op
  uses -S^+. Aut(X) preserves D as an orientation but cannot send D->D^op.
  This D is the EXACT analog of a tournament under its converse Z_2.""")
for (n, pos, val, naut) in half_arc_circulants:
    # build S and find S^+ : the half that the orientation uses.
    S=set()
    for p in pos: S.add(p%n); S.add((-p)%n)
    # An invariant orientation: x -> x+s for s in S^+ where S^+ chosen so the
    # multiplier group fixes S^+ setwise. The multipliers fixing S setwise:
    mults=[a for a in range(1,n) if gcd(a,n)==1 and {(a*s)%n for s in S}==S]
    # does any multiplier send S^+ -> -S^+ (=reverse)? if NONE, orientation rigid.
    # pick S^+ = pos (the positive diffs)
    Splus=set()
    for p in pos: Splus.add(p%n)
    reversers=[a for a in mults if {(a*s)%n for s in Splus}=={(-t)%n for t in Splus}]
    preservers=[a for a in mults if {(a*s)%n for s in Splus}==Splus]
    print("  n=%d pos=%s val=%d: |mult group fixing S|=%d, "
          "preserve orientation=%d, reverse orientation=%d -> %s"
          % (n, pos, val, len(mults), len(preservers), len(reversers),
             "ORIENTATION RIGID (tournament-like, converse NOT realized)"
             if len(reversers)==0 else "converse realized (symmetric)"))

print("\nDone. kind-pasteur-2026-06-21-kpswf5")
