#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Alcuin number of the OCF conflict graph Ω(T), Hamiltonicity of Ω, planarity (Kuratowski),
and the graph→tournament tiling map G↦T_G.  kind-pasteur-2026-06-16-S1.

Alcuin(G) ∈ {τ(G), τ(G)+1} (Csorba–Hurkens–Woeginger, SIAM J Discrete Math 24(3) 2010):
 small-boat ⟺ Alcuin=τ ; large-boat ⟺ Alcuin=τ+1.
 Lemma 4.3: TWO distinct maximum stable sets ⟹ small-boat.
 Structure thm 3.1: feasible at b ⟺ ∃ stable X=X1∪X2∪X3, nonempty Y1,Y2⊆Y=V−X, |Y|≤b,
   X1∪Y1 & X2∪Y2 stable, |Y1|+|Y2|≥|X3|.
For the conflict graph: V(Ω)=directed odd cycles of T, edge = share ≥1 T-vertex.
 α(Ω)=max vertex-disjoint odd-cycle packing ν_odd ≤ ⌊n/3⌋ (each cycle ≥3 vtx) — SMALL.
 τ(Ω)=|V(Ω)|−α(Ω)=(#odd cycles)−ν_odd (Gallai). H(T)=I(Ω,2)=Σ α_k 2^k.
"""
import sys, itertools
from math import comb
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
import networkx as nx
from networkx.algorithms.planarity import check_planarity

def directed_odd_cycles(A, n):
    """all directed odd cycles (3,5,7,...) as (frozenset(verts), canonical rotation tuple)."""
    cyc=[]
    for L in range(3, n+1, 2):
        for verts in itertools.combinations(range(n), L):
            vs=list(verts)
            for perm in itertools.permutations(vs[1:]):
                seq=[vs[0]]+list(perm)
                if all(A[seq[i]][seq[(i+1)%L]] for i in range(L)):
                    # canonical: rotate to smallest-first to dedupe rotations of same directed cycle
                    rots=[tuple(seq[k:]+seq[:k]) for k in range(L)]
                    cyc.append((frozenset(verts), min(rots)))
    # dedupe identical directed cycles
    seen=set(); out=[]
    for fs,rot in cyc:
        if rot not in seen: seen.add(rot); out.append((fs,rot))
    return out

def build_omega(cycles):
    G=nx.Graph(); G.add_nodes_from(range(len(cycles)))
    for i in range(len(cycles)):
        for j in range(i+1,len(cycles)):
            if cycles[i][0] & cycles[j][0]: G.add_edge(i,j)
    return G

def max_indep_packings(cycles):
    """all MAXIMUM vertex-disjoint odd-cycle packings (independent sets in Ω). Returns (alpha, count)."""
    m=len(cycles)
    if m==0: return 0,1
    best=1; sets=[frozenset([i]) for i in range(m)]
    # greedily grow by size since alpha<=floor(n/3) small
    cur=[[i] for i in range(m)]
    allmax=[[i] for i in range(m)]; alpha=1
    k=1
    while True:
        nxt=[]
        for S in cur:
            mx=max(S)
            for j in range(mx+1,m):
                if all(not (cycles[j][0]&cycles[i][0]) for i in S):
                    nxt.append(S+[j])
        if not nxt: break
        alpha=len(nxt[0]); allmax=nxt; cur=nxt; k+=1
    return alpha, len(allmax)

def H_count(A,n):
    c=0
    for p in itertools.permutations(range(n)):
        if all(A[p[t]][p[t+1]] for t in range(n-1)): c+=1
    return c

def stable(G, S):
    S=list(S)
    return all(not G.has_edge(S[i],S[j]) for i in range(len(S)) for j in range(i+1,len(S)))

def is_small_boat(G):
    """decide small-boat (Alcuin=τ) via Lemma 4.3 (>=2 max stable sets) then structure thm at b=τ.
       Only called when |V| small enough; uses α(G)≤⌊n/3⌋ smallness for Ω."""
    n=G.number_of_nodes()
    if n==0: return True, 0, 0   # trivial
    # alpha via complement clique (small alpha)
    comp=nx.complement(G)
    # max clique in complement = max independent set in G
    alpha=0; maxsets=[]
    try:
        for clq in nx.find_cliques(comp):
            if len(clq)>alpha: alpha=len(clq); maxsets=[frozenset(clq)]
            elif len(clq)==alpha: maxsets.append(frozenset(clq))
    except Exception:
        return None, None, None
    # maxsets from find_cliques are maximal; filter to maximum & dedupe
    maxsets=[s for s in set(maxsets) if len(s)==alpha]
    tau=n-alpha
    if len(maxsets)>=2:
        return True, tau, alpha   # Lemma 4.3
    # unique max stable set X; test structure theorem at b=τ (X=that set, Y=V−X)
    X=set(maxsets[0]); Y=set(G.nodes())-X
    Xl=list(X)
    # need X=X1∪X2∪X3 (disjoint), Y1,Y2⊆Y nonempty, X1∪Y1 & X2∪Y2 stable, |Y1|+|Y2|≥|X3|
    # since |X|=alpha small, brute over 3-colorings of X into (X1,X2,X3)
    found=False
    for coloring in itertools.product(range(3), repeat=len(Xl)):
        X1={Xl[i] for i in range(len(Xl)) if coloring[i]==0}
        X2={Xl[i] for i in range(len(Xl)) if coloring[i]==1}
        X3={Xl[i] for i in range(len(Xl)) if coloring[i]==2}
        if not stable(G,X1) or not stable(G,X2): continue
        # need Y1: X1∪Y1 stable, Y2: X2∪Y2 stable, nonempty, |Y1|+|Y2|>=|X3|
        # candidate Y-vertices for Yi = those in Y nonadjacent to all of Xi
        c1=[y for y in Y if all(not G.has_edge(y,x) for x in X1)]
        c2=[y for y in Y if all(not G.has_edge(y,x) for x in X2)]
        # also Y1 itself must be stable with X1 (Y1 internal stable too since X1∪Y1 stable)
        # we just need existence of nonempty stable Y1⊆c1, Y2⊆c2 with |Y1|+|Y2|>=|X3|.
        # max stable subset sizes:
        def max_stable_subset(cands):
            if not cands: return 0
            # small: max independent set within G[cands]
            sub=G.subgraph(cands); cc=nx.complement(sub)
            best=0
            for clq in nx.find_cliques(cc): best=max(best,len(clq))
            return best
        m1=max_stable_subset(c1); m2=max_stable_subset(c2)
        if m1>=1 and m2>=1 and (m1+m2)>=len(X3):
            found=True; break
    return found, tau, alpha

def hamiltonicity(G):
    """return (connected, has_ham_path, has_ham_cycle) — brute for small, Dirac/Ore sufficient for large."""
    n=G.number_of_nodes()
    conn = (n<=1) or nx.is_connected(G)
    if n==0: return True, True, True
    if n==1: return True, True, False
    if not conn: return conn, False, False
    if n<=9:
        nodes=list(G.nodes()); hp=False; hc=False
        for perm in itertools.permutations(nodes[1:]):
            seq=[nodes[0]]+list(perm)
            if all(G.has_edge(seq[i],seq[i+1]) for i in range(n-1)):
                hp=True
                if G.has_edge(seq[-1],seq[0]): hc=True; break
        return conn, hp, hc
    # large: Dirac (min deg >= n/2 ⟹ Ham cycle ⟹ Ham path)
    mind=min(dict(G.degree()).values())
    if mind*2>=n: return conn, True, True
    return conn, None, None  # unknown

def analyze_tournament(A,n,label):
    cyc=directed_odd_cycles(A,n); Om=build_omega(cyc)
    m=Om.number_of_nodes()
    alpha,nmax = max_indep_packings(cyc)
    tau=m-alpha
    H=H_count(A,n)
    sb=is_small_boat(Om) if m<=120 else (True if nmax>=2 else (None,)*3)
    small,tau2,alpha2 = sb if isinstance(sb,tuple) else (None,None,None)
    conn,hp,hc = hamiltonicity(Om)
    planar,_=check_planarity(Om) if m<=2000 else (None,None)
    alc = tau if small else (tau+1 if small is False else None)
    print(f"{label}: |Ω|={m} (odd cycles), α(Ω)=ν_odd={alpha} (#maxpackings={nmax}), τ(Ω)={tau}, "
          f"H={H} | Alcuin={alc} ({'small' if small else 'large' if small is False else '?'}-boat) | "
          f"Ω: connected={conn}, HamPath={hp}, HamCycle={hc}, planar={planar}")
    return dict(m=m,alpha=alpha,tau=tau,H=H,small=small,conn=conn,hp=hp,hc=hc,planar=planar,nmax=nmax)

def iso_reps(n):
    """one representative per iso class via canonical bit-string (small n)."""
    pairs=list(itertools.combinations(range(n),2)); seen={}
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        best=None
        for perm in itertools.permutations(range(n)):
            key=tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if best is None or key<best: best=key
        seen.setdefault(best,A)
    return list(seen.values())

if __name__=="__main__":
    print("="*70); print("CONFLICT-GRAPH ALCUIN ANALYSIS (per iso class, n=3..6; + special)"); print("="*70)
    for n in (3,4,5,6):
        reps=iso_reps(n)
        print(f"\n--- n={n}: {len(reps)} iso classes ---")
        large=0; nonplanar=0; hamcyc=0
        for A in reps:
            r=analyze_tournament(A,n,f"  T")
            if r['small'] is False: large+=1
            if r['planar'] is False: nonplanar+=1
            if r['hc']: hamcyc+=1
        print(f"  SUMMARY n={n}: large-boat(+1)={large}/{len(reps)}, Ω-nonplanar={nonplanar}, Ω-HamCycle={hamcyc}")
    # special tournaments
    print("\n--- special ---")
    paley7=[[1 if (j-i)%7 in(1,2,4) else 0 for j in range(7)] for i in range(7)]
    for i in range(7): paley7[i][i]=0
    analyze_tournament(paley7,7,"Paley T7")
    reg5=[[1 if (j-i)%5 in(1,2) else 0 for j in range(5)] for i in range(5)]
    for i in range(5): reg5[i][i]=0
    analyze_tournament(reg5,5,"regular T5")
