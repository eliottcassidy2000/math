#!/usr/bin/env python3
"""
S635 — the ALTERNATING GROUP GRAPH AG_n = Cay(A_n, {(1 2 i),(1 i 2): 3<=i<=n}), and how it changes
with n, in light of the chi-bridge (S634): both LRC and unit distance are chi(Cay(G, U)); AG_n is the
NON-ABELIAN, PERMANENTLY-ODD rung.
KEY: a Cayley graph Cay(G,S) is BIPARTITE (chi=2) iff there is a homomorphism G->Z/2 with S->1 (a
'sign'/parity grading). S_n has the sign map => Cay(S_n, transpositions) is bipartite (chi=2). A_n is
the KERNEL of sign => no such map => AG_n is non-bipartite (chi>=3). The alternating group is exactly
where parity runs out — the 'odd sector' (canon S587: worry = odd sector).
"""
import itertools
from functools import reduce

def compose(p, q):  # p after q  (p,q tuples; apply q then p) -- use right action consistently
    return tuple(p[q[i]] for i in range(len(p)))

def inv(p):
    r=[0]*len(p)
    for i,pi in enumerate(p): r[pi]=i
    return tuple(r)

def parity(p):
    seen=[False]*len(p); s=0
    for i in range(len(p)):
        if not seen[i]:
            j=i; c=0
            while not seen[j]: seen[j]=True; j=p[j]; c+=1
            s+=c-1
    return s%2

def three_cycle(n,a,b,c):  # the 3-cycle (a b c) as a permutation tuple on {0..n-1}
    p=list(range(n)); p[a],p[b],p[c]=p[b],p[c],p[a]; return tuple(p)

def alternating_group(n):
    return [p for p in itertools.permutations(range(n)) if parity(p)==0]

def AG_generators(n):
    # (1 2 i) and (1 i 2) in 1-indexed -> 0-indexed (0 1 i-1)
    gens=[]
    for i in range(2,n):       # i-1 from 2..n-1 (1-indexed i=3..n)
        gens.append(three_cycle(n,0,1,i))
        gens.append(three_cycle(n,0,i,1))
    return gens

def cayley_graph(group, gens):
    idx={g:k for k,g in enumerate(group)}; N=len(group)
    adj=[set() for _ in range(N)]
    for g in group:
        for s in gens:
            h=compose(s,g)            # left multiplication by generator
            adj[idx[g]].add(idx[h])
    return [sorted(a) for a in adj]

def chromatic_number(adj, ub=8):
    N=len(adj)
    order=sorted(range(N), key=lambda v:-len(adj[v]))
    def colorable(k):
        col=[-1]*N
        def bt(t):
            if t==N: return True
            v=order[t]
            for c in range(k):
                if all(col[u]!=c for u in adj[v]):
                    col[v]=c
                    if bt(t+1): return True
                    col[v]=-1
            return False
        return bt(0)
    for k in range(2,ub+1):
        if colorable(k): return k
    return ub+1

def is_bipartite(adj):
    N=len(adj); col=[-1]*N
    for s in range(N):
        if col[s]==-1:
            col[s]=0; st=[s]
            while st:
                v=st.pop()
                for u in adj[v]:
                    if col[u]==-1: col[u]=col[v]^1; st.append(u)
                    elif col[u]==col[v]: return False
    return True

def odd_girth(adj, cap=9):
    # BFS shortest odd cycle (approx via shortest cycle through each vertex - small graphs)
    from collections import deque
    N=len(adj); best=cap+1
    for s in range(N):
        dist={s:0}; par={s:-1}; dq=deque([s])
        while dq:
            v=dq.popleft()
            for u in adj[v]:
                if u not in dist: dist[u]=dist[v]+1; par[u]=v; dq.append(u)
                elif (dist[v]+dist[u]+1)%2==1: best=min(best,dist[v]+dist[u]+1)
        if best<=3: break
    return best

if __name__=="__main__":
    print("ALTERNATING GROUP GRAPH AG_n vs the bipartite transposition graph of S_n")
    print("="*72)
    for n in range(3,7):
        A=alternating_group(n); gens=AG_generators(n); adj=cayley_graph(A,gens)
        N=len(A); deg=len(adj[0]); bip=is_bipartite(adj)
        if N<=120:
            chi=chromatic_number(adj); og=odd_girth(adj)
        else:
            chi="(skip)"; og="(skip)"
        print(f"  AG_{n}: |A_n|={N} vertices, (2n-4)={deg}-regular, bipartite={bip}, odd-girth={og}, chi={chi}")
    print("\n  contrast: Cay(S_n, transpositions) is bipartite (sign map) => chi=2.")
    # verify S_n transposition graph bipartite for n=4
    Sn=list(itertools.permutations(range(4)))
    transp=[three_cycle.__self__ if False else None]  # placeholder
    def transposition(n,a,b):
        p=list(range(n)); p[a],p[b]=p[b],p[a]; return tuple(p)
    tg=[transposition(4,a,b) for a in range(4) for b in range(a+1,4)]
    adjS=cayley_graph(Sn,tg)
    print(f"  S_4 transposition graph: |S_4|=24, bipartite={is_bipartite(adjS)} (chi=2, the sign 2-coloring)")
    print("\n  => AG_n is the PERMANENTLY-ODD (non-bipartite) rung: A_n = ker(sign), no parity 2-coloring,")
    print("     so chi(AG_n)>=3 for all n. The alternating group is where parity runs out (S587 odd sector).")
