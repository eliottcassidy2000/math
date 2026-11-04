# sanity check of the disjoint-orbit lemma on all tournaments of order <= 8: in Aut(T) (odd order), supports of
# prime-length cycles are pairwise equal or disjoint, and supports of cycles of different prime lengths are disjoint.
import itertools, pynauty, subprocess, sys
def closure(gens, n):
    ident=tuple(range(n)); seen={ident}; fr=[ident]
    while fr:
        nxt=[]
        for g in fr:
            for h in gens:
                x=tuple(h[g[i]] for i in range(n))
                if x not in seen: seen.add(x); nxt.append(x)
        fr=nxt
    return seen
def cycles(g):
    n=len(g); seen=[False]*n; out=[]
    for i in range(n):
        if seen[i]: continue
        c=[]; j=i
        while not seen[j]: seen[j]=True; c.append(j); j=g[j]
        if len(c)>1: out.append(frozenset(c))
    return out
def isprime(m): return m>1 and all(m%d for d in range(2,int(m**0.5)+1))
bad=0; tested=0; overlaps=0
for n in range(3,9):
    lines=subprocess.run(["gentourng","-q",str(n)],capture_output=True,text=True).stdout.split()
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for line in lines:
        adj={v:[] for v in range(n)}
        for (i,j),b in zip(pairs,line):
            if b=='1': adj[i].append(j)
            else: adj[j].append(i)
        gens=pynauty.autgrp(pynauty.Graph(n,directed=True,adjacency_dict=adj))[0]
        G=closure([tuple(g) for g in gens],n) if gens else {tuple(range(n))}
        assert len(G)%2==1
        single=[(len(c[0]),c[0]) for g in G for c in [cycles(g)] if len(c)==1 and isprime(len(c[0]))]
        tested+=len(single)
        for (r,O),(s,O2) in itertools.combinations(set(single),2):
            if O&O2:
                overlaps+=1
                if not (O==O2 and r==s): bad+=1; print("VIOLATION n=%d"%n, r, sorted(O), s, sorted(O2))
print("orders 3-8 checked: prime-cycle elements tested =", tested, "overlapping pairs =", overlaps, "violations =", bad)
