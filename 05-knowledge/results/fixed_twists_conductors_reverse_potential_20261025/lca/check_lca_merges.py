"""Two small exact facts for the LCA/Petersen merge note (stdlib only).
(1) The clone-fibre polynomial of C3[TT2,TT2,TT2] (the Gale tournament of a linear K6 with three linked triangle pairs,
    THM-4529 Theorem L; the canonical clone hostile of CLONE_INVOLUTION) is 1+t+t^2+t^3 with depth triple (1 | mod3:1 | mod2:3).
(2) The crossing graph of the ten 2-subsets of a 5-set (two subsets cross iff they meet in exactly one point) is L(K5), the
    complement of the Petersen graph, with chromatic number 5; by Session56's hierarchy-cover theorem the ten pairs need
    exactly five laminar families.
"""
import json, itertools
from math import comb
def need(ok, m):
    if not ok: raise RuntimeError(m)
# (1) regular C3 on three gaps
perms=[tuple((i+c)%3 for i in range(3)) for c in range(3)]
def cycles(g):
    seen=set(); out=[]
    for i in range(3):
        if i in seen: continue
        j=i; c=0
        while j not in seen: seen.add(j); c+=1; j=g[j]
        out.append(c)
    return out
tot=[0]*4
for g in perms:
    f=[1]
    for L in cycles(g):
        fac=[0]*(L+1); fac[0]=fac[-1]=1
        nf=[0]*(len(f)+L)
        for i,x in enumerate(f):
            for j,y in enumerate(fac): nf[i+j]+=x*y
        f=nf
    tot=[a+b for a,b in zip(tot,f)]
P=[c//3 for c in tot]; need(all(c%3==0 for c in tot) and P==[1,1,1,1],'P=1+t+t^2+t^3')
Q=[sum(P[k]*comb(k,j)*(-1)**(k-j) for k in range(j,4)) for j in range(4)]
def depth(Q,l=None):
    for j,c in enumerate(Q):
        if (c!=0) if l is None else (c%l!=0): return j
d=(depth(Q),depth(Q,3),depth(Q,2)); need(d==(1,1,3),'depth triple')
# (2) crossing graph of 2-subsets of [5]
pairs=list(itertools.combinations(range(5),2))
adj={p:[q for q in pairs if q!=p and len(set(p)&set(q))==1] for p in pairs}
need(all(len(adj[p])==6 for p in pairs),'L(K5) is 6-regular')
petersen_adj={p:[q for q in pairs if q!=p and not set(p)&set(q)] for p in pairs}
need(all(len(petersen_adj[p])==3 for p in pairs),'complement is Petersen (3-regular)')
def colourable(k):
    order=pairs; col={}
    def bt(i):
        if i==len(order): return True
        p=order[i]
        for c in range(k):
            if all(col.get(q)!=c for q in adj[p]):
                col[p]=c
                if bt(i+1): return True
                del col[p]
        return False
    return bt(0)
need(not colourable(4) and colourable(5),'chi(L(K5))=5')
print(json.dumps({'clone_polynomial_C3_TT2_cubed':P,'depth_triple_(Q,mod3,mod2)':d,'chi_crossing_graph_of_2subsets_of_5':5,'crossing_graph':'L(K5) = complement of Petersen'}))
