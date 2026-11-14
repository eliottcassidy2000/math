"""Literal certificates for the induced-hull construction."""
from random import Random

def need(cond, *w):
    if not cond: raise AssertionError(w)

def degree(n, es):
    return max([sum(x==u for x,y in es) for u in range(n)] +
               [sum(y==u for x,y in es) for u in range(n)])

def color(n, es, k):
    # Complete the tail/head incidence graph to a regular bipartite multigraph.
    W=[[int((x,y) in es) for y in range(n)] for x in range(n)]
    d=[k-sum(row) for row in W]
    e=[k-sum(W[x][y] for x in range(n)) for y in range(n)]
    I=[x for x in range(n) for _ in range(d[x])]
    J=[y for y in range(n) for _ in range(e[y])]
    for x,y in zip(I,J,strict=True):W[x][y]+=1
    uncolored=set(es)
    for _ in range(k):
        into={}
        def push(x,S):
            for y in range(n):
                if W[x][y] and y not in S:
                    S.add(y)
                    if y not in into or push(into[y],S):
                        into[y]=x
                        return True
            return False
        for x in range(n):need(push(x,set()),W,x)
        F=set()
        for y,x in into.items():
            W[x][y]-=1
            if (x,y) in uncolored:
                uncolored.remove((x,y));F.add((x,y))
        yield F
    need(not uncolored,uncolored)

def cover(n, es):
    k=degree(n,es)
    chains=[]
    private=[]
    counter=n
    for F in color(n,es,k):
        succ=dict(F)
        heads={y for x,y in F}
        C=[]
        B=[]
        for u in range(n):
            if u in heads:continue
            if C:
                B.append(counter); C.append(counter);counter+=1
            while True:
                C.append(u)
                if u not in succ:break
                u=succ[u]
        chains.append(C);private.append(B)
    z=counter
    rank=list(range(n,z+1))+list(range(n))
    common=set(range(n,z))
    rank_idx={x:i for i,x in enumerate(rank)}
    paths=[sorted(common.difference(B))+[z]+C for C,B in zip(chains,private)]
    edge={(min(i,j),max(i,j)):(j,i) if (j,i) in es else (i,j)
          for i in range(n) for j in range(i+1,n)}
    used={}
    for p in paths:
        need(sorted(p)==list(range(z+1)),n,z,p)
        for x,y in zip(p,p[1:]):
            key=tuple(sorted((x,y)))
            need(used.get(key,(x,y))==(x,y),key)
            need(edge.get(key,(x,y))==(x,y),key)
            used[key]=(x,y)
    for x,y in es:need(used.get((y,x))==(x,y),x,y)
    augmented={(min(x,y),max(x,y)):used.get((min(x,y),max(x,y)),(x,y))
               for i,x in enumerate(rank) for y in rank[i+1:]}
    need(all(augmented[p]==a for p,a in edge.items()),edge)
    need(z+1==n+1+k*(n-1)-len(es),n,z,k)
    need(all(used.get(p)==(x,y) for p,(x,y) in augmented.items()
             if rank_idx[x]>rank_idx[y]),augmented)

def main():
    rng=Random(652)
    runs=0
    for n in range(2,19):
        for rate in [0.04,0.2,0.5,0.8,1]:
            for _ in range(3):
                es={(j,i) for j in range(n) for i in range(j) if rng.random()<rate}
                if es:cover(n,es);runs+=1
    print(f'{runs} literal induced hull embeddings: PASS')

if __name__=='__main__':main()
