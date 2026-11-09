"""Brute-force validation of the construction lemmas on explicit graphs.
Graph: list of adjacency bitmasks.  d(G) by subset enumeration; evenness (von Bromssen) by enumerating
Aut(G) with a backtracking search and counting reversed edges for each automorphism."""
import random, sys, itertools
sys.setrecursionlimit(10000)
def dcount(adj):
    n=len(adj); cl=[adj[v]|(1<<v) for v in range(n)]; full=(1<<n)-1; c=0
    # incremental: dom[S] = OR of cl over S
    dom=[0]*(1<<n)
    for S in range(1,1<<n):
        b=(S&-S).bit_length()-1; dom[S]=dom[S&(S-1)]|cl[b]
        if dom[S]==full: c+=1
    return c
def automorphisms(adj):
    n=len(adj); deg=[bin(a).count('1') for a in adj]
    order=sorted(range(n),key=lambda v:(-deg[v],v))
    res=[]; img=[-1]*n; used=[False]*n
    def bt(i):
        if i==n: res.append(img[:]); return
        v=order[i]
        for w in range(n):
            if used[w] or deg[w]!=deg[v]: continue
            ok=True
            for k in range(i):
                u=order[k]
                if ((adj[v]>>u)&1)!=((adj[w]>>img[u])&1): ok=False;break
            if ok:
                img[v]=w; used[w]=True; bt(i+1); used[w]=False; img[v]=-1
    bt(0); return res
def is_even(adj):
    n=len(adj); edges=[(u,v) for u in range(n) for v in range(u+1,n) if (adj[u]>>v)&1]
    for g in automorphisms(adj):
        rev=sum(1 for (u,v) in edges if g[u]>g[v])
        if rev%2: return False
    return True
def pad(adj,k=1): return adj+[0]*k
def cone(adj):
    n=len(adj); new=[a|(1<<n) for a in adj]; new.append((1<<n)-1); return new
def e2join(adj):
    n=len(adj); full=(1<<n)-1
    new=[a|(1<<n)|(1<<(n+1)) for a in adj]; new.append(full); new.append(full); return new
def ph(adj,u,l):
    n=len(adj); h=n; new=[a for a in adj]
    hmask=0
    for v in range(n):
        if v!=u: new[v]|=1<<h; hmask|=1<<v
    for i in range(l): hmask|=1<<(n+1+i)
    new.append(hmask)
    for i in range(l): new.append(1<<h)
    return new
def flags(adj):
    n=len(adj); deg=[bin(a).count('1') for a in adj]
    univ=[deg[v]==n-1 for v in range(n)]
    pend=[v for v in range(n) if deg[v]==1]
    pnu=[v for v in pend if not univ[(adj[v]&-adj[v]).bit_length()-1]]
    r3=[v for v in range(n) if deg[v]==2]
    return dict(U=any(univ),pnu=pnu,r3=r3,good=(len(pnu)>=2 and len(r3)>=1))
def rand_even(n,p=0.4,tries=200):
    for _ in range(tries):
        adj=[0]*n
        for u in range(n):
            for v in range(u+1,n):
                if random.random()<p: adj[u]|=1<<v; adj[v]|=1<<u
        if is_even(adj): return adj
    return None
if __name__=="__main__":
    random.seed(1)
    fails=0; tests=0
    for trial in range(60):
        n=random.randint(3,7); G=rand_even(n)
        if G is None: continue
        dG=dcount(G); f=flags(G)
        if not f['U'] and n>=2:                       # Lemma 2 (cone)
            C=cone(G); tests+=1
            if dcount(C)!=2**n+dG or not is_even(C): fails+=1; print('cone fail',G)
        J=e2join(G); tests+=1                         # Lemma 3 (E_2-join)
        if dcount(J)!=3*2**n-2+dG or not is_even(J): fails+=1; print('e2 fail',G)
        deg=[bin(a).count('1') for a in G]            # Lemma 4 (hub move), under the unique-max-degree hypothesis
        for u in range(n):
            r=deg[u]+1
            for l in (1,2):
                hdeg=n-1+l; others=max([deg[v]+1 for v in range(n) if v!=u]+[deg[u],1])
                if not hdeg>others: continue
                H=ph(G,u,l); tests+=1
                ok_d = dcount(H)==dG+(2**r-1)*2**(n+l-r)
                ok_e = is_even(H)
                if not(ok_d and ok_e): fails+=1; print('ph fail',G,u,l,ok_d,ok_e)
    print('lemma tests',tests,'fails',fails)
