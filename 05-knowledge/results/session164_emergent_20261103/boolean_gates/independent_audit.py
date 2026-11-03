#!/usr/bin/env python3
"""Independent audit. No imports from the submitted core or check modules.
Run: python3 /tmp/s164_audit/independent.py [also valid with -O].
"""
from itertools import permutations, combinations, product
from functools import lru_cache
from collections import Counter, deque
from math import comb

def need(ok, what):
    if not ok: raise RuntimeError(what)

def tour(n, code):
    a=[[False]*n for _ in range(n)]
    for i,(u,v) in enumerate(combinations(range(n),2)):
        a[u][v]=bool(code>>i&1); a[v][u]=not a[u][v]
    return tuple(map(tuple,a))

def fromrows(rows):
    return tuple(tuple(bool(r>>v&1) for v in range(len(rows))) for r in rows)

def apex(a, mask):
    n=len(a)
    return tuple(tuple(a[u])+(not bool(mask>>u&1),) for u in range(n)) + (tuple(bool(mask>>u&1) for u in range(n))+(False,),)

def sub(a, vs): return tuple(tuple(a[u][v] for v in vs) for u in vs)

def H(a):
    n=len(a)
    if not n:return 1
    @lru_cache(None)
    def extend(last, remaining):
        if not remaining:return 1
        return sum(extend(v,remaining^(1<<v)) for v in range(n) if remaining>>v&1 and a[last][v])
    return sum(extend(v,((1<<n)-1)^(1<<v)) for v in range(n))

def perm_h(a):
    return sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in permutations(range(len(a))))

def direct_covers(a):
    n=len(a);q=[[0]*n for _ in range(n)]
    # Each ordered nonempty two-path cover is one permutation plus one split.
    for p in permutations(range(n)):
        bad=[not a[u][v] for u,v in zip(p,p[1:])]
        if sum(bad)>1:continue
        for split in range(1,n):
            if sum(bad)==0 or bad[split-1]:q[p[split-1]][p[split]]+=1
    return q

def sccs(a):
    n=len(a);r=[list(x) for x in a]
    for i in range(n):r[i][i]=True
    for k in range(n):
        for i in range(n):
            for j in range(n):r[i][j] = r[i][j] or r[i][k] and r[k][j]
    cs=set(tuple(j for j in range(n) if r[i][j] and r[j][i]) for i in range(n))
    return sorted(cs,key=lambda c:-sum(a[c[0]][d[0]] for d in cs if d!=c))

def cut(c,D):return sum(c[u][v] for u in range(len(c)) for v in range(len(c)) if D>>u&1 and not D>>v&1)

def mincut(c,L,R):
    n=len(c);s=n;t=n+1; cap=[[0]*(n+2) for _ in range(n+2)]
    for u in range(n):
        for v in range(n):cap[u][v]=c[u][v]
    M=sum(map(sum,c))+1
    for v in range(n):
        if L>>v&1:cap[s][v]=M
        if R>>v&1:cap[v][t]=M
    ans=0
    while True:
        prev=[None]*(n+2);prev[s]=s;q=deque([s])
        while q and prev[t] is None:
            u=q.popleft()
            for v in range(n+2):
                if prev[v] is None and cap[u][v]>0:prev[v]=u;q.append(v)
        if prev[t] is None:return ans
        f=M;v=t
        while v!=s:u=prev[v];f=min(f,cap[u][v]);v=u
        v=t
        while v!=s:u=prev[v];cap[u][v]-=f;cap[v][u]+=f;v=u
        ans+=f

def test_carrier(a, test_cuts=False, brute_children=False):
    n=len(a);old=perm_h(a);need(H(a)==old,('H direct',a))
    q=direct_covers(a);c=[[0]*n for _ in range(n)]
    for u,v in combinations(range(n),2):
        k=q[u][v]+q[v][u]
        need(k%2==0,('odd K',a,u,v))
        need(k//2>=perm_h(sub(a,[i for i in range(n) if i not in [u,v]])),('floor',a,u,v))
        c[u][v]=k//2 if a[u][v] else 0;c[v][u]=k//2 if a[v][u] else 0
    children=[];strong_children=[]
    for D in range(1<<n):
        child=apex(a,D);got=H(child);children.append(got)
        if brute_children:need(got==perm_h(child),('child DP',a,D))
        need(got==old+2*cut(c,D),('carrier',a,D))
        need((got==old)==not_any_out(a,D),('equality geometry',a,D))
        if len(sccs(child))==1:strong_children.append(got)
    cs=sccs(a)
    need(sum(h==old for h in children)==len(cs)+1,('k+1',a))
    predicted=[]
    for u in cs[0]:
        for v in cs[-1]:
            if u!=v:predicted.append(old+2*mincut(c,1<<u,1<<v))
    need(min(predicted,default=None)==min(strong_children,default=None),('strong mincut',a))
    if test_cuts:
        for terminal in product(range(3),repeat=n):
            L=sum(1<<i for i,z in enumerate(terminal) if z==1)
            R=sum(1<<i for i,z in enumerate(terminal) if z==2)
            want=min(children[D] for D in range(1<<n) if D&L==L and D&R==0)
            need(old+2*mincut(c,L,R)==want,('forced mincut',a,L,R))
    return children

def not_any_out(a,D):return not any(a[u][v] for u in range(len(a)) for v in range(len(a)) if D>>u&1 and not D>>v&1)

def mobius(vals):
    c=list(vals);N=(len(c)-1).bit_length()
    for i in range(N):
        for S in range(1<<N):
            if S>>i&1:c[S]-=c[S^(1<<i)]
    return c

def forest(edges):
    degree=Counter();components=[]
    for u,v in edges:
        degree[u]+=1;degree[v]+=1
        if max(degree.values())>2:return False
        U=next((C for C in components if u in C),{u})
        V=next((C for C in components if v in C),{v})
        if U is V:return False
        components=[C for C in components if C is not U and C is not V]+[U|V]
    return True

def variable_values(n, edges):
    base=[list(row) for row in tour(n,0)];vals=[]
    for D in range(1<<len(edges)):
        a=[list(row) for row in base]
        for i,(u,v) in enumerate(edges):a[u][v]=bool(D>>i&1);a[v][u]=not a[u][v]
        vals.append(perm_h(a))
    return vals

def check_support(n,edges,S=None):
    vals=variable_values(n,edges); coeff=mobius(vals)
    lam=max(sum(mask>>i&1 for i in range(len(edges))) for mask in range(1<<len(edges)) if forest([e for i,e in enumerate(edges) if mask>>i&1]))
    for mask,c in enumerate(coeff):
        chosen=[e for i,e in enumerate(edges) if mask>>i&1]
        need(c==0 or forest(chosen),('support',chosen,c))
        if S is not None:
            r=sum(u in S and v in S for u,v in chosen)
            need(c==0 or len(chosen)<=2*len(S)-r,('internal edges refinement',chosen,c))
    for m in [3,4,6,9,15]:
        for c in range(m):
            support=sum((x-c)%m!=0 for x in vals)
            need(not support or support>=2**(len(edges)-lam),('scalar group',m,c,support))
    # A genuinely vector-valued control: (H(d)-a,H(complement(d))-b).
    for m in [3,4,6]:
        for a,b in product(range(m),repeat=2):
            support=sum((vals[D]-a)%m!=0 or (vals[(len(vals)-1)^D]-b)%m!=0 for D in range(len(vals)))
            need(not support or support>=2**(len(edges)-lam),('vector group',m,a,b,support))
    return len(edges),lam,sum(c!=0 for c in coeff)

def ordinal(blocks):
    n=sum(map(len,blocks));a=[[False]*n for _ in range(n)];p=0
    for block in blocks:
        for u,row in enumerate(block):
            for v,val in enumerate(row):a[p+u][p+v]=val
            for v in range(p+len(block),n):a[p+u][v]=True
        p+=len(block)
    return tuple(map(tuple,a))

def wrap(a):
    n=len(a)
    return tuple(tuple(row)+(True,False) for row in a)+(tuple([False]*n+[False,True]),tuple([True]*n+[False,False]))

def v3(x):
    k=0
    while x%3==0:x//=3;k+=1
    return k

def main():
    for n in range(1,6):
        for code in range(1<<comb(n,2)):
            test_carrier(tour(n,code),test_cuts=n<=4,brute_children=n<=4)
        print('Independent permutation covers/carrier/floor/SCC/strong mincut: all',1<<comb(n,2),'order',n,'tournaments',flush=True)
    for code in [0,37,4097,32767]:test_carrier(tour(6,code),test_cuts=True)
    print('Forced terminal mincuts: every L/R assignment for all n<=4 and four n=6 controls',flush=True)
    for n,edges,S in [(4,list(combinations(range(4),2)),None),(5,[(u,v) for u,v in combinations(range(5),2) if u in (0,1) or v in (0,1)],{0,1})]:
        print('Independent support / composite / vector controls:',n,check_support(n,edges,S),flush=True)
    gates={}
    for x,y in product(range(2),repeat=2):
        g=apex(fromrows((8,1,3,6)),2+x+4*y);gates[x,y]=g;t=y*(1-x)
        need(perm_h(g)==9+2*t,('gate H',x,y))
        q=direct_covers(g);need(2*perm_h(g)+sum(map(sum,q))==84+16*t,('gate F2',x,y))
        h0=perm_h(g);g1=apex(g,1<<4);h1=perm_h(g1);h2=perm_h(apex(g1,1<<5))
        need((h0,(h1-h0)//2,(h2-h1)//2)==(9+2*t,4+3*x+5*y,7+2*x+4*y),('observer',x,y))
    for s in (1,2):
        hist=Counter()
        for bits in product(range(2),repeat=2*s):
            blocks=[gates[bits[2*i],bits[2*i+1]] for i in range(s)]
            if s%2:blocks.append(((False,),))
            z=wrap(ordinal(blocks));r=s-sum(bits[2*i+1]*(1-bits[2*i]) for i in range(s));k=s-r
            val=H(z);need(len(sccs(z))==1,('strong wrapper',s,bits))
            need(val==9**r*11**k+2**(s%2)*84**r*100**k,('count wrapper',s,bits,val))
            need(v3(val)==r,('v3',s,bits,val));hist[val]+=1
        print('Independent literal wrapped family:',s,dict(sorted(hist.items())),flush=True)
    for s in range(1,41):
        for k in range(s+1):
            need(v3(9**(s-k)*11**k+2**(s%2)*84**(s-k)*100**k)==s-k,('arithmetic valuation',s,k))
    print('PASS: all independent tests; all-order assertions still rest on proofs, not finite controls.',flush=True)

if __name__=='__main__':main()
