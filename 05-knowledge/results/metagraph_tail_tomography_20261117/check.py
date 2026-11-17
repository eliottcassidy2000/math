from itertools import combinations, permutations, product
from functools import lru_cache
from math import isqrt
import json

def scores(rows):
    return tuple(r.bit_count() for r in rows)

def labelled(n, word):
    rows=[0]*n
    for (i,j),bit in zip(combinations(range(n),2),word):
        u,v=(i,j) if bit else (j,i)
        rows[u]|=1<<v
    return tuple(rows)

def z(m):
    return labelled(m,(j==i+1 for i,j in combinations(range(m),2)))

def reader(U,v,m):
    n=len(U)
    rows=list(U)+[(r<<n)|((1<<n)-1) for r in z(m)]
    rows[n]^=1<<v
    rows[v]|=1<<n
    return tuple(rows)

def pathfun(rows):
    @lru_cache(None)
    def end(S,v):
        R=S^(1<<v)
        if not R:
            return 1
        return sum(end(R,u) for u in range(len(rows)) if R>>u&1 and rows[u]>>v&1)
    @lru_cache(None)
    def H(S):
        if not S:
            return 1
        return sum(end(S,v) for v in range(len(rows)) if S>>v&1)
    return H,end

def profile(rows,v):
    h,end=pathfun(rows)
    n=len(rows)
    S=(1<<n)-1
    ending=end(S,v)
    b=sum(end(A,v)*h(S^A) for A in range(1,S+1) if A>>v&1)
    rev=tuple(sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n))
    start=pathfun(rev)[1](S,v)
    return ending,(h(S),b,start)

def aut(rows,mark=None):
    groups={}
    for i,s in enumerate(scores(rows)):
        groups.setdefault((s,i==mark),[]).append(i)
    cnt=0
    for blocks in product(*(permutations(g) for g in groups.values())):
        p=dict(zip((i for g in groups.values() for i in g),(j for g in blocks for j in g)))
        cnt+=all((rows[i]>>j&1)==(rows[p[i]]>>p[j]&1) for i,j in combinations(range(len(rows)),2))
    return cnt

def check(ok):
    if not ok:
        raise RuntimeError('failed check')

def main():
    h=[1,1,1,3]
    for m in range(4,32):
        h.append(sum(h[-3:]))
    e=[None,1]
    for m in range(2,32):
        e.append(h[m-1]-e[-1])
    t=[sum(h[:m]) for m in range(32)]
    def det(M):
        (a,b,c),(d,e,f),(g,h,i)=M
        return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    for m in range(1,30):
        check(det([(h[j],t[j],-e[j]) for j in range(m,m+3)])==2)
    for n in range(1,16):
        for m in range(1,16):
            check(h[n+m]==h[n]*h[m]+t[n]*t[m]-e[n]*e[m])
    checks=0
    witness={}
    for n in range(1,5):
        for word in product((0,1),repeat=n*(n-1)//2):
            U=labelled(n,word)
            for v in range(n):
                en,(a,b,c)=profile(U,v)
                if not en: continue
                for m in range(1,7):
                    W=reader(U,v,m)
                    N=a*h[m]+b*t[m]-c*e[m]
                    check(N==pathfun(W)[0]((1<<(m+n))-1))
                    check((b-c)%2==0 and b>=a)
                    if m>=4: check(aut(W)==aut(U,v))
                    checks+=1
                    if N<=127 and (N not in witness or (n+m,n,m)<tuple(witness[N][:3])):
                        witness[N]=[n+m,n,m,list(U),v,[a,b,c],list(W)]
    U=(20,13,16,5,10,31)
    check(profile(U,4)[1]==(11,26,0) and profile(U,4)[0]>0)
    W=reader(U,4,2)
    check(pathfun(W)[0](255)==63)
    witness[63]=[8,6,2,list(U),4,[11,26,0],list(W)]
    G={x:[y for y in range(1,33) if y!=x and isqrt(x+y)**2==x+y] for x in range(1,33)}
    triangles=[p for p in combinations(G,3) if p[1] in G[p[0]] and p[2] in G[p[0]] and p[2] in G[p[1]]]
    cycles=[]
    def dfs(P,S):
        if not S:
            if P[0] in G[P[-1]] and P[1]<P[-1]: cycles.append(P)
            return
        for x in G[P[-1]]:
            if x in S: dfs(P+[x],S-{x})
    dfs([1],set(G)-{1})
    check(triangles==[(6,19,30)] and len(cycles)==1)
    print(json.dumps(dict(checked_padded_reads=checks,example_63=witness.get(63),sample_counts_through_63=sorted(k for k in witness if k<=63),square32_triangle=triangles[0],square32_cycle=cycles[0]),indent=2))
if __name__=='__main__': main()
