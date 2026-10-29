"""Independent stdlib audit. No imports from the session's computational lanes."""
from itertools import permutations, product, combinations
from collections import Counter
from pathlib import Path
import json

def need(p, message):
    if not p: raise RuntimeError(message)

def covers(A,k):
    # Grow ordered paths vertex by vertex, rather than permuting and cutting.
    n=len(A); out=Counter()
    def grow(done, path, remaining):
        if len(done)==k-1 and not remaining:
            paths=done+(path,)
            out[tuple(p[0] for p in paths),tuple(p[-1] for p in paths)]+=1
        for v in remaining:
            rest=remaining-{v}
            if A[path[-1]][v]: grow(done,path+(v,),rest)
            if len(done)<k-1: grow(done+(path,),(v,),rest)
    if 1<=k<=n:
        for v in range(n): grow((),(v,),frozenset(range(n))-{v})
    return out

def update(A,a,b):
    n=len(A); x=n; old={k:covers(A,k) for k in (1,2,3)}; out={}
    for s in permutations(range(n+1),2):
      for t in permutations(range(n+1),2):
        if x in s and x in t:
            j=s.index(x); h=t.index(x)
            v=old[1][((s[1-j],),(t[1-j],))] if j==h else 0
        elif x in s:
            j=s.index(x); v=0
            for u in range(n):
                ss=list(s); ss[j]=u
                v+=b[u]*old[2][tuple(ss),t]
        elif x in t:
            j=t.index(x); v=0
            for u in range(n):
                tt=list(t); tt[j]=u
                v+=a[u]*old[2][s,tuple(tt)]
        else:
            v=sum(a[u]*b[w]*(old[3][(s[0],w,s[1]),(u,t[0],t[1])]+old[3][(s[0],s[1],w),(t[0],u,t[1])]) for u in range(n) for w in range(n))
        if v: out[s,t]=v
    return out

def digraph_checks():
    tally={}
    # Exhaustive loopless directed graphs including opposite directed arcs.
    for n in (1,2,3):
        edges=list(permutations(range(n),2)); count=0
        for bits in range(1<<len(edges)):
            A=[[0]*n for _ in range(n)]
            for q,(u,v) in enumerate(edges): A[u][v]=(bits>>q)&1
            for z in range(1<<(2*n)):
                a=[(z>>i)&1 for i in range(n)]; b=[(z>>(n+i))&1 for i in range(n)]
                B=[row+[a[u]] for u,row in enumerate(A)]+[b+[0]]
                need(update(A,a,b)==covers(B,2),('insertion',n,bits,z))
                count+=1
        tally[n]=count
    return tally

def language_checks():
    safe=[]; cyclic=[]; necklaces=[]
    x,y,z=2,3,5; weighted=[]
    for m in range(10):
        words=list(product('abe',repeat=m))
        good=lambda w: all(w[i-1]+w[i] not in ('be','ee') for i in range(1,len(w)))
        s=[w for w in words if good(w)]; safe.append(len(s))
        weighted.append(sum(x**w.count('a')*y**w.count('b')*z**w.count('e') for w in s))
        # First-hit partition, including overlapping eee.
        hits=Counter(next((i+1 for i in range(1,m) if w[i-1]+w[i] in ('be','ee')),0) for w in words)
        need(hits[0]==len(s) and sum(hits.values())==3**m,('firsthit',m))
        if m:
            c=[w for w in words if all(w[i-1]+w[i] not in ('be','ee') for i in range(m))]
            cyclic.append(len(c))
            orbits={min(w[i:]+w[:i] for i in range(m)) for w in c}; necklaces.append(len(orbits))
            # Burnside independently via rotations' fixed words, not divisor formula.
            fixed=sum(sum(w==w[r:]+w[:r] for w in c) for r in range(m))
            need(fixed==m*len(orbits),('Burnside',m))
    need(safe==[1,3,7,17,41,99,239,577,1393,3363],'safe coefficients')
    need(weighted[0]==1 and weighted[1]==x+y+z,'weighted boundary')
    for m in range(2,10): need(weighted[m]==(x+y)*weighted[m-1]+x*z*weighted[m-2],('weighted GF',m))
    for m in range(2,9): need(cyclic[m]==2*cyclic[m-1]+cyclic[m-2],('cyclic recurrence',m))
    return dict(safe=safe,cyclic=cyclic,necklaces=necklaces,weighted_at_2_3_5=weighted)

if __name__=='__main__':
    result={'status':'FINITE-EXACT independent audit; no tournament word criterion assumed', 'loopless_digraph_extensions':digraph_checks(),'formal_language':language_checks()}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
