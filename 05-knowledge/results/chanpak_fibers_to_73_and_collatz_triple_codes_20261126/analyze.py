from collections import Counter, defaultdict
from functools import cache
from itertools import combinations
from math import prod
from pathlib import Path
from fractions import Fraction
import sys
here=Path(__file__).resolve().parent
rows=[tuple(map(int,line.split())) for line in Path(sys.argv[1] if len(sys.argv)>1 else here/'census.out').read_text().splitlines()]
G=defaultdict(Counter); Gc=defaultdict(Counter); tags=defaultdict(list)
for n,h,m,a,s,*p in rows:
    G[h][n]+=1; Gc[h][n]+=s; tags[n,h].append((m,a,s,tuple(p)))
@cache
def W(h):
    if h==1:return {(0,0):1}
    c=Counter()
    for d, v in G.items():
        if h%d:continue
        for j,x in v.items():
            for (k,n),y in W(h//d).items():c[k+1,n+j]+=x*y
    return dict(c)
def term(k,n,c):
    return (str(c)+'*' if c!=1 else '')+f't^{n}/R^{k+1}'
for h in range(27,74,2):
    terms=sorted(W(h).items())
    print(h,' + '.join(term(k,n,c) for (k,n),c in terms) or '0')
    assert terms
print('census by size',sorted(Counter(r[0] for r in rows).items()))
def d(p):
    a,b,c,f,g=p
    return a,b+2*f,c+2*g
for mode in ['same-D','same-profile']:
    ans=[]
    for (n,h),ar in tags.items():
        if h<=25:continue
        for p,q in combinations(ar,2):
            if mode=='same-D' and d(p[-1])==d(q[-1]) and p[-1]!=q[-1] or mode=='same-profile' and p[-1]==q[-1] and p[1:3]!=q[1:3]:
                ans.append((h,n,p,q))
    print(mode,sorted(ans)[:7])
def triples(a):
    n=len(a)
    return {i for i,(u,v,w) in enumerate(combinations(range(n),3))
             if a[u][v]==a[v][w]==a[w][u]}
def analyze_patterns(n):
    E=list(combinations(range(n),2));T=list(combinations(range(n),3));K=len(E)
    pattern=Counter()
    for word in range(1<<K):
        a=[[0]*n for _ in range(n)]
        for i,(u,v) in enumerate(E):
            a[u][v]=word>>i&1;a[v][u]=1-a[u][v]
        pattern[tuple(sorted(triples(a)))]+=1
    def eq_for(B):
        D=[(i,0) for i in range(K)]
        def find(i):
            a,b=D[i]
            if a==i:return a,b
            r,x=find(a);D[i]=r,x^b;return D[i]
        for j in B:
            u,v,w=T[j];b=E.index((u,v))
            for a,z in [(E.index((v,w)),0),(E.index((u,w)),1)]:
                p,x=find(a);q,y=find(b)
                if p==q:
                    if x^y!=z:return 0
                else:D[p]=q,x^y^z
        return 2**sum(i==x for i,(x,z) in enumerate(D))
    scope=[{i for i,t in enumerate(T) if set(t)<=set(q)} for q in combinations(range(n),4)]
    local=0;witness={}; hist=Counter()
    for mask in range(1<<len(T)):
        R={j for j in range(len(T)) if mask>>j&1}; B=set(range(len(T)))-R
        ans=0
        for z in range(1<<len(B)):
            J={j for i,j in enumerate(sorted(B)) if z>>i&1}
            ans+=(-1)**len(J)*eq_for(R|J)
        assert ans==pattern[tuple(sorted(R))]
        if all(len(q&R)<=2 for q in scope):
            local+=1
            if not ans:
                forced=tuple(T[j] for j in sorted(B) if eq_for(R|{j})==eq_for(R))
                kind='required inconsistent' if not eq_for(R) else ('forced forbidden' if forced else 'jointly forbidden')
                hist[kind]+=1
                w=(len(R),tuple(T[j] for j in sorted(R)),forced)
                witness[kind]=min(w,witness.get(kind,w))
    print('triples',n,len(pattern),local,dict(hist),witness)
for n in range(3,6):analyze_patterns(n)

def accelerated(n):
    return (3*n+1)//2 if n%2 else n//2
def word(x,k):
    a=[]
    for i in range(k):
        a.append(x%2);x=accelerated(x)
    return tuple(a)
assert [x for x in range(8) if len(set(word(x,3)))==1]==[0,7]
for k in range(1,11):
    assert len({word(x,k) for x in range(1<<k)})==1<<k
print('Collatz parity cylinder tests passed')
def odd_step(x):
    y=3*x+1
    return y//(y&-y)
for root in [7,9,21,27,37,61,63,65]:
    a=[root]
    while a[-1]!=1 and len(a)<7:a.append(odd_step(a[-1]))
    print('3x+1',a)
