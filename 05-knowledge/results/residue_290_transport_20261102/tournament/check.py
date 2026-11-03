#!/usr/bin/env python3
"""Fresh exact scalar-kernel and native black-word checks; no repo imports."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
import hashlib
import json

checks=0
def require(ok,why):
    global checks
    checks+=1
    if not ok: raise RuntimeError(why)

def odd_parts(n,lo=1):
    if n==0:
        yield ()
    for x in range(lo,n+1,2):
        for rest in odd_parts(n-x,x): yield (x,)+rest

def exponent(t):
    return sum((x-1)//2 for x in t)+sum(gcd(x,y) for x,y in combinations(t,2))
def centralizer(t):
    ans=1
    for x,c in Counter(t).items(): ans*=x**c*factorial(c)
    return ans
def all_graph(t):
    return 2**sum(gcd(x,y) for x,y in combinations(t,2))

@lru_cache(None)
def connected(t):
    if len(t)==1: return 1
    if not t: return 0
    # Root is one distinguished occurrence of t[0]; group choices of all others.
    vals=sorted(Counter(t[1:]).items())
    ans=all_graph(t)
    for counts in product(*(range(c+1) for x,c in vals)):
        left=(t[0],); right=(); ways=1
        for (x,c),take in zip(vals,counts):
            left+=(x,)*take; right+=(x,)*(c-take); ways*=comb(c,take)
        if not right: continue
        ans-=ways*connected(tuple(sorted(left)))*all_graph(right)
    return ans

def set_partitions(n):
    if n==0:
        yield ()
        return
    for p in set_partitions(n-1):
        yield p+((n-1,),)
        for i in range(len(p)):
            yield p[:i]+(p[i]+(n-1,),)+p[i+1:]

def connected_mobius(t):
    return sum((-1)**(len(p)-1)*factorial(len(p)-1)*
               2**sum(gcd(t[i],t[j]) for b in p for i,j in combinations(b,2))
               for p in set_partitions(len(t)))

def connected_edges(t):
    n=len(t); es=list(combinations(range(n),2)); ans=0
    for mask in range(1<<len(es)):
        seen={0}; changed=True
        while changed:
            changed=False
            for j,(u,v) in enumerate(es):
                if mask>>j&1 and ((u in seen)!=(v in seen)):
                    seen|={u,v}; changed=True
        if len(seen)!=n: continue
        weight=1
        for j,(u,v) in enumerate(es):
            if mask>>j&1: weight*=2**gcd(t[u],t[v])-1
        ans+=weight
    return ans

def literal_fixed(p):
    n=len(p); visited=set(); orbits=0
    for i,j in combinations(range(n),2):
        if (i,j) in visited: continue
        a,b=i,j
        while True:
            visited.add(tuple(sorted((a,b))))
            a,b=p[a],p[b]
            if {a,b}=={i,j}:
                if (a,b)!=(i,j): return 0
                break
        orbits+=1
    return 2**orbits

def matrix_rank(a):
    a=[[F(x) for x in row] for row in a]; row=0
    for col in range(len(a[0]) if a else 0):
        pivot=next((i for i in range(row,len(a)) if a[i][col]),None)
        if pivot is None: continue
        a[row],a[pivot]=a[pivot],a[row]
        v=a[row][col]; a[row]=[x/v for x in a[row]]
        for i in range(len(a)):
            if i!=row:
                q=a[i][col]; a[i]=[x-q*y for x,y in zip(a[i],a[row])]
        row+=1
    return row

# Native tournament controls use immutable rows of 0/1 entries, independent
# of the producer's prior bit-row implementation.
def from_masks(rows):
    return tuple(tuple((r>>j)&1 for j in range(len(rows))) for r in rows)
def hp(a,end=None):
    n=len(a)
    for p in permutations(range(n)):
        if end is not None and p[-1]!=end: continue
        if all(a[u][v] for u,v in zip(p,p[1:])): yield p

def is_blue(a,p):
    n=len(a)
    return all(a[p[i]][p[j]]==a[p[n-1-j]][p[n-1-i]]
               for i,j in combinations(range(n),2))
def opposite(a):
    return tuple(tuple(0 if i==j else 1-a[i][j] for j in range(len(a))) for i in range(len(a)))
def path_reverse(a,p):
    b=[list(row) for row in a]
    for u,v in zip(p,p[1:]): b[u][v]^=1; b[v][u]^=1
    return tuple(map(tuple,b))
def nonpath_complement(a,p):
    keep={frozenset((u,v)) for u,v in zip(p,p[1:])}
    return tuple(tuple(0 if i==j else a[i][j] if frozenset((i,j)) in keep
                       else 1-a[i][j] for j in range(len(a))) for i in range(len(a)))
def strong(a,vertices=None):
    vs=set(range(len(a))) if vertices is None else set(vertices)
    if not vs: return False
    for s in vs:
        reached={s}
        while True:
            nxt=reached|{v for u in reached for v in vs if a[u][v]}
            if nxt==reached: break
            reached=nxt
        if reached!=vs: return False
    return True

def verify_route(initial,first):
    n=len(initial); m=(n-1)//2; v=first[-1]
    a=initial; native=initial; starts=Counter(); words=[]
    for j in range(m):
        if j==0: p=first
        elif j==1: p=next(p for p in hp(a,v) if p[0]!=first[0])
        else: p=next(hp(a,v))
        require(not is_blue(a,p),'route contains blue path')
        starts[p[0]]+=1
        nxt=path_reverse(a,p)
        q=p if j%2==0 else p[::-1]
        require(all(native[u][w] for u,w in zip(q,q[1:])),'illegal native word')
        require(not is_blue(native,q),'native gauge changes color')
        native=nonpath_complement(native,q)
        require(native==(opposite(nxt) if (j+1)%2 else nxt),'native gauge mismatch')
        for u in range(n):
            require(sum(nxt[u])==(m+j+1 if u==v else m-starts[u]),'lost named current')
        for mask in range(1,1<<n):
            if mask==(1<<n)-1: continue
            shore=[u for u in range(n) if mask>>u&1]
            d=sum(sum(nxt[u])-m for u in shore)
            require(abs(d)<=j+1,'cut current bound')
        words.append(q); a=nxt
    require(sum(a[v])==n-1,'not source in reversal gauge')
    require(strong(a,set(range(n))-{v}),'not strong terminal core')
    require(not strong(a),'terminal still strong')
    return words

def main():
    N=30; parts=[list(odd_parts(n)) for n in range(N+1)]
    beta=[F(0)]*(N+1); burnside=[]
    kernel_controls=0
    for n,pp in enumerate(parts):
        burnside.append(sum((F(2**exponent(t),centralizer(t)) for t in pp),F(0)))
        require(burnside[-1].denominator==1,'nonintegral Burnside count')
        if n:
            beta[n]=n*sum((F(2**((n-len(t))//2)*connected(t),centralizer(t)) for t in pp),F(0))
            require(beta[n].denominator==1 and beta[n]>0,'bad coefficient')
        for t in pp:
            if not t or n>15: continue
            if len(t)<=6:
                require(connected(t)==connected_mobius(t),'partition interaction mismatch'); kernel_controls+=1
            if len(t)<=5:
                require(connected(t)==connected_edges(t),'edge interaction mismatch'); kernel_controls+=1
    seq=[F(1)]
    for n in range(1,N+1):
        seq.append(sum(beta[k]*seq[n-k] for k in range(1,n+1))/n)
        require(seq[n]==burnside[n],'scalar recurrence differs from Burnside')
    require(seq[:11]==list(map(F,[1,1,1,2,4,12,56,456,6880,191536,9733056])),'known small prefix')
    literal=[]
    for n in range(8):
        count=F(sum(literal_fixed(p) for p in permutations(range(n))),factorial(n))
        require(count==seq[n],'literal signed pair action differs'); literal.append(int(count))
    ranks={}
    for m in range(3,31):
        weights={}
        for t in parts[m]:
            if 1 in t: continue
            weights[len(t)]=weights.get(len(t),F(0))+F(2**exponent(t),centralizer(t))
        dim=len(weights); expected=m//6 if m%2==0 else (m+3)//6
        require(dim==expected,'core support size')
        values=[sum(w*2**(t*r) for t,w in weights.items()) for r in range(2*dim+1)]
        require(matrix_rank([[values[i+j] for j in range(dim)] for i in range(dim)])==dim,'actual Hankel rank')
        poly=[1]
        for t in weights:
            nxt=[0]*(len(poly)+1)
            for i,c in enumerate(poly): nxt[i]-=2**t*c; nxt[i+1]+=c
            poly=nxt
        for r in range(dim+1):
            require(sum(c*values[r+i] for i,c in enumerate(poly))==0,'minimal response recurrence')
        ranks[m]=dim
    # First altered-weight loss and first erased-gcd error.
    require(2**1==2 and 2**3==8,'mass-nine response hostile')
    erased=[]
    for n,pp in enumerate(parts[:7]):
        erased.append(sum((F(2**(sum((x-1)//2 for x in t)+comb(len(t),2)),centralizer(t)) for t in pp),F(0)))
    require(erased[:6]==seq[:6] and erased[6]==F(164,3) and seq[6]-erased[6]==F(4,3),'first erased-gcd failure')
    seeds=[(112,97,67,7,14,28,56),(112,97,19,7,74,28,44),(112,73,35,21,70,26,44)]
    black=[]; traces=[]; dig=hashlib.sha256()
    for raw in seeds:
        a=from_masks(raw); require(all(sum(row)==3 for row in a),'seed not regular')
        count=0
        for p in hp(a):
            if is_blue(a,p): continue
            route=verify_route(a,p); count+=1
            dig.update(json.dumps([raw,route],separators=(',',':')).encode())
            if count==1: traces.append({'source_rows':raw,'native_words':route})
        black.append(count)
    require(sum(black)==382,'black first-path universe changed')
    pure=[]
    for n in (3,5):
        a=tuple(tuple(int(i!=j and 1<=(j-i)%n<=(n-1)//2) for j in range(n)) for i in range(n))
        paths=list(hp(a)); require(all(is_blue(a,p) for p in paths),'pure-blue boundary changed');pure.append(len(paths))
    a=from_masks((112,37,73,67,14,28,50)); colors=[]
    for p in ((6,1,5,4,2,3,0),(6,4,5,3,1,2,0),(6,5,4,3,2,1,0)):
        require(all(a[u][v] for u,v in zip(p,p[1:])),'bad hostile path')
        colors.append('blue' if is_blue(a,p) else 'black'); a=path_reverse(a,p)
    require(colors==['black','blue','blue'] and sum(a[0])==6 and sum(a[6])==0,'same-start hostile changed')
    print(json.dumps({'status':'PASS','checks':checks,'scalar_through':N,'a_values':[int(v) for v in seq],
        'g_values':[int(v) for v in beta[1:]],'kernel_controls':kernel_controls,
        'literal_permutation_counts_through7':literal,'unit_tail_ranks':ranks,
        'regular7_black_first_paths':black,'regular7_route_digest':dig.hexdigest(),'sample_native_traces':traces,
        'pure_blue_H_counts_C3_R5':pure,'same_start_colors':colors,
        'scope':'fresh finite checks, not reuse of unavailable October26 receipts'},indent=2,sort_keys=True))

if __name__=='__main__': main()
