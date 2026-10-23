#!/usr/bin/env python3
"""Exact forward-step/cyclotomic scanner; standard library only."""
import itertools,math,json,random,collections,sys
from pathlib import Path

def require(test,msg):
    if not test: raise RuntimeError(msg)
def digit_sum(m,p):
    s=0
    while m:s,m=s+m%p,m//p
    return s
def bound(n,p):return n-max(digit_sum(m,p) for m in range(1,n+1))
def arcs(n,code):
    arc=[[0]*n for _ in range(n)]
    for k,(a,b) in enumerate(itertools.combinations(range(n),2)):
        arc[a][b]=(code>>k)&1;arc[b][a]=1-arc[a][b]
    return arc
def polynomial(n,code):
    pairs=list(itertools.combinations(range(n),2)); arc=[[0]*n for _ in range(n)]
    for k,(a,b) in enumerate(pairs):arc[a][b]=(code>>k)&1;arc[b][a]=1-arc[a][b]
    dp={1<<i:{i:[1]} for i in range(n)}
    for mask in range(1,1<<n):
        if mask not in dp:continue
        for i,coeff in dp[mask].items():
            for j in range(n):
                if (mask>>j)&1:continue
                target=dp.setdefault(mask|(1<<j),{}).setdefault(j,[0]*(mask.bit_count()+1))
                step=arc[i][j]
                for k,c in enumerate(coeff):target[k+step]+=c
    f=[sum(dp[(1<<n)-1][i][k] for i in range(n)) for k in range(n)]
    require(sum(f)==math.factorial(n),'factorial');require(f==f[::-1],'palindrome')
    for j in range(n):
        cj=sum(math.comb(k,j)*f[k] for k in range(j,n))
        require(cj%math.factorial(n-j)==0,'path-cover integrality')
    return f

def reduce_cyclotomic(f,p):
    sums=[sum(f[j] for j in range(r,len(f),p)) for r in range(p)]
    return [sums[j]-sums[-1] for j in range(p-1)]
def valuation(a,p):
    if not any(a):return None
    v=0
    while sum(a)%p==0:
        t=-sum(a)//p;run=0;b=[]
        for j,x in enumerate(a):run+=x;b.append(-run-(j+1)*t)
        a=b;v+=1
    return v

def binary_bound(n):
    if n==1:return 0
    require(n%2==1,'binary bound odd only')
    return min(n-((n+1)&-(n+1)).bit_length()+1,n+1-((n-1).bit_length()-1))
def eulerian(n):
    a=[1]
    for j in range(2,n+1):
        a=[(k+1)*(a[k] if k<len(a) else 0)+(j-k)*(a[k-1] if k else 0) for k in range(j)]
    return a
def odd_cycle_expansion(n,code):
    arc=arcs(n,code); cycles=collections.Counter()
    for length in range(3,n+1,2):
        for verts in itertools.combinations(range(n),length):
            root=verts[0];mask=sum(1<<v for v in verts)
            for tail in itertools.permutations(verts[1:]):
                c=(root,)+tail
                if all(arc[c[j]][c[(j+1)%length]] for j in range(length)):cycles[mask]+=1
    from functools import lru_cache
    @lru_cache(None)
    def collections_for(mask):
        if not mask:return {(0,0):1}
        bit=mask&-mask;out=collections.Counter(collections_for(mask^bit))
        for c,count in cycles.items():
            if c&bit and c&mask==c:
                for (vertices,r),v in collections_for(mask^c).items():out[(vertices+c.bit_count(),r+1)]+=count*v
        return dict(out)
    out=[0]*n
    for (vertices,r),count in collections_for((1<<n)-1).items():
        d=n-vertices+r;gap=n-d
        for k,coeff in enumerate(eulerian(d)):
            for j in range(gap+1):out[k+j]+=count*2**r*coeff*math.comb(gap,j)*(-1)**(gap-j)
    return out
def run():
    out={'description':'All labelled tournaments n<=5; deterministic samples n=6..10. Not an isomorphism census.','orders':[]}
    rng=random.Random(4304756)
    for n in range(1,11):
        total=1<<(n*(n-1)//2)
        codes=range(total) if n<=5 else sorted({0,total-1}|{rng.randrange(total) for _ in range(100)})
        rows={p:{'bound':bound(n,p),'min_valuation':None,'zero_count':0,'min_witness':None,'valuation_histogram':collections.Counter()} for p in (2,3,5,7)}
        count=0
        for code in codes:
            f=polynomial(n,code);count+=1
            if n<=6:require(f==odd_cycle_expansion(n,code),f'odd-cycle expansion n{n} code{code}')
            for p,r in rows.items():
                a=reduce_cyclotomic(f,p);v=valuation(a,p)
                if v is None:r['zero_count']+=1;continue
                require(v>=r['bound'],f'valuation bound n{n} p{p} code{code}')
                if p%2:require(v%2==0,'real-subfield parity');require(v>=2*((r['bound']+1)//2),'rounded odd-prime bound')
                if p==2 and n%2:require(v>=binary_bound(n),'binary odd-cycle bound')
                if p==2 and n%2==0:raise RuntimeError('even order nonzero at minus one')
                r['valuation_histogram'][v]+=1
                if r['min_valuation'] is None or v<r['min_valuation']:r['min_valuation']=v;r['min_witness']={'code':code,'polynomial':f,'cyclotomic_coefficients':a}
        for r in rows.values():r['valuation_histogram']=dict(sorted(r['valuation_histogram'].items()))
        out['orders'].append({'n':n,'tested':count,'exhaustive':n<=5,'binary_cycle_bound':binary_bound(n) if n%2 else None,'primes':rows})
        print('completed',n,count,flush=True,file=sys.stderr)
    # Verify division formula by re-multiplication on arbitrary ring elements.
    for p in (2,3,5,7,11):
        for _ in range(100):
            b=[rng.randrange(-1000,1001) for j in range(p-1)]
            a=[-b[0]-b[-1]]+[b[j-1]-b[j]-b[-1] for j in range(1,p-1)]
            t=-sum(a)//p;run=0;got=[]
            for j,x in enumerate(a):run+=x;got.append(-run-(j+1)*t)
            require(got==b,'division formula')
    # Explicit all-order extremizer family; scalar evaluation uses proved (6).
    sharp=[]
    aval={n:sum((-1)**k*c for k,c in enumerate(eulerian(n))) for n in range(1,202,2)}
    for n in range(3,202,2):
        a=(n-1).bit_length()-1
        empty=n-((n+1)&-(n+1)).bit_length()+1
        if empty<=n+1-a:q=0;value=aval[n]
        else:
            q=n-(1<<a)+2
            value=aval[n]+sum(2*math.comb(q-2,l-2)*2**(l-1)*aval[n-l+1] for l in range(3,q+1,2))
        require(valuation([value],2)==binary_bound(n),'all-order explicit extremizer scalar')
        m=n+1
        power_two=m&(m-1)==0
        power_three_two=m%3==0 and (m//3)&((m//3)-1)==0
        require((valuation([aval[n]],2)==binary_bound(n))==(power_two or power_three_two),'transitive extremizer classification')
        if n<=9:
            pairs=list(itertools.combinations(range(n),2));code=(1<<len(pairs))-1
            if q:code^=1<<pairs.index((0,q-1))
            f=polynomial(n,code)
            require(sum((-1)**k*c for k,c in enumerate(f))==value,'actual extremizer tournament')
        sharp.append({'n':n,'q':q,'binary_valuation':binary_bound(n)})
    out['explicit_binary_extremizers_through_order_201']=sharp
    print(json.dumps(out,indent=2))
if __name__=='__main__':run()
