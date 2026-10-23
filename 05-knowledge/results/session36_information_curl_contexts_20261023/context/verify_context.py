#!/usr/bin/env python3
"""Independent exact checks of module-context tomography. No repository imports."""
import itertools, math, json, random
from functools import lru_cache
from fractions import Fraction
from pathlib import Path

OUT=Path(__file__).resolve().parent

def check(ok, label):
    if not ok: raise RuntimeError(label)

def tournament(n, mask):
    rows=[0]*n
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a,b=(j,i) if mask>>bit&1 else (i,j)
        rows[a]|=1<<b
    return rows

def transitive(n): return tournament(n,0)
def reversal(n,j):
    a=transitive(n)
    if j: a[0]^=1<<(2*j); a[2*j]^=1
    return a

def h_subsets(a):
    n=len(a); size=1<<n
    end=[[0]*n for _ in range(size)]
    h=[0]*size;h[0]=1
    for v in range(n): end[1<<v][v]=1
    for s in range(1,size):
        for v in range(n):
            val=end[s][v]
            if not val: continue
            avail=a[v] & ~s
            while avail:
                bit=avail&-avail; avail-=bit
                end[s|bit][bit.bit_length()-1]+=val
        h[s]=sum(end[s])
    return h

def profile(a):
    h=h_subsets(a); n=len(a)
    @lru_cache(None)
    def p(s):
        if not s: return (1,)
        root=s&-s; other=s^root
        out=[0]*(s.bit_count()+1)
        z=other
        while True:
            block=z|root
            val=h[block]
            for d,v in enumerate(p(s^block)): out[d+1]+=val*v
            if not z:break
            z=(z-1)&other
        return tuple(out)
    return p((1<<n)-1)

def permutation_profile(a):
    n=len(a); f=[0]*(n+1)
    for perm in itertools.permutations(range(n)):
        b=sum(not(a[u]>>v&1) for u,v in zip(perm,perm[1:]))
        for cuts in range(n-b): f[1+b+cuts]+=math.comb(n-1-b,cuts)
    return tuple(f[d]//math.factorial(d) for d in range(n+1))

def mulpoly(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def q_from_profile(p):
    n=len(p)-1; q=[0]*(n+1); ff=[1]
    for d in range(1,n+1):
        ff=mulpoly(ff,[-(d-1),1])
        for i,c in enumerate(ff):q[i]+=p[d]*c
    return q

def q_from_cycles(a):
    n=len(a); cycles=[]
    for size in range(3,n+1,2):
        for vertices in itertools.combinations(range(n),size):
            v=vertices[0]; mask=sum(1<<i for i in vertices)
            for rest in itertools.permutations(vertices[1:]):
                word=(v,)+rest+(v,)
                if all(a[u]>>w&1 for u,w in zip(word,word[1:])): cycles.append(mask)
    @lru_cache(None)
    def rec(s):
        if not s:return (1,)
        v=s&-s; out=[0]*(s.bit_count()+1)
        for k,c in enumerate(rec(s^v)):out[k+1]+=c
        for c in cycles:
            if c&v and c&s==c:
                for k,w in enumerate(rec(s^c)):out[k+1]+=2*w
        return tuple(out)
    return list(rec((1<<n)-1))

@lru_cache(None)
def stirling(n,d):
    if n==0:return int(d==0)
    if d<1 or d>n:return 0
    return d*stirling(n-1,d)+stirling(n-1,d-1)

def weights(m,n):
    if m==0:return [0,1]+[0]*(n-1)
    g=[math.factorial(k)*stirling(m,k) for k in range(m+2)]
    def G(k):return g[k] if 0<=k<len(g) else 0
    return [0]+[3*G(k)**2+2*G(k)*G(k+1)+G(k+1)**2+G(k-1)**2+2*G(k-1)*G(k) for k in range(1,n+1)]

def response(p,m):
    return sum(w*math.factorial(d)*p[d] for d,w in enumerate(weights(m,len(p)-1)))

def cyclic_context(a,m):
    if m==0:return a[:]
    n=len(a); result=[0]*(n+2*m)
    for i,row in enumerate(a):result[i]=row|(((1<<m)-1)<<n)
    for i in range(m):
        result[n+i]=sum(1<<(n+j) for j in range(i+1,m))|(((1<<m)-1)<<(n+m))
        result[n+m+i]=sum(1<<(n+m+j) for j in range(i+1,m))|((1<<n)-1)
    return result

def solve(a,b):
    k=len(a)
    if not k:return []
    a=[[Fraction(v) for v in row]+[Fraction(b[i])] for i,row in enumerate(a)]
    for col in range(k):
        pivot=next(i for i in range(col,k) if a[i][col])
        a[col],a[pivot]=a[pivot],a[col]
        z=a[col][col];a[col]=[x/z for x in a[col]]
        for i in range(k):
            if i!=col:
                z=a[i][col];a[i]=[x-z*y for x,y in zip(a[i],a[col])]
    return [row[-1] for row in a]

def reconstruct(n,observations):
    r=(n+1)//2; f=[0]*r
    for m,y in enumerate(observations):
        w=weights(m,n)
        f[m+1]=Fraction(y-sum(w[d]*f[d] for d in range(1,m+1)),w[m+1])
    exponents=list(range(n-2,0,-2))
    # F_d = Delta^d Q(0) = d! sum q_k S(k,d); q_n=1.
    a=[[math.factorial(d)*stirling(k,d) for k in exponents] for d in range(1,r)]
    b=[f[d]-math.factorial(d)*stirling(n,d) for d in range(1,r)]
    vals=solve(a,b)
    q=[Fraction(0)]*(n+1);q[n]=1
    for k,v in zip(exponents,vals):q[k]=v
    return q

def det(a):
    a=[list(map(Fraction,row)) for row in a];result=Fraction(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];result=-result
        z=a[k][k];result*=z
        for i in range(k+1,len(a)):
            z2=a[i][k]/z
            for j in range(k+1,len(a)):a[i][j]-=z2*a[k][j]
    check(result.denominator==1,'integral determinant');return result.numerator

def main():
    stats={'exhaustive':{},'larger_controls':{},'rank_witnesses':[]}
    for n in range(1,6):
        count=0; direct=0
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);p=profile(a);q=q_from_profile(p)
            check(p==permutation_profile(a),('profile',n,mask))
            check(q==q_from_cycles(a),('cycle polynomial',n,mask))
            obs=[response(p,m) for m in range(max(0,(n+1)//2-1))]
            check(reconstruct(n,obs)==q,('reconstruction',n,mask))
            for m,y in enumerate(obs):
                check(h_subsets(cyclic_context(a,m))[-1]==y,('direct context',n,mask,m));direct+=1
            count+=1
        stats['exhaustive'][n]={'tournaments':count,'direct_contexts':direct}
        print('exhaustive',n,stats['exhaustive'][n],flush=True)
    rng=random.Random(36003121)
    for n,count in [(6,128),(7,64),(8,32),(9,16)]:
        for _ in range(count):
            mask=rng.getrandbits(n*(n-1)//2);a=tournament(n,mask);p=profile(a);q=q_from_profile(p)
            obs=[response(p,m) for m in range((n+1)//2-1)]
            check(reconstruct(n,obs)==q,('larger reconstruction',n,mask))
            # Last probe is the largest; independent vertex-level DP.
            m=len(obs)-1
            check(h_subsets(cyclic_context(a,m))[-1]==obs[m],('larger direct context',n,mask))
        stats['larger_controls'][n]=count
        print('larger controls',n,count,flush=True)
    for n in range(1,13):
        r=(n+1)//2; ps=[]
        for j in range(r):
            p=profile(reversal(n,j));q=q_from_profile(p)
            expected=[0]*(n+1);expected[n]=1
            for ell in range(1,j+1):expected[n-2*ell]=2*math.comb(2*j-1,2*ell-1)
            check(q==expected,('basis cycle formula',n,j));ps.append(p)
        matrix=[[response(p,m) for m in range(r)] for p in ps]
        determinant=det(matrix)
        check(determinant!=0,('rank',n))
        stats['rank_witnesses'].append({'n':n,'rank':r,'determinant':str(determinant)})
        print('rank witness',n,r,determinant,flush=True)
    # Strong positive/hostile control inherited from THM-3134, rebuilt independently.
    pp=[]
    for mask in (40,76):
        p=profile(tournament(5,mask));pp.append(p)
    check(pp[0][1]==pp[1][1]==15,'same H hostile')
    check([response(p,1) for p in pp]==[123,135],'separating context hostile')
    stats['same_H_hostile']={'masks':[40,76],'H':15,'C1':[123,135], 'profiles':pp}
    (OUT/'verification.json').write_text(json.dumps(stats,indent=2)+'\n')
    print('ALL CHECKS PASSED',flush=True)

if __name__=='__main__':main()
