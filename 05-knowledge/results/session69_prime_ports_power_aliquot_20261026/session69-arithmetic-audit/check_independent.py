#!/usr/bin/env python3
"""Independent Session69 arithmetic audit; imports no producer modules.

Source sigma is computed multiplicatively from a smallest-prime-factor sieve;
power origins are independently recovered by exact integer-root tests.
"""
import argparse
from collections import Counter
from math import gcd, isqrt, lcm
import json
import hashlib

checks = 0
def need(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)

def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:
            out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=1
    return out

def radical(n):
    r=1
    for p in factors(n):r*=p
    return r

def exact_root(n,r):
    lo,hi=1,n
    while lo<=hi:
        z=(lo+hi)//2;v=z**r
        if v==n:return z
        if v<n:lo=z+1
        else:hi=z-1
    return None

def coordinates(n,square=False):
    a=0
    while n%2==0:n//=2;a+=1
    if n==1:return a,1,0
    if square:
        t=0
        while isqrt(n)**2==n:n=isqrt(n);t+=1
        return a,n,t
    best=n;e=1
    for r in range(2,n.bit_length()):
        b=exact_root(n,r)
        if b is not None and b<best:best=b;e=r
    return a,best,e

def support(n,cap,square=False):
    out=set();r=1
    while n**r<=cap:
        value=n**r
        while value<=cap:out.add(value);value*=2
        if n==1:break
        r=2*r if square else r+1
    return out

def sigma(n):
    out=1
    for p,e in factors(n).items():out*=(p**(e+1)-1)//(p-1)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--negative-control',action='store_true');args=ap.parse_args()
    N=20000
    origins={};counts={}
    for square in (False,True):
        first=[None]*(N+1);ct=Counter()
        for base in range(1,N+1):
            r=1
            while base**r<=N:
                v=base**r
                while v<=N:
                    if first[v] is None:first[v]=base
                    if base>1:ct[v]+=1
                    v*=2
                if base==1:break
                r=2*r if square else r+1
        for x in range(1,N+1):
            a,b,e=coordinates(x,square)
            need(first[x]==b,('origin',square,x))
            if b==1:
                predicted=sum(a//(2**j) for j in range(a.bit_length())) if square else sum(a//r for r in range(1,a+1))
            elif square:predicted=sum(a//(2**j)+1 for j in range(e+1))
            else:predicted=sum(a//r+1 for r in range(1,e+1) if e%r==0)
            need(ct[x]==predicted,('addresses',square,x))
        origins[square]=first;counts[square]=ct
    need(next(x for x in range(1,N+1) if origins[False][x]!=origins[True][x])==27,'first origins differ')
    need(next(x for x in range(1,N+1) if counts[False][x]!=counts[True][x])==8,'first counts differ')
    cap=2_000_000;bases=180;paircases=0
    for square in (False,True):
        rows={n:support(n,cap,square) for n in range(1,bases+1)}
        data={n:coordinates(n,square) for n in rows}
        for n in rows:
            for m in range(n,bases+1):
                paircases+=1
                a,b,e=data[n];c,d,f=data[m]
                if b!=d:L=None
                elif b==1:L=2**max(a,c)
                elif square:
                    t=max(e,f);L=2**max(a*2**(t-e),c*2**(t-f))*b**(2**t)
                else:
                    g=gcd(e,f);L=2**max(a*f//g,c*e//g)*b**lcm(e,f)
                need(rows[n]&rows[m]==(set() if L is None else support(L,cap,square)),('intersection',square,n,m))
    X=1_000_000;later=[];density=[]
    for n in range(1,X+1):
        z=n*n
        while z%2==0:z//=2
        u=n
        while u%2==0:u//=2
        criterion=n%2==0 and u>n//u
        need((n%2==0 and z>n)==criterion,('square-later',n))
        if criterion:later.append(n)
        if n in (100,1000,10000,100000,1000000):density.append([n,len(later),n//2-len(later)])
    need(later[0]==(4 if args.negative_control else 6),'negative-control canonical first-six')
    need(4**2==2*8 and 8>4,'arbitrary-root first-four')
    for r in range(2,10):
        hit=[]
        for n in range(1,101):
            z=n**r;k=0
            while z%2==0:z//=2;k+=1
            if k and z>n:hit.append(n)
        need(hit[0]==6,('higher-power first-six',r))
    triples={};hits=[];B=2
    while B*(B+1)<=100000:
        for C in range(B+1,100000//B+1,2):
            A=C-B
            need(gcd(A,B)==gcd(A,C)==gcd(B,C)==1,'binary primitive')
            triples[B*C]=(A,B,C)
            if C>radical(A*B*C):hits.append((B*C,A,B,C,radical(A*B*C)))
        B*=2
    need(set(triples)=={n for n in later if n<=100000},'binary triple bijection')
    need(min(hits)==(72,1,8,9,6),'first binary quality')
    M=1000;bound=(M-1)**2
    spf=list(range(bound+1))
    for p in range(2,isqrt(bound)+1):
        if spf[p]==p:
            for n in range(p*p,bound+1,p):
                if spf[n]==n:spf[n]=p
    sig=[0]*(bound+1);sig[1]=1
    fibers={m:[] for m in range(M+1)}
    inverse_bound_checks=0
    for n in range(2,bound+1):
        p=spf[n];r=n;power=1;geometric=1
        while r%p==0:r//=p;power*=p;geometric+=power
        sig[n]=sig[r]*geometric;m=sig[n]-n
        if m<=M:fibers[m].append(n)
        if m>1:
            need(n<=(m-1)**2,('inverse global bound',n));inverse_bound_checks+=1
            if not (isqrt(n)**2==n and spf[isqrt(n)]==isqrt(n)):
                need(4*n<=(m-1)**2,('inverse non-square bound',n))
    expected={2:[],3:[4],4:[9],5:[],6:[6,25],7:[8],8:[10,49],9:[15],10:[14],11:[21],12:[121],13:[27,35],14:[22,169],21:[18,51,91]}
    for m,F in expected.items():need(fibers[m]==F,('complete fiber',m))
    need(next(m for m in range(2,M+1) if len(fibers[m])>1)==6,'first finite collision')
    need([n for n in range(2,10000) if sig[n]==2*n]==[6,28,496,8128],'perfect normalization controls')
    for n in range(2,20001):
        m=sig[n]-n;g=gcd(n,m);A=n//g;B=m//g;C=(n+m)//g
        need(A+B==C and gcd(A,B)==gcd(A,C)==gcd(B,C)==1,('source ABC',n))
        need((A*g,B*g)==(n,m),('gcd source reconstruction',n))
        need(((A,B,C)==(1,1,2))==(m==n),('perfect ABC iff',n))
    primes=[p for p in range(2,500) if spf[p]==p]
    branches=overlap=cancel=twins=0
    for a in range(2,301):
        S=sig[a]-a;T=sig[a];d=gcd(S,T);U=S//d;V=T//d
        for p in primes:
            if a%p==0:continue
            branches+=1;h=gcd(p,V);A=U*p//h;B=V//h;C=(U*p+V)//h
            need(sigma(a*p)-a*p==S*p+T,'guarded branch')
            need(gcd(S*p,T)==d*h and gcd(A,B)==1 and A+B==C,'branch normalization')
            R=radical(A*B*C)
            need(C<=(U+V)*R,'fixed-core bound')
            if h==1:
                need(R==p*radical(U*V*C)//gcd(p,U),'support overlap formula')
                overlap+=int(U%p==0)
            else:cancel+=1
    for p in primes:
        if p>2 and spf[p+2]==p+2:
            twins+=1;need(sigma(2*p)-2*p==p+3==sigma((p+2)**2)-(p+2)**2,'twins')
            need((2*p==p+3)==(p==3),'twin self')
    need(sigma(8)-8==7 and sigma(27)-27==13 and sigma(64)-64==63,'Phi3 prime/composite')
    need(sigma(4)-4==3!=(sigma(2)-2)*2+sigma(2),'prime guard hostile')
    need((1+4+16)==(1+2+4)*(1-2+4)==21,'cyclotomic identity')
    # Independently repeat Booker's finite equal-prime supplement only.
    supplement=[]
    for p in range(5,9782):
        if spf[p]!=p:continue
        witnesses=[q for q in range(3,p,2) if spf[q]==q and spf[2*p-q]==2*p-q]
        need(bool(witnesses),('Booker distinct supplement',p))
        supplement.append([2*p,witnesses[0],2*p-witnesses[0]])
    def paths(adjacency):
        n=len(adjacency);dp=[[0]*n for _ in range(1<<n)]
        for v in range(n):dp[1<<v][v]=1
        for mask in range(1<<n):
            for last in range(n):
                if not dp[mask][last]:continue
                for nxt in adjacency[last]:
                    if not mask>>nxt&1:dp[mask|1<<nxt][nxt]+=dp[mask][last]
        return sum(dp[-1])
    paley=[{j for j in range(7) if (j-i)%7 in (1,2,4)} for i in range(7)]
    need(paths(paley)==189,'independent DP Paley H')
    import itertools
    edges=list(itertools.combinations(range(4),2));h5=0
    for mask in range(64):
        adj=[set() for _ in range(4)]
        for i,(u,v) in enumerate(edges):
            if mask>>i&1:adj[u].add(v)
            else:adj[v].add(u)
        h5+=int(paths(adj)==5)
    need(h5==24,'independent DP four-vertex H5')
    fiber_digest=hashlib.sha256(json.dumps({str(m):fibers[m] for m in range(2,M+1)},sort_keys=True,separators=(',',':')).encode()).hexdigest()
    print(json.dumps({'status':'PASS','checks':checks,'origin_and_count_cap':N,
      'row_bases':bases,'row_cap':cap,'row_pairs':paircases,'square_later_cap':X,
      'density_samples_X_qualifying_even_failures':density,'first_binary_quality':min(hits),
      'complete_aliquot_target_cap':M,'complete_aliquot_fibers_sha256':fiber_digest,'source_bound':bound,'inverse_bound_checks':inverse_bound_checks,
      'small_complete_fibers':expected,'source_ABC_cases':19999,
      'guarded_branch_core_cap':300,'guarded_branch_prime_cap':499,'guarded_branches':branches,
      'overlap_cases':overlap,'cancellation_cases':cancel,'twin_pairs':twins,
      'Booker_equal_prime_supplement_cases':len(supplement),'Booker_large_run_repeated':False,
      'independent_DP_Paley_H':189,'independent_DP_order4_tournaments_H5':h5},indent=2,sort_keys=True))

if __name__=='__main__':main()
