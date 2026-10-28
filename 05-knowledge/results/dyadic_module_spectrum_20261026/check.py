#!/usr/bin/env python3
"""Exact checks: unchanged ROOT sources, chart ensembles, lattices and grammar."""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import product, permutations
from pathlib import Path
from math import isqrt
import hashlib
import json

HERE=Path(__file__).resolve().parent

def check(p,msg):
    if not p: raise RuntimeError(msg)

def Q(v):
    a,b=v
    return a*a-a*b+b*b

def add(a,b):return (a[0]+b[0],a[1]+b[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]-a[1]*b[1])

def lattice_checks():
    V=((1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),(-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1))
    bad=None;count=0;strict=0
    for vs in product(V,repeat=4):
        S=sum(map(Q,vs));E=Q(tuple(map(sum,zip(*vs))))
        mids=[Q(add(vs[0],vs[j]))+Q(add(vs[k],vs[l])) for j,k,l in ((1,2,3),(2,1,3),(3,1,2))]
        check(sum(mids)==E+2*S,'three pairing interpolation')
        mean=Fraction(sum(mids),3)
        check(max(E-mean,0)+max(mean-S,0)==max(E-S,0),'ensemble has no intermediate debt')
        if min(max(E-m,0)+max(m-S,0) for m in mids)>max(E-S,0):
            strict+=1
            if bad is None:bad=dict(vectors=vs,source=S,parent=E,middle=mids)
        count+=1
    # Integral Haar image has the two coordinate parity guards.
    haar=0
    for a,b,c,d in product(range(-2,3),repeat=4):
        u=(a+c,b+d);v=(a-c,b-d)
        check(all((x-y)%2==0 for x,y in zip(u,v)),'Haar lattice image')
        check(Q(u)+Q(v)==2*(Q((a,b))+Q((c,d))),'Eisenstein Haar norm')
        haar+=1
    alpha=(3,1);cube=mul(mul(alpha,alpha),alpha)
    check(cube==(19,18),'norm-seven cube')
    orders=[]
    for r in range(1,13):
        modulus=2**r;z=(1,0);order=0
        while True:
            z=tuple(x%modulus for x in mul(z,alpha));order+=1
            if z==(1,0):break
            check(order<=3*2**r,'finite twist order search')
        check(order==3*2**(r-1),'exact dyadic twist order')
        orders.append([r,order])
    return dict(four_vector_cases=count,strict_single_chart_hostiles=strict,hostile=bad,
                integral_haar_controls=haar,alpha_cubed=cube,twist_orders=orders)

def edge(s):
    a=2 if s%8==1 else 1
    y=(3*s+1)>>a;t=a
    while y%8==5:y=(y-1)//4;t+=2
    return y,t

def charge(s,N):
    lower=((1,2),(-2,-1),(1,-1));upper=((0,1),(-1,-1),(1,0))
    a=b=0;m=s
    while m<4*N:
        x,y=(lower if m<2*N else upper)[m%3];a+=x;b+=y;m=4*m+1
    check(Q((a,b)) in (0,1,3),'unchanged source alphabet')
    return a,b

def budgets(E):
    K=len(E)-1
    def fixed(w):
        levels=list(range(0,K,w))+[K]
        return sum(max(E[a]-E[b],0) for a,b in zip(levels,levels[1:]))
    dp=[0]+[10**100]*K;word=[[] for _ in range(K+1)]
    for j in range(1,K+1):
        for w in range(1,min(3,j)+1):
            v=dp[j-w]+max(E[j-w]-E[j],0)
            if v<dp[j]:dp[j]=v;word[j]=word[j-w]+[w]
    return dict(binary=fixed(1),fourary=fixed(2),eightary=fixed(3),adaptive_1to3=dp[K],adaptive_word=word[K])

def energies(rows):
    # Finite sources are (s,u,t,depth,port,a,b); zeros remain in the records.
    maxindex=max((s-1)//2 for s,*_ in rows)
    mod=1;out=[]
    while True:
        bins=defaultdict(lambda:[0,0])
        for s,u,t,d,p,a,b in rows:
            v=bins[u,((s-1)//2)%mod];sign=1 if t%2==0 else -1
            v[0]+=sign*a;v[1]+=sign*b
        out.append(sum(map(Q,bins.values())))
        if mod>maxindex:break
        mod*=2
    return out

def root_checks():
    data=HERE.joinpath('frontier.txt').read_bytes()
    V={int(l.split()[0]):int(l.split()[1]) for l in data.decode().splitlines() if l and not l.startswith('#')}
    check(len(V)==627 and all(edge(s)[0] in V for s in V),'frontier fixed and closed')
    check(all(s==1 or V[edge(s)[0]]<r for s,r in V.items()),'frontier grounding')
    memo={s:(s,0,0,0) for s in V}
    def route(s):
        trail=[];seen=set();x=s
        while x not in memo:
            check(x not in seen and len(trail)<10000,'unresolved source is not dropped')
            seen.add(x);y,t=edge(x);trail.append((x,t));x=y
        u,t,d,p=memo[x]
        for x,dt in reversed(trail):
            if d==0:p=x
            t+=dt;d+=1;memo[x]=(u,t,d,p)
        return memo[s]
    out=[]
    for N in (4096,11002,18275,65536,131072,262144):
        rows=[];S=0
        for s in range(1,4*N,2):
            if s%8==5:continue
            a,b=charge(s,N);S+=Q((a,b));u,t,d,p=route(s);rows.append((s,u,t,d,p,a,b))
        E=energies(rows);check(E[-1]==S,'full resolution retains source norm')
        B=budgets(E)
        check(E[0]<=S+B['adaptive_1to3'],'blocked telescope')
        check(B['fourary']<=B['binary'] and B['eightary']<=B['binary'],'blocking removes intermediate debt')
        # Three-chart mean at every complete two-bit block.
        means=[]
        for k in range(0,len(E)-2,2):
            means.append(str(Fraction(E[k]+2*E[k+2],3)))
        prefixes={}
        if N in (4096,11002,18275):
            for name,pos in (('depth',3),('clock',2)):
                # Streaming all physical arrival cutoffs, retaining the full denominator S.
                levels=sorted({r[pos] for r in rows});K=len(E)-1
                bins=[defaultdict(lambda:[0,0]) for _ in range(K+1)];en=[0]*(K+1)
                batches=defaultdict(list)
                for r in rows:batches[r[pos]].append(r)
                peaks={k:[0,0] for k in ('binary','fourary','eightary','adaptive_1to3')}
                for level in levels:
                    for s,u,t,d,p,a,b in batches[level]:
                        sign=1 if t%2==0 else -1
                        for k in range(K+1):
                            v=bins[k][u,((s-1)//2)%(1<<k)];old=Q(v)
                            v[0]+=sign*a;v[1]+=sign*b;en[k]+=Q(v)-old
                    for key,value in budgets(en).items():
                        if key in peaks and value>peaks[key][0]:peaks[key]=[value,level]
                check(en==E,'prefix and full formulations agree')
                prefixes[name]=dict(cutoffs=len(levels),peaks=peaks)
        digest=hashlib.sha256(''.join(','.join(map(str,r))+'\n' for r in rows).encode()).hexdigest()
        out.append(dict(N=N,sources=len(rows),S=S,E=E,budgets=B,three_chart_means=means,
                        prefix=prefixes,max_clock=max(r[2] for r in rows),max_depth=max(r[3] for r in rows),
                        source_port_clock_digest=digest))
    return dict(frontier_sha256=hashlib.sha256(data).hexdigest(),scales=out)

def paley(p):
    check(p%4==3 and all(p%d for d in range(2,isqrt(p)+1)),'Paley prime scope')
    squares={x*x%p for x in range(1,p)}
    return [sum(1<<j for j in range(p) if (j-i)%p in squares) for i in range(p)]

def matrix(T):
    return [[0 if i==j else 1 if T[i]>>j&1 else -1 for j in range(len(T))] for i in range(len(T))]

def parse(T):
    n=len(T);FULL=(1<<n)-1
    @lru_cache(None)
    def back(A,B):
        total=0
        while A:
            v=A&-A;A-=v;total+=(B&~T[v.bit_length()-1]).bit_count()
        return total
    costs=[0]*(FULL+1);trees=[None]*(FULL+1)
    for S in range(1,FULL+1):
        if S&(S-1)==0:trees[S]=S.bit_length()-1;continue
        best=10**9;tree=None
        A=(S-1)&S
        while A:
            B=S^A;c=costs[A]+costs[B]+back(A,B)
            if c<best:best=c;tree=('J',A,B)
            A=(A-1)&S
        # Every unordered triple once, then both directed cycles.
        bit=S&-S;rest=S^bit;Asub=rest
        while True:
            A=Asub|bit;R=S^A
            if R and R&(R-1):
                first=R&-R;B=(R-1)&R
                while B:
                    C=R^B
                    if B&first:
                        internal=costs[A]+costs[B]+costs[C]
                        v=back(A,B)+back(B,C)+back(C,A)
                        allcross=A.bit_count()*B.bit_count()+B.bit_count()*C.bit_count()+C.bit_count()*A.bit_count()
                        for c,children in ((internal+v,(A,B,C)),(internal+allcross-v,(A,C,B))):
                            if c<best:best=c;tree=('C',*children)
                    B=(B-1)&R
            if Asub==0:break
            Asub=(Asub-1)&rest
        costs[S]=best;trees[S]=tree
    B=[0]*n
    def emit(S):
        tr=trees[S]
        if isinstance(tr,int):return [tr]
        children=[emit(s) for s in tr[1:]]
        pairs=((0,1),) if tr[0]=='J' else ((0,1),(1,2),(2,0))
        for i,j in pairs:
            for x in children[i]:
                for y in children[j]:B[x]|=1<<y
        return sum(children,[])
    labels=emit(FULL)
    check(sorted(labels)==list(range(n)),'each named port occurs exactly once')
    actual=sum((T[i]^B[i]).bit_count() for i in range(n))//2
    check(actual==costs[FULL],'reparsed baseline exact owned exceptions')
    return costs[FULL],B

def hpath(T,shift=0):
    P=[];n=len(T)
    for v in [(i+shift)%n for i in range(n)]:
        j=0
        while j<len(P) and not(T[v]>>P[j]&1):j+=1
        P.insert(j,v)
    check(all(T[x]>>y&1 for x,y in zip(P,P[1:])),'literal path guard')
    return P

def complement(T,P):
    n=len(T);U=[((1<<n)-1)^T[i]^(1<<i) for i in range(n)]
    for x,y in zip(P,P[1:]):U[x]^=1<<y;U[y]^=1<<x
    return U

def grammar_checks():
    allfour=[]
    for mask in range(64):
        T=[0]*4;k=0
        for i in range(4):
            for j in range(i+1,4):
                if mask>>k&1:T[i]|=1<<j
                else:T[j]|=1<<i
                k+=1
        c,_=parse(T);allfour.append(c)
    records=[]
    for n in (7,9,11):
        T=paley(n) if n!=9 else [sum(1<<j for j in range(9) if (j//3-i//3)%3==1 or (j//3==i//3 and (j-i)%3==1)) for i in range(9)]
        count=7 if n<11 else 3
        for step in range(count):
            cost,B=parse(T);P=hpath(T,step);U=complement(T,P);K=matrix(T);L=matrix(U)
            D=[[0]*n for _ in range(n)]
            for x,y in zip(P,P[1:]):D[x][y]=1;D[y][x]=-1
            check(all(L[i][j]==-K[i][j]+2*D[i][j] for i in range(n) for j in range(n)),'skew path update')
            for i in range(n):
                expected=n-1-T[i].bit_count()+(i==P[0])-(i==P[-1])
                check(U[i].bit_count()==expected,'score endpoints retain path addresses')
            if step:check(abs(cost-records[-1]['cost'])<=n-1,'one-step reparse tariff')
            records.append(dict(n=n,step=step,cost=cost,path=P,tournament=T,baseline=B))
            T=U
    paley_checks=[]
    for p in (3,7,11,19,31,43):
        K=matrix(paley(p))
        check(all(sum(K[i][k]*K[j][k] for k in range(p))==(p-1 if i==j else -1) for i in range(p) for j in range(p)),'Paley flat Gram')
        paley_checks.append(p)
    q=65519;paley_prime=all(q%d for d in range(2,isqrt(q)+1))
    check(paley_prime and q%4==3,'large explicit prime')
    L=(q-1).bit_length();root=isqrt(q)+1
    lower=(q*(q-1-3*root*L)+3)//4
    return dict(order_four_repair_costs=allfour,native_walks=records,paley_gram_controls=paley_checks,
                large_bound=dict(q=q,ceil_sqrt=root,ceil_log2=L,cyclic_module_exception_lower_bound=lower))

if __name__=='__main__':
    print(json.dumps(dict(status='PASS',lattice=lattice_checks(),root=root_checks(),grammar=grammar_checks()),indent=2,sort_keys=True))
