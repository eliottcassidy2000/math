#!/usr/bin/env python3
"""Exact response-bucket census of tournaments with prime quotients only C3.
Standard library; no graph-isomorphism calls in production census.
"""
import argparse, itertools, json, math, time
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path


def require(ok, message):
    if not ok: raise RuntimeError(message)

@lru_cache(None)
def join(a,b):
    """Ordered path-cover profile of the ordinal sum."""
    out=[0]*(len(a)+len(b)+1)
    for i,ai in enumerate(a,1):
        for j,bj in enumerate(b,1):
            for ell in range(min(i,j)+1):
                k=i+j-ell
                coef=math.factorial(k)//(math.factorial(i-ell)*math.factorial(j-ell)*math.factorial(ell))
                out[k]+=ai*bj*coef
    return tuple(out[1:])


def build_kernel(N):
    """Coefficients of z R/(1-zR) for the cyclic three-state skeleton.
    K[r][k] counts colored words of multiplicity r with k ordered path slots.
    A same-color transition MUST be cut; a forward C3 transition may be cut.
    """
    dp={}; K={}
    for total in range(1,N+1):
        for r0 in range(total+1):
            for r1 in range(total-r0+1):
                r=(r0,r1,total-r0-r1)
                ends=[]
                for j in range(3):
                    out=[0]*(total+1)
                    if r[j]:
                        prev=list(r);prev[j]-=1;prev=tuple(prev)
                        if total==1: out[1]=1
                        else:
                            for i in range(3):
                                for k,val in enumerate(dp[(prev,i)]):
                                    if val:
                                        out[k+1]+=val
                                        if (j-i)%3==1:out[k]+=val
                    dp[r,j]=tuple(out);ends.append(out)
                if min(r)>0:K[r]=tuple(map(sum,zip(*ends)))
    return K


def make_cycle(K):
    @lru_cache(None)
    def cycle(a,b,c):
        out=[0]*(len(a)+len(b)+len(c)+1)
        for i,ai in enumerate(a,1):
            for j,bj in enumerate(b,1):
                for k,ck in enumerate(c,1):
                    scale=ai*bj*ck
                    for ell,v in enumerate(K[i,j,k]):out[ell]+=scale*v
        return tuple(out[1:])
    return cycle


def add(d,q,count,mass):
    old=d.get(q,(0,Fraction(0)))
    d[q]=(old[0]+count,old[1]+mass)


def census(N):
    K=build_kernel(N);cycle=make_cycle(K)
    A=[{} for _ in range(N+1)];S=[{} for _ in range(N+1)]
    rows=[]
    for n in range(1,N+1):
        if n==1:add(S[n],(1,),1,Fraction(1))
        # Multisets of RESPONSE types, not graph types. The diagonal needs
        # +2*m in class counts, whereas reciprocal-Aut mass has no such term.
        pool=[(len(q),q,m,u) for sz in range(1,n-1) for q,(m,u) in A[sz].items()]
        pool.sort(key=lambda v:(v[0],v[1]))
        for ia,(sa,a,ma,ua) in enumerate(pool):
            if 3*sa>n:break
            for ib in range(ia,len(pool)):
                sb,b,mb,ub=pool[ib]
                if sa+2*sb>n:break
                sc=n-sa-sb
                for ic in range(ib,len(pool)):
                    sz,c,mc,uc=pool[ic]
                    if sz<sc:continue
                    if sz>sc:break
                    q=cycle(a,b,c)
                    if ia==ic: count=(ma**3+2*ma)//3;mass=ua**3/3
                    elif ia==ib:count=ma*ma*mc;mass=ua*ua*uc
                    elif ib==ic:count=ma*mb*mb;mass=ua*ub*ub
                    else:count=2*ma*mb*mc;mass=2*ua*ub*uc
                    add(S[n],q,count,mass)
        A[n]=dict(S[n])
        for i in range(1,n):
            for a,(ma,ua) in S[i].items():
                for b,(mb,ub) in A[n-i].items():
                    q=join(a,b) if a<=b else join(b,a)
                    add(A[n],q,ma*mb,ua*ub)
        mass=sum((u for m,u in A[n].values()),Fraction(0))*math.factorial(n)
        require(mass.denominator==1,'nonintegral labeled count')
        rows.append(dict(n=n,classes=sum(m for m,u in A[n].values()),strong_classes=sum(m for m,u in S[n].values()),response_types=len(A[n]),strong_response_types=len(S[n]),H_values=len(set(q[0] for q in A[n])),labeled_tournaments=int(mass),largest_response_fiber=max(m for m,u in A[n].values()),max_H=max(q[0] for q in A[n])))
    return A,S,rows,cycle


def independent_scalar_counts(N):
    # A=S/(1-S), S=x+(A(x)^3+2A(x^3))/3. Ordinary coefficients.
    a=[0]*(N+1);s=[0]*(N+1)
    # Reciprocal-aut mass is the labeled EGF: s=x+a^3/3.
    am=[Fraction(0)]*(N+1);sm=am.copy()
    for n in range(1,N+1):
        cubic=sum(a[i]*a[j]*a[n-i-j] for i in range(1,n) for j in range(1,n-i))
        require((cubic+(2*a[n//3] if n%3==0 else 0))%3==0,'Burnside division')
        s[n]=(n==1)+(cubic+(2*a[n//3] if n%3==0 else 0))//3
        sm[n]=Fraction(n==1)+sum((am[i]*am[j]*am[n-i-j] for i in range(1,n) for j in range(1,n-i)),Fraction(0))/3
        a[n]=s[n]+sum(s[i]*a[n-i] for i in range(1,n))
        am[n]=sm[n]+sum(sm[i]*am[n-i] for i in range(1,n))
    return a,s,[int(am[n]*math.factorial(n)) for n in range(N+1)]


def substitute(adj,children):
    sizes=[len(t) for t in children];starts=[sum(sizes[:i]) for i in range(len(sizes))]
    n=sum(sizes);out=[[0]*n for _ in range(n)]
    for i,t in enumerate(children):
        for a in range(len(t)):
            for b in range(len(t)):out[starts[i]+a][starts[i]+b]=t[a][b]
        for j,u in enumerate(children):
            if i!=j:
                for a in range(len(t)):
                    for b in range(len(u)):out[starts[i]+a][starts[j]+b]=adj[i][j]
    return out


def direct_profile(t):
    n=len(t);q=[0]*n
    for perm in itertools.permutations(range(n)):
        bad=sum(not t[a][b] for a,b in zip(perm,perm[1:]))
        good=n-1-bad
        for k in range(bad+1,n+1):q[k-1]+=math.comb(good,k-1-bad)
    return tuple(q)


def controls(cycle):
    one=[[0]];C=[[0,1,0],[0,0,1],[1,0,0]];J=[[0,1],[0,0]]
    chain2=substitute(J,[one,one]);A=substitute(J,[C,one]);B=substitute(J,[one,C])
    cases=[(C,[one,one,one]),(C,[chain2,one,one]),(C,[C,one,one]),(C,[chain2,C,one]),(C,[A,one,one]),(J,[C,C]),(J,[A,chain2])]
    checked=[]
    for adj,kids in cases:
        qs=[direct_profile(t) for t in kids]
        predicted=cycle(*qs) if len(qs)==3 else join(*qs)
        actual=direct_profile(substitute(adj,kids))
        require(predicted==actual,'independent permutation profile mismatch')
        checked.append(dict(n=sum(map(len,kids)),q=list(actual)))
    require(direct_profile(A)==direct_profile(B),'A/B response control')
    # The 12-vertex 81/27 boundary is an inherited fact, not re-audited here.
    # m=2 equal-response child classes -> (8+4)/3=4 ordinary root classes,
    # versus naive response-orbit count 1 and naive labeled-tuple/3=8/3.
    return dict(permutation_checks=checked,diagonal_Burnside_count_for_two_child_classes=(2**3+2*2)//3,A_B_profile=list(direct_profile(A)))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-n',type=int,default=12);parser.add_argument('--out',default='receipt.json');parser.add_argument('--bucket-max-n',type=int,default=12);args=parser.parse_args()
    t=time.perf_counter();A,S,rows,cycle=census(args.max_n)
    a,s,labeled=independent_scalar_counts(args.max_n)
    for row in rows:
        n=row['n'];require((row['classes'],row['strong_classes'],row['labeled_tournaments'])==(a[n],s[n],labeled[n]),'scalar census mismatch')
    hist={}
    for n in range(1,args.max_n+1):
        counts=defaultdict(int)
        for q,(m,u) in A[n].items():counts[q[0]]+=m
        hist[str(n)]={str(h):counts[h] for h in sorted(counts)}
    out=dict(status='FINITE-EXACT',family='Canonical prime quotients restricted to C3; arbitrary ordered strong-component sums',max_n=args.max_n,rows=rows,controls=controls(cycle),cache=dict(join=join.cache_info()._asdict(),C3=cycle.cache_info()._asdict()),elapsed_seconds=round(time.perf_counter()-t,3),H_distribution=hist,response_buckets={str(n):[dict(q=list(q),classes=m,reciprocal_aut_mass=str(u)) for q,(m,u) in sorted(A[n].items())] for n in range(1,min(args.max_n,args.bucket_max_n)+1)})
    Path(args.out).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('H_distribution','response_buckets')},indent=2))
if __name__=='__main__':main()
