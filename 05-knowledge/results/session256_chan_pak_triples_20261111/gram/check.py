#!/usr/bin/env python3
"""New Gram claims only; exact profiles and independent vertex subset DP."""
from functools import lru_cache
from itertools import product
from math import factorial, comb
import json, random
from fractions import Fraction
from pathlib import Path

def tournament(n, mask):
    rows=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if (mask>>k)&1 else (j,i)
            rows[u]|=1<<v; k+=1
    return tuple(rows)

def h_all(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    ans=[0]*(1<<n); ans[0]=1
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for j in range(n):
            value=dp[mask][j]
            if not value: continue
            todo=rows[j]&~mask
            while todo:
                bit=todo&-todo; todo-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
        ans[mask]=sum(dp[mask])
    return ans

@lru_cache(None)
def profile(rows):
    h=h_all(rows)
    @lru_cache(None)
    def pc(mask):
        if not mask:return (1,)
        n=mask.bit_count(); first=mask&-mask; rem=mask^first
        out=[0]*(n+1); sub=rem
        while True:
            part=sub|first
            for k,v in enumerate(pc(mask^part)):out[k+1]+=h[part]*v
            if sub==0:break
            sub=(sub-1)&rem
        return tuple(out)
    return (0,)+tuple(factorial(r)*v for r,v in enumerate(pc((1<<len(rows))-1)) if r)

def c3(S,T,U):
    blocks=(S,T,U); offsets=(0,len(S),len(S)+len(T)); rows=[]
    for b,B in enumerate(blocks):
        nxt=(b+1)%3
        cross=((1<<len(blocks[nxt]))-1)<<offsets[nxt]
        rows.extend((r<<offsets[b])|cross for r in B)
    return tuple(rows)

def entry(f,i):return f[i] if 0<=i<len(f) else 0

def gram(s,t,u):
    ans=0
    for r in range(1,len(u)):
        sm,sr,sp=(entry(s,r+d) for d in (-1,0,1))
        tm,tr,tp=(entry(t,r+d) for d in (-1,0,1))
        ans+=u[r]*((sm+sr)*(tm+tr)+(sr+sp)*(tr+tp)+sr*tr)
    return ans

def tt(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

@lru_cache(None)
def stirling(n,k):
    if n==0:return int(k==0)
    if k<=0:return 0
    return stirling(n-1,k-1)+k*stirling(n-1,k)

def transitive_profile(n):return (0,)+tuple(factorial(r)*stirling(n,r) for r in range(1,n+1))

def det(M):
    M=[list(map(Fraction,row)) for row in M]; out=Fraction(1)
    for i in range(len(M)):
        pivot=next((j for j in range(i,len(M)) if M[j][i]),None)
        if pivot is None:return 0
        if pivot!=i:M[i],M[pivot]=M[pivot],M[i];out=-out
        d=M[i][i];out*=d
        for j in range(i+1,len(M)):
            q=M[j][i]/d
            for k in range(i+1,len(M)):M[j][k]-=q*M[i][k]
            M[j][i]=0
    return int(out) if out.denominator==1 else str(out)

def main():
    rng=random.Random(256); counts={}; witnesses=[]
    # Every labelled source through order four, singleton observer and
    # a seeded partner of order at most four. Actual child orders <=9.
    cases=[]
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            S=tournament(n,mask); m=1+mask%4
            T=tournament(m,rng.randrange(1<<(m*(m-1)//2)))
            cases.append((S,T,tt(1)))
    # Nontransitive observers and larger profile widths.
    for k in range(50):
        sizes=(1+rng.randrange(4),1+rng.randrange(4),1+rng.randrange(3))
        cases.append(tuple(tournament(n,rng.randrange(1<<(n*(n-1)//2))) for n in sizes))
    for S,T,U in cases:
        s,t,u=map(profile,(S,T,U)); actual=h_all(c3(S,T,U))[-1]
        assert actual==gram(s,t,u),(S,T,U,actual,gram(s,t,u))
        assert gram(s,t,u)**2<=gram(s,s,u)*gram(t,t,u)
        if len(U)==1:
            h,p=entry(s,1),entry(s,2); k,q=entry(t,1),entry(t,2)
            assert gram(s,s,u)*gram(t,t,u)-gram(s,t,u)**2==2*(h*q-k*p)**2
    counts['independent_actual_graph_counts']=len(cases)
    ranks=[]
    for c in range(1,7):
        u=transitive_profile(c)
        M=[[gram(transitive_profile(a),transitive_profile(b),u) for b in range(1,c+2)] for a in range(1,c+2)]
        d=det(M);assert d>0
        enlarged=[[gram(transitive_profile(a),transitive_profile(b),u) for b in range(1,c+3)] for a in range(1,c+3)]
        assert det(enlarged)==0
        ranks.append({'observer_order':c,'rank':c+1,'positive_gram_determinant':d})
    counts['exact_rank_controls']=len(ranks)
    for a,b in product(range(1,15),repeat=2):
        P=lambda a:2**a-1
        R=lambda a:3**a-2*2**a+1
        V=lambda a:2**a-2
        expected=2+3*P(a)*P(b)+2*R(a)*R(b)+2*V(a)*V(b)
        assert expected==gram(transitive_profile(a),transitive_profile(b),transitive_profile(2))
    counts['packing_two_formula_controls']=14**2
    recurrences=[]
    for c in range(1,7):
        coeff=[1]
        for j in range(1,c+2):
            nxt=[0]*(len(coeff)+1)
            for k,v in enumerate(coeff):nxt[k]-=j*v;nxt[k+1]+=v
            coeff=nxt
        for b in range(1,8):
            vals=[gram(transitive_profile(a),transitive_profile(b),transitive_profile(c)) for a in range(1,20)]
            for a in range(len(vals)-len(coeff)+1):assert sum(q*vals[a+k] for k,q in enumerate(coeff))==0
        recurrences.append({'c':c,'coefficients_low_to_high':coeff})
    counts['recurrence_families']=42
    result={'status':'PASS','scope':'new claims only; proofs carry the unbounded statements','counts':counts,'rank_controls':ranks,'recurrences':recurrences,'hostile_control':{'graph':'C3[TT2,TT2,TT2]','actual_H':h_all(c3(tt(2),tt(2),tt(2)))[-1],'illegal_intercyclic_value':29},'boundaries':['Gram vectors must be realized by tournaments.','Fixed observer rank does not bound moving-observer rank.','Singleton observer norm has discriminant 2, not a Pythagorean norm.']}
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
