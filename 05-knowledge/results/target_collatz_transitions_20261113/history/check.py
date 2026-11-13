#!/usr/bin/env python3
"""Exact small controls for the norm chamber and physical Farey partner."""
from itertools import product
from math import gcd,isqrt
from pathlib import Path
import json

def need(ok,why):
    if not ok: raise AssertionError(why)

def inventory(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            rem=T[v]&~mask
            while c and rem:
                b=rem&-rem;rem-=b
                dp[mask|b][b.bit_length()-1]+=c
    H=[1]+[sum(d)for d in dp[1:]]
    return H,dp

def observations(T):
    H,_=inventory(T);full=len(H)-1
    return H[-1],sum(H[s]*H[full^s]for s in range(full+1))

def opposite(T):return tuple(((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T))

def state(T,r):
    H,E=inventory(T);_,S=inventory(opposite(T));full=len(H)-1
    return H[-1],S[-1][r],H[full^(1<<r)],sum(E[s][r]*H[full^s]for s in range(full+1))

def tt(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))

def construct(M,N,prepared):
    blocks=[(0,),(0,),M,(0,),N,(0,),(0,)]
    starts=[];n=0
    for block in blocks:starts.append(n);n+=len(block)
    rows=[]
    for block,start in zip(blocks,starts):
        tail=((1<<n)-1)^((1<<(start+len(block)))-1)
        rows.extend((r<<start)|tail for r in block)
    wins={0,1,4,5}if prepared else{0,3,4}
    a=sum(((1<<len(blocks[i]))-1)<<starts[i]for i in wins)
    return tuple(r|((1<<n) if not(a>>i&1) else 0) for i,r in enumerate(rows))+(a,),starts

def birth(T,r,retain):
    n=len(T);U=list(T);U[r]|=1<<n;U.append(((1<<n)-1)^(1<<r))
    return tuple(U),r if retain else n

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def valuation(n,l):
    k=0
    while n%l==0:n//=l;k+=1
    return k

def padding_selector(h,p,j,q,e):
    if (e-h*j)%6:return []
    C=(e-h*j)//6
    if C==0:raise ValueError('trivial infinite family handled separately')
    P,Q=p*j,h*q;out=set()
    for m in range(valuation(abs(C),2)+1):
        c=C//(1<<m)
        for direction,numerator,denominator in ((0,P-c,Q),(1,Q+c,P)):
            if numerator<=0 or numerator%denominator:continue
            z=numerator//denominator
            if z&(z-1):continue
            d=z.bit_length()-1
            out.add((m,m+d)if direction==0 else(m+d,m))
    return sorted(out)

def main():
    units=[]
    for v in range(1,5001):
        u=isqrt(5*v*v+1)
        if u*u-5*v*v==1:units.append((u,v))
    unit_checks=0;valuation_checks=0
    inert=[l for l in range(3,50,2)if prime(l)and pow(5,(l-1)//2,l)==l-1]
    for h in range(1,26,2):
        for p in range(2*h,10*h+1,2):
            x=6*p+2*h;Y=x*x-5*h*h
            for u,v in units:
                for sign in (-1,1):
                    xx=u*x+5*sign*v*h;hh=sign*v*x+u*h
                    need(abs(xx)<3*abs(hh),'Pell cone exclusion');unit_checks+=1
            for l in inert:
                need(valuation(Y,l)==2*min(valuation(x,l),valuation(h,l)),'inert valuation');valuation_checks+=1
    candidates=[];prime_targets=0
    for Y in range(1,50001,2):
        hits=[]
        for h in range(1,isqrt(Y//191)+1,2):
            x=isqrt(Y+5*h*h)
            if x*x==Y+5*h*h and (x-2*h)%12==0:
                p=(x-2*h)//6
                if p>=2*h:hits.append((h,p))
        if prime(Y):
            need(len(hits)<=1,'prime candidate uniqueness')
            if hits:need(Y%120 in(71,119),'prime residue guard');prime_targets+=1
        if len(hits)>1 and len(candidates)<4:candidates.append((Y,hits))
    composite_hostile=[(11,126),(37,118)]
    need(all((6*p+2*h)**2-5*h*h==604679 and p>=2*h and p%2==0 for h,p in composite_hostile),'composite uniqueness hostile')
    modules=[tt(1),tt(2),(2,4,1)];pairs=[]
    for M,N in product(modules,repeat=2):
        h,p=observations(M);j,q=observations(N)
        T,starts=construct(M,N,False);U,_=construct(M,N,True)
        before,after=state(T,0),state(U,0)
        A=12*p*q+2*p*j+6*h*q-h*j;Z=12*p*q+6*p*j+2*h*q-h*j
        need(before[0]==A and after[0]==Z,'two independent module counts')
        need(Z-A==4*(p*j-h*q),'determinant drift')
        need(before[2]-before[1]==0 and after[2]-after[1]==h*j,'physical defects')
        need(T[0]==U[0] and all((t&1)==(u&1)for t,u in zip(T,U)),'root star')
        changes=[(i,k)for i in range(len(T))for k in range(i+1,len(T))if (T[i]^U[i])>>k&1]
        need(changes==[(starts[i],len(T)-1)for i in(1,3,5)],'three apex changes')
        pairs.append(dict(M=M,N=N,before=before,after=after,drift=Z-A))
    T,_=construct(tt(2),tt(1),True);V,r=birth(T,0,False);V,r=birth(V,r,True)
    need(state(V,r)==(371,123,123,371),'polarized371')
    M=(4,9,82,21,35,143,187,31);N=(6,12,48,21,131,219,159,15)
    h,p=observations(M);j,q=observations(N)
    need((h,p,j,q)==(181,3424,179,3416),'root profile witness')
    D=p*j-h*q;e=6*D+h*j
    A=12*p*q+2*p*j+6*h*q-h*j;Z=A+4*D;Y=3*Z+2*h*j
    need((D,e,A,Z,Y)==(-5400,-1,145259777,145238177,435779329),'Farey witness')
    need(gcd(Z,Y)==1,'primitive witness')
    bound=(6*p-5*h+1)//(6*h)
    need(bound==18,'finite partner bound')
    negative_padding=padding_selector(h,p,j,q,-1)
    positive_padding=padding_selector(h,p,j,q,1)
    need(negative_padding==[(0,0)]and positive_padding==[],'actual witness padding closure')
    padding_checks=0
    for hh,jj in product((1,3,5),repeat=2):
        if hh*jj==1:continue
        for pp,qq in product((2,4),repeat=2):
            pp*=hh;qq*=jj
            for sign in(-1,1):
                selected=padding_selector(hh,pp,jj,qq,sign)
                brute=[(t,s)for t in range(13)for s in range(13)if hh*jj+6*((1<<t)*pp*jj-(1<<s)*hh*qq)==sign]
                need([z for z in selected if max(z)<=12]==brute,'finite padding selector')
                need(all(hh*jj+6*((1<<t)*pp*jj-(1<<s)*hh*qq)==sign for t,s in selected),'padding substitution')
                padding_checks+=1
    out=dict(status='PASS',universe='Pell v1..5000, hodd1..25 and even2h<=p<=10h; targetoddY<=50000; nine literal small module pairs; two inherited8vertex profiles',
             units=units,unit_checks=unit_checks,inert_primes=inert,valuation_checks=valuation_checks,
             prime_targets_with_candidate=prime_targets,composite_multiple_candidates=candidates,
             composite_hostile=dict(Y=604679,pairs=composite_hostile,second_order_bound=2),
             literal_pairs=pairs,polarized371=dict(rows=V,root=r,state=state(V,r)),
             farey=dict(M=M,N=N,h=h,p=p,j=j,q=q,D=D,e=e,A=A,Z=Z,Y=Y,partner_order_bound=bound,
                        negative_padding=negative_padding,positive_padding=positive_padding),
             padding_selector_controls=padding_checks,
             two_unit_steps_from191=dict(initial=[14,1],first=[146,65],second=[2614,1169],second_p=46))
    Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS',unit_checks,'unit tests;',valuation_checks,'valuations;',prime_targets,'prime candidates; nine literal module pairs; actual Farey witness')

if __name__=='__main__':main()
