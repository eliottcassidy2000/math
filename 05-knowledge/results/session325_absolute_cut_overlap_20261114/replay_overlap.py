#!/usr/bin/env python3
"""Standalone exact S325 overlap certificate, standard library only.
Use --full to independently count both order-19 witnesses of 7*3625619.
Every check survives python -O. No exploration file is read.
"""
from array import array
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
import json, sys

R13=(24,9,35,4,14,27)
R15=(108,97,2,6,79,92,12)
U571=(478,412,448,100,460,407,290,72,136)
W571=(138,212,105,66,45,139,49,92)
U401=(320,117,473,211,1,477,144,275,90)
W401=(130,28,193,53,5,151,59,90)

def need(c,why):
    if not c: raise RuntimeError(why)

def bits(s):
    while s:
        b=s&-s; yield b.bit_length()-1; s-=b

def valid(r):
    n=len(r)
    for i,row in enumerate(r):
        need(0<=row<1<<n and not(row>>i&1),'invalid row')
        for j in range(i):need((row>>j&1)+(r[j]>>i&1)==1,'invalid arc')

def strong(r):
    n=len(r);F=(1<<n)-1
    for rows in (r,tuple(F^row^(1<<i) for i,row in enumerate(r))):
        reached=1
        while True:
            after=reached
            for i in bits(reached):after|=rows[i]
            if after==reached:break
            reached=after
        if reached!=F:return False
    return True

def hp_all(r):
    n=len(r);size=1<<n;E=[[0]*n for _ in range(size)];H=[0]*size;H[0]=1
    for i in range(n):E[1<<i][i]=1
    for s in range(1,size):
        for i in bits(s):
            value=E[s][i];H[s]+=value
            if value:
                for j in bits(r[i]&~s):E[s|1<<j][j]+=value
    return H

def hp_large(r):
    # Independent endpoint recurrence: compute predecessors, not outgoing updates.
    n=len(r);size=1<<n;E=array('Q',[0])*(size*n)
    inc=[sum(1<<j for j in range(n) if r[j]>>i&1) for i in range(n)]
    for i in range(n):E[(1<<i)*n+i]=1
    for s in range(1,size):
        if s&(s-1)==0:continue
        for i in bits(s):
            rem=s^(1<<i);base=rem*n
            E[s*n+i]=sum(E[base+j] for j in bits(inc[i]&rem))
    return sum(E[(size-1)*n+i] for i in range(n))

def profile(r):
    valid(r);H=hp_all(r);F=len(H)-1;Q=[0]*len(H)
    for s in range(len(H)):
        a=s
        while True:
            Q[s]+=H[a]*H[s^a]
            if not a:break
            a=(a-1)&s
    return H[F],Q[F],sum(H[a]*Q[F^a] for a in range(F+1))

def ordinal(blocks):
    total=sum(map(len,blocks));off=0;out=[]
    for B in blocks:
        k=len(B);future=((1<<total)-1)^((1<<(off+k))-1)
        out.extend((row<<off)|future for row in B);off+=k
    return tuple(out)

def tt(m):return ordinal([(0,)]*m)

def inflate(Q,U,V):
    B=[U]+[(0,)]*(len(Q)-2)+[V];offs=[];n=0
    for r in B:offs.append(n);n+=len(r)
    out=[]
    for i,r in enumerate(B):
        cross=sum(((1<<len(B[j]))-1)<<offs[j] for j in bits(Q[i]))
        out.extend((row<<offs[i])|cross for row in r)
    valid(out);return tuple(out)

def itinerary(Q,i,j):
    rem=(i,)+(1,)*(len(Q)-2)+(j,)
    @lru_cache(None)
    def walk(last,left):
        if not any(left):return 1
        value=0
        for a,k in enumerate(left):
            if k and (last<0 or Q[last]>>a&1):
                new=list(left);new[a]-=1
                value+=walk(a,tuple(new))
        return value
    return walk(-1,rem)

def power_index(n):
    return n.bit_length()-1 if n>0 and n&(n-1)==0 else None

def v2(n):
    need(n!=0,'v2 zero');n=abs(n);return (n&-n).bit_length()-1

def overlap(C,d,E,f):
    """Return ('tail',odd_slope,first_power,intercept) or finite (t,u,N) list."""
    need(C>0 and E>0,'positive rail slopes')
    r=v2(C);s=v2(E);a=C>>r;b=E>>s;delta=d-f
    if not delta:return ('tail',a,max(r,s),d) if a==b else []
    ans=[]
    def put(x,y):
        if x>=r and y>=s:
            t=x-r;u=y-s;N=C*(1<<t)-d
            need(N==E*(1<<u)-f,'decoder equality')
            ans.append((t,u,N))
    if a!=b and delta%(a-b)==0:
        z=power_index(delta//(a-b))
        if z is not None:put(z,z)
    v=v2(delta);D=delta//(1<<v)
    if D<a-b and (a-D)%b==0:
        k=power_index((a-D)//b)
        if k is not None and k>0:put(v,v+k)
    elif D>a-b and (D+b)%a==0:
        k=power_index((D+b)//a)
        if k is not None and k>0:put(v+k,v)
    need(len(ans)<=2,'two-intersection bound')
    return sorted(ans)

def isprime(n):
    return n>=2 and all(n%p for p in range(2,isqrt(n)+1))

def main():
    tables=[];table_checks=0
    for Q,wanted in ((R13,((29,15),(24,12),(4,2))),(R15,((53,27),(48,24),(12,6)))):
        valid(Q);need(strong(Q),'quotient strong');bound=len(Q)-1;table=[]
        for i in range(1,bound+1):
            row=[]
            for j in range(1,bound+1):
                x=itinerary(Q,i,j);table_checks+=1
                need(x==(wanted[i-1][j-1] if i<=3 and j<=2 else 0),'complete itinerary table')
                row.append(x)
            table.append(row)
        tables.append(table)
    P=[]
    for U,W,h,A,p,q in ((U571,W571,571,396561,(571,14336,156129),(571,8032,58347)),
                       (U401,W401,401,317427,(401,10960,126435),(401,6080,47025))):
        need(profile(U)==p and profile(W)==q,'literal profiles')
        need(strong(U) and strong(W),'strong source')
        need(-3*h+6*p[1]+2*p[2]==A==-3*h+6*q[1]+6*q[2],'common absolute state')
        need(gcd(h,A)==1,'primitive state');P.append(dict(h=h,A=A,U=p,W=q))
    literal=[]
    tri=(2,4,1)
    for V in ((0,),tt(2),tri):
        d,b,_=profile(V);expected=396561*b-571*d
        for Q,S in ((R13,U571),(R15,W571)):
            T=inflate(Q,S,V);value=hp_all(T)[-1]
            need(value==expected and strong(T),'literal full-chart count')
            deg=[r.bit_count() for r in T]
            need((2 not in deg) if Q==R13 else (2 in deg),'nonisomorphic degree invariant')
            literal.append(dict(quotient=len(Q),reservoir_order=len(V),order=len(T),H=value))
    # Test the strict source-side defect on actual ordinal source profiles.
    growth=0
    for Z in ((0,),tt(2),tri):
        z=profile(Z);p=profile(ordinal((U571,Z)));q=profile(ordinal((W571,Z)))
        diff=(-3*p[0]+6*p[1]+2*p[2])-(-3*q[0]+6*q[1]+6*q[2])
        need(p[0]==q[0] and diff==6*6304*(z[1]-z[2])<0,'source growth obstruction')
        left=-3*p[0]+6*p[1]+2*p[2];right=-3*q[0]+6*q[1]+6*q[2]
        for V0 in ((0,),tri):
            vd,vb,_=profile(V0)
            need(overlap(left*vb,p[0]*vd,right*vb,q[0]*vd)==[],'asynchronous padding cannot repair source growth')
        growth+=1
    need(isprime(792551) and isprime(3625619),'finite atom primality')
    need(396561*64-571==7*3625619,'7p common target')
    need(all(((396561*(1<<m)-571)%7==0)==(m%3==0) for m in range(1,301)),'mod7 phase')
    # Exhaustive finite tests of the complete unbounded arithmetic decoder.
    checks=0;tails=0;max_finite=0;cap=9
    for C,E,d,f in product(range(1,17),range(1,17),range(-7,8),range(-7,8)):
        decoded=overlap(C,d,E,f)
        brute={(t,u,C*(1<<t)-d) for t in range(cap+1) for u in range(cap+1) if C*(1<<t)-d==E*(1<<u)-f}
        if isinstance(decoded,tuple):
            _,a,k,e=decoded;tails+=1;r=v2(C);s=v2(E)
            got={(j-r,j-s,a*(1<<j)-e) for j in range(k,max(r,s)+cap+1) if j-r<=cap and j-s<=cap}
        else:
            max_finite=max(max_finite,len(decoded));got={x for x in decoded if x[0]<=cap and x[1]<=cap}
        need(got==brute,'complete bounded replay of unbounded rail decoder');checks+=1
    need(overlap(3,5,1,1)==[(1,1,1),(2,3,7)],'sharp two-hit example')
    full=[]
    if '--full' in sys.argv:
        for Q,S in ((R13,U571),(R15,W571)):
            T=inflate(Q,S,tt(6));value=hp_large(T)
            need(value==25379333 and strong(T),'independent order-19 path count')
            full.append(dict(order=len(T),H=value,quotient=len(Q)))
    print(json.dumps(dict(status='PASS',itinerary_entries=table_checks,tables=tables,profiles=P,literal_counts=literal,
        source_growth_checks=growth,rail_decoder_checks=checks,coincident_tail_cases=tails,max_finite_hits=max_finite,
        prime_target=792551,seven_prime_target=25379333,prime_cofactor=3625619,full_order19=full),sort_keys=True,indent=2))

if __name__=='__main__':main()
