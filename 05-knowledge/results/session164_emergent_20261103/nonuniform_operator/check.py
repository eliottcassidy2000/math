#!/usr/bin/env python3
"""Independent literal controls for the S164 interface packet; assertions stay active with -O."""
from itertools import combinations, product
from math import gcd, comb
from functools import reduce
import json

def require(ok, msg):
    if not ok: raise RuntimeError(msg)

def validate(rows):
    n=len(rows)
    require(all(not(r>>i&1) and r>>n==0 for i,r in enumerate(rows)), 'invalid rows')
    require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)), 'not tournament')

def ordinal(A,B):
    n=len(A); m=len(B)
    return tuple(r|(((1<<m)-1)<<n) for r in A)+tuple(r<<n for r in B)

def extend(A,D):
    n=len(A)
    return tuple(r|((1<<n) if not(D>>i&1) else 0) for i,r in enumerate(A))+(D,)

def subset_paths(A):
    n=len(A); N=1<<n
    dp=[[0]*n for _ in range(N)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,N):
        for u,c in enumerate(dp[mask]):
            if not c:continue
            rest=A[u]&~mask
            while rest:
                bit=rest&-rest;rest-=bit
                dp[mask|bit][bit.bit_length()-1]+=c
    h=[sum(row) for row in dp];h[0]=1
    return h,dp

def H(A):return subset_paths(A)[0][-1]

def covers_from_h(h,rmax):
    N=len(h); fs=[[0]*N];fs[0][0]=1
    for r in range(1,rmax+1):
        f=[0]*N
        for mask in range(N):
            sub=mask;v=0
            while True:
                v+=h[sub]*fs[-1][mask^sub]
                if sub==0:break
                sub=(sub-1)&mask
            f[mask]=v
        fs.append(f)
    return fs

def Fs(A,rmax):return [f[-1] for f in covers_from_h(subset_paths(A)[0],rmax)]

def strong(A):
    if not A:return False
    n=len(A);seen=1;todo=1
    while todo:
        b=todo&-todo;todo-=b;v=b.bit_length()-1;new=A[v]&~seen;seen|=new;todo|=new
    if seen!=(1<<n)-1:return False
    Ar=tuple(((1<<n)-1)^r^(1<<i) for i,r in enumerate(A))
    seen=1;todo=1
    while todo:
        b=todo&-todo;todo-=b;v=b.bit_length()-1;new=Ar[v]&~seen;seen|=new;todo|=new
    return seen==(1<<n)-1

def observer(A,p):
    h=H(A); Aa=extend(A,1<<p);ha=H(Aa);haa=H(extend(Aa,1<<len(A)))
    require((ha-h)%2==0 and (haa-ha)%2==0,'observer parity')
    return h,(ha-h)//2,(haa-ha)//2

def attach(A,p,Q):
    out=list(ordinal(A,Q));n=len(A)
    out[p]^=1<<n;out[n]|=1<<p
    return tuple(out)

def mul(v,M):return tuple(sum(v[i]*M[i][j] for i in range(3)) for j in range(3))
def content(v):return reduce(gcd,v)
def coeffs(vals):
    c=list(vals);k=(len(vals)-1).bit_length()
    for i in range(k):
        for mask in range(1<<k):
            if mask>>i&1:c[mask]-=c[mask^(1<<i)]
    return c

def general_coefficients(W,D,r):
    h,ends=subset_paths(W);n=len(W);N=1<<n;full=N-1
    Wr=tuple(full^row^(1<<i) for i,row in enumerate(W));_,starts=subset_paths(Wr)
    e=[sum(row[i] for i in range(n) if D>>i&1) for row in ends]
    s=[sum(row[i] for i in range(n) if D>>i&1) for row in starts]
    fs=covers_from_h(h,r)
    E=sum(e[m]*fs[r-1][full^m] for m in range(N))
    B=sum(s[m]*fs[r][full^m] for m in range(N))
    C=0
    for left in range(1,N):
        rem=full^left;right=rem
        while right:
            C+=e[left]*s[right]*fs[r-1][rem^right]
            right=(right-1)&rem
    return B,fs[r][full]-C,E

TRI=(2,4,1);W=ordinal(TRI,TRI);D=9
P=((5,3,4),(8,7,8),(6,6,7));Bmat=((1,2,2),(4,3,4),(4,2,3))
Q=(8,1,3,6); Batom=(0,1,3)
def G(Z):return extend(ordinal(W,Z),D)
def f(r):return r**3+2*r
def q(r):return r*r-r+1
def op_coeff(r):return q(r+1)*(f(r+1)+f(r)),f(r)**2-q(r)**2-q(r+1)**2,q(r)*(f(r)+f(r-1))

def main():
    report={}
    # General one-apex operator, not merely the special two-triangle case.
    generic=0
    heads=[(),(0,),TRI,(2,4,1,7)]
    tails=[(),(0,),(2,0),TRI]
    for head in heads:
        validate(head)
        for d in range(1<<len(head)):
            for tail in tails:
                out=extend(ordinal(head,tail),d)
                outF=Fs(out,3);zF=Fs(tail,4)
                for r in range(1,4):
                    b,a,e=general_coefficients(head,d,r)
                    require(outF[r]==r*(b*zF[r+1]+a*zF[r]-e*zF[r-1]),'general tail operator')
                    generic+=1
    report['general_operator_equalities']=generic
    for r in range(1,6):require(general_coefficients(W,D,r)==op_coeff(r),'head coefficient polynomial')
    # Independent direct subset convolution controls through fourteen vertices.
    g1=G(());g2=G(g1)
    validate(g1);validate(g2)
    g1F=Fs(g1,5);g2F=Fs(g2,3)
    require(g1F[1:4]==[41,712,5763],'G(empty) controls')
    require(g2F[1]==31999 and g2F[2]==3749464,'G twice controls')
    require(strong(g1) and strong(g2),'strong iteration controls')
    for r in range(1,4):
        b,a,e=op_coeff(r)
        require(g2F[r]==r*(b*g1F[r+1]+a*g1F[r]-e*g1F[r-1]),'iterated operator')
    report['iterates']={'H_G1':g1F[1],'F_G1_2':g1F[2],'H_G2':g2F[1],'F_G2_2':g2F[2],
                        'H_G3_formula_only':45*g2F[2]-g2F[1]}
    # Sharp strong witnesses with exactly 3t+s vertices.
    sharp=[]
    for s,d in [(1,0),(1,1),(1,2),(2,0)]:
        Z=()
        for _ in range(s-1):Z=ordinal(Z,g1)
        for _ in range(d):Z=ordinal(Z,TRI)
        out=G(Z);actual=H(out)
        expected=41 if not Z else 45*712**(s-1)*12**d-41**(s-1)*3**d
        require(actual==expected and strong(out),'sharp strong count')
        require(actual%3**d==0 and actual%3**(d+1)!=0,'sharp valuation')
        require(len(out)==3*(2*s+d)+s,'sharp vertex count')
        sharp.append({'s':s,'t':2*s+d,'n':len(out),'H':actual,'v3':d})
    report['sharp_witnesses']=sharp
    # Full 64-interface response on the two triangles.
    vals=[H(extend(W,d)) for d in range(64)];co=coeffs(vals)
    for mask,c in enumerate(co):
        require(not c or mask.bit_count()<=2,'one-apex degree')
        if mask:require(c%2==0,'coefficient parity')
        k=(mask&7).bit_count();l=(mask>>3).bit_count()
        expected=9+36*k+6*l-6*(comb(k,2)+comb(l,2))-10*k*l
        require(vals[mask]==expected,'complete triangle interface polynomial')
    report['two_triangle_interface_assignments']=64
    # Two exceptional vertices: support <=2 per exception and no cycles/old degree >2.
    vals=[]
    for mask in range(64):
        A=extend(TRI,mask&7)
        A=extend(A,(mask>>3)&7) # The first new vertex beats the second.
        vals.append(H(A))
    co=coeffs(vals)
    for mask,c in enumerate(co):
        if not c:continue
        require((mask&7).bit_count()<=2 and (mask>>3).bit_count()<=2,'multi-apex bidegree')
        # The sole possible cycle is a four-cycle using both apices and two old vertices.
        common=(mask&7)&(mask>>3)
        require(common.bit_count()<2,'linear forest support')
    report['two_apex_interface_assignments']=64
    # Exact observer content and transfers, on all labeled tournaments of orders 1..4.
    observer_controls=0
    for n in range(1,5):
        edges=list(combinations(range(n),2))
        for mask in range(1<<len(edges)):
            A=[0]*n
            for bit,(u,v) in enumerate(edges):
                if mask>>bit&1:A[u]|=1<<v
                else:A[v]|=1<<u
            A=tuple(A)
            for p in range(n):
                v=observer(A,p)
                vq=observer(attach(A,p,Q),n+3)
                vb=observer(attach(A,p,Batom),n+2)
                require(vq==mul(v,P) and vb==mul(v,Bmat),'literal transfer')
                require(content((v[0],vq[0],vb[0]))==content(v),'count observer content')
                require(content(vb)==content(v),'unimodular b content')
                require(content(vq)//content(v) in (1,5),'q content jump')
                observer_controls+=1
    report['literal_observer_controls']=observer_controls
    R4=(996,993,994,999,975,912,800,576,128,256)
    one=extend(R4,1);v=observer(one,9)
    require(v==(1965,1103,1564),'R4 one-edge observer')
    qout=attach(one,9,Q);actual=H(qout)
    require(actual==28033 and strong(qout) and gcd(actual,15)==1,'R4 delayed simultaneous escape')
    report['R4_delayed_unit_mod15']=actual
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
