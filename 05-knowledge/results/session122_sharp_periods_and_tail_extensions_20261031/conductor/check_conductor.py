#!/usr/bin/env python3
"""Exact rational envelope and unshifted-integer hostile controls."""
from fractions import Fraction as F
from itertools import product
import json

def require(v,msg):
    if not v: raise RuntimeError(msg)
def ceil(x): return -(-x.numerator//x.denominator)
def env(x,a):
    n=x.numerator//x.denominator
    return (a*n+min(x-n,a))/x

def bad_measure(lo,hi,v,r):
    ans=F(0)
    left=(v*lo-r).numerator//(v*lo-r).denominator-1
    right=ceil(v*hi+r)+1
    for j in range(left,right+1):
        ans+=max(F(0),min(hi,(F(j)+r)/v)-max(lo,(F(j)-r)/v))
    return ans

def main():
    count=peaks=literal=hostiles=0
    for k in range(2,9):
        for den in range(2*k+1,4*k+8):
            a=F(2,den)
            if k*a>=1: continue
            N=ceil(a*(k-1)/(1-k*a)); c=k*a*N
            require(N-1+a<c<=N,'crossing position')
            for xn in range(1,400):
                x=F(xn,17)
                if x>=c:
                    require(env(x,a)<=F(1,k),'envelope gate');count+=1
            for n in range(1,N+4):
                require((env(n+a,a)<=F(1,k))==(n>=N),'peak iff');peaks+=1
            require(env(c,a)==F(1,k),'cutoff equality')
            eps=min(F(1,1000),(c-(N-1+a))/2)
            require(env(c-eps,a)>F(1,k),'minimal envelope cutoff')
    for v,r,lo,W in product(range(1,8),(F(1,14),F(1,6),F(1,4)),(F(0),F(1,11),F(-3,7)),(F(1,19),F(3,11),F(8,7))):
        require(bad_measure(lo,lo+W,F(v),r)<=W*env(v*W,2*r),'literal envelope');literal+=1
    examples=[]
    for k in range(2,7):
        r=F(1,4*k-2); c0=2*k*r-F(1,100);s=(c0+2*k*r)/2
        b=ceil(2/(r-s/(2*k)))
        vs=[]
        for i in range(1,k+1):
            z=F(2*i-1,2*k); raw=-s*z-s*b
            delta=raw-(raw.numerator//raw.denominator)
            vs.append(s+delta/b)
        import math
        D=math.lcm(*(v.denominator for v in vs));speeds=[int(D*v) for v in vs]
        require(min(speeds)*F(1,D)>c0,'hostile product')
        # Exact covering certificate: the i-th assigned subinterval stays strictly bad.
        for i,v in enumerate(vs):
            z=F(2*i+1,2*k)
            for x in (F(i,k),F(i+1,k)):
                val=v*(b+x);nearest=round(val)
                require(abs(val-nearest)<r,'hostile endpoint bad')
            # both endpoints lie in the same bad tooth, so the full subinterval is bad.
            require(round(v*(b+F(i,k)))==round(v*(b+F(i+1,k))),'hostile tooth')
            hostiles+=1
        examples.append({'k':k,'r':str(r),'speeds':speeds,'interval':[str(F(b,D)),str(F(b+1,D))],'product':str(min(vs)),'rejected_constant':str(c0)})
    return {'status':'PASS','rational_envelope_checks':count,'peak_iff_checks':peaks,'literal_measure_checks':literal,'hostile_partition_cells':hostiles,'integer_hostiles':examples,'radius_1_over_14_constants':[str(F(1,7))]+[str(F(k,7)*ceil(F(k-1,7-k))) for k in range(2,7)]}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
