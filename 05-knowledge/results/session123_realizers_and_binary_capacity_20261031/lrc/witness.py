"""Exact rational witness extraction from a common-clock bad-measure budget."""
from fractions import Fraction as F
from math import lcm

def need(ok,message):
    if not ok: raise ValueError(message)
def floor(x):return x.numerator//x.denominator
def norm(x):
    f=x-floor(x);return min(f,1-f)
def antiderivative(x,h):
    n=floor(x);return h*n+min(x-n,h)
def bad_measure(lo,hi,v,phi,r):
    return (antiderivative(v*hi+phi+r,2*r)-antiderivative(v*lo+phi+r,2*r))/v

def choose(interval,speeds,r,phases=None):
    lo,hi=map(F,interval);r=F(r);vs=tuple(speeds)
    need(lo<hi and 0<r<F(1,2),'input interval/radius')
    need(vs and all(type(v) is int and v>0 for v in vs),'positive integer speeds')
    ps=tuple(map(F,phases)) if phases is not None else (F(0),)*len(vs)
    need(len(ps)==len(vs),'phase arity')
    stats={'measure_calls':0,'safety_calls':0,'bisections':0,'branch':None}
    def safe(t):
        stats['safety_calls']+=len(vs)
        return all(norm(v*t+p)>=r for v,p in zip(vs,ps))
    def measure(l,u,i):
        stats['measure_calls']+=1
        return bad_measure(l,u,vs[i],ps[i],r)
    W=hi-lo;mus=[measure(lo,hi,i) for i in range(len(vs))]
    need(all(x<W for x in mus),'each bad set must be proper in measure')
    delta=W-sum(mus);need(delta>=0,'bad-measure budget not certified')
    Q=lcm(r.denominator,*(p.denominator for p in ps));V=max(vs);eps=F(1,Q*Q*V*V)
    depth=0;temp=W
    while temp>=eps:temp/=2;depth+=1
    stats['depth_bound']=depth;stats['event_separation']=str(eps)
    def result(t,branch=None):
        need(lo<=t<=hi and safe(t),'actual common-time witness')
        if branch is not None:stats['branch']=branch
        return t,stats
    if safe(lo):return result(lo,'left_endpoint')
    overlap=None
    if delta==0:
        i=next(i for i,(v,p) in enumerate(zip(vs,ps)) if norm(v*lo+p)<r)
        v,p=vs[i],ps[i];n=floor(v*lo+p+F(1,2));li=(n-r-p)/v;e=(n+r-p)/v
        need(lo<e<hi,'proper first bad component')
        if safe(e):return result(e,'equality_boundary')
        j=next(j for j,(v,p) in enumerate(zip(vs,ps)) if norm(v*e+p)<r)
        need(j!=i,'different bad component')
        v,p=vs[j],ps[j];n=floor(v*e+p+F(1,2));lj=(n-r-p)/v;hj=(n+r-p)/v
        overlap=(max(lo,li,lj),min(hi,e,hj));need(overlap[0]<overlap[1],'strict overlap')
        stats['branch']='overlap_repair';stats['overlap']=list(map(str,overlap))
    else:stats['branch']='strict_budget'
    def integral(l,u):
        val=u-l-sum(measure(l,u,i) for i in range(len(vs)))
        if overlap:val+=max(F(0),min(u,overlap[1])-max(l,overlap[0]))
        return val
    a,b=lo,hi;total=integral(a,b);need(total>0,'positive initial corrected density')
    while True:
        if safe(a):return result(a)
        if safe(b):return result(b)
        need(b-a>=eps,'boundary-separation termination contradicted')
        mid=(a+b)/2;left=integral(a,mid)
        if left>0:b=mid;total=left
        else:a=mid;total-=left
        need(total>0,'positive density invariant')
        stats['bisections']+=1
        need(stats['bisections']<=depth,'declared logarithmic depth bound')

def conductor(interval,speeds,r,phases=None):
    r=F(r);k=len(speeds);need(k>=2 and 0<2*k*r<1,'conductor domain')
    x=2*r*(k-1)/(1-2*k*r);N=-(-x.numerator//x.denominator);c=2*k*r*N
    W=F(interval[1])-F(interval[0]);need(min(speeds)*W>=c,'conductor clock cutoff')
    return choose(interval,speeds,r,phases)
