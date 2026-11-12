"""Exact rational backbone resummation and first marked-edge transfer."""
from fractions import Fraction as F
from math import gcd, factorial, comb
from collections import Counter
import sys,json,time
from pathlib import Path
OLD=Path(__file__).resolve().parents[2]/'session278_signed_corrections_20261112'/'enumeration'
sys.path.insert(0,str(OLD))
from gas import parts,elementary

def validate(n,backbone):
    if not isinstance(n,int) or n<0:raise ValueError("n must be a nonnegative integer")
    backbone=tuple(sorted(backbone))
    if len(set(backbone))!=len(backbone) or any(not isinstance(a,int) or a<3 or a%2==0 for a in backbone):
        raise ValueError("backbone must contain distinct odd integers at least three")
    return backbone

def configs(n,backbone):
    def rec(i,left,ms):
        if i==len(backbone):
            yield tuple(ms);return
        for j in range(left//backbone[i]+1):
            yield from rec(i+1,left-j*backbone[i],ms+[j])
    yield from rec(0,n,[])

def bweight(backbone,ms):
    e=sum(m*(a-1)//2+a*comb(m,2) for a,m in zip(backbone,ms))
    e+=sum(ms[i]*ms[j]*gcd(a,b) for i,a in enumerate(backbone) for j,b in enumerate(backbone) if i<j)
    z=1
    for a,m in zip(backbone,ms):z*=a**m*factorial(m)
    return F(1<<e,z)

def envelope(n,backbone=(3,),marked=False):
    """Compute U and optionally C1 without enumerating unrestricted partitions.

    Each conditioned backbone is followed by a size DP. The moments track
    multiplicity of each previous off-backbone length, not full partitions.
    """
    backbone=validate(n,backbone)
    high=[a for a in range(3,n+1,2) if a not in backbone]
    U=F(0);C=F(0);stats={'backbone_configs':0,'transitions':0}
    for ms in configs(n,backbone):
        stats['backbone_configs']+=1
        used=sum(a*m for a,m in zip(backbone,ms)); N=n-used;r=sum(ms)
        v=[F(1<<(comb(f,2)+f*r),factorial(f)) for f in range(N+1)]
        c=[F(0)]*(N+1);mom=[]
        previous=[]
        for b in high:
            if b>N:break
            R=sum(m*gcd(a,b) for a,m in zip(backbone,ms))
            out=[F(0)]*(N+1); oc=[F(0)]*(N+1)
            om=[[F(0)]*(N+1) for _ in range(len(previous)+1)] if marked else []
            qs=[1-F(1,1<<(min(a,b)-gcd(a,b))) for a in previous]
            for m in range(N+1):
                if not v[m]:continue
                Q=sum((q*arr[m] for q,arr in zip(qs,mom)),F(0)) if marked else 0
                for j in range((N-m)//b+1):
                    stats['transitions']+=1
                    t=m+b*j
                    w=F(1<<(j*m+j*(b-1)//2+b*comb(j,2)+j*R),b**j*factorial(j))
                    out[t]+=w*v[m]
                    if marked:
                        oc[t]+=w*(c[m]+j*Q)
                        for z,arr in enumerate(mom):om[z][t]+=w*arr[m]
                        om[-1][t]+=w*j*v[m]
            v,c,mom=out,oc,om
            previous.append(b)
        w=bweight(backbone,ms)
        U+=w*v[N];C+=w*c[N]
    return U,C,stats

def sector(p,backbone):
    e=sum((a-1)//2 for a in p);qs=[]
    for i,a in enumerate(p):
        for b in p[i+1:]:
            base=gcd(a,b) if a in backbone or b in backbone else min(a,b)
            e+=base;d=base-gcd(a,b)
            if d:qs.append(1-F(1,1<<d))
    z=1
    for a,m in Counter(p).items():z*=a**m*factorial(m)
    return F(1<<e,z),qs

def signed(n,backbone,k=3):
    backbone=validate(n,backbone)
    if not isinstance(k,int) or k<0:raise ValueError("k must be nonnegative integer")
    cs=[F(0)]*(k+2); low=[F(0)]*(k+1);exact=F(0);pcount=0
    for p in parts(n):
        pcount+=1;w,qs=sector(p,backbone);es=elementary(qs,k+1)
        for j in range(k+2):cs[j]+=w*es[j]
        for j in range(k+1):
            if len(qs)>j:low[j]+=w*F(j+1,len(qs))*es[j+1]
        exact+=w*__import__('functools').reduce(lambda x,y:x*(1-y),qs,F(1))
    ss=[sum(((-1)**j*cs[j] for j in range(t+1)),F(0)) for t in range(k+1)]
    lo=max(ss[t]+low[t] if t%2 else ss[t]-cs[t+1] for t in range(k+1))
    hi=min(ss[t]+cs[t+1] if t%2 else ss[t]-low[t] for t in range(k+1))
    ilo=-(-lo.numerator//lo.denominator);ihi=hi.numerator//hi.denominator
    assert ilo<=exact<=ihi and exact.denominator==1
    return {'backbone':list(backbone),'n':n,'partitions':pcount,'residue_bits':(ihi-ilo).bit_length(),'low':str(ilo),'high':str(ihi),'exact':str(exact),'C0':str(cs[0]),'C1':str(cs[1]),'errors':[str(x-exact) for x in ss]}

if __name__=='__main__':
    st=time.time();controls=[]
    for bad in [(3,3),(2,),(1,),(-3,)]:
        try:envelope(6,bad)
        except ValueError:pass
        else:raise AssertionError('bad backbone accepted')
    for n in range(12):
        u,c,_=envelope(n,tuple(range(3,n+1,2)),True)
        assert str(u)==signed(n,())['exact'] and c==0
        assert envelope(n,(101,),True)[:2]==envelope(n,(),True)[:2]
    for n in range(21):
        for B in [(),(3,),(3,5),(3,7)]:
            u,c,stats=envelope(n,B,True);s=signed(n,B)
            assert str(u)==s['C0'] and str(c)==s['C1']
            controls.append({'n':n,'backbone':B,**stats})
    rows=[]
    for n in [24,32,40,48,56,64]:
        for B in [(3,),(3,5),(3,5,7),(3,5,7,9),(3,5,7,9,11)]:
            row=signed(n,B);rows.append(row)
    print(json.dumps({'status':'PASS','seconds':time.time()-st,'dp_controls':controls,'rows':rows},indent=2))
