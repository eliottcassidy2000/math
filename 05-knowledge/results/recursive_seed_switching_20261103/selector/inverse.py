"""Exact binary-section inverse; no word-length cap or interval enumeration.

The parameter state is unbounded. Monotone sections are searched by exact
integer bisection, not assumed globally monotone.
"""
from math import gcd

def need(ok,message):
    if not ok:raise RuntimeError(message)

def excess(r):
    """ordinary one-apex Hamilton count minus1 for decode_r(r), including r0."""
    need(r>=0,'Nonnegative binary mask')
    total=0
    while r:
        u,b=divmod(r,2)
        if not b:total+=u+(u&1)
        r=u
    return total

def evaluate(K,q,d,r):return K*r+q*excess(r)+d*(r&1)

def sections(K,q,d,D):
    """All nonnegative integer solutions Kr+q*excess(r)+d*parity(r)=D."""
    need(K>0 and q>0 and 0<=d<=q,'Positive section coefficients')
    stats=dict(states=0,bisections=0,evaluations=0,negative_prunes=0,gcd_prunes=0,max_depth=0)
    def solve(k,delta,residual,depth):
        stats['states']+=1;stats['max_depth']=max(stats['max_depth'],depth)
        if residual<0:stats['negative_prunes']+=1;return []
        if residual==0:return [0]
        common=gcd(gcd(k,q),delta)
        if residual%common:stats['gcd_prunes']+=1;return []
        upper=residual//k
        if upper==0:return []
        if k>q*(upper-1)+delta:
            stats['bisections']+=1;lo=0;hi=upper
            while lo<=hi:
                mid=(lo+hi)//2;value=evaluate(k,q,delta,mid);stats['evaluations']+=1
                if value==residual:return [mid]
                if value<residual:lo=mid+1
                else:hi=mid-1
            return []
        left=solve(2*k+q,q,residual,depth+1)
        right=solve(2*k,0,residual-k-delta,depth+1)
        return [2*u for u in left]+[2*u+1 for u in right]
    answer=sorted(solve(K,d,D,0))
    need(all(evaluate(K,q,d,r)==D for r in answer),'Exact section acceptance')
    return answer,stats

def decode(r):
    need(r>0,'Positive nonempty reader mask')
    return ''.join('0' if r>>i&1 else '1' for i in range(r.bit_length()))

def inverse(seed,target,minimum_clock=1,include_empty=True):
    h,e,t,F,p,q,*_=seed
    need(0<=e<=q and t>=q and p>0 and F>0,'Actual seed inequalities')
    hits=[];clocks=[];M=1<<minimum_clock;s=minimum_clock
    while h+M*F<=target:
        D=target-h-M*F;K=t-q+M*p;d=q-e
        roots,stats=sections(K,q,d,D)
        for r in roots:
            if r==0:
                if include_empty:hits.append(dict(word='',prefix=s,r=0))
            else:hits.append(dict(word=decode(r),prefix=s,r=r))
        clocks.append(dict(prefix=s,K=K,d=d,D=D,**stats))
        M*=2;s+=1
    return hits,clocks
