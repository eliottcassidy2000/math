"""Exact five-adic cyclic-subgroup depth and periodic content rates.
All verification gates use explicit exceptions, including under python -O.
"""
from fractions import Fraction as F
from itertools import product
from math import gcd
import json
import random
INF = float('inf')
def require(x, msg):
    if not x: raise ValueError(msg)
def val(x):
    if x == 0: return INF
    x=abs(x); n=0
    while x%5==0: x//=5; n+=1
    return n

def coeff(w):
    P=Q=1; B=0
    for a in w: P,Q,B=3*P,Q*2**a,3*B+Q
    return P,Q,B

def compose(a,b): # chronological a then b
    P,Q,B=a; R,S,C=b
    return P*R,Q*S,R*B+C*Q

def params(W):
    P,Q,B=W; m=P*pow(Q,-1,5)%5
    t=next(t for t in (1,2,4) if pow(m,t,5)==1)
    if t==1: r=min(val(P-Q),val(B))
    else: r=val(P**t-Q**t)
    return t,r

def cutoff(W,H):
    P,Q,B=W; R,S,C=H
    t,r=params(W)
    c=val((P-Q)*C-B*(R-S))
    if t>1:
        x=R*pow(S,-1,5)%5
        if pow(x,t,5)!=1: return 0
        beta=val(R**t-S**t)
        return min(c,beta if beta<r else INF)
    s,h=val(P-Q),val(B)
    alpha=val(R-S) if s<=h else val(C)
    return max(0,min(c-r,alpha if alpha<r else INF))

def mass(t,r,H):
    if H==INF: return F(4*r+1,4*t)
    H=int(H)
    if H<=r: return F(H,t)
    return F(r,t)+F(1,4*t)*(1-F(1,5**(H-r)))

def phase_rate(w,u,v):
    W=coeff(w); H=coeff(v+u)
    t,r=params(W); depth=cutoff(W,H)
    return mass(t,r,depth),depth

def rate(w,u=()):
    return sum((phase_rate(w,u,w[:j])[0] for j in range(1,len(w)+1)),F(0))

def pair(C,q):
    P,Q,B=C; inv=pow(Q,-1,q)
    return P*inv%q,B*inv%q

def orbit(W,q):
    m,b=pair(W,q); x,y=1,0; out=set()
    while (x,y) not in out:
        out.add((x,y)); x,y=m*x%q,(m*y+b)%q
    require((x,y)==(1,0),'orbit failed to close')
    return out

def words(budget,prefix=()):
    if prefix: yield prefix
    for a in range(1,budget+1): yield from words(budget-a,prefix+(a,))

def primitive_root(w):
    return next(w[:d] for d in range(1,len(w)+1) if len(w)%d==0 and w[:d]*(len(w)//d)==w)

def run():
    ws=list(words(7)); hs=[()]+list(words(5)); gates=0
    for w in ws:
        W=coeff(w); t,r=params(W)
        for e in range(1,4):
            q=5**e; orb=orbit(W,q)
            require(len(orb)==t*5**max(0,e-r),('order',w,e))
            for hword in hs:
                H=coeff(hword); depth=cutoff(W,H)
                require((pair(H,q) in orb)==(e<=depth),('membership',w,hword,e,depth))
                gates+=1
    rng=random.Random(114)
    affine_gates=0
    cases=[(26,1,25),(126,1,125),(6,1,25),(26,1,5),(126,1,5),(6,1,1),(26,1,0)]
    cases += [(rng.choice([m for m in range(2,250) if m%5]),1,rng.randrange(250)) for _ in range(60)]
    for W in cases:
        t,r=params(W)
        targets=[(1,1,0),W]+[(rng.choice([m for m in range(1,250) if m%5]),1,rng.randrange(250)) for _ in range(50)]
        for e in (1,2,3,4):
            q=5**e; orb=orbit(W,q)
            require(len(orb)==t*5**max(0,e-r),('synthetic order',W,e))
            for H in targets:
                require((pair(H,q) in orb)==(e<=cutoff(W,H)),('synthetic membership',W,H,e))
                affine_gates+=1
    # Exact return averages at finite precision use a different direct state walk.
    finite_checks=0
    for w in list(words(5)):
        W=coeff(w)
        for u in [(),(1,),(2,),(1,2)]:
            for E in (1,2,3):
                q=5**E; T=len(orbit(W,q)); state=pair(coeff(u),q); total=0
                for k in range(T):
                    for a in w:
                        inv=pow(2**a,-1,q); state=(3*state[0]*inv%q,(3*state[1]+1)*inv%q)
                        total+=min(val(state[0]-1),val(state[1]),E)
                expected=sum((mass(*params(W),min(E,phase_rate(w,u,w[:j])[1])) for j in range(1,len(w)+1)),F(0))
                require(F(total,T)==expected,('mean',w,u,E,total,T,expected))
                finite_checks+=1
    all_rates=[]; maximal=[]; accepted=0; unaccepted=[]
    for w in words(12):
        rho=rate(w); ratio=rho/sum(w)
        require(ratio<=F(5,8),('short-controller excess',w,rho))
        all_rates.append(ratio)
        root=primitive_root(w)
        if params(coeff(root))[1]<=1: accepted+=1
        else: unaccepted.append(w)
        if ratio==F(5,8): maximal.append(w)
    # Nontrivial marked-entry examples, including an erased return rate.
    examples=[]
    for w,u in [((1,),()),((1,2),()),((1,2),(1,)),((1,2),(2,)),((3,1,3),()),((5,)*2,()),((1,)*2,()),((1,2,1,2),())]:
        examples.append({'controller':w,'entry':u,'rate_per_period':str(rate(w,u)),'rate_per_clock':str(rate(w,u)/sum(w)), 'phase_depths':[str(phase_rate(w,u,w[:j])[1]) for j in range(1,len(w)+1)]})
    hostile=(1,)*10+(2,)
    low=sum((mass(*params(coeff(hostile)),min(1,phase_rate(hostile,(),hostile[:j])[1])) for j in range(1,len(hostile)+1)),F(0))
    require(params(coeff(hostile))==(4,2) and low==F(3,2) and rate(hostile)==F(33,16), 'order-stall hostile')
    require(rate(hostile)>F(5,4)*low, 'false universal tail factor')
    return {'membership_gates' :gates,'synthetic_affine_membership_gates':affine_gates,'finite_mean_checks':finite_checks,'controllers_clock_at_most_12':len(all_rates),'root_order_growth_accepted':accepted,'root_order_growth_uncovered':len(unaccepted),'first_uncovered':unaccepted[:5],'maximum_clock_rate':str(max(all_rates)),'equality_controllers':len(maximal),'equality_only_all_ones':all(set(w)=={1} for w in maximal),'examples':examples}

if __name__=='__main__': print(json.dumps(run(),indent=2))
