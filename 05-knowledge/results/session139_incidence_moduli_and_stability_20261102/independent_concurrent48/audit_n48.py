#!/usr/bin/env python3
"""Independent exact audit of the pinned n>=48 two-regime certificate.

Continuum verification uses polynomial identities, not interpolation or a
grid. Explicit checks survive -O; negative controls: --threshold 47,
--bad-dual-ratio (the scalar ratio 24/5 hostile, not a tournament).
"""
import argparse
from fractions import Fraction as F
from itertools import combinations,product
from math import factorial
import json
import sys


class CheckFailure(Exception):pass


def add(*polys):
    z={}
    for p in polys:
        for k,v in p.items():z[k]=z.get(k,F(0))+v
    return {k:v for k,v in z.items() if v}


def scale(p,a):return {k:a*v for k,v in p.items() if a*v}


def mul(p,q):
    z={}
    for e,x in p.items():
        for f,y in q.items():
            k=(e[0]+f[0],e[1]+f[1]);z[k]=z.get(k,F(0))+x*y
    return {k:v for k,v in z.items() if v}


def graph_controls(t,d,lam):
    count=0
    for n in range(1,5):
        edges=list(combinations(range(n),2))
        for values in product((F(0),F(1,6),F(1,4)),repeat=len(edges)):
            q=dict(zip(edges,values));remaining={e for e in edges if q[e]};groups=[]
            while remaining:
                e=max(remaining,key=lambda e:(q[e],e))
                members={f for f in remaining if set(e)&set(f)}
                remaining-=members
                w=q[e];S=sum((q[f] for f in members if f!=e),F(0))
                energy=sum((q[f]*q[g] for f,g in combinations(members,2) if set(f)&set(g)),F(0))
                G=max(w*S,w*S/2+S*S/4)
                if energy<G:raise CheckFailure(f'graph_max_invoice_{count}')
                if w<=t and G<t*S-lam*(d-4*w):raise CheckFailure(f'graph_light_dual_{count}')
                groups.append((w,S,energy,G))
            W=sum((q[e]*q[f] for e,f in combinations(edges,2) if set(e)&set(f)),F(0))
            D=sum((S for w,S,E,G in groups),F(0))
            L=sum(w<=t for w,S,E,G in groups)
            U=sum((w for w,S,E,G in groups if w<=t),F(0))
            if W<sum((G for w,S,E,G in groups),F(0)) or W<t*D-lam*(d*L-4*U):
                raise CheckFailure(f'graph_disjoint_dual_{count}')
            count+=1
    return count


def run(n,bad_ratio=False):
    out={}
    def check(name,condition,value=None):
        if not condition:raise CheckFailure(f'{name}: {value}')
        out[name]=True if value is None else str(value)
    R=n//2
    check('domain',n>=8 and R>12)
    check('g1_floor',sum(F(4)**j/factorial(j) for j in range(8))>50)
    check('small_variance_excluded',F(49,200)>F(4,R))
    penalty=F(1,10000)
    check('coarse_perfect_penalty',F(1,4**R)<penalty)
    levels=[F(49,200),F(7,20),F(13,20),F(47,50),F(99,100)]
    boot=[]
    for i,(c,target,k) in enumerate(zip(levels,levels[1:],(2,3,6,11))):
        m=c*R/4-1-F(1,2*R)
        T=sum(m**j/factorial(j) for j in range(k+1))
        check(f'bootstrap_mass_{i}',m>0,m)
        check(f'bootstrap_uniform_slope_{i}',c/4>0)
        check(f'bootstrap_remaining_probability_{i}',1-target-penalty>0)
        slack=(1-target-penalty)*T-1
        check(f'bootstrap_advance_{i}',slack>0,slack)
        boot.append({'from':str(c),'to':str(target),'mass':str(m),
                     'degree':k,'coarse_slack':str(slack)})
    c,t=levels[-1],F(1,6)
    ratio=F(24,5) if bad_ratio else F(5)
    d=ratio*t;lam=t*t/d;beta=4*t/d;h=t-4*lam
    check('slot_positive_divisor',1-4*t>0,1-4*t)
    g=(c-4*t)*R-4-F(2,R)
    check('slot_forces_twelve',g>11*(1-4*t),g-11*(1-4*t))
    check('slot_uniform_slope',c-4*t>0,c-4*t)
    Q=F(R)+F(33,2)+F(8,R)
    check('variance_even_envelope',R+16<Q)
    check('variance_odd_envelope',F(2*R+1,2)+2*(2*R+1)*F(4,R)==Q)
    check('tail_left',Q*F(3,4)**(R-1)<c-d,c-d-Q*F(3,4)**(R-1))
    check('tail_right',Q*F(3,4)**12/(R-12)<c-d,c-d-Q*F(3,4)**12/(R-12))
    check('tail_left_uniform_ratio',F(3,4)*(1+F(1,R))<1)
    check('tail_right_uniform_form',Q/(R-12)==1+(F(57,2)+F(8,R))/(R-12))
    check('light_deficit_positive',d-4*t>0,d-4*t)
    # u,w/t is the first variable; v,S/t the second.
    one={(0,0):F(1)};u={(1,0):F(1)};v={(0,1):F(1)}
    linear=add(mul(u,v),scale(v,-1),one,scale(u,-beta))
    quad=add(scale(mul(u,v),F(1,2)),scale(mul(v,v),F(1,4)),scale(v,-1),one,scale(u,-beta))
    p1=add(scale(u,1-beta),scale(mul(u,u),-F(1,4)))
    p2=add(scale(mul(u,u),2),scale(u,-(2+beta)),one)
    check('linear_boundary_identity',add(linear,scale(p2,-1))==mul(add(one,scale(u,-1)),add(scale(u,2),scale(v,-1))))
    qvertex=add(v,scale(one,-2),u)
    check('quadratic_vertex_identity',add(quad,scale(p1,-1))==scale(mul(qvertex,qvertex),F(1,4)))
    boundarygap=add(scale(u,3),scale(one,-2))
    check('regime_minimum_gap_identity',add(p2,scale(p1,-1))==scale(mul(boundarygap,boundarygap),F(1,4)))
    check('quadratic_boundary_identity',add(quad,scale(p2,-1))==mul(add(v,scale(u,-2)),add(scale(v,F(1,4)),u,scale(one,-1))))
    center=(2+beta)/4
    remainder=1-(2+beta)**2/8
    squarebase=add(u,scale(one,-center))
    check('second_minimum_square_identity',p2==add(scale(mul(squarebase,squarebase),2),scale(one,remainder)))
    check('first_minimum_nonnegative',F(5,6)-beta>=0,F(5,6)-beta)
    check('second_minimum_nonnegative',remainder>=0,remainder)
    # Formal hostile in the two-invoice relaxation, never labelled a tournament.
    hu,hv=F(17,24),F(17,12)
    hostile=max(hu*hv,hu*hv/2+hv*hv/4)-hv+1-F(5,6)*hu
    check('scalar_ratio24over5_hostile',hostile==F(-1,288),hostile)
    check('dual_total_positive',h>0,h)
    check('uniform_Rt',R*t>=4,R*t-4)
    check('closing_coefficient',F(4,R)-4*lam<=h)
    check('no_light_closure',t>=F(4,R))
    check('exact_graph_controls',True,graph_controls(t,d,lam))
    check('slot47_negative', (c-4*t)*23-4-F(2,23)<11*(1-4*t))
    check('Rt47_negative',23*t<4)
    for order in (n,n+1):
        r=order//2;k=F(4,r)*(1-F(3,4)**r-(F(1,4)**r if order%2==0 else 0))
        check(f'parity_coefficient_{order}',0<k<=F(4,r))
    return {'status':'PASS','threshold':n,'R_min':R,'coarse_penalty':str(penalty),
            'bootstrap':boot,'constants':{'t':str(t),'d':str(d),'lambda':str(lam),'h':str(h)},'active_checks':out}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--threshold',type=int,default=48)
    p.add_argument('--bad-dual-ratio',action='store_true')
    args=p.parse_args()
    try:result=run(args.threshold,args.bad_dual_ratio)
    except CheckFailure as e:
        print(json.dumps({'status':'REJECTED','reason':str(e)},sort_keys=True));return 1
    print(json.dumps(result,sort_keys=True,indent=2));return 0


if __name__=='__main__':sys.exit(main())
