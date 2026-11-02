#!/usr/bin/env python3
"""Independent exact endpoint audit for the pinned concurrent n>=56 theorem.

Uses Taylor degrees 4,5,5 and the weaker penalty 1/10000, independently of
the upstream degree-12 certificate. All checks stay active under -O.
"""
import argparse
from fractions import Fraction as F
from math import factorial
import json
import sys


class CheckFailure(Exception):
    pass


def run(n,bad_dual=False):
    out={}
    def check(name,condition,value=None):
        if not condition:raise CheckFailure(f'{name}: {value}')
        out[name]=True if value is None else str(value)
    R=n//2
    check('domain',n>=8 and R>8)
    check('g1_lower',sum(F(4)**j/factorial(j) for j in range(8))>50)
    check('small_variance_excluded',F(49,200)>F(4,R))
    penalty=F(1,10000)
    check('coarse_penalty',F(1,4**R)<penalty)
    for order in (2*R,2*R+1):
        for r in (1,R):
            correction=F(1,4**r) if order==2*r else F(0)
            check(f'perfect_penalty_n{order}_r{r}',correction<=F(1,4**R))
    levels=[F(49,200),F(1,2),F(91,100),F(99,100)]
    advances=[]
    for i,(c,target,degree) in enumerate(zip(levels,levels[1:],(4,5,5))):
        m=c*R/4-1-F(1,2*R)
        T=sum(m**j/factorial(j) for j in range(degree+1))
        check(f'bootstrap_mass_{i}',m>0,m)
        check(f'bootstrap_uniform_mass_slope_{i}',c/4>0)
        check(f'bootstrap_remaining_probability_{i}',1-target-penalty>0)
        slack=T*(1-target-penalty)-1
        check(f'bootstrap_advance_{i}',slack>0,slack)
        advances.append({'from':str(c),'to':str(target),'mass':str(m),
                         'Taylor_degree':degree,'coarse_slack':str(slack)})
    c,t=levels[-1],F(1,7)
    d=F(3,4) if bad_dual else F(381,500)
    lam=t*t/d
    heavy=(c-4*t)*R-4-F(2,R)
    check('heavy_integer_floor',heavy>7,heavy-7)
    check('heavy_uniform_slope',c-4*t>0,c-4*t)
    Q=F(R)+F(33,2)+F(8,R)
    check('variance_even_envelope',R+16<Q)
    check('variance_odd_envelope',F(2*R+1,2)+2*(2*R+1)*F(4,R)==Q)
    check('tail_interval',R-8>=1)
    check('left_tail',Q*F(3,4)**(R-1)<c-d,c-d-Q*F(3,4)**(R-1))
    check('right_tail',Q*F(3,4)**8/(R-8)<c-d,c-d-Q*F(3,4)**8/(R-8))
    check('left_tail_uniform_ratio',F(3,4)*(1+F(1,R))<1)
    check('right_tail_uniform_form',Q/(R-8)==1+(F(49,2)+F(8,R))/(R-8))
    check('light_deficit_positive',d-4*t>0,d-4*t)
    check('dual_constant_identity',lam*d==t*t)
    # S,w exponents; multiply out the square rather than sampling S,w.
    linear={(1,0):F(1,2),(0,0):-t,(0,1):F(1,2)}
    rhs={}
    for e,xe in linear.items():
        for f,xf in linear.items():
            key=(e[0]+f[0],e[1]+f[1])
            rhs[key]=rhs.get(key,F(0))+xe*xf
    rhs[(0,1)]+=t-4*lam
    rhs[(0,2)]-=F(1,4)
    rhs={e:x for e,x in rhs.items() if x}
    lhs={(2,0):F(1,4),(1,1):F(1,2),(1,0):-t,(0,0):lam*d,(0,1):-4*lam}
    check('dual_polynomial_identity',lhs==rhs)
    check('dual_remainder_uniform_nonnegative',3*t/4-4*lam>=0,3*t/4-4*lam)
    check('dual_at_square_zero_boundary',t*(3*t/4-4*lam)>=0,t*(3*t/4-4*lam))
    check('total_dual_positive',t-4*lam>0,t-4*lam)
    check('uniform_Rt_endpoint',R*t>=4,R*t-4)
    check('closing_upper_bound',F(4,R)-4*lam<=t-4*lam)
    check('no_light_closure',t>=F(4,R))
    for order in (n,n+1):
        r=order//2
        k=F(4,r)*(1-F(3,4)**r-(F(1,4)**r if order%2==0 else 0))
        check(f'parity_coefficient_{order}',0<k<=F(4,r))
    return {'status':'PASS','threshold':n,'R_min':R,
            'coarse_penalty':str(penalty),'bootstrap':advances,
            'constants':{'c':str(c),'t':str(t),'d':str(d),'lambda':str(lam)},
            'active_checks':out}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--threshold',type=int,default=56)
    p.add_argument('--bad-dual',action='store_true')
    args=p.parse_args()
    try:result=run(args.threshold,args.bad_dual)
    except CheckFailure as e:
        print(json.dumps({'status':'REJECTED','reason':str(e)},sort_keys=True));return 1
    print(json.dumps(result,sort_keys=True,indent=2));return 0


if __name__=='__main__':sys.exit(main())
