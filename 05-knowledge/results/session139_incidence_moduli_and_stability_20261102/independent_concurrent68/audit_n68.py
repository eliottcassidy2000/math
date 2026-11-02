#!/usr/bin/env python3
"""Independent exact audit of the concurrent-input n >= 68 combination.

All comparisons are active under -O. Negative controls are --threshold 67
and --omit-selected-payment. Uniform proof and source credit are in AUDIT.md.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations, product
import json
from math import factorial
import sys


class CheckFailure(Exception):
    pass


def graph_controls():
    count = 0
    for n in range(1,5):
        edges=list(combinations(range(n),2))
        for values in product((F(0),F(1,8),F(1,4)),repeat=len(edges)):
            weights=dict(zip(edges,values))
            remaining={e:w for e,w in weights.items() if w>0}
            discarded={}
            selected={}
            assignments={}
            while remaining:
                e=max(remaining,key=lambda e:(remaining[e],e))
                w=remaining.pop(e); selected[e]=w
                removed=[f for f in remaining if set(f)&set(e)]
                assignments[e]=sum((remaining[f] for f in removed),F(0))
                for f in removed:
                    discarded[f]=(remaining.pop(f),e)
            W=sum((weights[e]*weights[f] for e,f in combinations(edges,2)
                   if set(e)&set(f)),F(0))
            D=sum((w for w,e in discarded.values()),F(0))
            deg=[sum((w for f,(w,e) in discarded.items() if v in f),F(0))
                 for v in range(n)]
            dd=sum((discarded[e][0]*discarded[f][0]
                    for e,f in combinations(discarded,2) if set(e)&set(f)),F(0))
            sq=sum((w*w for w,e in discarded.values()),F(0))
            paid=sum((w*selected[e] for w,e in discarded.values()),F(0))
            halfdeg=sum((v*v for v in deg),F(0))/2
            eta=F(2,17)
            L=sum(w<=eta for w in selected.values())
            if dd != halfdeg-sq or paid < sq or W < dd+paid or W < 2*D*D/n:
                raise CheckFailure(f'graph_cancellation_control_{count}')
            if W < eta*D-L*eta*eta:
                raise CheckFailure(f'graph_joint_budget_control_{count}')
            count+=1
    return count


def run(n,omit_selected_payment=False):
    checks={}
    def check(name,condition,value=None):
        if not condition: raise CheckFailure(f'{name}: {value}')
        checks[name]=True if value is None else str(value)
    R=n//2
    check('order_domain',n>=8 and R>14)
    check('initial_small_variance_excluded',F(49,200)>F(4,R))
    check('exp4_partial_sum',sum(F(4)**j/factorial(j) for j in range(8))==F(16319,315))
    check('exp4_exceeds51',F(16319,315)>51)
    check('initial_c_below_g1_lower',F(6,25)<F(49,200))
    check('perfect_matching_penalty',F(1,4**R)<F(1,10000))
    # Unit-weight path on three vertices: D=1, W_dd=0, W_sd=1.
    payment=F(0) if omit_selected_payment else F(1)
    check('selected_payment_essential',payment>=F(2,3),payment-F(2,3))
    check('exact_graph_controls',True,graph_controls())
    nk=n*F(4,R)
    left=(F(0),F(-1),F(2)/nk)
    right=((F(2)/nk)*(nk/4)**2-nk/8,-(F(2)/nk)*nk/2,F(2)/nk)
    check('factor_four_square_identity',left==right)
    bootstrap=[]
    for i,(c,b,exp_floor,target) in enumerate([
            (F(6,25),1,F(65,24),F(63,100)),
            (F(63,100),4,F(51),F(49,50)),
            (F(49,50),7,F(110),F(99,100))]):
        m=c*R/4-1-F(1,2*R)
        check(f'bootstrap_mass_{i}',m>b,m-b)
        check(f'bootstrap_mass_slope_{i}',c/4>0)
        slack=1-1/exp_floor-F(1,10000)-target
        check(f'bootstrap_advance_{i}',slack>0,slack)
        bootstrap.append({'from':str(c),'to':str(target),'mass':str(m),
                          'integer_mass_floor':b,'reserve_slack':str(slack)})
    check('e_taylor_floor',sum(F(1,factorial(j)) for j in range(5))==F(65,24))
    check('e7_taylor_floor',sum(F(7)**j/factorial(j) for j in range(5))>110)
    c,eta,delta,a=F(99,100),F(2,17),F(1,21),F(1979,4200)
    heavy=4*(R*(c/4-eta)-1-F(1,2*R))
    check('heavy_count_exceeds13',heavy>13,heavy-13)
    check('heavy_count_uniform_slope',c-4*eta>0,c-4*eta)
    Q=F(R)+F(33,2)+F(8,R)
    rho=F(3,4)
    check('variance_even_envelope',R+16<Q)
    check('variance_odd_envelope',F(2*R+1,2)+2*(2*R+1)*F(4,R)==Q)
    check('tail_interval_nonempty',R-14>=1)
    check('tail_left_endpoint',Q*rho**(R-1)<delta,delta-Q*rho**(R-1))
    check('tail_right_endpoint',Q*rho**14/(R-14)<delta,delta-Q*rho**14/(R-14))
    check('tail_left_uniform_ratio',rho*(1+F(1,R))<1)
    check('tail_right_decreasing_form',Q/(R-14)==1+(F(61,2)+F(8,R))/(R-14))
    check('positive_a',a>0,a)
    check('light_budget_absorption',c-4*eta-delta>a,c-4*eta-delta-a)
    check('strict_joint_coefficient',a-4*eta>0,a-4*eta)
    check('uniform_R_eta',R*eta>=4,R*eta-4)
    check('final_nonpositive_upper_difference',4*a/R-4*eta*eta<=eta*a-4*eta*eta)
    # Multiply the joint-budget identity by 4L and compare coefficients
    # in the two independent variables (D_l,L), not sample values.
    linear={(1,0):F(1),(0,1):-2*eta}
    square={}
    for p,cp in linear.items():
        for q,cq in linear.items():
            exponent=(p[0]+q[0],p[1]+q[1])
            square[exponent]=square.get(exponent,F(0))+cp*cq
    check('joint_budget_square_identity',
          {(2,0):F(1),(1,1):-4*eta,(0,2):4*eta*eta} == square)
    check('closure_eta_at_least_kappa_bound',eta>=F(4,R))
    for order in (n,n+1):
        r=order//2
        k=F(4,r)*(1-F(3,4)**r-(F(1,4)**r if order%2==0 else 0))
        check(f'parity_coefficient_{order}',0<k<=F(4,r))
    return {'status':'PASS','threshold':n,'R_min':R,
            'bootstrap':bootstrap,'active_checks':checks}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--threshold',type=int,default=68)
    p.add_argument('--omit-selected-payment',action='store_true')
    args=p.parse_args()
    try: result=run(args.threshold,args.omit_selected_payment)
    except CheckFailure as e:
        print(json.dumps({'status':'REJECTED','reason':str(e)},sort_keys=True));return 1
    print(json.dumps(result,sort_keys=True,indent=2));return 0


if __name__=='__main__':sys.exit(main())
