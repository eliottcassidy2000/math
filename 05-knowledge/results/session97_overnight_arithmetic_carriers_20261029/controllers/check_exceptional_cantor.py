#!/usr/bin/env python3
from fractions import Fraction as F
from check_integer_termination import require
from check_exceptional_cylinders import WORDS,chart
import json

INITIAL=[1,2,3,4,6,9,13,18,24,31,40]
def characteristic(x):return x**11-x**10-x**5-3

def cells(level,charts):
    # Any seed has the same length-level inverse-address residues because
    # every inverse branch contracts by at least one normalized dyadic level.
    mod=1<<(level+1);s={1}
    inverse=[(1<<c['D'],c['C'],pow(3**c['J'],-1,mod)) for c in charts]
    for _ in range(level):
        s={(power*y-carry)*inv%mod for y in s for power,carry,inv in inverse}
    return s

def main():
    charts=[chart(w) for w in WORDS]
    lower=F(1357705,1000000);upper=F(1357706,1000000)
    require(characteristic(lower)<0<characteristic(upper),'exact lambda bracket')
    require(lower**1000>2**441 and upper**1000<2**442,'exact s bracket')
    for i,n in enumerate(INITIAL):
        require(F(n)<=F(17,8)*lower**i,'inductive upper constant')
        require(F(n)>=upper**i,'inductive lower constant')
    N=INITIAL[:]
    for level in range(11,1001):N.append(N[-1]+N[level-6]+3*N[level-11])
    receipts=[]
    for level in range(25):
        s=cells(level,charts)
        require(len(s)==N[level],'inverse IFS vs count recurrence')
        receipts.append(len(s))
        # The final residue set is genuinely invariant, not seed dependent.
        mod=1<<(level+1)
        nextset={( (1<<c['D'])*y-c['C'])*pow(3**c['J'],-1,mod)%mod for y in s for c in charts}
        require(nextset==s,'finite precision attractor fixedpoint')
    # Completed-address counting differs from occupied depth-level cells.
    completed=[1]+[0]*100
    for n in range(1,101):
        completed[n]=sum(mult*completed[n-d] for d,mult in [(1,1),(6,1),(11,3)] if n>=d)
    numerator=[1,1,1,1,2,3,3,3,3,3,3]
    require(all(N[n]==sum(numerator[i]*completed[n-i] for i in range(min(n,10)+1)) for n in range(101)),'rational generating function')
    print(json.dumps({'lambda_interval':[str(lower),str(upper)],'s_interval':['441/1000','442/1000'],
        'cell_counts_through_24':receipts,'initial_counts':INITIAL,'generating_numerator':numerator,
        'cell_count_level_1000':str(N[1000]),'completed_counts_through_20':completed[:21],
        'integer_count_bound_constant':'17/8'},sort_keys=True))
if __name__=='__main__':main()
