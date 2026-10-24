#!/usr/bin/env python3
"""Independent decoder. Reads ONLY {n, moments}; no atlas, support cap, or count."""
import json,sys
from fractions import Fraction as F


def require(ok,message):
    if not ok: raise ValueError(message)


def add(a,b):
    return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]


def mul(a,b):
    ans=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): ans[i+j]+=x*y
    return ans


def solve(A,b):
    """Exact Gauss-Jordan elimination; no producer code or numerical solver."""
    d=len(b); rows=[[F(x) for x in row]+[F(y)] for row,y in zip(A,b)]
    for col in range(d):
        piv=next((r for r in range(col,d) if rows[r][col]),None)
        require(piv is not None,'singular moment matrix')
        rows[col],rows[piv]=rows[piv],rows[col]
        lead=rows[col][col];rows[col]=[v/lead for v in rows[col]]
        for r in range(d):
            if r!=col:
                factor=rows[r][col]
                rows[r]=[v-factor*w for v,w in zip(rows[r],rows[col])]
    return [row[-1] for row in rows]


def expectation(poly,moments):
    require(len(poly)<=len(moments),'insufficient supplied moments')
    return sum((F(a)*b for a,b in zip(poly,moments)),F(0))


def congruence(moments):
    k=len(moments);Q=[1]
    for i in range(k): Q=mul(Q,[-2*i-1,1])
    exponent=2*k-k.bit_count(); modulus=1<<exponent
    residue=(-sum(Q[j+1]*moments[j] for j in range(k))*pow(Q[0],-1,modulus))%modulus
    return Q,exponent,modulus,residue


def polynomial_bounds(moments):
    """Bounds global on h>=1, with exact coefficient identity certificates."""
    require(len(moments)>=8,'eight moments required by this bound certificate')
    b=moments[:4]
    C=[[moments[i+j+1] for j in range(4)] for i in range(4)]
    lower_q=solve(C,b)
    lower_p=add([2*c for c in lower_q],[0]+[-c for c in mul(lower_q,lower_q)])
    # 1 - x P_L = (1-x q_L)^2, coefficientwise, not finite sampling.
    error=add([1],[0]+[-c for c in lower_q])
    require(add([1],[0]+[-c for c in lower_p])==mul(error,error),'lower square identity')
    lower=expectation(lower_p,moments)
    require(lower==sum(q*x for q,x in zip(lower_q,b)),'lower normal-equation value')
    b=[moments[i+1]-moments[i] for i in range(3)]
    C=[[moments[i+j+2]-moments[i+j+1] for j in range(3)] for i in range(3)]
    upper_q=solve(C,b)
    upper_p=add(add([1],mul([2,-2],upper_q)),mul([0,-1,1],mul(upper_q,upper_q)))
    # x P_U - 1 = (x-1)(1-x q_U)^2.
    error=add([1],[0]+[-c for c in upper_q])
    require(add([-1],[0]+upper_p)==mul([-1,1],mul(error,error)),'upper square identity')
    upper=expectation(upper_p,moments)
    require(upper==moments[0]-sum(q*x for q,x in zip(upper_q,b)),'upper normal-equation value')
    require(lower<=upper,'inconsistent bounds')
    return lower,upper,{'lower_q':[str(x) for x in lower_q],'lower_polynomial':[str(x) for x in lower_p],'upper_q':[str(x) for x in upper_q],'upper_polynomial':[str(x) for x in upper_p],'square_identities_verified':True}


def decode(payload):
    require(set(payload)=={'n','moments'},'decoder accepts only n and moments')
    n=payload['n'];moments=payload['moments']
    require(type(n) is int and n>=1,'n must be positive')
    require(type(moments) is list and moments and all(type(v) is int and v>0 for v in moments),'positive integer moments required')
    m=(n-1)*(n-2)//2;N=1<<m
    require(moments[0]==N,'K0 disagrees with fixed-path word count')
    Q,e,M,residue=congruence(moments)
    lower,upper=F(1),F(N);certificate=None
    if len(moments)>=8:
        lower,upper,certificate=polynomial_bounds(moments)
        lower=max(F(1),lower);upper=min(F(N),upper)
    low=-(-lower.numerator//lower.denominator);high=upper.numerator//upper.denominator
    first=low+(residue-low)%M
    count=max(0,1+(high-first)//M)
    return {'n':n,'moments_used':len(moments),'Q_coefficients':Q,'exponent':e,'modulus':M,'residue':residue,'lower':str(lower),'upper':str(upper),'integer_interval':[low,high],'first_candidate':first if count else None,'candidate_count':count,'prediction':first if count==1 else None,'certificate':certificate}

if __name__=='__main__':
    if len(sys.argv)!=2: raise SystemExit('usage: reconstruct.py INPUT.json')
    print(json.dumps(decode(json.load(open(sys.argv[1]))),indent=2,sort_keys=True))
