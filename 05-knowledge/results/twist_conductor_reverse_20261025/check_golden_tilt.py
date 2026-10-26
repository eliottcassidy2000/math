#!/usr/bin/env python3
"""Exact Q(sqrt(5)) Doob-transform controls, with q^2=1-q."""
from fractions import Fraction as F
import json

def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
def inv(x):
    a,b=x;n=a*a-a*b-b*b
    return ((a-b)/n,-b/n)
def power(x,n):
    if n<0:return power(inv(x),-n)
    out=(F(1),F(0))
    while n:
        if n&1:out=mul(out,x)
        x=mul(x,x);n//=2
    return out
def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    zero=(F(0),F(0));one=(F(1),F(0));q=(F(0),F(1))
    r=mul(power(q,2),(F(1,2),F(0)));one_minus_r=add(one,(-r[0],-r[1]))
    need(mul(power(q,-1),(F(1,2),F(0)))==one_minus_r,'geometric leading probability')
    need(inv(one_minus_r)==(F(0),F(2)),'mean letter 2q')
    need(add(power(q,2),q)==one,'golden field identity')
    comparisons=0;samples=[];max_t=32
    ordinary=[F(int(a==0)) for a in range(max_t+3)]
    tilted=[one if a==0 else zero for a in range(max_t+3)]
    for t in range(1,max_t+1):
        old,oldt=ordinary,tilted
        ordinary=[F(1)]+[F(0)]*(max_t+2)
        tilted=[one]+[zero]*(max_t+2)
        for a in range(1,t+1):
            for k in range(1,(t-a+2)//2+1):
                target=a+2*k-3
                p=F(1,2**k);pk=mul((p,F(0)),power(q,2*k-3))
                need(pk==mul(one_minus_r,power(r,k-1)),'letter law')
                ordinary[a]+=p*old[target]
                tilted[a]=add(tilted[a],mul(pk,oldt[target]))
            need(tilted[a]==mul((ordinary[a],F(0)),power(q,-a)),'independent horizon recursions')
            comparisons+=1
        if t in (6,12,24,32):samples.append({'horizon':t,'a6_ordinary':str(ordinary[6]),'a6_conditional':[str(x) for x in tilted[6]]})
    need(power(q,6)==(F(5),F(-8)),'resource-six probability in basis (1,q)')
    print(json.dumps({'status':'FINITE-EXACT PASS','horizon_range':[1,max_t],
        'counter_range':'1..horizon','independent_horizon_comparisons':comparisons,
        'basis':'(1,q), q=(sqrt(5)-1)/2','q6':[5,-8],
        'q6_in_basis_1_sqrt5':[9,-4],'conditional_mean_letter':[0,2],'conditional_counter_drift':[-3,4],
        'samples':samples,'tail_handling':'letters landing above remaining horizon contribute exactly zero; never renormalized'},indent=2,sort_keys=True))
if __name__=='__main__':main()
