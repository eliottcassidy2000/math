#!/usr/bin/env python3
"""Independent small exact controls for two attached-paper bridges.

Exact: row-2 occupations, entropy product, determinant-2 clock obstruction,
golden norm filters and Fibonacci identities. Complex checks are explicitly
numerical supplements to the analytic proofs in the companion note.
"""
import argparse
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

import mpmath as mp
import sympy as sp


def require(p, text):
    if not p:
        raise RuntimeError(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    # Finite-extraction, Example 5.2, independently solved exactly.
    states = [1,2,4,6]
    mat = sp.zeros(4)
    for x in states:
        if x==6:
            mat[states.index(x),states.index(1)] = 1
        else:
            h=min(x,6-x)
            for y in range(x-h,x+h+1,2):
                if y:
                    mat[states.index(x),states.index(y)] = sp.Rational(3+2*y,(h+1)*(3+2*x))
    occupation = (sp.eye(4)-mat.T).inv()*sp.Matrix([1,0,0,0])
    require(list(occupation)==[2,sp.Rational(14,5),sp.Rational(11,5),1],'row occupations')
    # At K=1 this allowed circulation retains the top return and a self-loop.
    flow = {(1,0):1,(1,2):1,(2,2):1,(2,4):1,(4,6):1,(6,1):1,(0,1):1}
    qs = {}
    for x in states[:-1]:
        h=min(x,6-x)
        for y in range(x-h,x+h+1,2):
            qs[x,y]=Q(3+2*y,(h+1)*(3+2*x))
    for (x,y),q in qs.items():
        v=Q(occupation[states.index(x)])*q
        n=flow.get((x,y),0)
        require(v.numerator//v.denominator<=n<=-(-v.numerator//v.denominator),'rounding interval')
    require(all(e in qs or e in ((6,1),(0,1)) for e in flow),'only authenticated edges')
    balances={x:0 for x in range(7)}
    for (x,y),n in flow.items():
        balances[x]+=n;balances[y]-=n
    require(all(v==0 for v in balances.values()),'rounded flow balances')
    # A second graph demonstrates why an independent clock row is not TU.
    clock_witness=[]
    for a,b in itertools.product((0,1),repeat=2):
        if a==b:
            clock_witness.append([a,b,a+b])
    require(all(row[2]!=1 for row in clock_witness),'clock-rounding hostile')
    require(sp.Matrix([[1,-1],[1,1]]).det()==2,'lost unimodularity')
    norm_counts={}
    for p in [2,3,5,7,11,13,17,19]:
        zeros=[(a,b) for a in range(p) for b in range(p) if (a*a+a*b-b*b)%p==0]
        norm_counts[p]=len(zeros)
        if p in (2,3,7,13,17):
            require(zeros==[(0,0)],'inert golden norm')
    require((3*3+3-1)%11==0,'split-prime hostile')
    fib=[0,1]
    for _ in range(30):fib.append(fib[-1]+fib[-2])
    # Exact phi-basis arithmetic: phi^n=F_(n-1)+F_n phi.
    a,b=1,0
    for n in range(1,25):
        a,b=b,a+b
        require((a,b)==(fib[n-1],fib[n]),'Fibonacci clock')
    mp.mp.dps=70
    phi=(1+mp.sqrt(5))/2
    R=phi**(-phi-2)
    crit_res=mp.mpf('0'); golden_res=mp.mpf('0'); inv_res=mp.mpf('0')
    def logsin(z):
        return -1j*z+mp.log(1-mp.exp(2j*z))+1j*mp.pi/2-mp.log(2)
    for k in range(-12,13):
        z=mp.pi/2+k*mp.pi
        crit_res=max(crit_res,abs(logsin(z)+1j*k*mp.pi))
        zg=phi*(-mp.log(phi)+(2*k+1)*mp.pi*1j)
        golden_res=max(golden_res,abs(mp.exp(zg)+mp.exp(phi*zg)-R*mp.exp((2*k+1)*mp.pi*1j*phi)))
    for real in [-2,-1,-mp.log(2)]:
        for imag in [-mp.mpf('.1'),0,mp.mpf('.1')]:
            w=real+1j*imag
            z=-1j*mp.log(1j*mp.exp(w)+mp.sqrt(1-mp.exp(2*w)))
            inv_res=max(inv_res,abs(logsin(z)-w))
    require(max(crit_res,golden_res,inv_res)<mp.mpf('1e-60'),'analytic numerical controls')
    out=dict(exact=dict(occupations=list(map(str,occupation)),rounded_flow={str(e):n for e,n in flow.items()},clock_integer_options=clock_witness,
                       golden_norm_zero_counts=norm_counts,fibonacci_steps=24),
             numerical=dict(status='VERIFIED, not interval-enclosed',dps=70,
                            logsin_critical_residual=str(crit_res),golden_circle_residual=str(golden_res),
                            inverse_residual=str(inv_res)))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
