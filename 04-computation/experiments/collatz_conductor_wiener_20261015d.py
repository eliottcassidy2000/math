#!/usr/bin/env python3
"""Exact grading/atom-cost controls; the complex-parameter norm theorem is proved in the note."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
def need(test,message):
    if not test:raise RuntimeError(message)
def conductor(k,Q):
    if not k:return 0
    while k%3==0:k//=3;Q-=1
    return Q

def main():
    levels=[]
    for Q in range(1,9):
        m=3**Q;i2=pow(2,-1,m);counts=[0]*(Q+1)
        for k in range(1,m):
            q=conductor(k,Q);counts[q]+=1
            need(conductor(i2*k%m,Q)==q,'halving preserves conductor')
            target=3*i2*k%m
            need(conductor(target,Q)==q-1,'odd branch lowers conductor exactly once')
            if q==1:
                need((i2*k-2*k)%m==0,'rank-one correction cancels constant odd image')
        need(counts[1:]==[2*3**(q-1) for q in range(1,Q+1)],'full conductor multiplicities')
        R=F(3,2)
        mass=sum(F(counts[q])*R**q for q in range(1,Q+1))
        formula=2*R*((3*R)**Q-1)/(3*R-1)
        need(mass==formula,'marked source atom exact Fourier cost')
        levels.append(dict(Q=Q,nonzero_modes=m-1,counts=counts[1:],atom_R_three_halves=[mass.numerator,mass.denominator]))
    # At sigma=1/2, R=sqrt(2), gamma^2=7/12+1/sqrt(6)<1 iff 24<25.
    need(F(1,6)<F(25,144),'strict weighted contraction radical certificate')
    # Boundary N=1: a mode with Re(rho)>0 has 0^rho=0, so the floored sum is 1.
    need(1//3==0 and (2*1)//3==0,'floor-annihilation hostile')
    # Literal floor defect vanishes at every multiple of 3.
    for N in range(3,3001,3):
        need(F(2*N,3)==(2*N)//3 and F(N,3)==N//3,'divisible-by-three exact cancellation boundary')
    report=dict(status='FINITE-EXACT grading and atom controls for the proved weighted-norm theorem',
                levels=levels,gamma_square_certificate=dict(constant=[7,12],radical='1/sqrt(6)',strict_comparison='24<25'),
                floor_hostile=dict(N=1,floored_mode_value=1,scope='any Re(rho)>0, including a D-zero'),
                divisible_by_three_controls=1000)
    data=json.dumps(report,indent=2,sort_keys=True)+'\n'
    path=ROOT/'05-knowledge/results/collatz_conductor_wiener_20261015d.json';path.write_text(data)
    print(json.dumps(dict(output=str(path),sha256=hashlib.sha256(data.encode()).hexdigest()),sort_keys=True))
if __name__=='__main__':main()
