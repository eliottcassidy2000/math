#!/usr/bin/env python3
"""Independent strict-clock checks of the signed transfer and cutoff identity.

Python stdlib only. Reproduces census transfer rows through M=2^18 by
single Terras steps, and checks finite Dirichlet identities with Fractions.
"""
import argparse
import json
from fractions import Fraction as Q


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--census',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args()
    N=3*2**18
    tau=[-1]*(N+1);tau[1]=0
    checks=0
    def require(ok,message):
        nonlocal checks
        checks+=1
        if not ok:raise RuntimeError(message)
    for n in range(2,N+1):
        m=n;path=[]
        while m>N or tau[m]<0:
            path.append(m);m=(3*m+1)//2 if m%2 else m//2
            require(len(path)<100000,'finite control cap')
        clock=tau[m]
        for m in reversed(path):
            clock+=1
            if m<=N:tau[m]=clock
    sign=[0]+[1 if t%2==0 else -1 for t in tau[1:]]
    S=[0]*(N+1);V=[0]*(N+1)
    for n in range(1,N+1):
        S[n]=S[n-1]+sign[n]
        V[n]=V[n-1]+(n>1 and sign[n]!=sign[n-1])
    matched=0
    for M,sm,s2,s3,cuts,numerator in json.load(open(args.census))['dyadic_transfers']:
        if 3*M>N:continue
        require([S[M],S[2*M],S[3*M],V[3*M],3*S[2*M]+3*S[M]+S[3*M]]==
                [sm,s2,s3,cuts,numerator],'independent transfer census')
        matched+=1
    # The exact first-hit clock, including the single terminal exception.
    for m in range(1,N//3+1):
        require(tau[2*m]==tau[m]+1,'even clock')
        require(tau[2*m-1]==(0 if m==1 else tau[3*m-1]+1),'odd clock defect')
    # Cutoff Dirichlet identity: keep each image endpoint, including +1.
    for exponent in range(1,5):
        F=[Q(0)]
        for n in range(1,97):F.append(F[-1]+Q(sign[n],n**exponent))
        for cutoff in range(1,64):
            endpoint=3*((cutoff+1)//2)-1
            residue=Q(0);affine=Q(0)
            for n in range(2,endpoint+1,3):
                term=Q(sign[n],n**exponent)
                residue+=term
                affine+=term*((1-Q(1,2*n))**(-exponent)-1)
            lhs=F[cutoff]+Q(1,2**exponent)*F[cutoff//2]+Q(3,2)**exponent*(residue+affine)
            require(lhs==0,'finite Dirichlet/affine identity')
    out={'status':'FINITE-EXACT','strict_clock_sources':N,
         'matched_transfer_rows':matched,'rational_dirichlet_cutoffs':252,
         'checks':checks,'terminal_defect':'1-zeta^2 at odd source1 only'}
    with open(args.output,'w') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
