#!/usr/bin/env python3
"""Separate coarse-cylinder proof search. Does not import the producer.

Stop at the first coefficient contraction; combine all greater terminal
valuations into one coarse cylinder. Independently use suffix sums for C.
"""
import argparse
import json
from fractions import Fraction

def need(b,msg):
    if not b:raise RuntimeError(msg)

def carry(word):
    B=C=0;K=0;m=len(word)
    for i,a in enumerate(word):
        B+=3**(m-i-1)*2**K
        K+=a
        if a==1:C+=3**(m-i-1)*2**K
    return B,C,K

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--depth',type=int,default=18)
    args=ap.parse_args();need(args.depth>=1,'depth')
    stack=[((),0,1)];counts=[0]*args.depth;leastmargin=None
    exceptions=[];maxsharp=(Fraction(0),0,0)
    for m in range(1,65):
        A=3**m;K=A.bit_length();d=2**K
        w=(K-m+1,)+(1,)*(m-1);B,C,kk=carry(w)
        need(kk==K,'sharp clock')
        ratio=(1+Fraction(C,d-A))/d
        if ratio>maxsharp[0]:maxsharp=ratio,m,K
    need(maxsharp==(Fraction(161,256),5,8),'finite sharp maximum')
    while stack:
        w,K,A=stack.pop();m=len(w)+1;AA=3*A
        a=AA.bit_length()-K;need(a>=2,'first crossing exponent')
        ww=w+(a,);B,C,kk=carry(ww);d=2**kk
        r=(-B*pow(AA,-1,d))%d
        n=r if r>1 else r+d
        need(n>1 and n%2==1,'positive odd source')
        margin=(d-AA)*(n-1)-C
        if margin<=0:exceptions.append({'word':ww,'source':n,'margin':margin})
        counts[m-1]+=1
        if leastmargin is None or margin<leastmargin:leastmargin=margin
        if m<args.depth:
            for b in range(1,a):stack.append((w+(b,),K+b,AA))
    print(json.dumps({'status':'FINITE-EXACT separate coarse-cylinder checker',
        'depth':args.depth,'exit_cylinders_by_length':counts,
        'total_exit_cylinders':sum(counts),'exceptions':exceptions,
        'minimum_integer_payment_margin':leastmargin,
        'sharp_threshold_check_m_1_to_64':{'ratio':str(maxsharp[0]),'m':maxsharp[1],'K':maxsharp[2]}},
        indent=2,sort_keys=True))
    need(not exceptions,'unpaid source')

if __name__=='__main__':main()
