#!/usr/bin/env python3
"""Exact affine-chart order contrasts and source guards."""
from fractions import Fraction as F
from itertools import product
import json

def need(b,m):
    if not b:raise RuntimeError(m)
def v2(n):return (n&-n).bit_length()-1
def carry(b,w):
    c=0;d=0
    for a in w:c=b*c+2**d;d+=a
    return c,d
def residue(b,w):
    c,d=carry(b,w);M=2**(d+1)
    return ((2**d-c)*pow(b**len(w),-1,M))%M,M

def main():
    checks=0;chart=[]
    W=[(1,1,2),(1,2,1),(2,1,1)]
    for b in range(3,32,2):
        cs=[carry(b,w)[0] for w in W]
        need(cs==[b*b+2*b+4,b*b+2*b+8,b*b+4*b+8],'triple carry formula');checks+=1
        need(F(cs[2]-cs[1],cs[1]-cs[0])==F(b,2),'multiplier fingerprint');checks+=1
        rows=[]
        for w,c in zip(W,cs):
            r,M=residue(b,w)
            for t in range(25):
                n=r+M*t;x=n
                for a in w:
                    need(v2(b*x+1)==a,'exact chart guard');x=(b*x+1)//2**a;checks+=1
                need(x==F(b**3*n+c,16),'native affine chart');checks+=1
            rows.append({'word':w,'carry':c,'source_residue':r,'modulus':M})
        if b in (3,5):chart.append({'b':b,'contrast':str(F(b,2)),'charts':rows})
        for a,c in product(range(1,8),repeat=2):
            need(carry(b,(a,c))[0]-carry(b,(c,a))[0]==2**a-2**c,'two-letter blindness');checks+=1
    # A source cannot lie in two different exact three-word guard cylinders.
    for b in (3,5):
        rs=[residue(b,w)[0] for w in W];need(len(set(rs))==3,'distinct guards');checks+=1
    # Two native branches force the coefficient constraints for any affine candidate.
    for b,c,A,B in product((3,5,7),(3,5,7),range(-3,4),range(-3,4)):
        if A==0:continue
        ok=all(F(A*(b*x+1),2**a)+B==F(c*(A*x+B)+1,2**a) for a in (1,2) for x in (1,3))
        need(ok==(b==c and A==1 and B==0),'affine rigidity finite hostile');checks+=1
    print(json.dumps({'status':'PASS','checks':checks,'charts':chart,'scope':'Finite controls accompany the rational-function proof; no convergence claim.'},indent=2))
if __name__=='__main__':main()
