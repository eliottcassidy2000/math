#!/usr/bin/env python3
"""Exact reduced forms and independent PARI class-group controls.

All-order claims are proved in the companion note; these are finite controls.
"""
import json
import subprocess
from math import gcd
from quadratic_context_20261017c import forms


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def main():
    rows=[]
    for k in range(7):
        D=-1848*4**k; fs=forms(D)
        fixed=[f for f in fs if f[1]==0 or abs(f[1])==f[0] or f[0]==f[2]]
        need(len(fs)==8*2**k,'conductor class-number formula')
        need(len(fixed)==(8 if k==0 else 16),'ambiguous class count')
        witness=None
        if k>=2:
            witness=(16,8,1+462*4**(k-2))
            a,b,c=witness
            need(b*b-4*a*c==D and gcd(gcd(a,b),c)==1 and 0<b<a<c,'strictly interior primitive reduced inverse pair')
            need(witness in fs and (a,-b,c) in fs,'both distinct inverse forms appear')
        rows.append(dict(k=k,idoneal_parameter=462*4**k,conductor=2**k,
                         discriminant=D,class_number=len(fs),inversion_fixed=len(fixed),
                         inversion_pairs=(len(fs)-len(fixed))//2,genus_number=len(fixed),
                         merged_inversion_orbits=(len(fs)+len(fixed))//2,
                         strict_witness=witness,forms=fs))
    gp='for(k=0,6,print(quadclassunit(-1848*4^k)[1..2]))\n'
    result=subprocess.run(['gp','-fq'],input=gp,text=True,capture_output=True,check=True)
    groups=[json.loads(line) for line in result.stdout.splitlines()]
    need(groups==[[8,[2,2,2]]]+[[8*2**k,[2**k,2,2,2]] for k in range(1,7)],'independent PARI group invariants')
    # Inversion orbits, the genus quotient, and class number are distinct.
    controls=[]
    for D in (-44,-76,-84):
        fs=forms(D);a=sum(b==0 or abs(b)==x or x==c for x,b,c in fs)
        controls.append(dict(D=D,h=len(fs),fixed=a,genus=a,orbits=(len(fs)+a)//2))
    return dict(status='FINITE-EXACT; all-order cutoff mechanism in companion proof',
                rows=rows,pari_groups=groups,odd_kernel_controls=controls)


if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
