#!/usr/bin/env python3
"""Exact arithmetic checks for the bounded actual-source profinite obstruction."""
import json

def need(v,msg):
    if not v:raise RuntimeError(msg)

def K(H):
    ans=0
    for b in range(3):
        v=5**b
        if v>H:continue
        a=0
        while v*3<=H:v*=3;a+=1
        ans=max(ans,3*a+4*b)
    return ans

def bruteK(H):
    ans=0;b=0;p=1
    while p<=H:
        a=0;q=p
        while q<=H:ans=max(ans,3*a+4*b);a+=1;q*=3
        b+=1;p*=5
    return ans

def budget(H,I):
    N=K(H)+(1 if I==0 else 3+I.bit_length()-1)
    C=sum(1<<(n*(n+1)//2) for n in range(1,N+1))
    r=1
    while C*(r+1)>=(1<<r)-1:r+=1
    return {'H':H,'I':I,'K':K(H),'N':N,'base_universe_bound':C,'r':r,
            'modulus':(1<<r)-1,'residue_bound':C*(r+1),
            'certified_missing_residues_at_least':(1<<r)-1-C*(r+1)}

def main():
    for H in range(1,10001):need(K(H)==bruteK(H),'core budget formula mismatch')
    samples=[budget(H,I) for H,I in ((1,0),(1,1),(3,0),(5,3),(13,16),(100,100),(1000,1000))]
    for a in samples:need(a['residue_bound']<a['modulus'],'no actual modulus obstruction')
    # Actual H1,L0 scalar language has exceptional1 and dyadic ray2^k+1,k>=1.
    residues=sorted({1}|{(pow(2,k,15)+1)%15 for k in range(1,5)})
    need(residues==[1,2,3,5,9],'unit height zero interior residue language')
    # Complete period checks for arbitrary signed affine coefficients/offsets.
    checks=0
    for r in range(1,12):
        M=(1<<r)-1
        for a in range(-5,6):
            for d in range(-5,6):
                vals={(a*pow(2,k,M)+d)%M for k in range(4*r+1)}
                need(len(vals)<=r,'affine dyadic period bound');checks+=1
    print(json.dumps({'status':'PASS','exact_core_budget_checks':10000,
        'signed_affine_ray_period_checks':checks,'H1_L0_mod15_residues':residues,
        'certified_budget_moduli':samples,
        'scope':'Corollary of the inherited actual-source normal form; not a full-spectrum modular omission.'},sort_keys=True,indent=2))

if __name__=='__main__':main()
