#!/usr/bin/env python3
"""Finite linear identities certify an infinite valuation-one seven family."""
from pathlib import Path
import json
from refill import AP, polymul, value, wrap, count, require

def step(v,q):
    return tuple((t*v[(t+1)%q]+t*t*v[t]-t*v[(t-1)%q])%q for t in range(q))

def main():
    f=polymul(AP,AP); cert=[]
    wf=f
    for _ in range(4): wf=wrap(wf)
    first=count(4,polymul(f,wf))
    for q,wanted in ((3,2),(5,3),(49,35)):
        v=tuple(value(f,t)%q for t in range(q)); saved={}
        for d in range(339):
            if d in (2,4,338): saved[d]=v
            v=step(v,q)
        require(saved[2]==saved[338],'linear return certificate failed')
        final=tuple(value(f,t)*x%q for t,x in enumerate(saved[4]))
        h,p,w=final[1],final[2%q],final[3%q]
        residue=(120*81*w+(240*81+1)*p+(1-120*81)*h)%q
        require(residue==wanted and first%q==wanted,'seven-phase readout failed')
        cert.append({'modulus':q,'state_at_2':saved[2],'state_at_338':saved[338],
                     'phase_4_target_residue':residue})
    receipt={'status':'FINITE-EXACT linear certificates imply the stated all-parameter family',
             'family':{'j0':1,'word':'(0,4+336t,0), t>=0',
                       'order':'49+672t','first_H':first,
                       'first_cofactor':first//7,
                       'scope':'exact v7=1; cofactors coprime105; no primality claim'},
             'initial_polynomial':f,'certificates':cert,
             'reproduce':'python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/phase_certificate.py'}
    Path(__file__).with_name('phase_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt['family'],indent=2))

if __name__=='__main__': main()
