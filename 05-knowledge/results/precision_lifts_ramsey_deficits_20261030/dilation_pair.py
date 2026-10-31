#!/usr/bin/env python3
"""Two-tail dilation descent: retain the one common lift label."""
from fractions import Fraction as F
from math import gcd,ceil
from pathlib import Path
import json
from joint_phase import old,need,pair_bank

def run():
    cases=0;exceptions=0
    for d in range(2,45):
        for u in range(1,2*d+1):
            h=d//gcd(d,u)
            for v in range(u,2*d+1):
                k=d//gcd(d,v)
                if h==1 or k==1:continue
                for s in (F(1,7),F(7,50),F(1,3),F(2,5),F(1,2)):
                    good=[j for j in range(d) if min(old.dist(u*(s+j)/d),old.dist(v*(s+j)/d))>=F(1,8)]
                    if (h,k)!=(2,2):
                        need(F((h+3)//4,h)+F((k+3)//4,k)<1,'grid bound')
                        need(good,('joint fibre lifting',d,u,v,s))
                    elif not good:exceptions+=1
                    cases+=1
    # At the working radius 1/14, even the exceptional point can genuinely fail.
    s=F(7,50);r=F(1,14)
    need(all(min(old.dist((s+j)/2),old.dist(7*(s+j)/2))<r for j in (0,1)), 'two-point hostile')
    source=json.loads(Path(__file__).with_name('joint_phase.json').read_text())
    banks=[]
    # It suffices to check odd tails on doubled cores, but the general compiler
    # certifies all pairs, giving a stronger finite positive control.
    for entry in source['all_positive_pair_banks']:
        core=tuple(entry['core']);b=pair_bank(tuple(2*x for x in core),F(1,14))
        need(b['all_pairs_accepted'],'double bank fails')
        banks.append(b);print('doubled core',core,'queries',b['finite_queries'],flush=True)
    out=dict(status='PROVED reduction with FINITE-EXACT grid and doubled-bank controls',
             grid_cases=cases,two_point_failures=exceptions,
             hostile=dict(d=2,u=1,v=7,s='7/50',radius='1/14'),doubled_banks=banks)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(grid_cases=cases,two_point_failures=exceptions,bank_queries=sum(b['finite_queries'] for b in banks)),indent=2))

if __name__=='__main__':run()
