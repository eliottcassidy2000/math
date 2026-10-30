#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
from epoch_semantics import closed_epoch,lam,need
import json

checks=0
for r in range(1,4):
    T=2*5**(r-1);TT=5*T
    for h in range(r+1,r+4):
        for d in range(T):
            for j in range(5):
                dd=d+j*T;m=pow(3*pow(2,-1,5**h),-dd,5**h)
                for L in sorted({0,max(0,d-1),d,max(0,dd-1),dd,dd+1,dd+TT-1,dd+TT,dd+2*TT}):
                    gap=closed_epoch(m,h,L,r)-closed_epoch(m,h,L,r+1)
                    if L<d:formula=F(0)
                    elif L<dd:formula=h-r+(lam(h)-lam(r))*(L-d)
                    else:formula=F(L-dd,TT)%1+F(j,5)+(lam(h)-lam(r+1))*j*T
                    need(gap==formula and gap>=0,(r,h,d,j,L,gap,formula))
                    checks+=1
out={'status':'PASS','exact_delay_credit_checks':checks,'r':'1..3','h':'r+1..r+3',
     'delays':'every d_r and j=0..4','clock_cases':'thresholds and two later periods'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
