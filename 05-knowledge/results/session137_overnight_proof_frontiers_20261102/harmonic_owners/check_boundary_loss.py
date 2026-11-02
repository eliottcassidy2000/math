#!/usr/bin/env python3
"""Sharp nonunit-owner sidecar hostile inside the requested 1-in-A family."""
from fractions import Fraction as F
import json
from check import setup,windows,owners,safe,require

q=8; A={1,2}; B={1,2,6}; r=F(1,11); t=F(4,11); eta=F(1,10000)
_,da,_,EA,CA=setup(q,A); _,db,_,EB,CB=setup(q,B)
require((da,EA)==(db,EB),'identical chamber coordinates')
for x in (r,(r+F(1,q))/2,F(1,q)):
    require(windows(CA,x)==windows(CB,x),'same complete closed atlas')
require(owners(CA,t,r)==[(3,'+')],'old unique owner')
require(owners(CB,t,r)==[(3,'+'),(14,'+')],'new invisible tie')
u=t-eta/3
require(safe(CA,u,r-eta) and not safe(CB,u,r-eta),'different continuation')
print(json.dumps({'status':'FINITE-EXACT / PASS','q':q,'A':sorted(A),'B':sorted(B),
 'common_d':da,'common_E':list(EA),'radius':str(r),'point':str(t),
 'A_owners':owners(CA,t,r),'B_owners':owners(CB,t,r),
 'lower_radius':str(r-eta),'continuation_point':str(u),'A_safe':True,'B_safe':False},
 indent=2,sort_keys=True))
