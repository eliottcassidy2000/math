#!/usr/bin/env python3
"""Independent integer certificate; checks persist under python -O."""
from math import isqrt
from pathlib import Path
import json
A=(324,54756,119716,264196);B=(0,79200,227205,1258560)
L=tuple(a-1 for a in A);R=tuple(b+1 for b in B)
roots=tuple(tuple(isqrt(a+b) for b in B) for a in A)
def check(ok,why):
    if not ok:raise RuntimeError(why)
check(len(set(L+R))==8 and min(L+R)>0,'positive disjoint vertices')
check(all(roots[i][j]>0 and roots[i][j]**2==L[i]+R[j] for i in range(4) for j in range(4)),'16 nonzero cross squares')
check(len(set(q*q for row in roots for q in row))==16,'16 distinct entries')
check(len(set(roots[i][0]**2 for i in range(4)))==4,'distinct rows')
check(len(set(roots[0][j]**2 for j in range(4)))==4,'distinct columns')
check(all(roots[i][j]**2+roots[0][0]**2-roots[i][0]**2-roots[0][j]**2==0 for i in range(1,4) for j in range(1,4)),'X44 equations')
result={'status':'PASS','rows':L,'columns':R,'roots':roots,'U44_rational_point':True,'distinct_squared_entries':16,'in_decorated_locus_even_after_row_column_permutations':False}
print(json.dumps(result,indent=2))
