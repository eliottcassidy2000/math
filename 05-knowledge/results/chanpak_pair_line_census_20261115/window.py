"""Does singleton padding of a complete arithmetic window yield its dilates?"""
from score import band,read
from collections import defaultdict
import json,sys
T=read(sys.argv[1])[10]
orig=[(h,q,x) for h,q,x,z in T if 36835<=h+q<=59999]
data=[]
for s in range(1,7):
    scale=(1<<s)
    lo=scale*36835+1;hi=scale*59999-1
    ex={h+scale*q for h,q,x,z in T if lo<=h+scale*q<=hi}
    strict={h+scale*q for h,q,x in orig if lo<=h+scale*q<=hi}
    data.append([s,lo,hi,band(ex,lo,hi),band(strict,lo,hi)])
print(json.dumps(data,indent=2))
