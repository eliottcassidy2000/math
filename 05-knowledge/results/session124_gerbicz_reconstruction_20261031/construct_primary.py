#!/usr/bin/env python3
"""Independent constructive expansion of mathematically checked integer tables."""
import sys,json,hashlib
sys.dont_write_bytecode=True
from functools import lru_cache
from verify_gerbicz_data import arrays,pathcheck,check
D,sha=arrays(sys.argv[1])
@lru_cache(None)
def pair(n):
    if 41<=n<=2032:return tuple(tuple(x) for x in D['S'][2*(n-41):2*(n-41)+2])
    if n<41:raise ValueError('pair index outside invariant')
    m=(n-24)//49;d=n-49*m;inp=pair(m);out=[]
    for ext in (0,1):
        row=D['data'][4*(d-24)+2*(m%2)+ext];w=[]
        check(row[0]==[49,d,m%2,ext],'state addressing')
        for c,orient,e in row[1:]:
            if orient==0:w.append(c)
            else:
                block=[49*x+(c if j%2==0 else -c) for j,x in enumerate(inp[e])]
                w.extend(block if orient==1 else reversed(block))
        out.append(tuple(w))
    return tuple(out)
rows=[]
for n in [2033,2034,2081,1000000]:
    a,b=pair(n)
    for e,w in enumerate((a,b)):pathcheck(list(w),n+e,True)
    pa={v:j%2 for j,v in enumerate(a)}
    check(all(pa[v]==j%2 for j,v in enumerate(b) if v<=n),'expanded pair parity')
    rows.append({'n':n,'sha256':hashlib.sha256((' '.join(map(str,a))+'\n').encode()).hexdigest(),'prefix':list(a[:8]),'end':a[-1]})
print(json.dumps({'status':'PASS','source_sha256':sha,'expanded_pairs':rows,'external_code_executed':False},sort_keys=True,indent=2))
