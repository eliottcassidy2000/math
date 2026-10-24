"""Exact graded charge sections on the twenty two-clock sources."""
from common import *
from collections import Counter,defaultdict
def add(a,b):return ((a>>4)+(b>>4))%128*16 | ((a^b)&15)
def neg(a):return (-(a>>4))%128*16 | (a&15)
def sub(a,b):return add(a,neg(b))
def code(r):return r['label'][1]*16+r['label'][0]
def cost(r):
    p={o['target']:o['path'] for o in r['options']}
    return sum(a!=b for a,b in zip(p[1965],p[2234]))
def charge(mask,labels):
    z=0
    while mask:
        b=mask&-mask;z=add(z,labels[b.bit_length()-1]);mask-=b
    return z
def sections(labels,weights=None,limit=None):
    """Minimal (weight,mask) witness at every grade and group element."""
    if weights is None:weights=[1]*len(labels)
    if limit is None:limit=len(labels)
    D=[{} for _ in range(limit+1)];D[0][0]=(0,0)
    for i,(a,w) in enumerate(zip(labels,weights)):
        for k in range(min(i+1,limit),0,-1):
            for g,(c,m) in list(D[k-1].items()):
                h=add(g,a);v=(c+w,m|(1<<i))
                if h not in D[k] or v<D[k][h]:D[k][h]=v
    return D
def main():
    rows=json.loads((HERE/'menus.json').read_text())['rows']
    R=[r for r in rows if r['role']=='core' and cost(r)==2]
    labels=list(map(code,R));D=sections(labels)
    need(len(R)==20,'twenty genuine two-clock rows')
    need(all(len(D[k])==2048 for k in range(5,16)),'complete entire band, no fixed helper occupancy')
    # Independent full subset enumeration: tuples, no packed addition or DP.
    bins=[Counter() for _ in range(21)];prev=0;x=y=k=0
    for i in range(1<<20):
        m=i^(i>>1)
        if i:
            bit=(m^prev).bit_length()-1;sgn=1 if m>>bit&1 else -1
            y^=R[bit]['label'][0];x=(x+sgn*R[bit]['label'][1])%128;k+=sgn
        bins[k][(y,x)]+=1;prev=m
    need([len(x) for x in bins]==[len(x) for x in D],'independent Gray enumeration agrees at EVERY grade')
    need([len(x) for x in D]==[len(x) for x in D][::-1],'complement control')
    for k,d in enumerate(D):
        for g,(_,m) in d.items():need(m.bit_count()==k and charge(m,labels)==g,'every DP witness evaluated')
    # Authentic W signatures; the complete support has only five eligible sources.
    W=[r for r in rows if any(o['target']==570 for o in r['options'])];fibres=defaultdict(list)
    for m in range(1<<len(W)):fibres[(m.bit_count(),charge(m,list(map(code,W))))].append(m)
    need(len(W)==5 and len(fibres)==30,'complete W subset classification')
    need(all(len({bool(m>>3&1) for m in v})==1 for v in fibres.values()),'2393 occupancy fixed by W-owned charge')
    return {'status':'FINITE-EXACT independent DP and full 2^20 subset enumeration',
            'reservoir_sources':[r['source'] for r in R],'coverage':[len(x) for x in D],
            'multiplicity_ranges':[[min(b.values()),max(b.values())] for b in bins],
            'sections':[{str(g):m for g,(_,m) in d.items()} for d in D],
            'W_sources':[r['source'] for r in W],
            'W_fibres':[{'k':k,'charge':g,'masks':v} for (k,g),v in sorted(fibres.items())]}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
