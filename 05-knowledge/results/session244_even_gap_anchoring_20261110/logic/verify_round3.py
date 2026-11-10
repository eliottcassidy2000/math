#!/usr/bin/env python3
"""Independent labelled census, forbidden-core orbit check, and literal paths."""
import itertools,json
from functools import lru_cache

def need(x,msg):
    if not x:raise RuntimeError(msg)

def labelled(n,code):
    a=[0]*n
    for i,j in itertools.combinations(range(n),2):
        if code&1:a[i]|=1<<j
        else:a[j]|=1<<i
        code>>=1
    return tuple(a)

def sub(a,vs):
    return tuple(sum(1<<j for j,w in enumerate(vs) if a[v]>>w&1) for v in vs)

def delete(a,v):return sub(a,[u for u in range(len(a)) if u!=v])

def reach(a,v):
    s=1<<v
    while True:
        t=s
        for i in range(len(a)):
            if s>>i&1:t|=a[i]
        if s==t:return s==(1<<len(a))-1
        s=t

def strong(a):return all(reach(a,v) for v in range(len(a)))

@lru_cache(None)
def degone(a):
    return len(a)<=1 or any(x.bit_count()<=1 and degone(delete(a,v)) for v,x in enumerate(a))

@lru_cache(None)
def guarded(a):
    if len(a)<=2:return True
    for v,x in enumerate(a):
        if x.bit_count()>1:continue
        b=delete(a,v)
        if x:
            r=x.bit_length()-1;r-=r>v
            if not reach(b,r):continue
        if guarded(b):return True
    return False

def hp(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            z=dp[s][v]
            if not z:continue
            for w in range(n):
                if a[v]>>w&1 and not s>>w&1:dp[s|1<<w][w]+=z
    return sum(dp[-1])

cores=[(6,12,24,17,3),(6,12,48,5,11,27),(6,12,24,33,11,23),(6,12,24,33,43,7),(6,24,26,33,41,7),(6,24,42,33,13,19)]
orbits=[]
for a in cores:
    orbit={sub(a,p) for p in itertools.permutations(range(len(a)))}
    need(all(not orbit.intersection(o) for o in orbits),'duplicate core orbit')
    orbits.append(orbit)
expected={n:set().union(*(o for a,o in zip(cores,orbits) if len(a)==n)) for n in (5,6)}
census=[]
for n in range(1,7):
    minimal=set();nd=ng=0;mismatch=strong_mismatch=None
    for code in range(1<<(n*(n-1)//2)):
        a=labelled(n,code)
        d=degone(a);g=guarded(a)
        nd+=d;ng+=g
        need(not g or d,'guard subset')
        if d and not g:
            if mismatch is None:mismatch=a
            if strong(a) and strong_mismatch is None:strong_mismatch=a
        if min(x.bit_count() for x in a)>=2 and all(min(y.bit_count() for y in delete(a,v))<=1 for v in range(n)):
            minimal.add(a)
    need(minimal==expected.get(n,set()),'complete core census')
    census.append({'n':n,'labelled':1<<(n*(n-1)//2),'degone':nd,'guarded':ng,'minimal_cores':len(minimal),'first_mismatch':mismatch,'first_strong_mismatch':strong_mismatch})

# Literal guarded words from edge 0->1, rerooting only at J.
def word(w):
    a=(2,0);r=0
    for c in w:
        n=len(a)
        need(reach(a,r),'root guard')
        b=list(a)+[0 if c=='O' else 1<<r]
        for v in range(n):
            if c=='O' or v!=r:b[v]|=1<<n
        if c=='J':r=n
        a=tuple(b)
    return a
witnesses={15:'OJOZ',33:'OOOOZ',37:'OZZOZ',45:'OZOZZ'}
cert=[]
for h,w in witnesses.items():
    a=word(w);need(hp(a)==h and guarded(a) and strong(a),'scalar core replacement')
    cert.append({'H':h,'word':w,'rows':a})
# Strong five-vertex guarded obstruction: source followed by C3, then r->source.
b=(14,4,8,2);a=tuple(x|16 for x in b[0:1])+tuple(x|16 for x in b[1:])+ (1,)
# The source is the unique old vertex beaten by r, so remove its arc to r.
a=(a[0]^16,)+a[1:]
need(strong(a) and degone(a) and not guarded(a) and hp(a)==15,'strong order5 obstruction')
print(json.dumps({'status':'PASS','census':census,'cores':[{'rows':a,'orbit_size':len(o),'H':hp(a)} for a,o in zip(cores,orbits)],'guarded_scalar_witnesses':cert,'strong_guard_obstruction':{'rows':a,'H':hp(a)},'scope':'Finite census through order6 plus infinite proof in ROUND3.md.'},indent=2))
