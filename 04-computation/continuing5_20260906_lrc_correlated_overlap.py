"""Actual decoder control with zero individual interval credits and positive located count.
Standalone exact arithmetic; inherited profile JSON is data, not an oracle.
"""
from math import gcd, prod, comb
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json, sys
sys.stdout.reconfigure(newline='\n')
N=0

def need(ok, label):
    global N
    N+=1
    if not ok: raise ArithmeticError(label)

Q=91**6
CENTER=83
LEAVES=(198,215,251,257,263)
P=prod(LEAVES)
V=tuple(sorted((P,)+tuple(CENTER*P//q for q in LEAVES)))
GROUPS=(107,149,167,173,179,191,197)
H=prod(GROUPS)
U=tuple(sorted(H//q for q in GROUPS))
SCALES=(968,)
PROFILE_SHA='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'

def inert_sum(n):
    if n>356 or n<2:return False
    p=2
    while p*p<=n:
        e=0
        while n%p==0:n//=p;e+=1
        if e and (p%3!=2 or e>2):return False
        p+=1
    return n==1 or n%3==2

def graph(row):
    return tuple((i,j) for i,j in combinations(range(len(row)),2)
                 if inert_sum((row[i]+row[j])//gcd(row[i],row[j])))

def components(n,edges):
    adj=[set() for _ in range(n)]
    for i,j in edges:adj[i].add(j);adj[j].add(i)
    unseen=set(range(n));out=[]
    while unseen:
        seen={min(unseen)};queue=list(seen)
        for i in queue:
            for j in sorted(adj[i]-seen):seen.add(j);queue.append(j)
        unseen-=seen;out.append(tuple(sorted(seen)))
    return tuple(out)

def rank(rows):
    a=[list(map(F,r)) for r in rows];at=0
    for j in range(len(a[0])):
        k=next((k for k in range(at,len(a)) if a[k][j]),None)
        if k is None:continue
        a[at],a[k]=a[k],a[at];v=a[at][j]
        a[at]=[x/v for x in a[at]]
        for k in range(len(a)):
            if k!=at and a[k][j]:
                v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[at])]
        at+=1
        if at==len(a):break
    return at

def packet(row,d):
    return tuple(j for j in range(d)
                 if all(14*min((v*j)%d,(-v*j)%d)>=d for v in row))

def clear(row,x):return min(min((v*x)%1,(-v*x)%1) for v in row)

def pair_lengths(p,q):
    p,q=sorted((p,q))
    need(gcd(p,q)==1,'primitive pair for interval formula')
    kmax=(p+q+13)//14-1
    return (F(1,7*q),)+tuple(F(min(2*p,p+q-14*k),14*p*q)
                               for k in range(1,kmax+1) for _ in range(2))

def open_grid_min(length,t):
    a=length*t
    return (a.numerator+a.denominator-1)//a.denominator-1

def literal_intersection_lengths(p,q):
    # Build literal circle arcs, intersect by a sorted two-pointer sweep,
    # and rejoin the two pieces of the origin component.
    def arcs(v):
        out=[]
        for k in range(v+1):
            lo=max(F(0),F(14*k-1,14*v));hi=min(F(1),F(14*k+1,14*v))
            if lo<hi:out.append((lo,hi))
        return out
    a,b=arcs(p),arcs(q);i=j=0;out=[]
    while i<len(a) and j<len(b):
        lo=max(a[i][0],b[j][0]);hi=min(a[i][1],b[j][1])
        if lo<hi:out.append((lo,hi))
        if a[i][1]<b[j][1]:i+=1
        elif b[j][1]<a[i][1]:j+=1
        else:i+=1;j+=1
    lengths=[hi-lo for lo,hi in out]
    if out[0][0]==0 and out[-1][1]==1:
        lengths=[lengths[0]+lengths[-1]]+lengths[1:-1]
    return tuple(sorted(lengths))


def located_intervals(p,q):
    def arcs(v):
        return [(max(F(0),F(14*j-1,14*v)),min(F(1),F(14*j+1,14*v))) for j in range(v+1)]
    a,b=arcs(p),arcs(q);i=j=0;out=[]
    while i<len(a) and j<len(b):
        lo=max(a[i][0],b[j][0]);hi=min(a[i][1],b[j][1])
        if lo<hi:out.append((lo,hi))
        if a[i][1]<b[j][1]:i+=1
        elif b[j][1]<a[i][1]:j+=1
        else:i+=1;j+=1
    need(out[0][0]==0 and out[-1][1]==1,'wrap component exists')
    return [(out[-1][0]-1,out[0][1])]+out[1:-1]

def ceil(x):return -((-x.numerator)//x.denominator)
def floor(x):return x.numerator//x.denominator

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--profiles',type=Path,default=Path(__file__).with_name('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'))
    raw=parser.parse_args().profiles.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==PROFILE_SHA,'full inherited profile data pin')
    profiles={int(k):{(c,tuple(w)) for c,w in v['profiles']} for k,v in json.loads(raw)['levels'].items()}
    need(V==(227923228170,233244393030,238819956210,278808413994,302746510145,722214566370),'literal six shape')
    need(U==(15747768383257,16242462678019,17331342857551,17932429893073,18576708811387,20820874976521,28993554873847),'literal seven shape')
    need(P==722214566370 and H==3102310371501629,'cofactor products')
    for xs in (LEAVES,GROUPS):
        for a,b in combinations(xs,2):need(gcd(a,b)==1,'cofactor groups coprime')
    need(gcd(CENTER,P)==1 and reduce(gcd,V)==reduce(gcd,U)==1,'primitive shapes')
    need(gcd(prod(V),prod(U))==1,'disjoint component prime supports')
    need({v%2 for v in V}=={0,1} and all(u%2 for u in U),'mixed primitive V, odd primitive U')
    need(clear(V,F(1,4))==F(1,4),'explicit V-quarter supplier')
    need(min(V)>3*Q//28 and min(U)>1 and max(U)>28,'unitless high-minimum context')
    dmin=min(gcd(a,b) for a,b in combinations(V,2))
    threshold=max(Q//dmin+1,Q*(GROUPS[-1]+GROUPS[-2])//min(V)+1)
    need(dmin==886860810 and threshold==967,'exact actual-entry sufficient threshold')
    need(1042*sum(V)+sum(U)<Q*Q,'whole bounded scale class fits box')
    need(components(6,graph(V))==(tuple(range(6)),) and len(graph(V))==5,'strict actual six-star')
    need(components(7,graph(U))==(tuple(range(7)),) and len(graph(U))==6,'strict actual seven-tree')
    for size in range(7,13):need((1,(1,)*(13-size)) in profiles[13-size],'all-one full profile allowed')
    for p,q in combinations(GROUPS,2):
        ls=pair_lengths(p,q)
        need(tuple(sorted(ls))==literal_intersection_lengths(p,q),'complete clipped geometry by literal arcs')
        need(all(1042*ell<1 for ell in ls),'every individual interval credit zero throughout class')
    t=968;row=tuple(t*v for v in V)+U
    need(gcd(t,H)==1 and 967<=t<=1042,'main class membership')
    need(len(set(row))==13 and reduce(gcd,row)==1 and sum(row)==2075281984219247<Q*Q,'literal primitive boxed physical row')
    ed=graph(row)
    need(components(13,ed)==(tuple(range(6)),tuple(range(6,13))) and len(ed)==11,'complete actual graph')
    relations=[]
    for i,j in ed:
        d=gcd(row[i],row[j]);r=[0]*13;r[i],r[j]=row[j]//d,-row[i]//d
        need(sum(v*c for v,c in zip(row,r))==0 and max(map(abs,r))<=355,'bounded weighted edge relation')
        relations.append(r)
    need(rank(relations)==11,'literal weighted rank')
    forward=[];reverse=[]
    for a,b in combinations(row[:6],2):
        d=gcd(a,b)
        for y in U:
            c=d//gcd(d,y);need(c>Q,'two V one U crossing impossible');forward.append(c)
    for a,b in combinations(U,2):
        d=gcd(a,b);p,q=a//d,b//d
        for y in row[:6]:
            ratio=F(y//gcd(d,y),Q*(p+q));need(ratio>1,'one V two U full signed amplitude impossible');reverse.append(ratio)
    need((len(forward),len(reverse))==(105,126),'complete mixed support universe')
    need(min(forward)==858481264080 and min(reverse)==F(55157421217140,55083317447977),'exact crossing margins')
    nprofiles=0
    for size in range(7,13):
        for ids in combinations(range(13),size):
            d=reduce(gcd,(row[i] for i in ids));word=tuple(sorted(gcd(d,row[i]) for i in range(13) if i not in ids))
            need(d==1 and (d,word) in profiles[13-size],'literal complete inherited profile');nprofiles+=1
    need(nprofiles==4095,'complete subset count')
    for d in range(2,7):need(not packet(row,d),'first physical denominator lower exclusions')
    need(packet(row,7)==tuple(range(1,7)) and clear(row,F(1,7))==F(1,7),'independent elementary safety; not unresolved row')
    excess=7*ceil(F(t,7))-t
    need(excess==5 and all(gcd(t,u)==1 for u in U),'exact danger-cap excess')
    ints=located_intervals(107,167)
    need(len(ints)==39 and sum(hi-lo for lo,hi in ints)==F(2553,125083),'selected located geometry')
    need(sum(open_grid_min(hi-lo,t) for lo,hi in ints)==0,'separate interval minima all vanish')
    walls=sorted({(t*x)%1 for arc in ints for x in arc})
    probes=walls+[((a+(walls[(i+1)%len(walls)]+(i==len(walls)-1)))/2)%1 for i,a in enumerate(walls)]
    counts=[]
    for alpha in probes:
        count=sum(ceil(t*hi-alpha)-floor(t*lo-alpha)-1 for lo,hi in ints)
        # Separate pointwise modular danger predicates, no interval membership oracle.
        den=t*alpha.denominator;num=alpha.numerator
        literal=sum(all(14*min((v*(num+j*alpha.denominator))%den,(-v*(num+j*alpha.denominator))%den)<den for v in (107,167)) for j in range(t))
        need(count==literal,'all walls and cells by literal modular counts');counts.append(count)
    need(len(walls)==78 and min(counts)==15 and max(counts)==23,'exact uniform translated count')
    need(sum(ceil(t*hi)-floor(t*lo)-1 for lo,hi in ints)==15,'minimum attained at alpha zero')
    d=H//(107*167)
    need(gcd(d,t)==1,'selected cofactor clock permutes all grid points')
    safe=[(j,clear(row,F(4*j+1,4*t))) for j in range(t) if clear(row,F(4*j+1,4*t))>F(1,14)]
    need(len(safe)==335 and safe[0]==(3,F(885,3872)),'literal V-quarter safe lifts')
    need(15-excess==10<=len(safe),'located overlap supplier has positive full physical consumer')
    need(t%7 and all(v%7 for v in row),'no weak endpoints on these rational lifts')
    print('STATUS PROVED actual bounded entry class; FINITE-EXACT correlated overlap supplier; row already safe')
    print('V',V,'U',U,'GROUPS',GROUPS,'P',P,'H',H)
    print('CLASS 967<=t<=1042, gcd(t,H)=1; g=1; actual rank11 equality, all4095 full profiles, all individual pair-interval credits0')
    print('ACTUAL t968 sum',sum(row),'rank',rank(relations),'edges',len(ed),'crossings',len(forward),len(reverse),'margins',min(forward),min(reverse))
    print('LOCATED_PAIR107_167 intervals39 measure2553/125083 walls78 probes156 uniform_min15 max23 individual_min_sum0 excess5 certified_safe10 literal_quarter_safe335')
    print('PHYSICAL first7 packet',packet(row,7),'clearance',clear(row,F(1,7)),'quarter_first',F(13,3872),'clearance',F(885,3872))
    print('SCOPE relative interval positions repair separate-rounding loss; no unresolved-entry or general LRC14 claim')
    print('PASS',N,'always-active exact gates; raw LF')

if __name__=='__main__':main()
