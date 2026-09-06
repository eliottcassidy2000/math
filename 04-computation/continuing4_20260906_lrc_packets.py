"""Actual small-scale decoder packet controls, and their inherited grid closure.
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
LEAVES=(179,191,251,257,263)
P=prod(LEAVES)
V=tuple(sorted((P,)+tuple(CENTER*P//q for q in LEAVES)))
GROUPS=(81,101,125,128,161,209,221)
H=prod(GROUPS)
U=tuple(sorted(H//q for q in GROUPS))
SCALES=(1369,1373,1583)
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

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--profiles',type=Path,default=Path(__file__).with_name('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'))
    args=parser.parse_args();raw=args.profiles.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==PROFILE_SHA,'exact inherited full-profile data pin')
    data=json.loads(raw)
    profiles={int(k):{(c,tuple(w)) for c,w in v['profiles']} for k,v in data['levels'].items()}
    need(V==(183050675309,187324231931,191802102017,252054071237,268951550873,580028043449),'literal primitive six-star')
    need(U==(4404519504000,4657410576000,6045955344000,7604678206125,7787190483072,9637611984000,12017269264000),'literal seven-cofactor shape')
    need(P==580028043449 and H==973398810384000,'exact cofactor products')
    for xs in (LEAVES,GROUPS):
        for a,b in combinations(xs,2):need(gcd(a,b)==1,'pairwise coprime cofactor groups')
    need(gcd(CENTER,P)==1 and reduce(gcd,V)==reduce(gcd,U)==1,'primitive shapes')
    need(gcd(prod(V),prod(U))==1,'disjoint component prime supports')
    need(1 not in V and 1 not in U and min(V)>3*Q//28 and max(U)>28,'unitless high-minimum reduced-scale context')
    need(all(v%2 for v in V) and {u%2 for u in U}=={0,1},'primitive parity scopes')
    dmin=min(gcd(a,b) for a,b in combinations(V,2))
    threshold=max(Q//dmin+1,Q*(GROUPS[-1]+GROUPS[-2])//min(V)+1)
    need(dmin==712259437 and threshold==1334,'exact sufficient actual-entry threshold')
    need(97096*sum(V)+sum(U)<Q*Q,'entire bounded-scale class fits physical box')
    need(components(6,graph(V))==(tuple(range(6)),),'actual six component connectivity')
    need(components(7,graph(U))==(tuple(range(7)),),'actual seven component connectivity')
    need(len(graph(V))==5 and len(graph(U))==8,'complete internal strict-atlas edge counts')
    for size in range(7,13):
        need((1,(1,)*(13-size)) in profiles[13-size],'all-one complement word is genuinely in inherited bank')
    divisor_witness=[]
    for d in range(2,29):
        i=next(i for i,u in enumerate(U) if u%d==0)
        need(U[i]%d==0,'small denominator killed by literal U divisor')
        need(not packet(U,d),'full numerator bank excluded at every denominator through28')
        divisor_witness.append(i)
    a29=packet(V,29);b29=packet(U,29)
    need(a29==(3,7,8,9,11,13,14,15,16,18,20,21,22,26),'complete V29 bank')
    need(b29==(1,2,4,5,24,25,27,28),'complete U29 first bank')
    need(not packet(U,30) and bool(packet(U,31)),'first two useful U clocks are exactly29 and31')
    joins29=tuple(sum((r*t)%29 in a29 for r in b29) for t in range(29))
    need(joins29==(0,0,2,4,6,2,0,4,8,8,6,8,0,6,2,2,6,0,8,6,8,8,4,0,2,6,4,2,0),'complete scale-residue compatibility profile')
    for p,q in ((1,2),(81,101)):
        need(tuple(sorted(pair_lengths(p,q)))==literal_intersection_lengths(p,q),'independent clipped overlap lengths')
    lengths=pair_lengths(81,101)
    need(len(lengths)==25 and sum(lengths)==F(389,19089),'exact selected overlap geometry')
    threshold_credit=sum(open_grid_min(ell,threshold) for ell in lengths)
    need(threshold_credit==13>6,'unpooled interval credit closes entire displayed class')
    results=[]
    for t,first in zip(SCALES,(31,29,37)):
        need(threshold<=t<=97096 and gcd(t,H)==1,'scale lies inside exact declared coprime class')
        row=tuple(t*v for v in V)+U
        need(len(set(row))==13 and reduce(gcd,row)==1 and sum(row)<Q*Q,'literal physical primitive distinct boxed row')
        ed=graph(row)
        need(components(13,ed)==(tuple(range(6)),tuple(range(6,13))),'full actual two-component graph')
        relations=[]
        for i,j in ed:
            d=gcd(row[i],row[j]);r=[0]*13;r[i],r[j]=row[j]//d,-row[i]//d
            need(sum(v*c for v,c in zip(row,r))==0 and max(map(abs,r))<=355,'literal actual bounded edge row')
            relations.append(r)
        need(rank(relations)==11,'literal rational weighted decoder rank')
        nforward=nreverse=0;fmin=None;rmin=None
        for a,b in combinations(row[:6],2):
            d=gcd(a,b)
            for y in U:
                c=d//gcd(d,y)
                need(c>Q,'every two-V one-U nonzero outside coefficient impossible')
                fmin=c if fmin is None else min(fmin,c);nforward+=1
        for a,b in combinations(U,2):
            d=gcd(a,b);p,q=a//d,b//d
            for y in row[:6]:
                delta=gcd(d,y);ratio=F(y//delta,Q*(p+q))
                need(ratio>1,'every one-V two-U nonzero outside multiple exceeds full signed amplitude')
                rmin=ratio if rmin is None else min(rmin,ratio);nreverse+=1
        need((nforward,nreverse)==(105,126),'complete mixed support orientations including zero inner coefficients')
        profile_count=0;maxima={}
        for size in range(7,13):
            maximum=0
            for ids in combinations(range(13),size):
                d=reduce(gcd,(row[i] for i in ids));word=tuple(sorted(gcd(d,row[i]) for i in range(13) if i not in ids))
                need(d==1 and (d,word) in profiles[13-size],'full hereditary profile membership, not scalar substitution')
                maximum=max(maximum,d);profile_count+=1
            maxima[size]=maximum
        need(profile_count==4095,'complete inherited subset universe')
        for d in range(2,first):need(not packet(row,d),'complete first physical denominator exclusion')
        bank=packet(row,first);need(bool(bank),'first physical denominator attained')
        for d in (29,31):
            need(bool(packet(V,d)) and bool(packet(U,d)) and gcd(t,d)==1,'both primitive packets nonempty and scale invertible at29 and31')
            need(packet(row,d)==tuple(r for r in packet(U,d) if (t*r)%d in packet(V,d)),'literal full row equals located CRT packet intersection')
        phase=F(bank[0],first);clearance=clear(row,phase)
        need(clearance>F(1,14),'literal physical strict positive phase')
        excess=7*((t+6)//7)-t
        need(all(gcd(t,u)==1 for u in U) and 0<=excess<=6,'exact seven-tail ceiling excess')
        credits=[];rounded=[]
        for p,q in combinations(GROUPS,2):
            ls=pair_lengths(p,q)
            credits.append(t*sum(ls)-len(ls))
            rounded.append(open_grid_min(sum(ls),t)+1-len(ls))
        if t in (1369,1373):
            need(sum(max(F(0),c) for c in credits)<excess,'even sum of every positive pooled pair credit fails')
            need(sum(max(0,c) for c in rounded)<=excess,'ceil-before-total-length pooling also fails')
        interval_credit=sum(open_grid_min(ell,t) for ell in lengths)
        need(interval_credit>excess,'retained individual intervals certify positive weak-safe grid count')
        safe=[]
        for j in range(t):
            x=F(2*j+1,2*t)
            if clear(row,x)>=F(1,14):safe.append((j,clear(row,x)))
        need(len(safe)>=interval_credit-excess>0,'literal actual V-half-phase lifts validate inherited grid consequence')
        need(all(c>F(1,14) for j,c in safe),'coprime7 endpoint excludes weak equality in this control')
        results.append(dict(t=t,sum=sum(row),edges=len(ed),rank=11,forward_min=fmin,reverse_min=str(rmin),profiles=profile_count,large_gcds=maxima,first_denominator=first,bank=bank,phase=str(phase),clearance=str(clearance),excess=excess,pooled_all_pairs_max=str(max(credits)),unpooled_credit=interval_credit,literal_safe_lifts=len(safe),first_lift=str(F(2*safe[0][0]+1,2*t)),first_lift_clearance=str(safe[0][1])))
    print('STATUS PROVED elementary actual-entry class plus FINITE-EXACT packet controls; inherited grid method already closes this class')
    print('V',V,'U',U,'COFACTOR_GROUPS',GROUPS)
    print('CLASS 1334<=t<=97096, gcd(t,973398810384000)=1, g=1; actual rank11 equality and all4095 profiles')
    print('DIVISOR_WITNESS_U_INDEX_BY_DENOMINATOR_2_TO_28',divisor_witness)
    print('FIRST_U29',b29,'V29',a29,'JOIN29',joins29)
    print('INHERITED_GRID selected normalized pair81,101; separate interval credit>=13>6 for every t>=1334')
    for r in results:print('ACTUAL_CONTROL',json.dumps(r,sort_keys=True))
    print('SCOPE first-clock/next-clock forcing fails in reduced scale domain; actual rows remain safely closed by the unpooled grid count')
    print('PASS',N,'always-active exact gates; actual LF and no producer imports')

if __name__=='__main__':main()
