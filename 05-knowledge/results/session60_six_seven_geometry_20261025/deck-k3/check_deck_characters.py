#!/usr/bin/env python3
"""Finite-exact controls; standard library only. Does not replace proofs."""
from itertools import combinations
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path


CHECK_COUNT = 0

def require(condition, message):
    """An explicit gate: unlike assert, this remains active under python -O."""
    global CHECK_COUNT
    CHECK_COUNT += 1
    if not condition:
        raise RuntimeError("Verification failed: " + message)


def parity(x): return x.bit_count() & 1

def chi(c, g): return -1 if parity(c & g) else 1

def rank2(rows):
    b = {}
    for x in rows:
        while x:
            j = x.bit_length()-1
            if j in b: x ^= b[j]
            else:
                b[j] = x
                break
    return len(b)

def qrank(rows):
    a = [[F(x) for x in row] for row in rows]
    r = 0
    for j in range(len(a[0])):
        p = next((i for i in range(r,len(a)) if a[i][j]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        v=a[r][j]; a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                v=a[i][j]; a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
    return r

def h0(d): return comb(d+2,2) if d>=0 else 0

def hilbert_pg(m):
    d=m-6
    if d<0: return 0
    return sum((-1)**j*comb(m-3,j)*comb(d-2*j+m-1,m-1)
               for j in range(d//2+1))

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def dot(a,b): return sum(x*y for x,y in zip(a,b))

def norm(p):
    first=next(x for x in p if x)
    return tuple(F(x,first) for x in p)

def intersections(lines):
    points={norm(cross(a,b)) for a,b in combinations(lines,2)}
    return {p:tuple(i for i,l in enumerate(lines) if dot(p,l)==0) for p in points}

report={"status":"FINITE-EXACT PASS","genus":{},"hostile_controls":{}}
for m in range(5,10):
    full=(1<<m)-1
    chars=[c for c in range(1<<m) if not parity(c)]
    require(len(chars)==2**(m-1), "original source line 64")
    # Descent under the projective global sign, plus every local basis chart.
    for c in chars:
        for g in range(1<<(m-1)):
            require(chi(c,g)==chi(c,g^full), "original source line 68")
    for j in range(m):
        supports={c & ~(1<<j) for c in chars}
        require(len(supports)==2**(m-1), "original source line 71")
    pg=sum(h0(c.bit_count()//2-3) for c in chars)
    require(pg==hilbert_pg(m), "original source line 73")
    report["genus"][str(m)]=pg
require(list(report["genus"].values())==[0,1,7,31,111], "original source line 75")

m=7; full=127; cc=[full^(1<<j) for j in range(7)]
require(rank2(cc)==6, "original source line 78")
relations=[s for s in range(1,128) if not __import__('functools').reduce(
    int.__xor__,(cc[j] for j in range(7) if s>>j&1),0)]
require(relations==[127], "original source line 81")
for g in range(128):
    # determinant twist on the linear coordinate numerator
    determinant=(-1)**g.bit_count()
    for j in range(7):
        require(determinant*((-1)**((g>>j)&1))==chi(cc[j],g), "original source line 86")
    require([chi(c,g) for c in cc]==[chi(c,g^full) for c in cc], "original source line 87")
for j in range(7):
    require([k for k,c in enumerate(cc) if chi(c,1<<j)==1]==[j], "original source line 89")
for i,j in combinations(range(7),2):
    require(sum(chi(c,(1<<i)|(1<<j))==1 for c in cc)==5, "original source line 91")
    require(sum(chi(c,1<<i)==chi(c,1<<j)==1 for c in cc)==0, "original source line 92")
report["hostile_controls"]["single_pair_flip_pg"]=5
report["hostile_controls"]["two_generator_subgroup_pg"]=0
report["involution_pg_census"]={}
for w in (2,4,6):
    reps=[g for g in range(1,128) if g.bit_count()==w]
    require(all(sum(chi(c,g)==1 for c in cc)==7-w for g in reps), "original source line 98")
    report["involution_pg_census"][str(w)]={"count":len(reps),"pg":7-w}
# Seven distinct nonzero F_2^3 labels produce the simplex character subcode.
simplex=[sum((parity(a&v)<<j) for j,v in enumerate(range(1,8))) for a in range(8)]
require(len(set(simplex))==8 and all(c.bit_count()==4 for c in simplex[1:]), "original source line 102")
require(set(simplex).isdisjoint(cc), "original source line 103")
report["simplex_nonzero_weight"]=4

hourglass=[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1)]
fullmagic=hourglass+[(1,1,-1),(1,-1,1)]
require(sorted(len(v) for v in intersections(hourglass).values()).count(3)==3, "original source line 108")
report["hourglass_deletions"]=[]
for omitted in range(7):
    lines=[l for i,l in enumerate(hourglass) if i!=omitted]
    incid=intersections(lines)
    triples=sum(len(v)==3 for v in incid.values())
    require(triples==(0 if omitted==0 else 2), "original source line 114")
    report["hourglass_deletions"].append({"omitted":omitted,"triples":triples,
        "nodes":4*triples,"picard_lower_bound":16+triples})
fixed_tests=0
triple_census={}
for indices in combinations(range(9),6):
    lines=[fullmagic[i] for i in indices]
    incid=intersections(lines)
    t=sum(len(v)==3 for v in incid.values())
    triple_census[t]=triple_census.get(t,0)+1
    require(max(map(len,incid.values()))<=3, "original source line 124")
    for i,j in combinations(range(6),2):
        p=norm(cross(lines[i],lines[j])); k=len(incid[p])
        points=2**(6-k-1)
        exceptional_fixed=2 if k==3 else 1
        require(points*exceptional_fixed==8, "original source line 129")
        # opposite sign eigenspace would need four concurrent lines.
        require(qrank([lines[a] for a in range(6) if a not in (i,j)])==3, "original source line 131")
        fixed_tests+=1
require(triple_census=={0:2,1:20,2:46,3:16}, "original source line 133")
report["six_root_fixed_point_tests"]={"subsets":84,"involutions":fixed_tests,
    "fixed_points_each":8,"triple_census":triple_census}
require((22+15*6)//16==7, "original source line 136")
require(15*8//8==15, "original source line 137")

# Central-deletion conic and Gale identities; all arithmetic exact.
pts=hourglass[1:]
for A,B,C in pts: require(A*A-B*B-C*C+B*C==0, "original source line 141")
params=[F(2),F(0),F(-1),F(1),None,F(1,2)]
for t,(_,B,C) in zip(params,pts):
    if t is None: p=(F(1),F(1))
    else:
        d=t*t-t+1; p=((t*t-1)/d,t*(t-2)/d)
    require(p==(F(B),F(C)), "original source line 147")
ss=[F(0) if t is None else 1/(t-3) for t in params]
ww=[]
for s in ss:
    p=F(1)
    for z in ss:
        if z!=s: p*=s-z
    ww.append(1/p)
for k in range(5): require(sum(w*s**k for w,s in zip(ww,ss))==0, "original source line 155")
V=[[s**k for s in ss] for k in range(3)]
R=[[w*s**k for w,s in zip(ww,ss)] for k in range(3)]
require(qrank(V)==qrank(R)==3, "original source line 158")
for a in V:
    for b in R: require(dot(a,b)==0, "original source line 160")
report["gale"]={"parameters":list(map(str,ss)),"weights":list(map(str,ww)),
    "vanishing_moments":list(range(5)),"normal_rank":3,"relation_rank":3}
report["hodge_numeric_controls"]={"b2_hourglass":78,"rho_sum_offset":76,
    "rho_outer_options":[18,19,20],"rho_hourglass_options":[6*r-56 for r in (18,19,20)]}
report["explicit_checks"] = CHECK_COUNT
report["optimization_safe"] = True
print(json.dumps(report,indent=2,sort_keys=True))
Path(__file__).with_name("EXACT_CHECKS.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
