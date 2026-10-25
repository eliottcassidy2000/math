#!/usr/bin/env python3
"""Exact structural checks; no assertions, so python -O retains all checks."""
from fractions import Fraction
from itertools import product, combinations
from pathlib import Path
import json

def require(p, msg):
    if not p: raise RuntimeError(msg)
H=((0,1,2),(6,7,8),(0,4,8),(1,4,7),(2,4,6))
V=((0,3,6),(2,5,8),(0,4,8),(3,4,5),(2,4,6))
LINES=tuple(dict.fromkeys(H+V))
HC={i for L in H for i in L}; VC={i for L in V for i in L}
def arr(E,x,z): return (E+x,E-x-z,E+z,E+z-x,E,E+x-z,E-z,E+x+z,E-x)
def rank(rows):
    a=[list(map(Fraction,r)) for r in rows]; i=0
    for j in range(len(a[0])):
        k=next((k for k in range(i,len(a)) if a[k][j]),None)
        if k is None: continue
        a[i],a[k]=a[k],a[i]; d=a[i][j]; a[i]=[v/d for v in a[i]]
        for k in range(len(a)):
            if k!=i:
                t=a[k][j]; a[k]=[u-t*v for u,v in zip(a[k],a[i])]
        i+=1
    return i
def rows(lines): return [[int(j in L) for j in range(9)]+[-1] for L in lines]
def is_sq(n):
    from math import isqrt
    return n>=0 and isqrt(n)**2==n

def main():
    require(len(HC&VC)==5 and len(HC|VC)==9,'cell incidence')
    require(len(set(H)&set(V))==2 and len(LINES)==8,'line incidence')
    pairs={tuple(sorted(e)) for L in H for e in combinations(L,2)}
    unused=set(combinations(sorted(HC),2))-pairs
    require(len(unused)==6,'unused pair count')
    require(all(sum(v in e for e in unused)==2 for v in HC-{4}),'unused pairs C6 degrees')
    require(all(not set(combinations(t,2))<=unused for t in combinations(sorted(HC),3)),'maximal linear triple system')
    edges=set(product((-1,0,1),repeat=2))
    eh=edges-{(1,-1),(-1,1)}; ev=edges-{(1,1),(-1,-1)}
    tree=eh&ev
    require(tree=={(i,j) for i,j in edges if i*j==0},'AP tree intersection')
    # Edges interpreted on disjoint shores. Connectivity + five edges proves tree.
    seen={('A',0)}
    while True:
        nxt=seen|{('B',j) for i,j in tree if ('A',i) in seen}|{('A',i) for i,j in tree if ('B',j) in seen}
        if nxt==seen: break
        seen=nxt
    require(len(seen)==6 and len(tree)==5,'spanning tree')
    for E,x,z in product(range(-2,3),repeat=3):
        require(sorted(E+i*x+j*z for i,j in eh)==sorted(arr(E,x,z)[i] for i in HC),'hourglass AP-grid map')
    ranks={'horizontal' :rank(rows(H)),'vertical':rank(rows(V)),'union':rank(rows(LINES))}
    require(ranks=={'horizontal':5,'vertical':5,'union':7},'rank check')
    # Independent linear nullspace dimension + three independent displayed vectors.
    basis=[arr(1,0,0)+(3,),arr(0,1,0)+(0,),arr(0,0,1)+(0,)]
    require(rank(basis)==3,'basis independence')
    for b in basis:
        for r in rows(LINES): require(sum(x*y for x,y in zip(b,r))==0,'basis kernel')
    grid=0
    for E,x,z in product(range(-5,6),repeat=3):
        a=arr(E,x,z); require(all(sum(a[i] for i in L)==3*E for L in LINES),'line sums')
        full_guard=x*z*(x*x-z*z)*(x*x-4*z*z)*(4*x*x-z*z)!=0
        hour_guard=x*z*(x+z)*(x-z)*(2*x+z)*(x+2*z)!=0
        require((len(set(a))==9)==full_guard,'full collision guard')
        require((len({a[i] for i in HC})==7 and len({a[i] for i in VC})==7)==(len(set(a))==9),'two-copy injectivity iff')
        require((len({a[i] for i in HC})==7)==hour_guard,'hourglass guard')
        grid+=1
    # Independent seven-cell enumeration: choose actual squared entries A,B,C,E
    # and derive G,H,I from opposite sums. No offset membership used here.
    census=[]
    for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59):
        Q={r*r%p for r in range(p)}; actual=set(); bad=None; full=0
        for x,z in product(range(p),repeat=2):
            a=tuple(v%p for v in arr(0,x,z))
            require((len({a[i] for i in HC})==7 and len({a[i] for i in VC})==7)==(len(set(a))==9),'finite field joint injectivity')
        for A,B,C,E in product(Q,repeat=4):
            if (A+B+C-3*E)%p: continue
            G,Hh,I=(2*E-C)%p,(2*E-B)%p,(2*E-A)%p
            if not {G,Hh,I}<=Q: continue
            h=(A,B,C,E,G,Hh,I)
            if len(set(h))!=7: continue
            key=(E,(A-E)%p,(C-E)%p); actual.add(key)
            a=tuple(v%p for v in arr(*key))
            if len(set(a))==9 and set(a)<=Q: full+=1
            elif bad is None: bad={'parameters':key,'array':a,'missing_values':[a[3],a[5]],'missing_are_squares':[a[3] in Q,a[5] in Q],'nine_distinct':len(set(a))==9}
        predicted=set()
        for E in Q:
            D={t for t in range(p) if (E+t)%p in Q and (E-t)%p in Q}
            for x,z in product(D,repeat=2):
                if (x+z)%p in D:
                    a=tuple(v%p for v in arr(E,x,z))
                    if len({a[i] for i in HC})==7: predicted.add((E,x,z))
        require(actual==predicted,'independent finite field census '+str(p))
        census.append({'prime':p,'distinct_hourglasses':len(actual),'completing_distinct_magic':full,'first_noncompletion':bad})
    a=arr(25,24,-24)
    require(all(is_sq(a[i]) for i in HC),'integer repeated-entry control')
    require([a[3],a[5]]==[-23,73] and not is_sq(a[3]) and not is_sq(a[5]),'integer noncompletion')
    source=(367**2,1337**2,1151**2,1105**2,1057**2,809**2,1519**2)
    require(source[0]+source[6]==source[1]+source[5]==source[2]+source[4]==2*source[3],'Cain AP control')
    require(sum(source[:3])!=3*source[3] and sum(source[4:])!=3*source[3],'Cain row defect')
    out={'status':'PASS','ranks':ranks,'parameter_grid_cases':grid,'finite_field_census':census,'integer_repetition_hostile':a,'Cain_AP_only_top_defect':sum(source[:3])-3*source[3]}
    Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__': main()
