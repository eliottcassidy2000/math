#!/usr/bin/env python3
"""New local controls only; not a Collatz convergence search."""
from itertools import permutations, combinations, product
from collections import Counter
from pathlib import Path
import json

def check(x):
    if not x: raise RuntimeError('failed control')

def raw(x): return x//2 if x%2==0 else 3*x+1

def odd(x,s=1):
    z=3*x+s;k=0
    while z%2==0: z//=2;k+=1
    return z,k

def canon(edges):
    return min(tuple(int((p[i],p[j]) in edges) for i in range(4) for j in range(i+1,4)) for p in permutations(range(4)))

def hp(edges):
    return sum(all((p[i],p[i+1]) in edges for i in range(3)) for p in permutations(range(4)))

def triangles(e):
    return sum(any(all((p[i],p[(i+1)%3]) in e for i in range(3)) for p in permutations(c)) for c in combinations(range(4),3))

classes={};representatives={};pair=list(combinations(range(4),2))
for bits in product((0,1),repeat=6):
    e={ (a,b) if bit else (b,a) for (a,b),bit in zip(pair,bits)}
    h=hp(e);t=triangles(e);check(h==1+2*t)
    key=canon(e)
    if key not in classes: classes[key]={'H':h,'triangles':t,'labeled':0,'score':sorted(sum((v,w) in e for w in range(4)) for v in range(4)), 'aut':sum({(p[a],p[b]) for a,b in e}==e for p in permutations(range(4)))}
    classes[key]['labeled']+=1
    representatives[key]=e
check(len(classes)==4)
completions=Counter()
for bits in product((0,1),repeat=3):
    e={(0,1),(1,2),(2,3)}
    for (a,b),bit in zip([(0,2),(0,3),(1,3)],bits): e.add((a,b) if bit else (b,a))
    completions[tuple(classes[canon(e)]['score'])]+=1
check(sorted(completions.values())==[1,1,1,5])

# Punctured directed path profiles; class fiber weights, not raw label counts.
punctures=[]
for lengths in [(1,3),(2,2),(3,1),(1,1,2),(1,1,1,1)]:
    starts=[];pos=0;forced=set()
    for size in lengths:
        starts.append((pos,pos+size))
        forced.update((i,i+1) for i in range(pos,pos+size-1));pos+=size
    gaps=[p for p in pair if p not in forced]
    counts=Counter()
    for bits in product((0,1),repeat=len(gaps)):
        e=forced|{(a,b) if bit else (b,a) for (a,b),bit in zip(gaps,bits)}
        counts[canon(e)]+=1
    ordered=sorted(classes,key=lambda k:(classes[k]['H'],classes[k]['score']))
    weights=[counts[k] for k in ordered]
    for key in ordered:
        e=representatives[key]
        C=sum(all((order[i],order[i+1]) in e for lo,hi in starts for i in range(lo,hi-1)) for order in permutations(range(4)))
        check(C==counts[key]*classes[key]['aut'])
        check(counts[key]%2==0)
    punctures.append({'composition':lengths,'class_weights':weights,'total':sum(weights)})
check(punctures[0]['class_weights']==[4,2,2,8])
check(punctures[1]['class_weights']==[6,2,2,6])

# Word coefficients and legal ordinary shortcut cycles, exact rational signs.
cycle_records=[]
for seed in [1,-1,-5,-17]:
    xs=[];x=seed;A=1;B=0;L=0;r=0
    while x not in xs:
        xs.append(x);bit=x%2
        if bit: A*=3;B=3*B+2**L;r+=1
        L+=1
        x=(3*x+1)//2 if bit else x//2
        check(len(xs)<100)
    check(x==seed);D=2**L-A;check(D*seed==B)
    cycle_records.append({'seed':seed,'shortcut_cycle':xs,'L':L,'r':r,'B':B,'denominator':D})
check(odd(5)==(1,4));check(odd(1,-1)==(1,1))
check(odd(21)==(1,6))
# Equal residue cannot recognize the literal hole.
for q in range(1,100): check((5+q)%q==5%q and 5+q!=5)
# Exact signed accelerated table transport.
for s in [-1,1]:
    for x in range(1,501,2):
        y,k=odd(x,s);m=4*x+s;next_m=4*y+s
        check(2**k*(next_m-s)==3*m+s)
# Raw arcs on positive/negative cycles, with prescribed punctures.
check(raw(5)==16 and raw(-2)==-1)
check(all((x,raw(x))!=(5,16) for x in [1,4,2]))
check(any((x,raw(x))==(-2,-1) for x in [-1,-2]))
# Exact monoid transport for the signed edge-code coordinates.
for x in range(1,21):
    for y in range(1,21):
        for s,t in product([-1,1],repeat=2):
            z=4*x*y+t*x+s*y
            check(z>=1 and (4*x+s)*(4*y+t)==4*z+s*t)
# Transported raw dynamics on odd table codes.
def phi(h):
    return {1:lambda:(h+1)//2,3:lambda:3*h-2,5:lambda:3*h+2,7:lambda:(h-1)//2}[h%8]()
for h in range(3,1001,2):
    sign=1 if h%4==1 else -1
    x=(h-sign)//4
    y=x//2 if x%2==0 else 3*x+sign
    check(phi(h)==4*y+sign and phi(h)%4==h%4)
check([phi(h) for h in [5,17,9]]==[17,9,5])
check([phi(h) for h in [3,7,21,65,33]]==[7,3,65,33,17])
out={'scope':'64 labeled four-vertex tournaments; eight directed-P4 completions; five punctured path profiles; four explicit cycles; 500 signed table transports; 1600 code products; 499 transported raw steps; literal-hole controls only','four_classes':sorted(classes.values(),key=lambda x:(x['H'],x['score'])),'path_completions':[{'score':list(k),'count':v} for k,v in sorted(completions.items())],'punctured_path_profiles':punctures,'cycles':cycle_records,'all_passed':True}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
