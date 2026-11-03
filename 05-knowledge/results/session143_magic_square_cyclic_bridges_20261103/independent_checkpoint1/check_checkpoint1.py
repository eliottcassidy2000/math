#!/usr/bin/env python3
"""Independent standard-library audit, deliberately not importing producer code.
Positive: python3 [-O] check_checkpoint1.py
Negative: add --mutate-spectrum.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
import argparse, json, math

p=argparse.ArgumentParser()
p.add_argument('--mutate-spectrum',action='store_true')
args=p.parse_args()
checks=0
def demand(v,msg):
    global checks
    checks+=1
    if not v: raise RuntimeError(msg)

# Graphs are sets of ordered arcs, relabeling acts directly on those pairs.
V=range(5); pairs=list(combinations(V,2)); perms=list(permutations(V))
def graph(m): return frozenset((a,b) if m>>i&1 else (b,a) for i,(a,b) in enumerate(pairs))
def encode(E): return sum(1<<i for i,(a,b) in enumerate(pairs) if (a,b) in E)
def moved(E,s): return frozenset((s[a],s[b]) for a,b in E)
def converse(E): return frozenset((b,a) for a,b in E)
remaining=set(range(1024)); classes=[]; canonical={}
while remaining:
    m=min(remaining); E=graph(m)
    orbit={encode(moved(E,s)) for s in perms}
    demand(min(orbit)==m,'canonical minimum')
    remaining-=orbit
    for k in orbit:canonical[k]=m
    classes.append((m,E,len(orbit)))
demand(len(classes)==12,'12 unlabeled tournaments')

def mm(A,B): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def charpoly(A):
    n=len(A); power=[[int(i==j) for j in range(n)] for i in range(n)]; traces=[]
    for k in range(1,n+1):
        power=mm(power,A); traces.append(sum(power[i][i] for i in range(n)))
    coeff=[F(1)]
    for k in range(1,n+1): coeff.append(-sum(coeff[k-i]*traces[i-1] for i in range(1,k+1))/k)
    demand(all(x.denominator==1 for x in coeff),'integral Newton coefficients')
    return tuple(map(int,coeff))
rows=[]
for m,E,orbit in classes:
    A=[[int((a,b) in E) for b in V] for a in V]
    aut=[s for s in perms if moved(E,s)==E]
    anti=[s for s in perms if moved(E,s)==converse(E)]
    cp=charpoly(A); cycles=[]
    for k in (3,4,5):
        count=sum(all((s[i],s[(i+1)%k]) in E for i in range(k)) for s in permutations(V,k))
        demand(count%k==0,'cycle rotation divisibility');cycles.append(count//k)
    demand(cp==(1,0,0,-cycles[0],-cycles[1],-cycles[2]),'cycle charpoly identity')
    demand(len(aut)*orbit==120,'orbit stabilizer')
    demand(len(aut)%2==1,'odd tournament automorphism group')
    demand(len(anti) in (0,len(aut)),'anti coset size')
    demand(all(all(s[s[i]]==i for i in V) for s in anti),'order-five anti involutions')
    rows.append({'mask':m,'dual':canonical[encode(converse(E))],'score':sorted(map(sum,A)),
                 'aut':len(aut),'anti':len(anti),'cp':list(cp),'cycles':cycles})
selfdual=sum(r['mask']==r['dual'] for r in rows)
spectra={tuple(r['cp']) for r in rows}
demand(selfdual==8 and (12+selfdual)//2==10,'self-converse and quotient counts')
demand(len(spectra)==(8 if args.mutate_spectrum else 9),'nine distinct adjacency spectra')
regular=[r for r in rows if r['score']==[2]*5]
demand(len(regular)==1 and regular[0]['aut']==5 and regular[0]['anti']==5,'unique regular class and D10 size')
demand(regular[0]['cp']==[1,0,0,-5,-5,-2],'regular characteristic polynomial')
demand(len({(tuple(r['cp']),tuple(r['score'])) for r in rows})==12,'score restores injectivity')
# Family action independently computed in arcs.
sets=[frozenset(s) for s in product((1,4),(2,3))]
family={frozenset((a,b) for a in V for b in V if (b-a)%5 in S) for S in sets}
stabilizer=[s for s in perms if {moved(E,s) for E in family}==family]
family_kernel=[s for s in stabilizer if all(moved(E,s)==E for E in family)]
demand(len(family)==4 and len(stabilizer)==20 and len(family_kernel)==5,'AGL family action orders')
demand(set(stabilizer)=={tuple((a*i+b)%5 for i in V) for a in range(1,5) for b in V},'family stabilizer exactly affine')

# Exact Q[z]/Phi5 representation using 4 coefficient vectors.
def reduce5(v):
    out=[F(0)]*4
    for n,c in enumerate(v):
        r=n%5
        if r==4:out=[x-c for x in out]
        else:out[r]+=c
    return tuple(out)
def ringmul(a,b):
    v=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):v[i+j]+=x*y
    return reduce5(v)
def ringadd(a,b):return tuple(x+y for x,y in zip(a,b))
q=(1,3,4,2,1)
for r in range(1,5):
    ar=reduce5([int(i==r)+int(i==2*r) for i in range(2*r+1)])
    value=(F(0),)*4
    for c in reversed(q):value=ringadd(ringmul(value,ar),(F(c),F(0),F(0),F(0)))
    demand(value==(0,0,0,0),'all four exact cyclotomic roots')
bs=[]
for name,d,c in [('5A',4371,21),('5B',4371,-4)]:
    ms=[]
    for r in range(5):
        v=[F(d)]+[F(0)]*4
        for j in range(1,5):v[(-r*j)%5]+=c
        v=reduce5(v)
        demand(v[1:]==(0,0,0) and v[0].denominator==1 and v[0]%5==0,'exact Fourier multiplicities')
        ms.append(int(v[0]/5))
    demand(sum(ms)==d and 2*ms[0]-2*ms[1]==2*c,'dimension and trace')
    bs.append({'class':name,'multiplicities':ms,'trace':2*c})
demand(bs[0]['multiplicities']==[891,870,870,870,870],'5A multiplicities')
demand(bs[1]['multiplicities']==[871,875,875,875,875],'5B multiplicities')
demand([sum(math.gcd(i,n)==1 for i in range(1,n+1)) for n in (5,17,85)]==[4,16,64],'coprime-conductor degree obstruction')


print(json.dumps({'status':'PASS','checks':checks,'tournament_classes':rows,'self_dual':selfdual,'spectra':len(spectra),'baby_monster_cited_specialization':bs,'scope':'Independent n=5 graph/C5 operator algebra audit; Baby Monster character and construction facts remain CITED, not locally rebuilt.'},indent=2,sort_keys=True))
