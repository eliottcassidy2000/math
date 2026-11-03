#!/usr/bin/env python3
"""Independent standard-library audit, deliberately not importing producer code.
Positive: python3 [-O] check_independent.py
Negative: add --mutate-spectrum, --mutate-moment, or --mutate-octonion.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
import argparse, json, math

p=argparse.ArgumentParser()
p.add_argument('--mutate-spectrum',action='store_true')
p.add_argument('--mutate-moment',action='store_true')
p.add_argument('--mutate-octonion',action='store_true')
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

# Coefficient-by-coefficient cubic moment in variables A,B (C cancels).
def polymul(f,g):
    h={}
    for (a,b),u in f.items():
        for (c,d),v in g.items():h[a+c,b+d]=h.get((a+c,b+d),0)+u*v
    return {m:c for m,c in h.items() if c}
moment={}
for k in range(-4,5):
    term={(1,0):F(k*k)-F(20,3),(0,1):F(k)}
    cube=polymul(polymul(term,term),term)
    for mon,c in cube.items():moment[mon]=moment.get(mon,0)+c
moment={m:c for m,c in moment.items() if c}
expected={(3,0):F(2860,3)+(1 if args.mutate_moment else 0),(1,2):F(924)}
demand(moment==expected,'quadratic cubic-centered-moment coefficients')
lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def magic(v):return len({sum(v[i] for i in line) for line in lines})==1
loshu=[3,-4,1,-2,0,2,-1,4,-3]
demand(magic(loshu),'linear sufficiency witness')
odd_cubic=[18,-22,4,-14,0,14,-4,22,-18]
demand(magic(odd_cubic),'nonquadratic hostile magic witness')
demand(sorted(odd_cubic)==sorted(k**3-15*k for k in range(-4,5)),'cubic assignment witness')
demand(not magic([k*k for k in range(2,11)]),'literal squares do not form magic grid')
demand(sum(k*k for k in range(2,11))==384,'literal square sum')
center_candidates=[]
for j in range(-4,5):
    if not j:continue
    a=F(20-3*j*j,6*j); vals=[(a+k)**2 for k in range(-4,5)];mu=sum(vals)/9
    demand(vals[j+4]==mu,'center equation candidate')
    demand(max(vals)>2*mu,'candidate violates nonnegative complement')
    center_candidates.append(str(a))
# Nonzero real square-entry magic squares CAN have a real-square total:
# all entries 1, total 3=(sqrt(3))^2. The arithmetic obstruction needs Q or Z.
demand(magic([1]*9) and sum([1]*3)==3,'real-square arithmetic-boundary witness')

# Every centered 3x3 magic grid has opposite-pair values equivalent to
# {u,v,u+v,u-v}; solve all signed bijections exactly, including degeneracies.
cubic_survivors=[]; cubic_reject=Counter()
for inds in permutations(range(1,5)):
    i,j,h,l=inds
    for signs in product((-1,1),repeat=4):
        si,sj,sh,sl=signs
        equations=[(si*i+sj*j-sh*h,si*i**3+sj*j**3-sh*h**3),
                   (si*i-sj*j-sl*l,si*i**3-sj*j**3-sl*l**3)]
        if any(a==0 and b!=0 for a,b in equations):
            cubic_reject['constant']+=1;continue
        candidates={F(-b,a) for a,b in equations if a}
        demand(bool(candidates),'no entirely vacuous cubic branch')
        if len(candidates)!=1:
            cubic_reject['inconsistent']+=1;continue
        t=next(iter(candidates))
        demand(t==-15,'unique odd-cubic real parameter')
        cubic_survivors.append({'indices':inds,'signs':signs,'t':str(t)})
demand(len(cubic_survivors)==8 and dict(cubic_reject)=={'constant':88,'inconsistent':288},'384 cubic branch census')
for m in range(1,21):
    n=2*m+1;q0=m*(m+1)
    for A,B in [(1,0),(1,1),(2,-3),(-3,4)]:
        vals=[F(A*k*k+B*k) for k in range(-m,m+1)];mean=sum(vals)/n
        lhs=sum((x-mean)**3 for x in vals)
        rhs=F(n*A*q0*(4*q0-3),945)*((4*q0-15)*A*A+63*B*B)
        demand(lhs==rhs,'general odd-packet formula exact samples')
for a,b in [(1,F(1,3)),(1,F(-1,3))]:
    vals=sorted(a*k*k+b*k for k in (-1,0,1))
    demand(vals[0]+vals[2]==2*vals[1],'sharp three-argument boundary')

# Enumerate affine maps of F3^2 and determine their actual permutation orders.
points=list(product(range(3),repeat=2));index={v:i for i,v in enumerate(points)}
def order(s):
    q=tuple(range(len(s)))
    for k in range(1,100):
        q=tuple(s[i] for i in q)
        if q==tuple(range(len(s))):return k
    raise RuntimeError('permutation order bound')
affine=[]
for a,b,c,d in product(range(3),repeat=4):
    if (a*d-b*c)%3==0:continue
    for u,v in points:
        affine.append(tuple(index[((a*x+b*y+u)%3,(c*x+d*y+v)%3)] for x,y in points))
orders=Counter(map(order,affine))
demand(len(set(affine))==432 and 9 not in orders,'no affine order-nine element')
cyc=[4,8,3,7,2,6,10,5,9]
demand(all(2+(x-2+4)%9==y for x,y in zip(cyc,cyc[1:]+cyc[:1])),'stated +4 cycle')
for n in range(9):
    u,v=n%3,n//3;n2=(n+4)%9
    demand((n2%3,n2//3)==((u+1)%3,(v+1+u+2*u*u)%3),'carry polynomial')

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

# Octonions independently built from Hamilton quaternion arithmetic, not the
# producer's recursive signed-basis code. Its doubling convention is
# (a,b)(c,d)=(ac-d conjugate(b), conjugate(a)d+cb).
def qa(a,b):return tuple(x+y for x,y in zip(a,b))
def qn(a):return tuple(-x for x in a)
def qc(a):return (a[0],-a[1],-a[2],-a[3])
def qm(a,b):
    w,x,y,z=a;W,X,Y,Z=b
    return (w*W-x*X-y*Y-z*Z,w*X+x*W+y*Z-z*Y,w*Y-x*Z+y*W+z*X,w*Z+x*Y-y*X+z*W)
def om(x,y):
    a,b=x[:4],x[4:];c,d=y[:4],y[4:]
    return qa(qm(a,c),qn(qm(d,qc(b))))+qa(qm(qc(a),d),qm(c,b))
# Index 3 is -k for this doubling convention.
Q=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,-1)]
zero=(0,0,0,0);E=[v+zero for v in Q]+[zero+v for v in Q]
def scale(s,x):return tuple(s*a for a in x)
def label(x):
    for s in (-1,1):
        for k,e in enumerate(E):
            if x==scale(s,e):return s,k
    raise RuntimeError('not a signed unit')
def f(a,b):
    x=[a>>i&1 for i in range(3)];y=[b>>i&1 for i in range(3)]
    return (sum(x[i]*y[j] for i in range(3) for j in range(i,3))+x[1]*x[2]*y[0]+x[0]*x[2]*y[1]+x[0]*x[1]*y[2])%2
for a,b in product(range(8),repeat=2):
    expected=((-1)**f(a,b),a^b)
    if args.mutate_octonion and (a,b)==(1,2):expected=(-expected[0],expected[1])
    demand(label(om(E[a],E[b]))==expected,'Hamilton doubling matches cochain')
def det(a,b,c):
    return sum((a>>i&1)*(b>>j&1)*(c>>k&1) for i,j,k in permutations(range(3)))%2
assoc=Counter()
for a,b,c in product(range(8),repeat=3):
    left=om(om(E[a],E[b]),E[c]);right=om(E[a],om(E[b],E[c]))
    diff=tuple(x-y for x,y in zip(left,right)); s,k=label(left)
    demand(label(right)[1]==k==a^b^c,'associator common label')
    coefficient=label(left)[0]-label(right)[0];assoc[coefficient]+=1
    demand((coefficient!=0)==bool(det(a,b,c)),'associator nonzero iff independent')
    demand(left==scale((-1)**det(a,b,c),right),'associator ratio is determinant sign')
demand(dict(assoc)=={0:344,2:84,-2:84},'associator coefficient census')
# A gauge changing e7 flips the numerical associator coefficient at (1,2,4),
# while leaving the ratio -1 invariant. Raw coefficients are not gauge invariant.
h=lambda a:int(a==7)
fg=lambda a,b:f(a,b)^h(a)^h(b)^h(a^b)
def coefficient(ff,a,b,c):return (-1)**(ff(a,b)^ff(a^b,c))-(-1)**(ff(b,c)^ff(a,b^c))
demand(coefficient(f,1,2,4)==2 and coefficient(fg,1,2,4)==-2,'gauge ratio-versus-difference firewall')
perimeter=[0,1,2,5,8,7,6,3]
for a in range(8):
    s,b=label(om(E[4],E[a]));t,c=label(om(E[4],E[b]))
    demand(c==a and s*t==-1,'left unit multiplication square -I')
    demand(perimeter[b]==8-perimeter[a],'perimeter lift halfturn')
demand(om(E[4],E[0])==E[4]!=E[0],'left multiplication is not unital automorphism')
subspaces=[]
for mask in range(256):
    S={a for a in range(8) if mask>>a&1}
    if 0 in S and all(a^b in S for a,b in product(S,repeat=2)):subspaces.append(S)
demand(Counter(map(len,subspaces))=={1:1,2:7,4:7,8:1},'16-subspace lattice is not eight-point cube')

print(json.dumps({'status':'PASS','checks':checks,'tournament_classes':rows,'self_dual':selfdual,'spectra':len(spectra),
    'odd_cubic_survivors':cubic_survivors,'odd_cubic_rejected':dict(cubic_reject),'affine_orders':dict(sorted(orders.items())),'center_candidates':center_candidates,
    'moment_coefficients':{str(k):str(v) for k,v in moment.items()},'baby_monster_cited_specialization':bs,
    'octonion_associator_census':dict(sorted(assoc.items())),
    'hostiles':'transitive nilpotent; cospectral collision; nonlinear cubic; real-square arithmetic boundary; gauge coefficient flip; nonunital left action',
    'scope':'Independent finite algebra and proof audit. External character/construction facts remain cited. No GAP reconstruction, explicit Baby Monster matrices, or open-problem solution.'},indent=2,sort_keys=True))
